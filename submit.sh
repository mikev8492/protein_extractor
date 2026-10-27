#!/usr/bin/env bash
# Usage: bash submit.sh QUERY.faa GENOME_FOLDER [run_name]
# Submit a SLURM array, extract sequences, then collect the results.
set -euo pipefail

# Sets Project directory and passes it to child processes (running sbatch below)
PROJECT=${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}
export PROJECT_ROOT="$PROJECT"

# These modes are called automatically by SLURM.
if [[ ${1:-} == --extract || ${1:-} == --collect ]]; then
    module load miniforge/25.9
    eval "$(conda shell.bash hook)"
    conda activate "${CONDA_ENV:-protein_extractor}"
fi

# =================================
# Protein Extraction
# =================================
if [[ ${1:-} == --extract ]]; then
    RUN_DIR=$2
    THREADS=${SLURM_CPUS_PER_TASK:-4}

    # Match this array task to one row in the genome manifest.
    IFS=$'\t' read -r GENOME LABEL INDEX_NAME < <(
        sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RUN_DIR/manifest.tsv"
    )
    PREFIX="$RUN_DIR/genomes/$LABEL"
    INDEX="$PROJECT/output/indexes/$INDEX_NAME"

    # Build a genome index once, then reuse it for future queries.
    (
        flock 9
        if [[ ! -s "$INDEX" ]]; then
            miniprot -t "$THREADS" -d "$INDEX.tmp" "$GENOME"
            mv "$INDEX.tmp" "$INDEX"
        fi
    ) 9>"$INDEX.lock"

    # -----------------------------------------------
    # Align the query and save detailed QC output.
    # -----------------------------------------------
    miniprot -t "$THREADS" -I --outn=3 --outs=0.8 \
        --gff --aln --trans \
        "$INDEX" "$RUN_DIR/query.faa" > "$PREFIX.miniprot.txt"

    # ----------------------------------------------------------
    # Remove alignment-specific records to create a clean GFF3.
    # ----------------------------------------------------------
    sed -E '/^##(PAF|ATN|ATA|AQA|STA)([[:space:]]|$)/d' \
    "$PREFIX.miniprot.txt" > "$PREFIX.gff3"

    # -----------------------------------------------
    # Extract GFF3 feature annnotations
    # -----------------------------------------------
    if grep -q $'\tmRNA\t' "$PREFIX.gff3"; then
        gffread "$PREFIX.gff3" -g "$GENOME" \
            -x "$PREFIX.cds.fna" -y "$PREFIX.protein.faa" -S
    else
        : > "$PREFIX.cds.fna"
        : > "$PREFIX.protein.faa"
    fi

    touch "$PREFIX.done"
    exit
fi

# ======================
# Collect Results
# ======================
if [[ ${1:-} == --collect ]]; then
    python3 "$PROJECT/src/collect.py" "$2"
    exit
fi

# Normal mode: prepare a run and submit the two SLURM jobs.
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: bash submit.sh QUERY.faa GENOME_FOLDER [run_name]" >&2
    exit 1
fi

QUERY=$(realpath "$1")
GENOMES=$(realpath "$2")
RUN_NAME=${3:-$(basename "${QUERY%.*}")_$(date +%Y%m%d_%H%M%S)}
RUN_DIR="$PROJECT/output/runs/$RUN_NAME"
LOG_DIR="$PROJECT/output/logs/$RUN_NAME"

mkdir -p "$PROJECT/output/runs" "$PROJECT/output/indexes" "$PROJECT/output/logs"
mkdir "$RUN_DIR" "$RUN_DIR/genomes" "$LOG_DIR"
cp "$QUERY" "$RUN_DIR/query.faa"
python3 "$PROJECT/src/prepare.py" "$GENOMES" "$RUN_DIR/manifest.tsv"
N=$(( $(wc -l < "$RUN_DIR/manifest.tsv") - 1 ))

PARTITION=${PARTITION:-Orion}
THREADS=${THREADS:-4}
MAX_CONCURRENT=${MAX_CONCURRENT:-4}

ARRAY_JOB=$(sbatch --parsable --job-name=protein_extract \
    --partition="$PARTITION" --cpus-per-task="$THREADS" \
    --mem-per-cpu=8G --time=08:00:00 \
    --array="1-${N}%${MAX_CONCURRENT}" \
    --output="$LOG_DIR/${RUN_NAME}_%A_%a.log" \
    "$PROJECT/submit.sh" --extract "$RUN_DIR")
ARRAY_JOB=${ARRAY_JOB%%;*}

COLLECT_JOB=$(sbatch --parsable --job-name=protein_collect \
    --partition="$PARTITION" --cpus-per-task=1 \
    --mem=1G --time=00:10:00 \
    --dependency="afterany:$ARRAY_JOB" \
    --output="$LOG_DIR/${RUN_NAME}_collect_%j.log" \
    "$PROJECT/submit.sh" --collect "$RUN_DIR")

printf 'Submitted %s genome tasks (job %s). Collection job: %s\nResults: %s\n' \
    "$N" "$ARRAY_JOB" "$COLLECT_JOB" "$RUN_DIR" "$LOG_DIR"
