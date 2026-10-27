# Protein Extractor -  Code Walkthrough

This guide explains the **actual control flow** of `submit.sh`, `src/prepare.py`, and
`src/collect.py`. The examples below are *pseudocode*: they describe the program's
logic rather than providing executable Bash or Python.

## 1. What the pipeline does

**Inputs**

- `QUERY.faa`: a FASTA intended to contain **one protein sequence** (not explicitly checked by Bash).
- `GENOME_FOLDER`: a directory containing genome FASTA files (`.fa`, `.fna`,
  `.fasta`), optionally nested in genus directories.
- `[run_name]`: an optional name for this extraction run.

**Outputs**

- A GFF3 file, CDS FASTA, and translated protein FASTA for each genome.
- `all_cds.fna` and `all_proteins.faa` combining results across genomes.
- `summary.tsv` reporting each genome's candidate count and completion status.

The pipeline **extracts candidate sequences only**. It does not align the
recovered proteins to each other, determine orthology, or assess gene function.

### Overall workflow

```text
USER runs: bash submit.sh QUERY.faa GENOME_FOLDER [run_name]
    |
    v
submit.sh (submission mode)
    |-- Read query and genome-folder arguments
    |-- Create output folder and copy query
    |-- Run prepare.py -> manifest.tsv
    |-- Submit one SLURM array task per genome
    `-- Submit collection job dependent on the array
                      |
                      v
          submit.sh --extract (each array task)
            |-- Look up one genome in manifest.tsv
            |-- Build/reuse a miniprot genome index
            |-- Align protein query to genome with miniprot
            |-- Extract CDS/protein with gffread (if hits exist)
            `-- Mark successful completion (.done)
                      |
                      v
          submit.sh --collect (after array finishes)
            `-- Run collect.py
                 |-- Read manifest.tsv
                 |-- Combine protein and CDS FASTAs
                 `-- Write summary.tsv
```

---

## 2. `submit.sh`: one Bash entry point with three modes

The Bash script has three paths. The user calls it with a protein query and
folder of genomes. SLURM runs it again with `--extract` for each genome, then
with `--collect` after the array finishes. `set -euo pipefail` stops execution
on most failed commands and unset variables; the only explicit argument check
is the basic usage check.

### A. Submission mode — called by the user

**Command:** `bash submit.sh QUERY.faa GENOME_FOLDER [run_name]`

```text
START submit.sh

FIND the project directory (where submit.sh is located)

IF called with --extract or --collect:
    LOAD Miniforge and activate the Conda environment
    CONTINUE with the requested SLURM mode instead

CHECK the user gave 2 or 3 arguments
RESOLVE query and genome folder to absolute paths
USE the supplied run name, or make one from query name + timestamp
CREATE the log, index, run, and per-genome output folders
    (a preexisting run directory causes mkdir to stop)
COPY the query as RUN_DIR/query.faa
RUN prepare.py to discover genomes and write manifest.tsv
COUNT genomes using the manifest's line count minus its header

SET SLURM defaults unless overridden:
    PARTITION = Orion
    THREADS = 4
    MAX_CONCURRENT = 4

SUBMIT a SLURM array with one task per genome:
    RUN submit.sh --extract RUN_DIR
    REQUEST 8 GB memory per CPU and up to 8 hours
    SAVE the array job ID

SUBMIT a collection job dependent on the finished array:
    RUN submit.sh --collect RUN_DIR
    USE afterany so incomplete tasks can still be reported
    REQUEST 1 CPU, 1 GB memory, and up to 10 minutes

PRINT the job IDs and results folder
END
```

**Why a job array?** The same extraction code runs independently on each
genome. **Why `afterany`?** Collection runs after all array tasks end, even if
some fail, and records those tasks as `incomplete`.

### B. Extraction mode — called for each genome by SLURM

**Internal command:** `submit.sh --extract RUN_DIR`

```text
LOAD Miniforge and activate the Conda environment
SET THREADS from SLURM_CPUS_PER_TASK

READ the manifest row for SLURM_ARRAY_TASK_ID
    (+1 skips the column-header row)
READ GENOME, LABEL, INDEX_NAME from that row
SET the output file prefix and cached index location

LOCK the genome index using flock
IF the index does not exist or is empty:
    BUILD the miniprot index to a temporary file
    RENAME the temporary file to the shared cache location
UNLOCK the index

RUN miniprot to align the query against the indexed genome:
    USE THREADS CPUs
    REPORT up to 3 candidate loci as GFF3
    WRITE PREFIX.gff3

IF GFF3 contains at least one mRNA feature:
    USE gffread and the genome FASTA to extract CDS + protein
    WRITE PREFIX.cds.fna and PREFIX.protein.faa
ELSE:
    WRITE empty CDS and protein FASTA files

CREATE PREFIX.done to mark a successfully finished task
END
```

The cached index is reused in later runs. The `flock` block is retained
because separate queries may run concurrently against the same genome.
A `.done` file distinguishes a completed no-hit result from a task that
failed or was interrupted. `miniprot` does not confirm that a candidate is a
functional ortholog.

### C. Collection mode — called after the array

**Internal command:** `submit.sh --collect RUN_DIR`

```text
LOAD Miniforge and activate the Conda environment
RUN python3 src/collect.py RUN_DIR
END
```

---


## 3. `src/prepare.py`: discover genomes and build the task manifest

The Python standard-library modules have specific roles: `pathlib` handles file
paths, `re` cleans labels, `hashlib` generates index-cache names, `csv` writes
TSV files, and `sys` reads command-line arguments or exits on errors.

### `find_genomes(folder) -> list[Path]`

**Purpose:** Find every uncompressed genome FASTA under the supplied folder,
including files inside genus subfolders.

```text
FUNCTION find_genomes(folder):
    IF folder is not a directory: STOP with an error
    SEARCH recursively for all files under folder
    KEEP files whose extension is .fa, .fna, or .fasta
         (extension matching ignores case)
    SORT paths for a consistent task order
    IF no files found: STOP with an error
    RETURN the sorted list of genome paths
```

A directory such as `Chiroptera/Myotis/` is searched, but its folder name is not included in the output label.

### `make_label(genome) -> str`

**Purpose:** Generate a filesystem-safe output prefix from just the genome filename, without the FASTA extension.

```text
FUNCTION make_label(genome):
    TAKE the genome filename without its extension
    REPLACE characters outside letters, digits, _, ., and - with "_"
    RETURN the label
```

Example:

```text
Genome: /data/Chiroptera/Myotis/Myotis_lucifugus_GCF_001.fna
Label:  Myotis_lucifugus_GCF_001
```

### `make_index_name(genome) -> str`

**Purpose:** Determine a reusable `.mpi` cache filename for this genome file.

```text
FUNCTION make_index_name(genome):
    READ filesystem metadata for genome
    BUILD a string from:
        absolute genome path + file size + modification time
    COMPUTE SHA-256 hash of that string
    TAKE the first 16 hexadecimal characters
    APPEND ".mpi"
    RETURN the index filename
```

**Note:** This fingerprints the **path and metadata**, not the actual genome
contents. Updating the file's size or modification time changes its index name.

### `write_manifest(genomes, output) -> None`

**Purpose:** Write one line per genome so each SLURM array task knows which
input to process and where to write its output.

```text
FUNCTION write_manifest(genomes, output):
    CREATE the manifest's parent directory if needed
    CREATE an empty set called used_labels
    OPEN output as a tab-separated file
    WRITE header: genome, label, index_name

    FOR each genome in genomes:
        label = make_label(genome)
        IF label is already in used_labels:
            STOP to prevent output filename collisions
        ADD label to used_labels
        index_name = make_index_name(genome)
        WRITE genome path, label, index_name as one TSV row
```

Example manifest (`manifest.tsv`):

```text
genome\tlabel\tindex_name
/data/Chiroptera/Myotis/Myotis_lucifugus_GCF_001.fna\tMyotis_lucifugus_GCF_001\tabc123...mpi
/data/Chiroptera/Pteropus/Pteropus_vampyrus_GCF_003.fna\tPteropus_vampyrus_GCF_003\tdef456...mpi
```

The example index hashes above are illustrative, not actual calculated hashes.

### `main() -> None`

**Purpose:** Handle the command-line interface and coordinate the helpers.

```text
FUNCTION main():
    REQUIRE exactly two arguments: GENOME_FOLDER and MANIFEST_TSV
    RESOLVE GENOME_FOLDER to an absolute path
    CONVERT MANIFEST_TSV into a Path
    genomes = find_genomes(GENOME_FOLDER)
    write_manifest(genomes, GENOME_FOLDER, MANIFEST_TSV)
    PRINT the number of discovered genomes

IF this file was started as a script:
    CALL main()
```

**Called by Bash:** `python3 src/prepare.py GENOME_FOLDER MANIFEST_TSV`

---

## 4. `src/collect.py`: combine FASTAs and summarize each genome

This script uses `csv` to read/write TSV, `pathlib` for filenames, `sys` for
arguments, and `TextIO` only as a type hint for writable output files.

### `copy_fasta(source, destination, label) -> int`

**Purpose:** Append all records from one genome's FASTA to an already-open
combined FASTA, prefixing each header with the genome label.

```text
FUNCTION copy_fasta(source, destination, label):
    IF source file does not exist: RETURN 0
    count = 0
    OPEN source for reading

    FOR each line in source:
        IF line starts with ">":
            WRITE ">" + label + "|" + original header text after ">"
            INCREMENT count
        ELSE:
            COPY sequence line exactly to destination

    RETURN count
```

For example, a record headed `>MP000001` in a genome-specific FASTA becomes
`>Myotis_lucifugus_GCF_001|MP000001` in `all_proteins.faa`. It preserves
sequence text and does **not** filter, translate, or align anything.

### `count_hits(gff) -> int`

**Purpose:** Count miniprot's reported candidate models using the GFF3 `mRNA`
feature lines.

```text
FUNCTION count_hits(gff):
    IF GFF3 file does not exist: RETURN 0
    OPEN GFF3
    COUNT non-comment lines containing a tab + "mRNA" + tab
    RETURN that count
```

This is a count of **candidate loci/features** in the annotation file, not a
count of experimentally confirmed genes.

### `summarize_genome(item, run_dir, proteins, cds) -> list`

**Purpose:** Process one manifest row, append its extracted FASTA sequences,
and return one summary row.

- `item`: dictionary with `genome`, `label`, and `index_name` from the manifest.
- `proteins` and `cds`: open output streams for the combined FASTA files.

```text
FUNCTION summarize_genome(item, run_dir, proteins, cds):
    label = item["label"]
    prefix = run_dir/genomes/label

    IF prefix.done does not exist:
        RETURN [genome_path, label, 0, 0, "incomplete"]
        (do not copy any files from this incomplete task)

    hits = count_hits(prefix.gff3)
    protein_count = copy_fasta(prefix.protein.faa, proteins, label)
    copy_fasta(prefix.cds.fna, cds, label)

    IF hits == 0:
        status = "no_hit"
    ELSE IF protein_count == 0:
        status = "no_translation"
    ELSE:
        status = "extracted"

    RETURN [genome_path, label, hits, protein_count, status]
```

### `collect_results(run_dir) -> None`

**Purpose:** Read the manifest, aggregate the completed per-genome results, and
write one summary row for every original genome.

```text
FUNCTION collect_results(run_dir):
    OPEN run_dir/manifest.tsv as a tab-separated file
    READ each genome row as a dictionary

    CREATE or overwrite all_proteins.faa and all_cds.fna
    summary = empty list
    FOR each manifest row:
        row = summarize_genome(row, run_dir, protein_output, cds_output)
        ADD row to summary
    CLOSE both combined FASTA files

    CREATE or overwrite summary.tsv
    WRITE column names:
        genome, label, candidate_loci, protein_sequences, status
    WRITE each row in summary

    SUM the protein_sequences column
    PRINT total recovered protein sequences and genomes processed
    PRINT the output directory
```

### `main() -> None`

**Purpose:** Read the run directory from the command line and start collection.

```text
FUNCTION main():
    REQUIRE exactly one argument: RUN_DIR
    collect_results(Path(RUN_DIR))

IF this file was started as a script:
    CALL main()
```

**Called by Bash:** `python3 src/collect.py RUN_DIR`

---

## 5. How the scripts pass information to each other

| Stage | Reads | Writes | Next consumer |
|---|---|---|---|
| `prepare.py` | Root genome directory | `manifest.tsv` | Each SLURM array task; `collect.py` |
| `submit.sh --extract` | Manifest row, one genome, `query.faa` | `.gff3`, `.cds.fna`, `.protein.faa`, `.done` | `collect.py` |
| `collect.py` | Manifest and completed genome outputs | `all_proteins.faa`, `all_cds.fna`, `summary.tsv` | Future gene-specific analysis |

### Example output directory

```text
output/
├── indexes/
│   └── <genome-cache-key>.mpi
├── logs/
│   ├── MUC2_<jobid>_<arrayid>.log
│   └── MUC2_collect_<jobid>.log
└── runs/
    └── MUC2/
        ├── query.faa
        ├── manifest.tsv
        ├── all_proteins.faa
        ├── all_cds.fna
        ├── summary.tsv
        └── genomes/
            ├── Myotis_lucifugus_GCF_001.gff3
            ├── Myotis_lucifugus_GCF_001.cds.fna
            ├── Myotis_lucifugus_GCF_001.protein.faa
            └── Myotis_lucifugus_GCF_001.done
```

### Interpreting `summary.tsv`

| Status | Meaning in the current code |
|---|---|
| `extracted` | Task finished and at least one protein sequence was written. |
| `no_hit` | Task finished but no GFF3 `mRNA` features were found. |
| `no_translation` | Task finished and had a GFF3 `mRNA` feature, but no protein FASTA records were recovered. |
| `incomplete` | No `.done` marker; the task may have failed or been interrupted. |

The summary is **operational QC**, not a biological classification. No output
sequence is automatically accepted as the correct ortholog.

---

## 6. Other project files

### `environment.yml`

```text
DEFINE conda environment named protein_extractor
USE conda-forge and bioconda channels
INSTALL Python >= 3.10, miniprot, and gffread
```

The Python scripts otherwise use only Python's standard library. The Bash
workflow additionally expects SLURM (`sbatch`) and common shell utilities,
including `flock` for index locking.

### `.gitignore`

```text
IGNORE Python bytecode caches (__pycache__, *.pyc)
IGNORE generated output/indexes/, output/logs/, output/runs/ content
```

This keeps potentially large genomes-derived results and cluster logs out of
version control by default.

---

## 7. One-run example

```bash
bash submit.sh input/MUC2.faa /data/Chiroptera MUC2
```

1. `submit.sh` validates the inputs and creates `output/runs/MUC2/`.
2. `prepare.py` discovers every eligible FASTA under `/data/Chiroptera` and
   creates one manifest row per genome.
3. SLURM schedules array task 1 for the first genome, task 2 for the second,
   and so on (with the configured concurrency limit).
4. Each array task reuses or constructs its `.mpi` genome index, aligns MUC2
   with miniprot, extracts candidate CDS and protein sequences using gffread,
   and writes a `.done` marker on success.
5. When the array is finished, a dependent collection job runs `collect.py`.
6. `collect.py` combines all completed sequences into `all_proteins.faa` and
   `all_cds.fna`, then records one status row per genome in `summary.tsv`.
7. To run another query (for example CMAH), supply its separate protein FASTA
   and a **new** run name. Cached genome indexes can be reused.

**Scope boundary:** Any multiple-sequence alignment, candidate validation,
paralog/ortholog selection, pseudogene inference, or species-level comparison
belongs in a separate downstream workflow.
