# Protein Extractor

A reusable SLURM pipeline to extract candidate CDS and translated protein sequences from a folder of genome assemblies using **miniprot** and **gffread**. Run it separately for any query gene. 

## Project structure

```text
protein_genome_extractor/
├── submit.sh           
├── src/
│   ├── prepare.py      
│   └── collect.py      
├── environment.yml
└── output/             
    ├── indexes/        
    ├── logs/           
    └── runs/           
```

## Code organization

`submit.sh` handles SLURM submission, array extraction, and result collection.
- `prepare.py`: `find_genomes()` locates assemblies, `make_label()` names them, `make_index_name()` identifies cached indexes, and `write_manifest()` creates the task list.
- `collect.py`: `count_hits()` reads candidate counts, `copy_fasta()` adds genome identifiers, `summarize_genome()` processes each genome, and `collect_results()` writes the combined outputs.



## Inputs

- **Protein FASTA:** supply one amino acid sequence (`query.faa`) per run (the script assumes, but does not check, this).
- **Genome folder:** uncompressed `.fa`, `.fna`, or `.fasta` files, found recursively, including in genus folders. Output labels use the genome filename.

For example:

```text
Chiroptera/
├── Myotis/
│   ├── Myotis_lucifugus_GCF_001.fna
│   └── Myotis_brandtii_GCF_002.fna
└── Pteropus/
    └── Pteropus_vampyrus_GCF_003.fna
```

## Setup

```bash
module load miniforge/25.9
conda env create -f environment.yml
conda activate protein_extractor
```

The defaults assume the `Orion` SLURM partition and `miniforge/25.9` module. If your HPC uses different names, edit the two settings in `submit.sh`. The script uses 4 CPUs per genome (8 GB per CPU) and runs up to 4 genomes concurrently by default.

## Usage


```bash
bash submit.sh query.faa /path/to/genomes [query_name]
```

Example:
```bash
bash submit.sh input/MUC2.faa /path/to/Chiroptera MUC2
```

The optional third argument names the output run. If omitted, a name based on the query filename and the current time is generated. Every run name must be new; `mkdir` stops rather than overwriting an existing run.

To change concurrency, CPUs, or SLURM partition:

```bash
MAX_CONCURRENT=2 THREADS=8 PARTITION=Orion bash submit.sh query.faa /path/to/genomes my_run
```

**You only run `submit.sh` once.** It submits a SLURM array with one task per genome, then a dependent collection job. SLURM invokes the *same* `submit.sh` using its internal `--extract` and `--collect` modes; no separate `.sbatch` files are needed. The collection job runs even when an array task fails so the summary can report incomplete genomes.

## Results

```text
output/runs/MUC2/
├── query.faa             # Copy of the exact reference sequence
├── manifest.tsv          # List of genomes and corresponding task IDs
├── summary.tsv           # One row per genome with result status
├── all_proteins.faa      # Combined predicted proteins
├── all_cds.fna           # Combined CDS nucleotide sequences
└── genomes/
    ├── Myotis_lucifugus_GCF_001.miniprot.txt
    ├── Myotis_lucifugus_GCF_001.gff3
    ├── Myotis_lucifugus_GCF_001.cds.fna
    ├── Myotis_lucifugus_GCF_001.protein.faa
    └── Myotis_lucifugus_GCF_001.done
```

The combined FASTA headers include the genome label, for example `>Myotis_lucifugus_GCF_001|MP000001`.

Possible `summary.tsv` statuses: `extracted` (at least one protein), `no_hit` (no candidate alignment), `no_translation` (hit but no translated protein), `incomplete` (task did not finish).

To inspect alignment details (example):

```bash
grep -E '^##(ATN|ATA|AQA|STA)' \
    output/runs/MUC2/genomes/Myotis_lucifugus_GCF_001.miniprot.txt \
    | less -S
```
## Notes

- `miniprot` reports up to three candidate loci with `--outn=3 --outs=0.8`; candidates are **not** automatically labeled as orthologs.
- Predicted proteins may be partial or disrupted, especially for complex genes. The pipeline does **not** determine whether a gene is functional or pseudogenized.
- Indexes are shared across runs. A short `flock` block prevents two jobs from rebuilding the same index simultaneously.
- Genome filenames (without extensions) must be unique across all genus folders to avoid output name collisions.
- Each SLURM job loads `miniforge/25.9` and activates `protein_extractor` (or `CONDA_ENV` if set).
- Check `output/logs/` and `summary.tsv` for task failures.

Software documentation: [miniprot](https://lh3.github.io/miniprot/miniprot.html) and [gffread](https://github.com/gpertea/gffread).
