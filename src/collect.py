#!/usr/bin/env python3
"""Combine extracted sequences and summarize one protein extraction run."""

import csv
import sys
from pathlib import Path
from typing import TextIO


def copy_fasta(source: Path, destination: TextIO, label: str) -> int:
    """Append FASTA sequences with a genome label and return their count."""
    if not source.exists():
        return 0

    count = 0
    with source.open(encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                destination.write(f">{label}|{line[1:]}")
                count += 1
            else:
                destination.write(line)
    return count


def count_hits(gff: Path) -> int:
    """Count mRNA features (candidate loci) in a miniprot GFF3 file."""
    if not gff.exists():
        return 0

    with gff.open(encoding="utf-8") as file:
        return sum("\tmRNA\t" in line for line in file if not line.startswith("#"))


def summarize_genome(item: dict, run_dir: Path, proteins: TextIO, cds: TextIO) -> list:
    """Collect one genome's sequences and return a row for the run summary."""
    label = item["label"]
    prefix = run_dir / "genomes" / label

    if not Path(f"{prefix}.done").exists():
        return [item["genome"], label, 0, 0, "incomplete"]

    hits = count_hits(Path(f"{prefix}.gff3"))
    protein_count = copy_fasta(Path(f"{prefix}.protein.faa"), proteins, label)
    copy_fasta(Path(f"{prefix}.cds.fna"), cds, label)

    if hits == 0:
        status = "no_hit"
    elif protein_count == 0:
        status = "no_translation"
    else:
        status = "extracted"

    return [item["genome"], label, hits, protein_count, status]


def collect_results(run_dir: Path) -> None:
    """Combine per-genome FASTAs and write a one-row-per-genome summary."""
    with (run_dir / "manifest.tsv").open(newline="", encoding="utf-8") as file:
        genomes = list(csv.DictReader(file, delimiter="\t"))

    with (run_dir / "all_proteins.faa").open("w", encoding="utf-8") as proteins, \
        (run_dir / "all_cds.fna").open("w", encoding="utf-8") as cds:
        summary = []
        for item in genomes:
            summary.append(summarize_genome(item, run_dir, proteins, cds))

    with (run_dir / "summary.tsv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "label", "candidate_loci", "protein_sequences", "status"])
        writer.writerows(summary)

    total = sum(row[3] for row in summary)
    print(f"Collected {total} protein sequences from {len(genomes)} genomes")
    print(f"Output: {run_dir}")


def main() -> None:
    """Read the run directory from the command line and collect its results."""
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 src/collect.py RUN_DIR")

    collect_results(Path(sys.argv[1]))


if __name__ == "__main__":
    main()
