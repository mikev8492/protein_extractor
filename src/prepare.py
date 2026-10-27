#!/usr/bin/env python3
"""Find genome FASTA files and prepare a SLURM task manifest."""

import csv
import hashlib
import re
import sys
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".fna", ".fasta"}


def find_genomes(folder: Path) -> list[Path]:
    """Find genome FASTA files in a folder and all its subfolders."""
    if not folder.is_dir():
        sys.exit(f"Genome folder does not exist: {folder}")

    genomes = sorted(
        path for path in folder.rglob("*")
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )
    if not genomes:
        sys.exit("No .fa, .fna, or .fasta genomes found.")

    return genomes


def make_label(genome: Path) -> str:
    """Create a filesystem-safe output label from the genome filename."""
    return re.sub(r"[^a-zA-Z0-9_.-]", "_", genome.stem)


def make_index_name(genome: Path) -> str:
    """Name a reusable miniprot index using its path, file size, and modification time."""
    stat = genome.stat()
    key = f"{genome}:{stat.st_size}:{stat.st_mtime_ns}"
    return hashlib.sha256(key.encode()).hexdigest()[:16] + ".mpi"


def write_manifest(genomes: list[Path], output: Path) -> None:
    """Write each genome's path, output label, and index name to a TSV."""
    output.parent.mkdir(parents=True, exist_ok=True)
    used_labels = set()

    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "label", "index_name"])

        for genome in genomes:
            label = make_label(genome)
            if label in used_labels:
                sys.exit(f"Two genomes produce the same output name: {label}")
            used_labels.add(label)
            writer.writerow([str(genome), label, make_index_name(genome)])


def main() -> None:
    """Read command-line inputs and create the genome task manifest."""
    if len(sys.argv) != 3:
        sys.exit("Usage: python3 src/prepare.py GENOME_FOLDER MANIFEST_TSV")

    folder = Path(sys.argv[1]).expanduser().resolve()
    output = Path(sys.argv[2])

    genomes = find_genomes(folder)
    write_manifest(genomes, output)
    print(f"Found {len(genomes)} genome files.")


if __name__ == "__main__":
    main()
