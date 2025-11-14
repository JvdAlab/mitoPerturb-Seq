#!/usr/bin/env python3

"""
Create working directories for pairwise DamID-seq samples and symlink FASTQ files

This script prepares the directory structure for DamID-seq pairwise peak calling
by creating all combinations of Dam-fusion vs Dam-only samples.

Usage:
    python 00-1_prepare_workdir.py -i samplesheet.csv -o /path/to/workdir

Samplesheet format (CSV):
    sample,sample_type,fastq1,fastq2
    ATF4_rep1,Dam-fusion,/path/to/ATF4_rep1_R1.fastq.gz,/path/to/ATF4_rep1_R2.fastq.gz
    DAM_rep1,Dam-only,/path/to/DAM_rep1_R1.fastq.gz,/path/to/DAM_rep1_R2.fastq.gz
"""

import os
import csv
import argparse
from itertools import product


def parse_args(args=None):
    Description = "Create working directories for pairwise samples and symlink files"
    Epilog = "Example usage: python3 00-1_prepare_workdir.py -i samplesheet.csv -o /path/to/workdir"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        dest="SAMPLESHEET_FILE",
        required=True,
        help="Path to input samplesheet file (CSV format)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        dest="WORK_DIR",
        required=True,
        help="Path to output working directory",
    )
    return parser.parse_args(args)


def create_Dam_pairs(samplesheet_file, work_dir):
    """Create combinations of DamID-fusion and DamID-only samples"""

    dam_fusion_samples = []
    dam_only_samples = []
    sample_info = {}

    # Read samplesheet
    with open(samplesheet_file, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            sample_info[row["sample"]] = (row["fastq1"], row["fastq2"])
            if row["sample_type"] == "Dam-fusion":
                dam_fusion_samples.append(row["sample"])
            elif row["sample_type"] == "Dam-only":
                dam_only_samples.append(row["sample"])

    print(f"Found {len(dam_fusion_samples)} Dam-fusion samples: {dam_fusion_samples}")
    print(f"Found {len(dam_only_samples)} Dam-only samples: {dam_only_samples}")

    # Create all pairwise combinations
    sample_pairs = list(product(dam_fusion_samples, dam_only_samples))
    print(f"\nCreating {len(sample_pairs)} pairwise directories...")

    # Create the working directory if it doesn't exist
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
        print(f"Created working directory: {work_dir}")

    # Create directory for each pair and symlink FASTQs
    for pair in sample_pairs:
        pair_dir = os.path.join(work_dir, f"{pair[0]}_vs_{pair[1]}")

        if not os.path.exists(pair_dir):
            os.makedirs(pair_dir)
            print(f"  Created: {pair[0]}_vs_{pair[1]}/")

        # Symlink the fastq files for each sample in the pair
        for sample in pair:
            fastq1, fastq2 = sample_info[sample]

            # Create symlinks
            link1 = os.path.join(pair_dir, os.path.basename(fastq1))
            link2 = os.path.join(pair_dir, os.path.basename(fastq2))

            if not os.path.exists(link1):
                os.symlink(fastq1, link1)
            if not os.path.exists(link2):
                os.symlink(fastq2, link2)

    print(f"\n✓ Successfully created {len(sample_pairs)} pairwise directories!")
    print(f"  Working directory: {work_dir}")


def main(args=None):
    args = parse_args(args)
    create_Dam_pairs(samplesheet_file=args.SAMPLESHEET_FILE, work_dir=args.WORK_DIR)


if __name__ == "__main__":
    main()
