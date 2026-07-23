#!/usr/bin/env python3
"""
mg-clust module 2: ORF prediction and coverage estimation.

Assumes execution inside the conda environment "mg-clust-module-2" (or equivalent)
where dependencies are available on PATH.

- Predicts ORFs from an assembled metagenome using FragGeneScanRs
- Converts ORF coordinates from the FragGeneScanRs .out file to BED format,
  converting 1-based start coordinates to 0-based as required by bedtools
- Builds a bedtools genome file from the BAM header and sorts the BED file to
  the same chromosome order, so bedtools coverage can run in -sorted mode
  (streams both inputs instead of loading the BAM into memory)
- Computes the number of reads per ORF using bedtools coverage
- Computes mean depth per ORF using bedtools coverage -mean
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import re
import shutil
import sys, os
import subprocess
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file

# os, subprocess, and sys are imported in utils.py, so they are available here as well

fraggenescan = "FragGeneScanRs"
bedtools = "bedtools"
samtools = "samtools"

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: ORF prediction and coverage estimation", 
        add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--assembly_file", dest="assembly_file", required=True,
        help="input assembled metagenome (fasta format)")

    parser.add_argument("--bam_file", dest="bam_file", required=True,
        help="bam input file (reads mapped to contigs)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--train_file_name", dest="train_file_name", default="illumina_1",
        help="train file name used to run FragGeneScan (default: illumina_1)")

    return parser.parse_args()

###############################################################################
# 2.2 Add sample name to coverage tables
###############################################################################

def add_sample(coverage_file: str, sample_name: str) -> None:
    try:
        with open(coverage_file, "r", encoding="utf-8") as fh_in, \
             open(f"{coverage_file}.tmp", "w", encoding="utf-8") as fh_out:
            for line in fh_in:
                fh_out.write(f"{line.strip()}\t{sample_name}\n")
        os.replace(f"{coverage_file}.tmp", coverage_file)
    except Exception as exc:
        print(f"Adding sample name to coverage file {coverage_file} failed: {exc}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    check_tools([fraggenescan, bedtools, samtools])
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(args.assembly_file, "input assembly file")
    check_file(args.bam_file, "input bam file")

    ###########################################################################
    # 3.2. Check output directory
    ###########################################################################

    if os.path.isdir(args.output_dir):
        if not args.overwrite:
            print(f"{args.output_dir} already exists; use --overwrite to overwrite")
            sys.exit(0)
        try:
            shutil.rmtree(args.output_dir)
        except Exception:
            print(f"rm -r output directory {args.output_dir} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.3. Create output directory
    ###########################################################################

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir output directory {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.4. Predict ORFs
    ###########################################################################

    try:
        run(
            [
                fraggenescan,
                "-s", args.assembly_file,
                "-o", os.path.join(args.output_dir, f"{args.sample_name}_orfs"),
                "-w", "0",
                "--unordered",  
                "-p", str(args.nslots),
                "-t", args.train_file_name
            ]
        )
    except subprocess.CalledProcessError:
        print("FragGeneScan failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Create BED file
    ###########################################################################

    out_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.out")
    bed_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.bed")

    # check that .out file was created
    check_file(out_file, "FragGeneScan output .out file")

    # Parse .out file and create BED file
    # .out format: contig_id    start    end    strand    ...
    # BED format: contig_id, start-1 (zero-based), end (one-based)

    try:
        with open(out_file, "r", encoding="utf-8") as fh_in, \
             open(bed_file, "w", encoding="utf-8") as fh_out:
            first_data_line = True
            for line in fh_in:
                if line.startswith("#") or not line.strip():
                    continue
                if first_data_line:
                    if not line.startswith(">"):
                        print(f"Unexpected format in {out_file}: \n"
                              f"first line does not start with '>'", file=sys.stderr)
                        sys.exit(1)
                    first_data_line = False
                parts = line.strip().split("\t")
                if re.match(r"^>(.+)$", parts[0]):
                    contig_id = parts[0][1:]  # remove leading '>'
                if not re.match(r"^>(.+)$", parts[0]) and len(parts) >= 3:
                    start = int(parts[0])
                    end = int(parts[1])
                    strand = parts[2]
                    orf_id = f"{contig_id}_{start}_{end}_{strand}"
                    fh_out.write(f"{contig_id}\t{start - 1}\t{end}\t{strand}\t{orf_id}\n")
    except Exception as exc:
        print(f"Creating bed file failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Build bedtools genome file from BAM header
    ###########################################################################

    # The genome file's chromosome order is derived directly from the same BAM
    # that bedtools coverage -sorted will scan below, so it is guaranteed by
    # construction to match that BAM's coordinate-sort order.
    genome_file = os.path.join(args.output_dir, f"{args.sample_name}_genome.txt")
    try:
        header = subprocess.run(
            [samtools, "view", "-H", args.bam_file],
            capture_output=True, text=True, check=False
        )
        if header.returncode != 0:
            print(f"samtools failed to read BAM header: {header.stderr}", file=sys.stderr)
            sys.exit(1)
        with open(genome_file, "w", encoding="utf-8") as fh_out:
            for line in header.stdout.splitlines():
                if line.startswith("@SQ"):
                    fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
                    fh_out.write(f"{fields['SN']}\t{fields['LN']}\n")
    except Exception as exc:
        print(f"Building bedtools genome file failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Sort BED file according to genome file order
    ###########################################################################

    # bedtools coverage -sorted requires -a and -b to follow the same chromosome
    # order as -g; this sorts the BED side, while the BAM side is already
    # guaranteed by construction since module 1's samtools sort produced it.
    sorted_bed_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs_sorted.bed")
    try:
        run([bedtools, "sort", "-i", bed_file, "-g", genome_file], stdout_path=sorted_bed_file)
    except subprocess.CalledProcessError:
        print("bedtools sort of bed file failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8. Get number of reads per ORF
    ###########################################################################

    # Run bedtools coverage -counts to get read counts per ORF
    # -sorted -g lets bedtools stream both inputs instead of loading the BAM into memory
    bedtool_reads = os.path.join(args.output_dir, f"{args.sample_name}_orfs_readscov.tsv")
    try:
        run(
            [bedtools, "coverage",
             "-a", sorted_bed_file,
             "-b", args.bam_file,
             "-counts",
             "-sorted",
             "-g", genome_file],
            stdout_path=bedtool_reads
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute read counts failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.9. Get mean coverage per ORF
    ###########################################################################

    # Run bedtools coverage -mean to get mean depth per ORF
    bedtools_mean = os.path.join(args.output_dir, f"{args.sample_name}_orfs_meancov.tsv")
    try:
        run(
            [bedtools, "coverage", "-a", sorted_bed_file, "-b", args.bam_file, "-mean",
             "-sorted", "-g", genome_file],
            stdout_path=bedtools_mean
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute mean coverage failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.10. Add sample names to mean coverage and read counts tables
    ###########################################################################

    add_sample(bedtool_reads, args.sample_name)
    add_sample(bedtools_mean, args.sample_name)

    ###########################################################################
    # 3.11. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###########################################################################
# 4. Execute main function
###########################################################################

if __name__ == "__main__":
    main()
