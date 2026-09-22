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

Alternatively, --precomputed_orfs_faa + --precomputed_orfs_bed can be given
together to skip FragGeneScanRs and the .out-to-BED conversion entirely, staging
the given files (plain or gzipped; gzipped input is decompressed while it is
staged) in place of those steps' outputs; --assembly_file is then unused
(coverage estimation only needs --bam_file and the ORF BED file either way). In
that case the BED's contig IDs are checked against the BAM header before bedtools
runs, since nothing else ties the two files together -- a mismatch there means the
ORFs were predicted from a different assembly, or from one prefixed with a
different module 1 --id_sep.
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
from utils import run, check_tools, check_file, gzip_file, stage_decompressed

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

    parser.add_argument("--assembly_file", dest="assembly_file", default=None,
        help="input assembled metagenome (fasta format); required unless both "
             "--precomputed_orfs_faa and --precomputed_orfs_bed are given")

    parser.add_argument("--bam_file", dest="bam_file", required=True,
        help="bam input file (reads mapped to contigs)")

    parser.add_argument("--precomputed_orfs_faa", dest="precomputed_orfs_faa", default=None,
        help="path to a precomputed ORF protein FASTA, plain or gzipped; skips "
             "FragGeneScanRs. Must be combined with --precomputed_orfs_bed (default: None)")

    parser.add_argument("--precomputed_orfs_bed", dest="precomputed_orfs_bed", default=None,
        help="path to a precomputed ORF BED file, plain or gzipped "
             "(contig_id, start(0-based), end, "
             "strand, orf_id) matching --precomputed_orfs_faa; skips the "
             ".out-to-BED conversion step. Must be combined with "
             "--precomputed_orfs_faa (default: None)")

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

    precomputed_orfs_faa = args.precomputed_orfs_faa
    precomputed_orfs_bed = args.precomputed_orfs_bed
    assembly_file = args.assembly_file
    bam_file = args.bam_file
    output_dir = args.output_dir
    overwrite = args.overwrite
    sample_name = args.sample_name
    nslots = args.nslots
    train_file_name = args.train_file_name

    have_precomputed_orfs = bool(precomputed_orfs_faa and precomputed_orfs_bed)

    if bool(precomputed_orfs_faa) != bool(precomputed_orfs_bed):
        print("--precomputed_orfs_faa and --precomputed_orfs_bed must be given together",
              file=sys.stderr)
        sys.exit(1)
    if not have_precomputed_orfs and not assembly_file:
        print("--assembly_file is required unless both --precomputed_orfs_faa and "
              "--precomputed_orfs_bed are given", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    if have_precomputed_orfs:
        check_file(precomputed_orfs_faa, "precomputed ORF protein FASTA")
        check_file(precomputed_orfs_bed, "precomputed ORF BED file")
    else:
        check_file(assembly_file, "input assembly file")
    check_file(bam_file, "input bam file")

    ###########################################################################
    # 3.2. Check output directory
    ###########################################################################

    if os.path.isdir(output_dir):
        if not overwrite:
            print(f"{output_dir} already exists; use --overwrite to overwrite")
            sys.exit(0)
        try:
            shutil.rmtree(output_dir)
        except Exception:
            print(f"rm -r output directory {output_dir} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.3. Create output directory
    ###########################################################################

    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir output directory {output_dir} failed", file=sys.stderr)
        sys.exit(1)

    faa_file = os.path.join(output_dir, f"{sample_name}_orfs.faa")
    bed_file = os.path.join(output_dir, f"{sample_name}_orfs.bed")

    if have_precomputed_orfs:

        ###########################################################################
        # 3.4. Stage precomputed ORFs (skip FragGeneScanRs + BED conversion)
        ###########################################################################

        # Decompressed on the way in: these may arrive gzipped (a prior run's
        # module 2 output is *_orfs.faa.gz / *_orfs.bed.gz), and every step below
        # -- the BED/BAM check, bedtools, and gzip_file() at the end -- expects the
        # same plain files the freshly-predicted path produces.
        stage_decompressed(precomputed_orfs_faa, faa_file, "--precomputed_orfs_faa")
        stage_decompressed(precomputed_orfs_bed, bed_file, "--precomputed_orfs_bed")

    else:

        ###########################################################################
        # 3.4. Predict ORFs
        ###########################################################################

        out_file = os.path.join(output_dir, f"{sample_name}_orfs.out")

        try:
            run(
                [
                    fraggenescan,
                    "-s", assembly_file,
                    "-a", faa_file,
                    "-m", out_file,
                    "-w", "0",
                    "--unordered",
                    "-p", str(nslots),
                    "-t", train_file_name
                ]
            )
        except subprocess.CalledProcessError:
            print("FragGeneScan failed", file=sys.stderr)
            sys.exit(1)

        ###########################################################################
        # 3.5. Create BED file
        ###########################################################################

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

        # The .out has now been fully parsed into the BED and has no other consumer,
        # so it gets deleted here instead of at the end of main() 
        os.remove(out_file)

    ###########################################################################
    # 3.6. Build bedtools genome file from BAM header
    ###########################################################################

    # The genome file's chromosome order is derived directly from the same BAM
    # that bedtools coverage -sorted will scan below, so it is guaranteed by
    # construction to match that BAM's coordinate-sort order.
    genome_file = os.path.join(output_dir, f"{sample_name}_genome.txt")
    try:
        header = subprocess.run(
            [samtools, "view", "-H", bam_file],
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
    # 3.7. Verify precomputed ORFs match this BAM
    ###########################################################################

    # bedtools sort -g and coverage -sorted both require every BED contig to be
    # present in the genome file. When the ORFs are predicted here that holds by
    # construction -- the BED is derived from the same assembly the BAM was built
    # against. Precomputed ORFs come from a separate file with nothing tying it
    # to this BAM, so check it up front: bedtools otherwise fails obscurely, or
    # (for a partial overlap) silently reports no coverage for the missing
    # contigs, which surfaces only as thin OPU tables in module 6.
    if have_precomputed_orfs:
        try:
            with open(bed_file, "r", encoding="utf-8") as fh_in:
                bed_chroms = {line.split("\t")[0] for line in fh_in if line.strip()}
            with open(genome_file, "r", encoding="utf-8") as fh_in:
                bam_chroms = {line.split("\t")[0] for line in fh_in if line.strip()}
        except Exception as exc:
            print(f"Comparing precomputed ORF contig IDs with the BAM header failed: {exc}",
                  file=sys.stderr)
            sys.exit(1)

        missing = bed_chroms - bam_chroms
        if missing:
            print(f"{len(missing)} of {len(bed_chroms)} contig IDs in "
                  f"--precomputed_orfs_bed are absent from the BAM header, e.g. "
                  f"{sorted(missing)[:3]}; the precomputed ORFs were not generated "
                  f"from this assembly/BAM. If they came from an earlier run, check "
                  f"that --id_sep matches the separator used then.", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.8. Sort BED file according to genome file order
    ###########################################################################

    # bedtools coverage -sorted requires -a and -b to follow the same chromosome
    # order as -g; this sorts the BED side, while the BAM side is already
    # guaranteed by construction since module 1's samtools sort produced it.
    sorted_bed_file = os.path.join(output_dir, f"{sample_name}_orfs_sorted.bed")
    try:
        run([bedtools, "sort", "-i", bed_file, "-g", genome_file], stdout_path=sorted_bed_file)
    except subprocess.CalledProcessError:
        print("bedtools sort of bed file failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.9. Get number of reads per ORF
    ###########################################################################

    # Run bedtools coverage -counts to get read counts per ORF
    # -sorted -g lets bedtools stream both inputs instead of loading the BAM into memory
    bedtool_reads = os.path.join(output_dir, f"{sample_name}_orfs_readscov.tsv")
    try:
        run(
            [bedtools, "coverage",
             "-a", sorted_bed_file,
             "-b", bam_file,
             "-counts",
             "-sorted",
             "-g", genome_file],
            stdout_path=bedtool_reads
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute read counts failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.10. Get mean coverage per ORF
    ###########################################################################

    # Run bedtools coverage -mean to get mean depth per ORF
    bedtools_mean = os.path.join(output_dir, f"{sample_name}_orfs_meancov.tsv")
    try:
        run(
            [bedtools, "coverage", "-a", sorted_bed_file, "-b", bam_file, "-mean",
             "-sorted", "-g", genome_file],
            stdout_path=bedtools_mean
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute mean coverage failed", file=sys.stderr)
        sys.exit(1)

    # Both bedtools coverage stages are done, so the sorted BED and the genome file
    # have no further use. They are removed here.
    for dead_path in (sorted_bed_file, genome_file):
        if os.path.isfile(dead_path):
            os.remove(dead_path)

    ###########################################################################
    # 3.11. Add sample names to mean coverage and read counts tables
    ###########################################################################

    add_sample(bedtool_reads, sample_name)
    add_sample(bedtools_mean, sample_name)

    ###########################################################################
    # 3.12. Compress the published outputs
    ###########################################################################

    # Done only now, after every bedtools stage has finished with the plain files.
    # Downstream consumers all read gzip directly: module 3 (byte-copy concat, bbduk,
    # mmseqs), module 5 (pyhmmer) and module 6 (DuckDB). The BED is accepted but
    # unused by module 4.
    for out_file in (faa_file, bed_file, bedtool_reads, bedtools_mean):
        gzip_file(out_file)

    ###########################################################################
    # 3.13. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###########################################################################
# 4. Execute main function
###########################################################################

if __name__ == "__main__":
    main()
