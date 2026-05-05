#! /usr/bin/env python3

"""
mg-clust module 1: de novo assembly and read mapping.

Assumes execution inside the conda environment "mg-clust-module-1" (or equivalent)
where dependencies are available on PATH.

- Assembles paired-end reads with MEGAHIT
- Maps reads back to the assembly with BWA-MEM
- Converts, filters, and sorts alignments with samtools
- Optionally marks and removes duplicates with Picard MarkDuplicates
- Removes intermediate files after completion
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import sys, os
import subprocess
import shutil
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file

megahit = "megahit"
bwa = "bwa"
samtools = "samtools"
picard = "picard"  # use picard wrapper instead of java -jar

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: de novo assembly and read mapping", 
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")
    
    parser.add_argument("--assem_preset", dest="assem_preset", default="meta-sensitive",
        help="MEGAHIT preset to generate assembly (default: meta-sensitive)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--min_contig_length", dest="min_contig_length", type=int, default=250,
        help="minimum length of contigs (smaller than this will be discarded; default: 250)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    parser.add_argument("--reads1", dest="reads1", required=True,
        help="input R1 metagenome data (as fasta/q file)")

    parser.add_argument("--reads2", dest="reads2", required=True,
        help="input R2 metagenome data (as fasta/q file)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--min_seq", dest="min_seq", type=int, default=5,
        help="minimum number of assembled sequences to continue (default: 5)")

    parser.add_argument("--markdup", dest="markdup", action="store_true", default=False,
        help="run Picard MarkDuplicates to remove duplicates (default: False)")
    
    return parser.parse_args()

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    check_tools([megahit, bwa, samtools, picard])
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################
    
    check_file(args.reads1, "read1")
    check_file(args.reads2, "read2")

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
        print(f"mkdir {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.4. De novo assembly (always run MEGAHIT; it will create output_dir)
    ###########################################################################

    try:
        megahit_ouput_dir = os.path.join(args.output_dir, "assembly")
        run(
            [
                megahit,
                "--num-cpu-threads",
                str(args.nslots),
                "-1",
                args.reads1,
                "-2",
                args.reads2,
                "--presets",
                args.assem_preset,
                "--min-contig-len",
                str(args.min_contig_length),
                "--out-prefix",
                args.sample_name,
                "--out-dir",
                megahit_ouput_dir,
            ]
        )
    except subprocess.CalledProcessError:
        print("megahit failed", file=sys.stderr)
        sys.exit(1)

    assembly_file = os.path.join(megahit_ouput_dir, f"{args.sample_name}.contigs.fa")
    # check if assembly_file exists
    check_file(assembly_file, "assembly_file")
    
    ###########################################################################
    # 3.5. Count number of sequences by counting lines starting with '>'
    ###########################################################################
    
    try:
        result = subprocess.run(
            ["grep", "-c", "^>", assembly_file],
            capture_output=True, text=True, check=False
        )
        assembly_file_nseq = int(result.stdout.strip())
    except Exception as e:
        print(f"Failed to count sequences in assembly file: {e}", file=sys.stderr)
        sys.exit(1)

    if assembly_file_nseq < args.min_seq:
        print(f"Not enough assembled sequences to continue: {assembly_file_nseq}")
        sys.exit(0)

    ###########################################################################
    # 3.6. Add the sample name to the contig headers to ensure uniqueness across samples (e.g. for downstream clustering)
    ###########################################################################

    tmp_assembly_file = assembly_file + ".tmp"
    try:
        with open(assembly_file, "r", encoding="utf-8") as in_fh, open(tmp_assembly_file, "w", encoding="utf-8") as out_fh:
            for line in in_fh:
                if line.startswith(">"):
                    header = line[1:].rstrip("\n")
                    if not header.startswith(f"{args.sample_name}|"):
                        header = f"{args.sample_name}|{header}"
                    out_fh.write(f">{header}\n")
                else:
                    out_fh.write(line)
        os.replace(tmp_assembly_file, assembly_file)
    except Exception as e:
        if os.path.exists(tmp_assembly_file):
            os.remove(tmp_assembly_file)
        print(f"Failed to rewrite contig headers with sample name: {e}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Map short reads
    ###########################################################################

    # bwa index
    try:
        run([bwa, "index", assembly_file])
    except subprocess.CalledProcessError:
        print("bwa index failed", file=sys.stderr)
        sys.exit(1)

    # bwa mem -> SAM
    sam_path = os.path.join(args.output_dir, f"{args.sample_name}.sam")
    try:
        rg = f"@RG\\tID:{args.sample_name}\\tSM:{args.sample_name}\\tLB:{args.sample_name}"
        run([bwa, "mem", "-M", "-t", str(args.nslots), "-R", rg,
             assembly_file, args.reads1, args.reads2], stdout_path=sam_path)
    except subprocess.CalledProcessError:
        print("bwa mem failed", file=sys.stderr)
        sys.exit(1)

    # samtools view -> BAM (filter)
    # convert to bam and filter high-quality primary alignments
    # flags taken from https://broadinstitute.github.io/picard/explain-flags.html
    # secondary alignments are moved; however this makes very little difference
    # ORFs coverage using F4 vs F260 flags had a MSE=0.00486 and Pearson cor 0.997 (same toydataset)
    # Here we use -F 2308 (4 + 256 + 2048 = 2308) to filter out unmapped, secondary alignments, and supplementary alignments

    bam_path = os.path.join(args.output_dir, f"{args.sample_name}.bam")
    try:
        run([samtools, "view", "-@", str(args.nslots), 
             "-q", "10", "-F", "2308", "-b", "-o", bam_path, sam_path])
    except subprocess.CalledProcessError:
        print("samtools convert to bam failed", file=sys.stderr)
        sys.exit(1)

    # samtools sort -> sorted BAM
    sorted_bam_path = os.path.join(args.output_dir, f"{args.sample_name}_sorted.bam")
    try:
        run([samtools, "sort", "-@", str(args.nslots), "-o", sorted_bam_path, bam_path])
    except subprocess.CalledProcessError:
        print("samtools sort failed", file=sys.stderr)
        sys.exit(1)

    # samtools index
    try:
        run([samtools, "index", "-@", str(args.nslots), sorted_bam_path])
    except subprocess.CalledProcessError:
        print("samtools inex failed", file=sys.stderr)
        sys.exit(1)

    # remove duplicates with Picard
    if args.markdup:
        tmp_dir = os.path.join(args.output_dir, "tmp")
        try:
            os.makedirs(tmp_dir, exist_ok=False)
        except Exception:
            print(f"mkdir {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

        markdup_bam = os.path.join(args.output_dir, f"{args.sample_name}_sorted_markdup.bam")
        metrics_file = os.path.join(args.output_dir, f"{args.sample_name}_sorted_markdup.metrics.txt")

        try:
            run(
                [
                picard,
                "MarkDuplicates",
                "-INPUT", sorted_bam_path,
                "-OUTPUT", markdup_bam,
                "-METRICS_FILE", metrics_file,
                "-REMOVE_DUPLICATES", "TRUE",
                "-ASSUME_SORTED", "TRUE",
                "-MAX_FILE_HANDLES_FOR_READ_ENDS_MAP", "900",
                "-TMP_DIR", tmp_dir,
                ]
            )
        except subprocess.CalledProcessError:
            print("picard failed", file=sys.stderr)
            sys.exit(1)

        #######################################################################
        # 3.8. Define cleanup paths including the temporary directory if markdup was run
        #######################################################################
          
        cleanup_paths = [sam_path, bam_path, tmp_dir]
    else:
        cleanup_paths = [sam_path, bam_path]

    ########################################################################### 
    # 3.9. Clean
    ###########################################################################

    try:
        for p in cleanup_paths:
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.remove(p)
    except Exception:
        print("removing intermediate mapping files failed", file=sys.stderr)
        sys.exit(1)

    interm_contigs = os.path.join(args.output_dir, "assembly", "intermediate_contigs")
    if os.path.isdir(interm_contigs):
        try:
            shutil.rmtree(interm_contigs)
        except Exception:
            print("removing assembly intermediate files failed", file=sys.stderr)
            sys.exit(1)

    ########################################################################### 
    # 3.10. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

########################################################################### 
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()


