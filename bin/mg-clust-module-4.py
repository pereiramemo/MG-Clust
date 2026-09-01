#!/usr/bin/env python3
"""
mg-clust module 4: Taxonomic annotation of contigs.

Assumes execution inside the conda environment "mg-clust-module-4" (or equivalent)
where dependencies are available on PATH.

- Creates an MMseqs2 nucleotide sequence database from the input contigs fasta
- Runs mmseqs taxonomy against a GTDB database to assign taxonomy to each contig
- Exports the per-contig taxonomy assignments as a TSV file
- Generates a Kraken-style taxonomy report

Note: --bed_file is optional and currently unused; taxonomy is reported at the
contig level only, no per-ORF join is produced.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import shutil
import sys, os
import subprocess
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file, gzip_file

mmseqs = "mmseqs"

GTDB_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "gtdb", "gtdb")

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Taxonomic annotation of contigs using MMseqs2 taxonomy", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--contigs", dest="contigs", required=True,
        help="input contigs fasta file (nucleotide sequences)")

    parser.add_argument("--bed_file", dest="bed_file", default=None,
        help="ORF BED file from module 2 (*_orfs.bed); accepted for compatibility but "
             "currently unused -- taxonomy is reported at the contig level only, no "
             "per-ORF join is produced (default: None)")

    parser.add_argument("--gtdb", dest="gtdb", default=GTDB_DEFAULT,
        help=f"path to the MMseqs2 GTDB taxonomy database (default: {GTDB_DEFAULT})")

    parser.add_argument("--lca_mode", dest="lca_mode", type=int, default=3,
        help="LCA mode passed to mmseqs taxonomy: 1=single-hit, 2=LCA of best hits, "
             "3=top-hit with LCA fallback (default: 3)")

    parser.add_argument("--sensitivity", dest="sensitivity", type=float, default=4.0,
        help="sensitivity passed as -s to mmseqs taxonomy (default: 4.0)")

    parser.add_argument("--tax_lineage", dest="tax_lineage", type=int, default=1,
        help="include full lineage in TSV output (1=yes, 0=no; default: 1)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    check_tools([mmseqs])
    args = parse_args()
    contigs = args.contigs
    bed_file = args.bed_file
    output_dir = args.output_dir
    overwrite = args.overwrite
    gtdb = args.gtdb
    sample_name = args.sample_name
    nslots = args.nslots
    sensitivity = args.sensitivity
    lca_mode = args.lca_mode
    tax_lineage = args.tax_lineage

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(contigs, "contigs fasta file")
    if bed_file:
        check_file(bed_file, "ORF BED file")

    ###########################################################################
    # 3.2. Check GTDB database; download if absent
    ###########################################################################

    # A valid mmseqs2 database always has a companion .dbtype file.
    # If it is missing, download the GTDB database via mmseqs databases.
    if not os.path.isfile(gtdb + ".dbtype"):
        print(f"GTDB database not found at {gtdb}; downloading now ...")

        gtdb_dir = os.path.dirname(gtdb)
        try:
            os.makedirs(gtdb_dir, exist_ok=True)
        except Exception:
            print(f"mkdir {gtdb_dir} failed", file=sys.stderr)
            sys.exit(1)

        download_tmp = os.path.join(gtdb_dir, "gtdb_download_tmp")
        try:
            run(
                [
                    mmseqs,
                    "databases",
                    "GTDB",
                    gtdb,
                    download_tmp,
                    "--threads", str(nslots)
                ]
            )
        except subprocess.CalledProcessError:
            print("mmseqs databases GTDB download failed", file=sys.stderr)
            sys.exit(1)

        if os.path.isdir(download_tmp):
            try:
                shutil.rmtree(download_tmp)
            except Exception:
                print(f"rm -r {download_tmp} failed", file=sys.stderr)
                sys.exit(1)

        print("GTDB database download complete.")

    ###########################################################################
    # 3.3. Check output directory
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
    # 3.4. Create output directory
    ###########################################################################

    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir {output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Create MMseqs2 nucleotide database from contigs
    ###########################################################################

    contigs_db = os.path.join(output_dir, "contigs_db")

    try:
        run(
            [
                mmseqs,
                "createdb",
                contigs,
                contigs_db,
                "--dbtype", "2"  # 2 = nucleotide
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createdb failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Run mmseqs taxonomy against GTDB
    ###########################################################################

    tax_db = os.path.join(output_dir, "contigs_tax_db")
    tmp_dir = os.path.join(output_dir, "tmp")

    try:
        run(
            [
                mmseqs,
                "taxonomy",
                contigs_db,
                gtdb,
                tax_db,
                tmp_dir,
                "--threads", str(nslots),
                "-s", str(sensitivity),
                "--lca-mode", str(lca_mode),
                "--tax-lineage", str(tax_lineage)
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs taxonomy failed", file=sys.stderr)
        sys.exit(1)

    # Remove tmp directory
    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            print(f"rm -r {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.7. Export taxonomy assignments to TSV
    ###########################################################################

    tax_tsv = os.path.join(output_dir, f"{sample_name}_contig_tax_annot.tsv")

    try:
        run(
            [
                mmseqs,
                "createtsv",
                contigs_db,
                tax_db,
                tax_tsv
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createtsv failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8. Generate Kraken-style taxonomy report
    ###########################################################################

    tax_report = os.path.join(output_dir, f"{sample_name}_contig_tax_report.txt")

    try:
        run(
            [
                mmseqs,
                "taxonomyreport",
                gtdb,
                tax_db,
                tax_report
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs taxonomyreport failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.9. Compress the taxonomy table
    ###########################################################################

    # Module 6 concatenates this with a raw byte copy and reads it via DuckDB; both
    # handle gzip.
    gzip_file(tax_tsv)

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
