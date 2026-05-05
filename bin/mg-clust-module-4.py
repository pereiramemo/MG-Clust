#!/usr/bin/env python3
"""
mg-clust module 4: Taxonomic annotation of contigs.

Assumes execution inside the conda environment "mg-clust-module-4" (or equivalent)
where dependencies are available on PATH.

- Creates an MMseqs2 nucleotide sequence database from the input contigs fasta
- Runs mmseqs taxonomy against a GTDB database to assign taxonomy to each contig
- Exports the per-contig taxonomy assignments as a TSV file
- Generates a Kraken-style taxonomy report
- Maps ORF IDs to contig-level taxonomy via a left join on the BED file
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import shutil
import sys, os
import subprocess
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file

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

    parser.add_argument("--bed_file", dest="bed_file", required=True,
        help="ORF BED file from module 2 (*_orfs.bed); used to map ORF IDs to contig taxonomy")

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

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(args.contigs, "contigs fasta file")
    check_file(args.bed_file, "ORF BED file")

    ###########################################################################
    # 3.2. Check GTDB database; download if absent
    ###########################################################################

    # A valid mmseqs2 database always has a companion .dbtype file.
    # If it is missing, download the GTDB database via mmseqs databases.
    if not os.path.isfile(args.gtdb + ".dbtype"):
        print(f"GTDB database not found at {args.gtdb}; downloading now ...")

        gtdb_dir = os.path.dirname(args.gtdb)
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
                    args.gtdb,
                    download_tmp,
                    "--threads", str(args.nslots)
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
    # 3.4. Create output directory
    ###########################################################################

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Create MMseqs2 nucleotide database from contigs
    ###########################################################################

    contigs_db = os.path.join(args.output_dir, "contigs_db")

    try:
        run(
            [
                mmseqs,
                "createdb",
                args.contigs,
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

    tax_db = os.path.join(args.output_dir, "contigs_tax_db")
    tmp_dir = os.path.join(args.output_dir, "tmp")

    try:
        run(
            [
                mmseqs,
                "taxonomy",
                contigs_db,
                args.gtdb,
                tax_db,
                tmp_dir,
                "--threads", str(args.nslots),
                "-s", str(args.sensitivity),
                "--lca-mode", str(args.lca_mode),
                "--tax-lineage", str(args.tax_lineage)
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

    tax_tsv = os.path.join(args.output_dir, f"{args.sample_name}_contig_tax_annot.tsv")

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

    tax_report = os.path.join(args.output_dir, f"{args.sample_name}_contig_tax_report.txt")

    try:
        run(
            [
                mmseqs,
                "taxonomyreport",
                args.gtdb,
                tax_db,
                tax_report
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs taxonomyreport failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.9. Map ORF IDs to contig-level taxonomy
    ###########################################################################

    # mmseqs createtsv output columns (0-indexed): contig_id, taxid, rank, name[, lineage]
    tax_map = {}
    try:
        with open(tax_tsv, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts:
                    tax_map[parts[0]] = parts[1:]  # contig_id -> [taxid, rank, name, ...]
    except Exception as exc:
        print(f"reading taxonomy TSV failed: {exc}", file=sys.stderr)
        sys.exit(1)

    orf_tax_tsv = os.path.join(args.output_dir, f"{args.sample_name}_orf_tax_annot.tsv")

    try:
        with open(args.bed_file, "r", encoding="utf-8") as bed_fh, \
             open(orf_tax_tsv, "w", encoding="utf-8") as out_fh:
            for line in bed_fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                contig_id = parts[0]
                orf_id    = parts[4]
                # Left join: fill with empty strings when contig has no taxonomy hit
                tax_fields = (tax_map.get(contig_id, "NA"))
                out_fh.write(f"{orf_id}\t" + "\t".join(tax_fields) + "\n")
    except Exception as exc:
        print(f"mapping ORFs to taxonomy failed: {exc}", file=sys.stderr)
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
