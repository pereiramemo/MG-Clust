#!/usr/bin/env python3
"""
mg-clust module 6: Build unified OPU-ORF coverage table.

Assumes execution inside the conda environment "mg-clust-module-6" (or equivalent).

- Concatenates per-sample mean coverage tables
- Concatenates per-sample read coverage tables
- LEFT JOINs both coverages onto the clustering TSV via DuckDB (clust_tsv is the
  driving table: an ORF with coverage but absent from clust_tsv, e.g. filtered out
  upstream by module 3, is silently excluded rather than kept with zero coverage) to
  produce a per-ORF table: sample_name, opu_id, orf_id, read_coverage, mean_coverage
- Collapses that table by summing coverage per (sample_name, opu_id), dropping orf_id,
  to produce the primary per-OPU abundance table
- If given, also concatenates per-sample taxonomy and/or functional annotation files
  into their own standalone tables (contigs_tax_annot.tsv,
  orfs-minlen<N>aa-fun_annot.tsv) -- these are NOT joined against the coverage/OPU
  tables above; they are separate outputs sharing sample-level provenance only
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import shutil
import sys
import duckdb

sys.path.insert(0, os.path.dirname(__file__))
from utils import check_file

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Build unified OPU-ORF coverage table",
        add_help=False,
    )

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--meancov_files", dest="meancov_files", required=True, nargs="+",
        help="list of per-sample mean coverage tables (*_orfs_meancov.tsv)",
    )

    parser.add_argument("--readscov_files", dest="readscov_files", required=True, nargs="+",
        help="list of per-sample read coverage tables (*_orfs_readscov.tsv)",
    )

    parser.add_argument("--tax_annot_files", dest="tax_annot_files", required=False, nargs="*",
        help="list of per-sample contig taxonomy annotation files from module 4 "
             "(*_contig_tax_annot.tsv); concatenated into a standalone output, "
             "not joined against the coverage/OPU tables",
    )

    parser.add_argument("--fun_annot_files", dest="fun_annot_files", required=False, nargs="*",
        help="list of per-sample functional annotation files from module 5 "
             "(*_orfs-minlen<N>aa-fun_annot.tsv); concatenated into a standalone "
             "output, not joined against the coverage/OPU tables",
    )

    parser.add_argument("--clust_tsv", dest="clust_tsv", required=True,
        help="MMseqs2 clustering TSV from module 3 (columns: opu_id, orf_id)",
    )

    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=0.7,
        help="clustering threshold used to name the output file (default: 0.7)",
    )

    parser.add_argument("--min_orf_len", dest="min_orf_len", type=int, default=100,
        help="minimum ORF length used to filter ORFs (default: 100)",
    )


    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data",
    )

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)",
    )

    return parser.parse_args()

###############################################################################
# 2.2 Concatenate table files
###############################################################################

def concat_tables(input_files: list[str], output_file: str) -> None:
    """Concatenate all files in input_files into output_file."""
    try:
        with open(output_file, "wb") as out_fh:
            for in_file in input_files:
                check_file(in_file, "coverage file")
                with open(in_file, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)
    except Exception as exc:
        print(f"concatenating files into {output_file} failed: {exc}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.3 Build unified OPU-ORF coverage table via DuckDB
###############################################################################

def build_unified_table(
    clust_tsv: str,
    meancov_concat: str,
    readscov_concat: str,
    output1_tsv: str,
) -> None:
    
    cov_cols = "{'contig_id':'VARCHAR','start':'BIGINT','end':'BIGINT','strand':'VARCHAR', \
                 'orf_id':'VARCHAR','coverage':'DOUBLE', 'sample_name':'VARCHAR'}"
    
    clust_cols = "{'opu_id':'VARCHAR','orf_id':'VARCHAR'}"

    sql = f"""
        COPY (
            SELECT
                r.sample_name,
                c.opu_id,
                c.orf_id,
                COALESCE(r.coverage, 0.0) AS read_coverage,
                COALESCE(m.coverage, 0.0) AS mean_coverage,
            FROM read_csv('{clust_tsv}', delim='\t', header=false,
                          columns={clust_cols}) AS c
            LEFT JOIN read_csv('{readscov_concat}', delim='\t', header=false,
                               columns={cov_cols}) AS r ON c.orf_id = r.orf_id
            LEFT JOIN read_csv('{meancov_concat}', delim='\t', header=false,
                               columns={cov_cols}) AS m ON c.orf_id = m.orf_id
        ) TO '{output1_tsv}' (HEADER true, DELIMITER '\t')
    """
    try:
        duckdb.execute(sql)
    except Exception as exc:
        print(f"building unified table {output1_tsv} failed: {exc}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.4 Collapse coverage tables OPU ID
###############################################################################

def build_collapsed_table(
    output1_tsv: str,
    output2_tsv: str,
) -> None:
    output1_cols = "{'sample_name':'VARCHAR','opu_id':'VARCHAR','orf_id':'VARCHAR', \
                    'read_coverage':'DOUBLE','mean_coverage':'DOUBLE'}"

    sql = f"""
        COPY (
            SELECT
                o1.sample_name,
                o1.opu_id,
                SUM(o1.read_coverage) AS read_coverage,
                SUM(o1.mean_coverage) AS mean_coverage
            FROM read_csv('{output1_tsv}', delim='\t', header=true,
                          columns={output1_cols}) AS o1
            GROUP BY o1.sample_name, o1.opu_id
        ) TO '{output2_tsv}' (HEADER true, DELIMITER '\t')
    """
    try:
        duckdb.execute(sql)
    except Exception as exc:
        print(f"building collapsed table {output2_tsv} failed: {exc}", file=sys.stderr)
        sys.exit(1)


###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(args.clust_tsv, "cluster TSV")
    for f in args.meancov_files:
        check_file(f, "mean coverage file")
    for f in args.readscov_files:
        check_file(f, "reads coverage file")

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
    # 3.4. Concatenate coverage files
    ###########################################################################

    meancov_concat = os.path.join(args.output_dir, "orfs_meancov.tsv")
    readscov_concat = os.path.join(args.output_dir, "orfs_readscov.tsv")

    concat_tables(args.meancov_files, meancov_concat)
    concat_tables(args.readscov_files, readscov_concat)

    ###########################################################################
    # 3.5. Concatenate annot files
    ###########################################################################
    
    if args.tax_annot_files:
        tax_annot_concat = os.path.join(args.output_dir, "contigs_tax_annot.tsv")
        concat_tables(args.tax_annot_files, tax_annot_concat)

    if args.fun_annot_files:
        fun_annot_concat = os.path.join(args.output_dir, f"orfs-minlen{args.min_orf_len}aa-fun_annot.tsv")
        concat_tables(args.fun_annot_files, fun_annot_concat)

    ###########################################################################
    # 3.5. Build unified OPU-ORF coverage table
    ###########################################################################

    clust_thres_str = str(args.clust_thres * 100).rstrip("0").rstrip(".")
    output1_tsv = os.path.join(
        args.output_dir,
        f"orfs_clust-minlen{args.min_orf_len}aa-id{clust_thres_str}perc-coverage.tsv",
    )

    build_unified_table(args.clust_tsv, meancov_concat, 
                        readscov_concat, output1_tsv)

    ###########################################################################
    # 3.6. Build collapsed by sample and OPU ID coverage table
    ###########################################################################

    output2_tsv = os.path.join(
        args.output_dir,
        f"orfs_clust-minlen{args.min_orf_len}aa-id{clust_thres_str}perc-collapsed_coverage.tsv",
    )

    build_collapsed_table(output1_tsv, output2_tsv)

    ########################################################################### 
    # 3.7. Write output log and exit
    ###########################################################################
    
    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###############################################################################
# 4. Run the main function
###############################################################################

if __name__ == "__main__":
    main()
