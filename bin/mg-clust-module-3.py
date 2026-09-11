#!/usr/bin/env python3
"""
mg-clust module 3: Concatenation, filtering, MMseqs2 DB creation, and ORF clustering.

Assumes execution inside the conda environment "mg-clust-module-3" (or equivalent)
where dependencies are available on PATH.

- Concatenates ORF protein sequences from all samples (raw byte copy; headers already
  carry the `<sample_name>|` prefix applied by module 1 to the contig names)
- Filters ORFs by minimum length using bbduk
- Creates an MMseqs2 sequence database from the filtered ORFs
- Clusters filtered ORFs using MMseqs2 at a given sequence identity threshold
- Exports clustering results as a TSV table
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import sys, os
import subprocess
import shutil
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, gzip_file

bbduk = "bbduk.sh"
mmseqs = "mmseqs"

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Concatenation, filtering, MMseqs2 DB creation, and ORF clustering", 
        add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--orf_files", dest="orf_files", required=True, nargs="+",
        help="list of per-sample ORF fasta files (*_orfs.faa)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--max_mem", dest="max_mem", type=int, default=0,
        help="memory budget for this task in GiB, passed on to bbduk as -Xmx")

    parser.add_argument("--min_orf_len", dest="min_orf_len", type=int, default=60,
        help="minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default: 60)")

    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=0.7,
        help="clustering threshold, passed as --min-seq-id to mmseqs cluster (default: 0.7)")

    parser.add_argument("--clust_cov_len", dest="clust_cov_len", type=float, default=0.85,
        help="minimum fraction of aligned residues for clustering, passed as -c to mmseqs cluster (default: 0.85)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()


# Dev only section
"""
orf_files = [
    "/home/epereira/workspace/repos/tools/MG-Clust/test/mg-clust-output/intermediate/module-2/P01-A01-1/P01-A01-1_orfs.faa.gz",
    "/home/epereira/workspace/repos/tools/MG-Clust/test/mg-clust-output/intermediate/module-2/P01-A02-9/P01-A02-9_orfs.faa.gz",
    "/home/epereira/workspace/repos/tools/MG-Clust/test/mg-clust-output/intermediate/module-2/P01-A03-17/P01-A03-17_orfs.faa.gz"]
nslots = 4
max_mem = 16
min_orf_len = 60
clust_thres = 0.7
clust_cov_len = 0.85
output_dir = "/home/epereira/workspace/repos/tools/MG-Clust/test/mg-clust-output/intermediate/module-3"
overwrite = False
"""

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    check_tools([bbduk, mmseqs])
    args = parse_args()
    orf_files = args.orf_files
    nslots = args.nslots
    max_mem = args.max_mem
    min_orf_len = args.min_orf_len
    clust_thres = args.clust_thres
    clust_cov_len = args.clust_cov_len
    output_dir = args.output_dir
    overwrite = args.overwrite

    ###########################################################################
    # 3.1. Check output directory
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
    # 3.2. Create output directory
    ###########################################################################

    try:
        os.makedirs(output_dir, exist_ok=False)
    except Exception:
        print(f"mkdir {output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.3. Concat ORF files
    ###########################################################################

    # Inputs arrive gzipped from module 2. The raw byte-copy concatenation below
    # needs no change: concatenated gzip members are themselves a valid gzip stream,
    # and bbduk/mmseqs both read .gz directly.
    concat_orfs = os.path.join(output_dir, "orfs.faa.gz")

    with open(concat_orfs, "wb") as out_fh:
        for orf_file in orf_files:
            try:
                with open(orf_file, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)
            except Exception as exc:
                print(f"concatenating {orf_file} failed: {exc}", file=sys.stderr)
                sys.exit(1)

    ###########################################################################
    # 3.4. Filter ORFs by length
    ###########################################################################

    concat_orfs_filt = os.path.join(output_dir, f"orfs_filt-minlen{min_orf_len}aa.faa.gz")

    try:
        run(
            [
                bbduk,
                # Must precede the bbduk arguments: the wrapper forwards a leading
                # -Xmx straight to the JVM instead of computing one from the node's
                # RAM. 0.75 leaves headroom for the JVM's own non-heap memory, which
                # is outside -Xmx but inside the cgroup.
                *([f"-Xmx{max(1, int(max_mem * 0.75))}g"] if max_mem else []),
                f"in={concat_orfs}",
                f"out={concat_orfs_filt}",
                "overwrite=t",
                f"threads={nslots}",
                f"minlength={min_orf_len}",
                "amino=t"
            ]
        )
    except subprocess.CalledProcessError:
        print("bbduk failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Create mmseqs database
    ###########################################################################

    orfs_filt_db_dir = os.path.join(output_dir, f"orfs_filt_db-minlen{min_orf_len}aa")
    orfs_filt_db = os.path.join(orfs_filt_db_dir, "orfs_filt_db")

    try:
        os.makedirs(orfs_filt_db_dir, exist_ok=False)
    except Exception:
        print(f"mkdir {orfs_filt_db_dir} failed", file=sys.stderr)
        sys.exit(1)

    try:
        run(
            [
                mmseqs,
                "createdb",
                concat_orfs_filt,
                orfs_filt_db,
                "--dbtype", "1"
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createdb failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Run ORF clustering
    ###########################################################################

    clust_thres_str = str(clust_thres * 100).rstrip("0").rstrip(".")
    clust_dir = os.path.join(output_dir, f"orfs_clust-minlen{min_orf_len}aa-id{clust_thres_str}perc")

    try:
        os.makedirs(clust_dir, exist_ok=False)
    except Exception:
        print(f"mkdir {clust_dir} failed", file=sys.stderr)
        sys.exit(1)

    tmp_dir = os.path.join(clust_dir, "tmp")
    clust_db = os.path.join(clust_dir, f"orfs_clust-minlen{min_orf_len}aa-id{clust_thres_str}perc")

    try:
        run(
            [
                mmseqs,
                "cluster",
                orfs_filt_db,
                clust_db,
                tmp_dir,
                "--min-seq-id", str(clust_thres),
                "--threads", str(nslots),
                "--cov-mode", "0",
                "-c", str(clust_cov_len)
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs cluster failed", file=sys.stderr)
        sys.exit(1)

    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            print(f"rm {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.7. Convert clustering results to TSV
    ###########################################################################

    mmseqs_clust_table = os.path.join(clust_dir, f"orfs_clust-minlen{min_orf_len}aa-id{clust_thres_str}perc.tsv")

    try:
        run(
            [
                mmseqs,
                "createtsv",
                orfs_filt_db,
                orfs_filt_db,
                clust_db,
                mmseqs_clust_table
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createtsv failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8. Compress the clustering table and drop the mmseqs databases
    ###########################################################################

    # mmseqs createtsv only writes plain text, so compress after the fact. Module 6
    # reads this through DuckDB read_csv, which handles .gz transparently.
    gzip_file(mmseqs_clust_table)

    # createtsv in 3.7 was the last consumer of both mmseqs databases, so they 
    # are deleted here.

    # The cluster DB shares clust_dir with the published table and differs from it
    # only by suffix (<clust_db> vs <clust_db>.tsv.gz), so the table is excluded by
    # name rather than trusted to fall outside the glob. glob.escape guards against
    # metacharacters in --output_dir.
    keep = {mmseqs_clust_table, mmseqs_clust_table + ".gz"}
    for db_path in glob.glob(glob.escape(clust_db) + "*"):
        if db_path in keep or not os.path.isfile(db_path):
            continue
        try:
            os.remove(db_path)
        except OSError as exc:
            print(f"rm {db_path} failed: {exc}", file=sys.stderr)
            sys.exit(1)

    # The sequence DB has a directory to itself, so it goes wholesale.
    if os.path.isdir(orfs_filt_db_dir):
        try:
            shutil.rmtree(orfs_filt_db_dir)
        except Exception:
            print(f"rm {orfs_filt_db_dir} failed", file=sys.stderr)
            sys.exit(1)

    ########################################################################### 
    # 3.9. Write output log and exit
    ###########################################################################
    
    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()
