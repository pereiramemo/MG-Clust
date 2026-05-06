#!/usr/bin/env python3
"""
mg-clust module 5: Functional annotation of ORFs using pyHMMER.

Assumes execution inside the conda environment "mg-clust-module-5" (or equivalent)
where pyhmmer is importable.

- Searches ORF protein sequences against a HMM database (default: KEGG KO profiles)
- Uses per-model gathering thresholds if available (--cut_ga), otherwise applies an e-value threshold
- Exports a domain table and a best-hit annotation table (one KO per ORF)
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import os
import shutil
import sys
import tarfile
import urllib.request
import pyhmmer

sys.path.insert(0, os.path.dirname(__file__))
from utils import check_file

KO_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "ko", "ko_profiles.hmm")
KO_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Functional annotation of ORFs using pyHMMER", 
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--orfs_faa", dest="orfs_faa", required=True,
        help="input ORF protein sequences (FASTA)")

    parser.add_argument("--hmm_db", dest="hmm_db", default=KO_DEFAULT,
        help=f"path to HMM database (default: KO profiles at {KO_DEFAULT})")

    parser.add_argument("--evalue_thres", dest="evalue_thres", type=float, default=1e-3,
        help="e-value threshold for hmmsearch (default: 1e-3; ignored when --cut_ga is set)")

    parser.add_argument("--cut_ga", dest="cut_ga",
        action=argparse.BooleanOptionalAction, default=False,
        help="use per-model gathering thresholds (default: False; use --cut_ga to enable)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Run pyhmmer hmmsearch and write domain table
###############################################################################

def run_hmmsearch(hmm_db, orfs_faa, domtblout, evalue_thres, cut_ga, nslots):
    with pyhmmer.plan7.HMMFile(hmm_db) as hmms:
        with pyhmmer.easel.SequenceFile(orfs_faa, digital=True) as seqs:
            with open(domtblout, "wb") as out:
                if cut_ga:
                    # domE acts as fallback for profiles without GA thresholds
                    for hits in pyhmmer.hmmer.hmmsearch(
                        hmms, seqs, domE=evalue_thres, bit_cutoffs="gathering", cpus=nslots
                    ):
                        hits.write(out, format="domains", header=False)
                else:
                    for hits in pyhmmer.hmmer.hmmsearch(
                        hmms, seqs, domE=evalue_thres, cpus=nslots
                    ):
                        hits.write(out, format="domains", header=False)

###############################################################################
# 2.3 Parse domain table and write best-hit annotation table
###############################################################################

def write_best_hits(domtblout, annotation_table):
    args = parse_args()
    # domain table columns (0-indexed after whitespace split):
    #  0: target (ORF) name, 3: query (KO profile) name,
    #  6: full-sequence E-value,  7: full-sequence score
    best = {}  # orf_id -> (ko_id, score, evalue)
    with open(domtblout, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"):
                continue
            cols = line.split()
            if len(cols) < 8:
                continue
            orf_id = cols[0].decode()
            ko_id  = cols[3].decode()
            evalue = float(cols[6])
            score  = float(cols[7])
            if orf_id not in best or score > best[orf_id][1]:
                best[orf_id] = (ko_id, score, evalue)

    with open(annotation_table, "w") as out:
        for orf_id, (ko_id, score, evalue) in sorted(best.items()):
            out.write(f"{orf_id}\t{ko_id}\t{score:.2f}\t{evalue:.2e}\n")

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(args.orfs_faa, "ORF protein FASTA file")

    ###########################################################################
    # 3.2. Check KO HMM database; download if absent
    ###########################################################################

    if not os.path.isfile(args.hmm_db):
        print(f"HMM database not found at {args.hmm_db}; downloading KO profiles ...")
        ko_dir = os.path.dirname(args.hmm_db)
        try:
            os.makedirs(ko_dir, exist_ok=True)
        except Exception:
            print(f"mkdir {ko_dir} failed", file=sys.stderr)
            sys.exit(1)

        archive = os.path.join(ko_dir, "profiles.tar.gz")
        try:
            print(f"Downloading KO profiles from {KO_PROFILES_URL} ...")
            urllib.request.urlretrieve(KO_PROFILES_URL, archive)
        except Exception as exc:
            print(f"Download of KO profiles failed: {exc}", file=sys.stderr)
            sys.exit(1)

        try:
            print("Extracting profiles ...")
            with tarfile.open(archive) as tar:
                tar.extractall(ko_dir)
        except Exception as exc:
            print(f"Extraction failed: {exc}", file=sys.stderr)
            sys.exit(1)

        profiles_dir = os.path.join(ko_dir, "profiles")
        hmm_files = sorted(glob.glob(os.path.join(profiles_dir, "*.hmm")))
        if not hmm_files:
            print("No .hmm files found after extraction", file=sys.stderr)
            sys.exit(1)

        try:
            print(f"Concatenating {len(hmm_files)} HMM profiles into {args.hmm_db} ...")
            with open(args.hmm_db, "wb") as out:
                for hmm_file in hmm_files:
                    with open(hmm_file, "rb") as f:
                        out.write(f.read())
        except Exception as exc:
            print(f"Concatenation of HMM profiles failed: {exc}", file=sys.stderr)
            sys.exit(1)

        os.remove(archive)
        shutil.rmtree(profiles_dir)
        print("KO HMM database ready.")

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
    # 3.5. Run hmmsearch
    ###########################################################################

    domtblout = os.path.join(args.output_dir, "orfs_ko_domtblout.txt")

    try:
        run_hmmsearch(
            hmm_db=args.hmm_db,
            orfs_faa=args.orfs_faa,
            domtblout=domtblout,
            evalue_thres=args.evalue_thres,
            cut_ga=args.cut_ga,
            nslots=args.nslots,
        )
    except Exception as exc:
        print(f"hmmsearch failed: {exc}", file=sys.stderr)
        sys.exit(1)
        
    ###########################################################################
    # 3.6. Write best-hit annotation table
    ###########################################################################

    annotation_table = os.path.join(args.output_dir, f"{args.sample_name}_orfs_fun_annot.tsv")

    try:
        write_best_hits(domtblout, annotation_table)
    except Exception as exc:
        print(f"Writing annotation table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ########################################################################### 
    # 3.7. Write output log and exit
    ###########################################################################
    
    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()
