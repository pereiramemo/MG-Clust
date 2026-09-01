#!/usr/bin/env python3
"""
mg-clust module 5: Functional annotation of ORFs using pyHMMER.

Assumes execution inside the conda environment "mg-clust-module-5" (or equivalent)
where pyhmmer is importable.

Updated for pyhmmer >= 0.11.0: Sequence.name/.accession/.description/.source are
str, not bytes, as of that release (breaking change, pyhmmer issue #88). Do not
reintroduce .decode() on these attributes if backporting to pyhmmer < 0.11.0.

- Searches ORF protein sequences against a HMM database (default: KEGG KO profiles)
- --db_mode selects the scoring strategy:
    ko:     KOfam adaptive per-KO thresholds (ko_list) with coverage-aware rescue
            of truncated ORFs and an e-value-only fallback for KOs with no
            cross-validated threshold (nt-KOs)
    ga:     per-model GA/TC/NC cutoffs embedded in the HMM file (bit_cutoffs="gathering")
    evalue: a single global --evalue_thres applied uniformly to all models
- Exports two tables: a collapsed one-KO-per-ORF summary (best hit per ORF) and a
  rich multi-KO-per-ORF table with every accepted (ORF, KO) pair
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import gzip
import os
import shutil
import sys
import tarfile
import urllib.request
import pyhmmer

sys.path.insert(0, os.path.dirname(__file__))
from utils import check_file, gzip_file

KO_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "ko", "ko_profiles.hmm")
KO_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"

KO_LIST_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "ko", "ko_list.tsv")
KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"

# ga mode: how many missing-model names to print when reporting a partial/total failure
GA_MISSING_SAMPLE_N = 10

# ga mode: if the first this-many models are ALL missing GA cutoffs, stop scanning and
# fail fast rather than churning through the rest of a large (e.g. KOfam-sized) database
# one model at a time only to reach the same "wrong mode" conclusion
GA_EARLY_FAIL_SAMPLE = 50

# best-hit ranking: lower rank wins; ties broken by each tier's own secondary key

TIER_RANK = {
    "score_confident": 3,
    "score_putative": 1,
    "evalue_confident": 2,
    "evalue_putative": 1,
    "ga_pass": 0,
    "evalue_pass": 0,
}

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

    parser.add_argument("--db_mode", dest="db_mode", required=True,
        choices=["ko", "ga", "evalue"],
        help="scoring strategy: ko (KOfam adaptive thresholds + truncation-aware "
             "rescue), ga (embedded per-model GA/TC/NC cutoffs), evalue (single "
             "global --evalue_thres). Required; no default")

    parser.add_argument("--ko_list", dest="ko_list", default=KO_LIST_DEFAULT,
        help=f"path to KOfam ko_list threshold file (default: {KO_LIST_DEFAULT}); "
             "only used, fetched, and required when --db_mode ko")

    parser.add_argument("--evalue_thres", dest="evalue_thres", type=float, default=1e-3,
        help="e-value threshold (default: 1e-3); sole criterion in --db_mode evalue, "
             "fallback bar for nt-KOs adjusted by --tighten_evalue_factor in "
             "--db_mode ko, unused in --db_mode ga")

    parser.add_argument("--edge_tol", dest="edge_tol", type=int, default=3,
        help="aa tolerance for alignment-flush check when detecting truncated ORFs "
             "(default: 3; --db_mode ko only)")

    parser.add_argument("--tighten_evalue_factor", dest="tighten_evalue_factor", type=float, default=100,
        help="factor by which the e-value bar is tightened for rescued truncated-ORF "
             "hits (default: 100; --db_mode ko only)")

    parser.add_argument("--min_orf_len", dest="min_orf_len", type=int, default=60,
        help="minimum ORF length (default: 60)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

# Dev only
"""
orfs_faa = "/home/epereira/workspace/repos/tools/MG-Clust/test/data/test_orfs.faa.gz"
hmm_db = "/home/epereira/.mg-clust/db/ko/ko_profiles.hmm"
db_mode = "ko"
ko_list = "/home/epereira/.mg-clust/db/ko/ko_list.tsv"
evalue_thres = 1e-3
edge_tol = 10
tighten_evalue_factor = 100
min_orf_len = 60
sample_name = "test"
nslots = 4
output_dir = "/home/epereira/workspace/repos/tools/MG-Clust/test/test_output/module-5/dev"
overwrite = True
"""

###############################################################################
# 2.2 Fetch/cache the KO HMM profiles (download/extract/concat; scratch cleaned
#     in a finally so a failed build strands nothing)
###############################################################################

def _fetch_ko_profiles(hmm_db: str) -> None:
    ko_dir = os.path.dirname(hmm_db)
    os.makedirs(ko_dir, exist_ok=True)
    """"""
    archive = os.path.join(ko_dir, "profiles.tar.gz")
    profiles_dir = os.path.join(ko_dir, "profiles")
    tmp_hmm_db = hmm_db + ".tmp"
    # The download is ~1.5 GB and the extracted profiles/ tree is larger still. Both
    # are scratch: only the concatenated hmm_db is kept. Clean them in a finally so an
    # interrupted or failed build cannot strand gigabytes in the cache directory.
    try:
        print(f"Downloading KO profiles from {KO_PROFILES_URL} ...")
        urllib.request.urlretrieve(KO_PROFILES_URL, archive)
        """"""
        print("Extracting profiles ...")
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(ko_dir)    
        """"""
        hmm_files = sorted(glob.glob(os.path.join(profiles_dir, "*.hmm")))
        if not hmm_files:
            raise RuntimeError("no .hmm files found after extraction")
        """"""
        print(f"Concatenating {len(hmm_files)} HMM profiles into {hmm_db} ...")
        with open(tmp_hmm_db, "wb") as out:
            for hmm_file in hmm_files:
                with open(hmm_file, "rb") as f:
                    shutil.copyfileobj(f, out)
        os.replace(tmp_hmm_db, hmm_db)
    finally:
        for path in (archive, tmp_hmm_db):
            if os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass
        if os.path.isdir(profiles_dir):
            shutil.rmtree(profiles_dir, ignore_errors=True)

###############################################################################
# 2.3 Fetch/cache the KOfam ko_list threshold file
###############################################################################

def _fetch_ko_list(ko_list: str) -> None:
    ko_dir = os.path.dirname(ko_list)
    os.makedirs(ko_dir, exist_ok=True)
    """"""
    archive = os.path.join(ko_dir, "ko_list.gz")
    print(f"Downloading ko_list from {KO_LIST_URL} ...")
    urllib.request.urlretrieve(KO_LIST_URL, archive)
    """"""
    print(f"Decompressing ko_list into {ko_list} ...")
    tmp_ko_list = ko_list + ".tmp"
    with gzip.open(archive, "rb") as gz_in, open(tmp_ko_list, "wb") as out:
        shutil.copyfileobj(gz_in, out)
    os.replace(tmp_ko_list, ko_list)
    """"""
    os.remove(archive)

###############################################################################
# 2.4 Ensure ko_profiles.hmm + ko_list are present as a matched pair (db_mode="ko")
###############################################################################

def ensure_ko_database(hmm_db: str, ko_list: str) -> None:
    if os.path.isfile(hmm_db) and os.path.isfile(ko_list):
        return
    """"""
    print("ko_profiles.hmm + ko_list cache incomplete; (re-)fetching both together ...")
    failures = []
    """"""
    try:
        _fetch_ko_profiles(hmm_db)
    except Exception as exc:
        failures.append(f"ko_profiles.hmm ({exc})")
    """"""
    try:
        _fetch_ko_list(ko_list)
    except Exception as exc:
        failures.append(f"ko_list ({exc})")
    """"""
    if failures or not (os.path.isfile(hmm_db) and os.path.isfile(ko_list)):
        print(f"KO database setup failed for: {'; '.join(failures) if failures else 'unknown reason'}",
              file=sys.stderr)
        sys.exit(1)
    """"""
    print("KO profiles + ko_list ready.")

###############################################################################
# 2.5 Parse ko_list into ko_id -> (threshold, score_type)
###############################################################################

def parse_ko_list(ko_list: str) -> dict:
    thresholds = {}
    with open(ko_list, "r", encoding="utf-8") as fh:
        next(fh, None)  # skip header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            ko_id, threshold_str, score_type = cols[0], cols[1], cols[2]
            if threshold_str == "-":
                thresholds[ko_id] = (None, None)
            else:
                thresholds[ko_id] = (float(threshold_str), score_type)
    return thresholds

###############################################################################
# 2.6 Load ORF sequences: lengths dict + digital sequence block for hmmsearch
###############################################################################

def load_orfs(orfs_faa: str):
    with pyhmmer.easel.SequenceFile(
        orfs_faa, format="fasta", digital=True, 
        alphabet=pyhmmer.easel.Alphabet.amino()
    ) as seq_file:
        # Prefetch the entire digitized set into a DigitalSequenceBlock in memory up front,
        # which according to pyhmmer documentation, improves performance at the cost of higher memory usage.
        block = seq_file.read_block() 
        # Compute the lengths of each ORF sequence in the block and save them in a dictionary
        orf_lengths = {seq.name: len(seq) for seq in block}
    return orf_lengths, block 

###############################################################################
# 2.7 Run hmmsearch: ga mode (per-model GA/TC/NC cutoffs)
###############################################################################

def run_ga_mode(hmm_db: str, block, nslots: int):
    # A single batched `hmmsearch(hmms, ...)` call aborts entirely (the generator
    # raises StopIteration on the next() call right after the first MissingCutoffs),
    # so per-model failures can't be skipped within one call -- each query HMM must
    # be searched individually to isolate failures from one another.
    missing_models = []
    total_models = 0
    all_hits = []
    """"""
    with pyhmmer.plan7.HMMFile(hmm_db) as hmms:
        for hmm in hmms:
            total_models += 1
            try:
                hits = list(pyhmmer.hmmer.hmmsearch(hmm, block, bit_cutoffs="gathering", cpus=nslots))
            except pyhmmer.errors.MissingCutoffs as exc:
                missing_models.append(getattr(exc, "model_name", None) or str(exc))
                if total_models >= GA_EARLY_FAIL_SAMPLE and len(missing_models) == total_models:
                    # The break applies when all 50 first hmm models are missing GA cutoffs.
                    break
                continue
            all_hits.extend(hits)
    """"""
    if missing_models:
        sample = missing_models[:GA_MISSING_SAMPLE_N]
        print(f"{len(missing_models)}/{total_models} models scanned are missing embedded "
              f"GA cutoffs (sample: {', '.join(sample)}"
              f"{', ...' if len(missing_models) > GA_MISSING_SAMPLE_N else ''})",
              file=sys.stderr)
        if len(missing_models) == total_models:
            print("No model in this HMM database has embedded GA cutoffs; --db_mode ga "
                  "is not usable against it. Use --db_mode ko or --db_mode evalue instead.",
                  file=sys.stderr)
            sys.exit(1)
    """"""
    return all_hits

###############################################################################
# 2.8 Run hmmsearch: evalue mode (E/domE bound at --evalue_thres; tiering happens in Python)
###############################################################################

def run_evalue_mode(hmm_db: str, block, evalue_thres: float, nslots: int):
    with pyhmmer.plan7.HMMFile(hmm_db) as hmms:
        # inclusion thresholds follow the reporting ones so `included` cannot mean
        # something looser than what this search actually accepted
        yield from pyhmmer.hmmer.hmmsearch(hmms, block,
                                           E=evalue_thres, domE=evalue_thres,
                                           incE=evalue_thres, incdomE=evalue_thres,
                                           cpus=nslots)

###############################################################################
# 2.9 Compute the fraction of the HMM model covered by a domain alignment
###############################################################################

def domain_coverage(aln, model_len: int) -> float:
    return (aln.hmm_to - aln.hmm_from + 1) / model_len

###############################################################################
# 2.10 Build rows for the simple pass/fail modes (ga, evalue)
###############################################################################

def rows_from_simple_hits(hits_iter, tier_label: str) -> list:
    rows = []
    for hits in hits_iter:
        ko_id = hits.query.name
        for hit in hits:
            dom = hit.best_domain
            aln = dom.alignment
            coverage = domain_coverage(aln, hits.query.M)
            rows.append({
                "orf_id": hit.name,
                "ko_id": ko_id,
                "score": hit.score,
                "evalue": hit.evalue,
                "coverage": coverage,
                "tier": tier_label,
            })
    return rows

###############################################################################
# 2.11 Classify a single (ORF, KO) hit against its ko_list threshold (db_mode="ko")
###############################################################################

def classify_ko_hit(orf_id, ko_id, hit, orf_len, model_len, threshold, score_type,
                     edge_tol, evalue_thres, tighten_evalue_factor):
    best_dom = hit.best_domain
    aln = best_dom.alignment
    """"""
    hmm_from, hmm_to = aln.hmm_from, aln.hmm_to
    target_from, target_to = aln.target_from, aln.target_to
    """"""
    coverage = domain_coverage(aln, model_len)
    missing_n = hmm_from > 1
    missing_c = hmm_to < model_len
    flush_n = target_from <= edge_tol
    flush_c = target_to >= orf_len - edge_tol
    truncation_explains_gap_case_1 = (missing_c and flush_c) and (not missing_n)
    truncation_explains_gap_case_2 = (missing_n and flush_n) and (not missing_c)
    truncation_explains_gap_case_3 = (missing_n and flush_n) and (missing_c and flush_c)
    truncation_explains_gap = truncation_explains_gap_case_1 or truncation_explains_gap_case_2 or truncation_explains_gap_case_3
    # The three cases above all say the same thing: the model is only partially
    # aligned, but every missing piece of the model falls off an end of the ORF,
    # so the gap is explained by a fragmentary ORF (contig edge) rather than by a
    # poor match. Legend: [===] aligned region, --- unaligned model remainder,
    # | a real sequence end.
    #
    #   case 1 - ORF truncated at its C-terminal end:
    #       model |[==========]------|   model C-term unaligned (missing_c) ...
    #       ORF   |[==========]|         ... but the alignment runs into the ORF
    #                                    end (flush_c), and the model N-term is
    #                                    fully covered (not missing_n)
    #
    #   case 2 - ORF truncated at its N-terminal end:
    #       model |------[==========]|   model N-term unaligned (missing_n) ...
    #       ORF         |[==========]|   ... but the alignment starts at the ORF
    #                                    start (flush_n), and the model C-term is
    #                                    fully covered (not missing_c)
    #
    #   case 3 - ORF is an internal fragment, truncated at both ends:
    #       model |----[========]----|   both model ends unaligned, and the
    #       ORF        |[========]|      alignment is flush with both ORF ends
    #
    # Anything else - a model end missing while the alignment stops well inside
    # the ORF (further than edge_tol from its boundary) - is a genuinely partial
    # match and is not rescued.
    """"""
    score = best_dom.score
    evalue = best_dom.i_evalue
    tier = "evalue_pass"
    """"""
    # Check score level confidence
    if threshold is not None:
        # Determine the right score according to calibration
        if score_type == "full":
            score = hit.score
        if truncation_explains_gap:
            tier = "score_putative"
        if evalue < evalue_thres / tighten_evalue_factor:
            tier = "evalue_confident"
        if score >= threshold:
            tier = "score_confident"    
    # Check evalue level confidence
    else:
        if truncation_explains_gap:
            tier = "evalue_putative"
        if evalue < evalue_thres / tighten_evalue_factor:
            tier = "evalue_confident"
    """"""
    # Return the classification result as a dictionary
    return {
        "orf_id": orf_id,
        "ko_id": ko_id,
        "score": score,
        "evalue": evalue,
        "coverage": coverage,
        "tier": tier,
        "calibration": f"{score_type}, {threshold}"
    }

###############################################################################
# 2.12 Build rows for db_mode="ko" over all (ORF, KO) hits
###############################################################################

def rows_from_ko_hits(hits_iter, orf_lengths, ko_thresholds, edge_tol,
                       evalue_thres, tighten_evalue_factor) -> list:
    rows = []
    for hits in hits_iter:
        ko_id = hits.query.name
        model_len = hits.query.M
        threshold, score_type = ko_thresholds.get(ko_id, (None, None))
        # The following loop iterates over all reported hits for the current KO.
        # This is juts a safeward, given that the same evalues are used in for the
        # reported and included hits. 
        for hit in hits.reported:
            dom = next(iter(hit.domains.reported), None)
            if dom is None:
                continue
            # The following filtering is added to condition selected hits based 
            # on the i-Evalue the domE parameter in hmmsearch is applied only on 
            # the c-Evalue. See https://manpages.ubuntu.com/manpages/focal/man1/phmmer.1.html.
            if hit.best_domain.i_evalue > evalue_thres:
                continue
            orf_id = hit.name
            orf_len = orf_lengths.get(orf_id)
            if orf_len is None:
                # flag an alarm
                print(f"Warning: ORF length not found for {orf_id}", 
                      file=sys.stderr)
                continue
            row = classify_ko_hit(orf_id, ko_id, hit, orf_len, model_len, threshold, score_type,
                                   edge_tol, evalue_thres, tighten_evalue_factor)
            if row is not None:
                rows.append(row)
    return rows

###############################################################################
# 2.13 Rank rows (best-hit selection across competing KOs per ORF)
###############################################################################

def rank_key(row: dict, ko_thresholds: dict) -> tuple:
    tier = row["tier"]
    rank = TIER_RANK[tier]
    if tier in ("score_confident", "score_putative"):
        threshold, _ = ko_thresholds[row["ko_id"]]
        margin = row["score"] - threshold
        return (rank, margin)
    # tie-break by e-value
    return (rank, -row["evalue"])

###############################################################################
# 2.14 Collapse to one best KO per ORF
###############################################################################

def collapse_best_per_orf(rows: list, ko_thresholds: dict) -> list:
    best = {}
    for row in rows:
        orf_id = row["orf_id"]
        if orf_id not in best:
            best[orf_id] = row
        else:
            if rank_key(row, ko_thresholds) > rank_key(best[orf_id], ko_thresholds):
                best[orf_id] = row
    return [best[orf_id] for orf_id in sorted(best)]

###############################################################################
# 2.15 Write a rows table to disk (headerless, matches module 6's raw concatenation)
###############################################################################

def write_table(rows: list, path: str) -> None:
    with open(path, "w") as out:
        for row in sorted(rows, key=lambda r: (r["orf_id"], r["ko_id"])):
            out.write(
                f"{row['orf_id']}\t{row['ko_id']}\t{row['evalue']:.2e}\t"
                f"{row['score']:.2f}\t{row['calibration']}\t"
                f"{row['coverage']:.4f}\t{row['tier']}\n"
            )

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    args = parse_args()
    orfs_faa = args.orfs_faa
    hmm_db = args.hmm_db
    db_mode = args.db_mode
    ko_list = args.ko_list
    output_dir = args.output_dir
    overwrite = args.overwrite
    evalue_thres = args.evalue_thres
    edge_tol = args.edge_tol
    tighten_evalue_factor = args.tighten_evalue_factor
    nslots = args.nslots
    sample_name = args.sample_name
    min_orf_len = args.min_orf_len

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    check_file(orfs_faa, "ORF protein FASTA file")

    ###########################################################################
    # 3.2. Check/fetch the HMM database (behavior depends on --db_mode)
    ###########################################################################

    if db_mode == "ko":
        ensure_ko_database(hmm_db, ko_list)
    else:
        check_file(hmm_db, "HMM database")

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
    # 3.5. Load ORF sequences
    ###########################################################################

    try:
        orf_lengths, block = load_orfs(orfs_faa)
    except Exception as exc:
        print(f"loading ORF sequences failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Run hmmsearch and classify hits, dispatching on --db_mode
    ###########################################################################

    ko_thresholds = {}

    try:
        if db_mode == "ko":
            ko_thresholds = parse_ko_list(ko_list)
            hits_iter = run_evalue_mode(hmm_db, block, evalue_thres, nslots)
            rows = rows_from_ko_hits(hits_iter, orf_lengths, ko_thresholds, edge_tol,
                                      evalue_thres, tighten_evalue_factor)
        elif db_mode == "ga":
            hits = run_ga_mode(hmm_db, block, nslots)
            rows = rows_from_simple_hits(hits, "ga_pass")
        else:
            hits_iter = run_evalue_mode(hmm_db, block, evalue_thres, nslots)
            rows = rows_from_simple_hits(hits_iter, "evalue_pass")
    except Exception as exc:
        print(f"hmmsearch failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Collapse to best hit per ORF
    ###########################################################################

    best_rows = collapse_best_per_orf(rows, ko_thresholds)

    ###########################################################################
    # 3.8. Write best hits (one-KO-per-ORF) and rich (multi-KO-per-ORF) tables
    ###########################################################################

    best_hits_table = os.path.join(
        output_dir, f"{sample_name}_orfs-minlen{min_orf_len}aa-fun_annot.tsv")
    # all_table = os.path.join(
    #    output_dir, f"{sample_name}_orfs-minlen{min_orf_len}aa-fun_annot_all.tsv")

    try:
        write_table(best_rows, best_hits_table)
     #   write_table(rows, all_table)
    except Exception as exc:
        print(f"writing annotation tables failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.9. Compress the annotation table
    ###########################################################################

    # Consumed by module 6 via raw byte-copy concatenation and DuckDB, both gzip-safe.
    gzip_file(best_hits_table)

    ###########################################################################
    # 3.10. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###############################################################################
# 4. Run the main function
###############################################################################

if __name__ == "__main__":
    main()
