#!/usr/bin/env python3

"""
mg-clust database pre-fetch: populate ~/.mg-clust/db before a pipeline run.

Modules 4 and 5 need GTDB (~1.1 GB) and the KEGG KO profiles + ko_list (~7.6 GB).
Fetching them from inside the pipeline means every per-sample task races for the
same download into the same fixed scratch paths, and many clusters block outbound
HTTPS from compute nodes entirely. Run this once on the submit/login node before
the first real run.

The KO half is pure stdlib and runs anywhere. The GTDB half shells out to mmseqs,
which lives only in the module-4 container, so there are three ways to run this:

    download_dbs.py --databases ko              # no mmseqs needed
    download_dbs.py --container                 # KO natively, GTDB in the image
    download_dbs.py                             # all native; needs mmseqs on PATH

--container requires the module-4 image to be cached already: run
bin/download_imgs.py first.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file
# Reused so the SIF filename is derived in exactly one place; see --container.
from download_imgs import APPTAINER_CACHE_DIR, IMAGE_TEMPLATE, cache_filename

# Define output paths and URLs
GTDB_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "gtdb", "gtdb")

KO_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "ko", "ko_profiles.hmm")
KO_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"

KO_LIST_DEFAULT = os.path.join(os.path.expanduser("~"), ".mg-clust", "db", "ko", "ko_list.tsv")
KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"

# Define tools
mmseqs = "mmseqs"

# The only module image carrying mmseqs, and so the one --container re-execs into.
MMSEQS_MODULE = 4

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: pre-fetch the MG-Clust reference databases",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--databases", dest="databases", default="both",
        choices=["gtdb", "ko", "both"],
        help="which databases to fetch; 'gtdb' requires mmseqs on PATH, 'ko' does "
             "not (default: both)")

    parser.add_argument("--gtdb", dest="gtdb", default=GTDB_DEFAULT,
        help=f"path prefix for the MMseqs2 GTDB taxonomy database (default: {GTDB_DEFAULT})")

    parser.add_argument("--ko_db", dest="ko_db", default=KO_DEFAULT,
        help=f"path to the concatenated KO HMM profiles (default: {KO_DEFAULT})")

    parser.add_argument("--ko_list", dest="ko_list", default=KO_LIST_DEFAULT,
        help=f"path to the KOfam ko_list threshold file (default: {KO_LIST_DEFAULT})")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads passed to mmseqs databases (default: 4)")

    parser.add_argument("--container", dest="container", action="store_true", default=False,
        help=f"run the GTDB fetch inside the cached module-{MMSEQS_MODULE} Apptainer "
             "image instead of requiring mmseqs on PATH (default: False)")

    parser.add_argument("--force", dest="force", action="store_true", default=False,
        help="re-download databases that are already cached (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Fetch/cache the KO HMM profiles (download/extract/concat; scratch cleaned
#     in a finally so a failed build strands nothing)
###############################################################################

def _fetch_ko_profiles(ko_db: str) -> None:
    ko_dir = os.path.dirname(ko_db)
    os.makedirs(ko_dir, exist_ok=True)

    archive = os.path.join(ko_dir, "profiles.tar.gz")
    profiles_dir = os.path.join(ko_dir, "profiles")
    tmp_hmm_db = ko_db + ".tmp"
    # The download is ~1.5 GB and the extracted profiles/ tree is larger still. Both
    # are scratch: only the concatenated hmm_db is kept. Clean them in a finally so an
    # interrupted or failed build cannot strand gigabytes in the cache directory.
    try:
        print(f"Downloading KO profiles from {KO_PROFILES_URL} ...")
        urllib.request.urlretrieve(KO_PROFILES_URL, archive)

        print("Extracting profiles ...")
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(ko_dir)

        hmm_files = sorted(glob.glob(os.path.join(profiles_dir, "*.hmm")))
        if not hmm_files:
            raise RuntimeError("no .hmm files found after extraction")

        print(f"Concatenating {len(hmm_files)} HMM profiles into {ko_db} ...")
        with open(tmp_hmm_db, "wb") as out:
            for hmm_file in hmm_files:
                with open(hmm_file, "rb") as f:
                    shutil.copyfileobj(f, out)
        os.replace(tmp_hmm_db, ko_db)
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

    archive = os.path.join(ko_dir, "ko_list.gz")
    tmp_ko_list = ko_list + ".tmp"
    # Same finally discipline as _fetch_ko_profiles above: the archive and the
    # partial decompression are scratch, and a failure between download and
    # os.replace would otherwise leave both behind.
    try:
        print(f"Downloading ko_list from {KO_LIST_URL} ...")
        urllib.request.urlretrieve(KO_LIST_URL, archive)

        print(f"Decompressing ko_list into {ko_list} ...")
        with gzip.open(archive, "rb") as gz_in, open(tmp_ko_list, "wb") as out:
            shutil.copyfileobj(gz_in, out)
        os.replace(tmp_ko_list, ko_list)
    finally:
        for path in (archive, tmp_ko_list):
            if os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

###############################################################################
# 2.4 Ensure ko_profiles.hmm + ko_list are present as a matched pair
###############################################################################

def ensure_ko_database(ko_db: str, ko_list: str, force: bool = False) -> None:
    # Checked and fetched as a pair.
    if os.path.isfile(ko_db) and os.path.isfile(ko_list) and not force:
        print(f"KO profiles + ko_list already cached in {os.path.dirname(ko_db)}; skipping")
        return

    print("Fetching KO profiles + ko_list together ...")
    failures = []

    try:
        _fetch_ko_profiles(ko_db)
    except Exception as exc:
        failures.append(f"ko_profiles.hmm ({exc})")

    try:
        _fetch_ko_list(ko_list)
    except Exception as exc:
        failures.append(f"ko_list ({exc})")

    if failures or not (os.path.isfile(ko_db) and os.path.isfile(ko_list)):
        print(f"KO database setup failed for: {'; '.join(failures) if failures else 'unknown reason'}",
              file=sys.stderr)
        sys.exit(1)

    print("KO profiles + ko_list ready.")

###############################################################################
# 2.5 Fetch/cache the GTDB MMseqs2 taxonomy database
###############################################################################

def fetch_gtdb_database(gtdb: str, nslots: int, force: bool = False) -> None:

    # A valid mmseqs2 database always has a companion .dbtype file; that is the
    # presence test module 4 uses.
    if os.path.isfile(gtdb + ".dbtype") and not force:
        print(f"GTDB database already cached at {gtdb}; skipping")
        return

    # Only needed to actually download. Checked after the presence test so that a
    # warm cache does not require mmseqs at all -- it lives in the module-4 image,
    # and the common case is re-running this on a login node that has neither.
    check_tools([mmseqs])

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
    finally:
        # Several GB of scratch. Drop it whether or not the download succeeded,
        # so an interrupted run does not strand it in the shared cache.
        if os.path.isdir(download_tmp):
            shutil.rmtree(download_tmp, ignore_errors=True)

    print("GTDB database download complete.")

###############################################################################
# 2.6 Re-run the GTDB fetch inside the module-4 image (--container)
###############################################################################

def exec_gtdb_in_container(gtdb: str, nslots: int, force: bool) -> None:
    """Replace this process with the same fetch running inside the module-4 image.

    mmseqs exists only in that image. os.execvp does not return, so the container's
    exit status becomes this script's directly -- which is also why the caller must
    already have finished the KO half before calling this.

    --container is deliberately not forwarded to the child, so the inner run takes
    the ordinary native path and cannot recurse.
    """
    check_tools(["apptainer"])

    img = IMAGE_TEMPLATE.format(n=MMSEQS_MODULE)
    sif = os.path.join(APPTAINER_CACHE_DIR, cache_filename(img))
    if not os.path.isfile(sif):
        print(f"module-{MMSEQS_MODULE} image not cached at {sif}; "
              "run bin/download_imgs.py first", file=sys.stderr)
        sys.exit(1)

    script = os.path.abspath(__file__)
    # Apptainer auto-mounts $HOME and the cwd, which is enough while the checkout
    # and the cache both live under $HOME. Bind them explicitly anyway so this keeps
    # working when either sits on /project or /scratch, as it usually does on a
    # cluster. Binding a path that is already visible is harmless.
    binds = sorted({os.path.dirname(os.path.dirname(script)),
                    os.path.dirname(gtdb)})
    for b in binds:
        os.makedirs(b, exist_ok=True)

    cmd = ["apptainer", "exec"]
    for b in binds:
        cmd += ["--bind", b]
    cmd += [sif, script,
            "--databases", "gtdb",
            "--gtdb", gtdb,
            "--nslots", str(nslots)]
    if force:
        cmd.append("--force")

    print(f"Running the GTDB fetch inside {os.path.basename(sif)} ...")
    # execvp replaces the process image without unwinding Python, so anything still
    # sitting in the stdout buffer is lost -- and stdout is block-buffered whenever
    # it is a pipe rather than a terminal. Without this flush every message printed
    # before the hand-off (the KO result included) silently disappears when the
    # output is piped or redirected.
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        os.execvp("apptainer", cmd)
    except OSError as exc:
        print(f"failed to exec apptainer: {exc}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    args = parse_args()
    databases = args.databases 
    ko_db = args.ko_db
    ko_list = args.ko_list
    gtdb = args.gtdb
    nslots = args.nslots
    force = args.force
    container = args.container

    want_ko = databases in ("ko", "both")
    want_gtdb = databases in ("gtdb", "both")

    ###########################################################################
    # 3.1. Download databases
    ###########################################################################

    if want_ko:
        ensure_ko_database(ko_db, ko_list, force)

    if want_gtdb:
        # Last, because exec_gtdb_in_container never returns: the KO half above
        # must already be done by the time we hand off.
        if container:
            exec_gtdb_in_container(gtdb, nslots, force)
        fetch_gtdb_database(gtdb, nslots, force)

    ###########################################################################
    # 3.2. Check files
    ###########################################################################

    if want_ko:
        check_file(ko_db, "KO HMM database")
        check_file(ko_list, "KO list")

    if want_gtdb:
        check_file(gtdb + ".dbtype", "GTDB database")

    ###########################################################################
    # 3.3. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###############################################################################
# 4. Run the main function
###############################################################################

if __name__ == "__main__":
    main()
