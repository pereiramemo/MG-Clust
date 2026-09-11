#!/usr/bin/env python3

"""
mg-clust image pre-fetch: populate the container image cache before a pipeline run.

Letting the pipeline pull its own images means every per-sample task races for the
same download on a cold cache, and many clusters block outbound HTTPS from compute
nodes entirely. Run this once on the submit/login node before the first real run.

Apptainer images are written into the directory the `slurm` profile declares as
`apptainer.cacheDir`, under the exact filenames Nextflow derives from each image
URI -- see cache_filename() -- because Nextflow looks the cache up by name and
will silently pull its own copies if they do not match.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import shutil
import signal
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools, check_file

# Number of pipeline modules. Must track docker/dockerbuild_commands.sh.
N_MODULES = 6
IMAGE_TEMPLATE = "ghcr.io/pereiramemo/mg-clust/module-{n}:latest"

# Must match apptainer.cacheDir in nextflow.config, which is the single source of
# truth: Nextflow's SingularityCache.getCacheDir() consults the config value first
# and only falls back to NXF_APPTAINER_CACHEDIR when it is unset. Since the slurm
# profile always sets it, that environment variable is dead for the pipeline --
# honouring it here would make this script pre-pull into a directory Nextflow
# never reads. Use --cache_dir if you change the config.
APPTAINER_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".mg-clust", "apptainer")

# Apptainer unpacks each layer into a rootfs under its temp dir before assembling
# the SIF, so it needs several GB of scratch per image. That defaults to /tmp,
# which is routinely a small separate partition (and on many clusters a tmpfs in
# RAM), giving "no space left on device" mid-unpack. Keep the scratch beside the
# cache instead, where the images are going to live anyway.
APPTAINER_TMP_DIR = os.path.join(os.path.expanduser("~"), ".mg-clust", "tmp")

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: pre-fetch MG-Clust container images",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--engine", dest="engine", default="apptainer",
        choices=["apptainer", "docker", "both"],
        help="container runtime to fetch images for (default: apptainer)")

    parser.add_argument("--cache_dir", dest="cache_dir", default=APPTAINER_CACHE_DIR,
        help="directory to write Apptainer images into; must match apptainer.cacheDir "
             f"in nextflow.config or Nextflow will not find them (default: {APPTAINER_CACHE_DIR})")

    parser.add_argument("--tmp_dir", dest="tmp_dir", default=APPTAINER_TMP_DIR,
        help="scratch directory Apptainer unpacks layers into; needs several GB free "
             f"and must not be a small /tmp partition (default: {APPTAINER_TMP_DIR})")

    parser.add_argument("--force", dest="force", action="store_true", default=False,
        help="re-pull images that are already cached (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Build the list of module images
###############################################################################

def module_images() -> list:
    return [IMAGE_TEMPLATE.format(n=n) for n in range(1, N_MODULES + 1)]

###############################################################################
# 2.3 Derive the filename Nextflow expects for a cached image
###############################################################################

def cache_filename(img: str) -> str:
    """Name a cached image the way Nextflow's SingularityCache.simpleName() does.

    Nextflow looks the cache up by filename, so an image pre-pulled under any 
    other name is invisible to it and gets pulled a second time. For the plain 
    registry references this pipeline uses -- see IMAGE_TEMPLATE -- the rule is
    just: replace ":" and "/" with "-", append ".img". ApptainerCache inherits
    it unchanged from SingularityCache.

    ".img" here is only a filename, not a format: the file Apptainer writes is a
    SIF whatever the target is called. 
    """
    return img.replace(":", "-").replace("/", "-") + ".img"

###############################################################################
# 2.4 Fetch Apptainer images
###############################################################################

def fetch_apptainer_imgs(cache_dir: str, tmp_dir: str, force: bool) -> None:
    for d in (cache_dir, tmp_dir):
        try:
            os.makedirs(d, exist_ok=True)
        except Exception:
            print(f"mkdir {d} failed", file=sys.stderr)
            sys.exit(1)

    work_tmp = tempfile.mkdtemp(prefix="apptainer-", dir=tmp_dir)
    os.environ["APPTAINER_TMPDIR"] = work_tmp
    os.environ["TMPDIR"] = work_tmp

    try:
        _pull_all(cache_dir, force)
    finally:
        # Reached on success, on sys.exit, and on SIGINT/SIGTERM alike.
        shutil.rmtree(work_tmp, ignore_errors=True)


def _pull_all(cache_dir: str, force: bool) -> None:
    for img in module_images():
        cache_path = os.path.join(cache_dir, cache_filename(img))

        # An Apptainer pull is a full OCI-to-SIF conversion with no up-to-date 
        # check, and it refuses outright when the target exists ("Image file
        # already exists - will not overwrite", exit 255), so skipping is what
        # makes this script re-runnable at all. The cost is that a rebuilt
        # :latest is not picked up until --force.
        if os.path.isfile(cache_path) and not force:
            print(f"{cache_path} already cached; skipping")
            continue

        print(f"Pulling {img} -> {cache_path} ...")
        cmd = ["apptainer", "pull"]
        if force:
            cmd.append("--force")
        cmd += [cache_path, f"docker://{img}"]
        try:
            run(cmd)
        except subprocess.CalledProcessError:
            print(f"apptainer pull {img} failed", file=sys.stderr)
            sys.exit(1)

        check_file(cache_path, f"Apptainer image {cache_path}")

###############################################################################
# 2.5 Fetch Docker images
###############################################################################

def fetch_docker_imgs() -> None:
    # No --force equivalent is needed. Unlike an Apptainer pull, `docker pull`
    # checks the registry manifest and returns almost immediately when the image
    # is current, so pulling unconditionally keeps re-runs cheap AND picks up a
    # rebuilt :latest -- there is nothing for --force to change.
    for img in module_images():
        print(f"Pulling {img} ...")
        try:
            run(["docker", "pull", img])
        except subprocess.CalledProcessError:
            print(f"docker pull {img} failed", file=sys.stderr)
            sys.exit(1)

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    args = parse_args()
    engine = args.engine
    cache_dir = args.cache_dir
    force = args.force
    tmp_dir = args.tmp_dir

    ###########################################################################
    # 3.1. Download images
    ###########################################################################

    if engine in ("apptainer", "both"):
        check_tools(["apptainer"])
        fetch_apptainer_imgs(cache_dir, tmp_dir, force)

    if engine in ("docker", "both"):
        check_tools(["docker"])
        fetch_docker_imgs()

    ###########################################################################
    # 3.2. Write output log and exit
    ###########################################################################

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###############################################################################
# 4. Run the main function
###############################################################################

if __name__ == "__main__":
    main()
