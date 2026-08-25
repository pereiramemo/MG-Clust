###############################################################################
# 1. Define utility functions
###############################################################################

import gzip
import os
import shutil
import subprocess
import sys
from typing import List, Optional

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Run a command and check return code; optionally redirect stdout to file
###############################################################################

def run(cmd: List[str], stdout_path: Optional[str] = None) -> None:
    """Run a command, stream output or redirect to file, and fail on non-zero return code."""
    try:
        if stdout_path:
            with open(stdout_path, "w") as out_f:
                proc = subprocess.run(cmd, stdout=out_f, check=False)
        else: 
            proc = subprocess.run(cmd, check=False)
    except FileNotFoundError as exc:
        print(f"Command not found: {cmd[0]} ({exc})", file=sys.stderr)
        sys.exit(1)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)

###############################################################################
# 2.2 Check that required tools are available on PATH
###############################################################################

def check_tools(tools: List[str]) -> None:
    missing = [t for t in tools if subprocess.run(["which", t], capture_output=True).returncode != 0]
    if missing:
        print(f"Missing tools: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.3 Check a file exists; exit with error if not
###############################################################################

def check_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"{label} is not a real file", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.4 Compress a file in place, returning the new path
###############################################################################

def gzip_file(path: str, level: int = 6) -> str:
    """Replace path with path + ".gz" and return the new path.

    Used on module outputs that are published and then read again downstream.
    Every consumer in this pipeline handles gzip transparently: DuckDB read_csv,
    bbduk in=/out=, mmseqs createdb, and pyhmmer SequenceFile. Concatenating the
    results with a raw byte copy also stays valid, because concatenated gzip
    members are themselves a well-formed gzip stream.

    Writes to a .tmp first and only then replaces, so an interrupted compression
    cannot leave a truncated .gz standing in for real output.
    """
    gz_path = path + ".gz"
    tmp_path = gz_path + ".tmp"
    try:
        with open(path, "rb") as f_in, gzip.open(tmp_path, "wb", compresslevel=level) as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.replace(tmp_path, gz_path)
        os.remove(path)
    except Exception as exc:
        if os.path.isfile(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        print(f"compressing {path} failed: {exc}", file=sys.stderr)
        sys.exit(1)
    return gz_path
