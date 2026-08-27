"""Locate the module and the default inputs without hard-coding anyone's disk.

Every script in this directory imports `cutsearch_design` from beside itself and
resolves default input FASTAs relative to the repository root, found by walking
up from this file until a `.git` directory appears.  Nothing here depends on
where the repo is checked out, or on the directory ever being called `final`.

Defaults are a convenience only: each script also takes explicit paths on the
command line, which is the portable way to run them on a machine that does not
have this repo's untracked data directories.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)          # so `import cutsearch_design` works


def repo_root(start=HERE):
    """Nearest ancestor containing .git; falls back to this directory."""
    d = start
    while True:
        if os.path.isdir(os.path.join(d, ".git")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return start
        d = parent


ROOT = repo_root()


def at(*parts):
    return os.path.join(ROOT, *parts)


# Known alignments in this repo, as (label, path).  Missing ones are skipped by
# `existing()` rather than crashing, since several live in untracked data dirs.
CLUSTER1 = at("resources", "260629_issue84_algoplanning", "stage0_cluster1",
              "out", "clusters", "1", "c1.core.aln.fasta")
CLUSTER1_OLD = at("ninetypidorfs", "cluster1.core.aln.fasta")
CLUSTER2_OLD = at("ninetypidorfs", "cluster2.core.aln.fasta")

DEFAULT_INPUTS = (("cluster2 (old)", CLUSTER2_OLD),
                  ("cluster1 (old)", CLUSTER1_OLD),
                  ("c1 stage0", CLUSTER1))


def existing(pairs=DEFAULT_INPUTS):
    return [(lab, p) for lab, p in pairs if os.path.exists(p)]


def inputs_from_argv(argv, pairs=DEFAULT_INPUTS):
    """Explicit paths if given, else whichever defaults are present."""
    if len(argv) > 1:
        return [(os.path.basename(a), a) for a in argv[1:]]
    got = existing(pairs)
    if not got:
        sys.exit("no input alignments found; pass one or more paths, e.g.\n"
                 f"  python {os.path.basename(argv[0])} path/to/cores.aln.fasta")
    return got
