#!/usr/bin/env python3
"""Run viamd_nci_reference over a tree of VeloxChem calculations, then pool the expected regions.

Finds every *.scf.results.h5 or plain name.h5 with an SCF wavefunction under ROOT (no *.scf.tensors.h5), keeps one calculation per molecule (by formula, or per
distinct geometry with --per-geometry), preferring *.scf.results.h5, skipping NTO files, and runs the
CLI on each. Resumable: a molecule whose -density.csv exists is skipped. When a file cannot be used the next
candidate of the same molecule is tried. Ends by running nci_reference_regions.py over all outputs.

Without --settings the CLI uses the uniform reference settings (s <= 1.0, rho <= 0.07, R 2.7 Bohr for every CP, each
sample to its nearest CP). With --settings TAG:S:RHO[:R] (repeatable) every CP uses s <= S, rho <= RHO and radius R
(default 2.7), one NAME-TAG-density.csv (and NAME-TAG-cp-density.csv) per set from a single critical point search,
pooled into OUT_DIR/regions-TAG.

Usage: run_batch.py ROOT OUT_DIR [--cli PATH] [--spacing BOHR] [--max-atoms N] [--timeout S] [--jobs N]
                    [--settings TAG:S:RHO[:R] ...] [--per-geometry] [--exclude SUBSTRING ...] [--dry-run]
"""
import argparse
import collections
import hashlib
import os
import re
import subprocess
import sys
import time
import warnings

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Ground-state files only: VeloxChem's *.scf.results.h5 or a plain name.h5 (no further suffix), not
# *.scf.tensors.h5, response solutions or natural transition orbitals
SKIP = re.compile(r"_nto_|_nto\.h5$")


def wanted(filename):
    f = filename.lower()
    if not f.endswith(".h5") or SKIP.search(f):
        return False
    return f.endswith(".scf.results.h5") or "." not in f[:-len(".h5")]


def inventory(root):
    """[(path, natoms, formula, geometry signature)] for files with an SCF wavefunction."""
    rows = []
    for dp, _, fn in os.walk(root):
        for f in sorted(fn):
            if not wanted(f):
                continue
            p = os.path.join(dp, f)
            try:
                with h5py.File(p, "r") as h:
                    keys = set()
                    h.visit(keys.add)
                    if not any(k.endswith("C_alpha") or k.endswith("orbital/alpha/coefficient") for k in keys):
                        continue
                    z = next((k for k in keys if k.endswith("nuclear_charges")), None)
                    xyz = next((k for k in keys if k.endswith("atom_coordinates")), None)
                    if not z or not xyz:
                        continue
                    Z = np.array(h[z]).astype(int).ravel()
                    X = np.array(h[xyz]).ravel()
            except Exception:
                continue
            formula = "".join(f"{k}_{v}." for k, v in sorted(collections.Counter(Z.tolist()).items()))
            geometry = hashlib.md5(Z.tobytes() + np.round(X, 3).tobytes()).hexdigest()[:10]
            rows.append((p, len(Z), formula, geometry))
    return rows


def suffixes(args):
    return [f"-{t}-density.csv" for t in args.tags] if args.tags else ["-density.csv"]


def is_done(args, root, path):
    return all(os.path.isfile(os.path.join(args.out, name_of(root, path) + x)) for x in suffixes(args))


def name_of(root, path):
    return os.path.relpath(path, root)[:-len(".h5")].replace(os.sep, "_").replace(" ", "_")


def run_one(args, root, candidates):
    """Runs the first usable candidate; returns (name, status, seconds)."""
    for path in candidates:
        name = name_of(root, path)
        prefix = os.path.join(args.out, name)
        log = os.path.join(args.out, "log", name + ".log")
        t = time.time()
        try:
            with open(log, "w") as f:
                settings = [a for spec in args.settings for a in ("--settings", spec)]
                r = subprocess.run([args.cli, "--spacing", str(args.spacing), *settings, path, prefix], stdout=f,
                                   stderr=subprocess.STDOUT, timeout=args.timeout)
        except subprocess.TimeoutExpired:
            return name, "timeout", time.time() - t
        if r.returncode == 0:
            return name, "ok", time.time() - t
    return name_of(root, candidates[0]), "unusable", 0.0


def main():
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("out")
    ap.add_argument("--cli", default=os.path.join(HERE, "..", "..", "build", "bin", "viamd_nci_reference"))
    ap.add_argument("--spacing", type=float, default=0.2)
    ap.add_argument("--max-atoms", type=int, default=0, help="skip larger molecules (0: no limit)")
    ap.add_argument("--timeout", type=float, default=0, help="seconds per molecule (0: none)")
    ap.add_argument("--jobs", type=int, default=1, help="molecules run at once (each run is multithreaded)")
    ap.add_argument("--settings", action="append", default=[], metavar="TAG:S:RHO[:R]",
                    help="uniform s, rho (and radius) for every category (repeatable)")
    ap.add_argument("--per-geometry", action="store_true", help="one run per distinct geometry, not per formula")
    ap.add_argument("--exclude", nargs="*", default=[], help="skip paths containing any of these")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    args.timeout = args.timeout or None
    args.tags = [spec.split(":")[0] for spec in args.settings]
    root = os.path.abspath(args.root)
    if not args.dry_run and not os.access(args.cli, os.X_OK):
        print(f"{args.cli}: not found; build the viamd_nci_reference target (VIAMD_ENABLE_QM=ON) or pass --cli")
        return 1
    os.makedirs(os.path.join(args.out, "log"), exist_ok=True)

    groups = collections.defaultdict(list)
    for path, natoms, formula, geometry in inventory(root):
        if any(e in path for e in args.exclude) or (args.max_atoms and natoms > args.max_atoms):
            continue
        groups[geometry if args.per_geometry else formula].append((path, natoms))
    todo, done = [], []
    for key, members in sorted(groups.items(), key=lambda kv: kv[1][0][1]):
        # *.scf.results.h5 first, then the smallest file
        members.sort(key=lambda m: (not m[0].lower().endswith(".scf.results.h5"), os.path.getsize(m[0])))
        candidates = [m[0] for m in members]
        if any(is_done(args, root, c) for c in candidates):
            done.append(candidates[0])
        else:
            todo.append((members[0][1], candidates))
    print(f"{len(groups)} molecules: {len(done)} already done, {len(todo)} to run")
    for natoms, candidates in todo:
        print(f"  {natoms:6d} atoms  {os.path.relpath(candidates[0], root)}"
              + (f"  (+{len(candidates) - 1} alternatives)" if len(candidates) > 1 else ""))
    if args.dry_run:
        return 0

    status = collections.Counter()
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for name, result, seconds in pool.map(lambda t: run_one(args, root, t[1]), todo):
            status[result] += 1
            print(f"  {result:8s} {seconds:8.1f} s  {name}", flush=True)
    print("runs: " + ", ".join(f"{k} {v}" for k, v in status.items()) if status else "nothing to run")

    code = 0
    for tag in args.tags or [""]:
        if tag:
            densities = [f for f in os.listdir(args.out) if f.endswith(f"-{tag}-density.csv")]
        else:
            densities = [f for f in os.listdir(args.out) if f.endswith("-density.csv") and not f.endswith("-cp-density.csv")]
        if not densities:
            code = 1
            continue
        print(f"== {tag or 'default settings'}: {len(densities)} files", flush=True)
        code |= subprocess.call([sys.executable, os.path.join(HERE, "nci_reference_regions.py"),
                                 os.path.join(args.out, f"regions-{tag}" if tag else "regions")]
                                + sorted(os.path.join(args.out, f) for f in densities))
    return code


if __name__ == "__main__":
    sys.exit(main())
