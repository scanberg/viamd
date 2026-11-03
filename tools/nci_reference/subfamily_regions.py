#!/usr/bin/env python3
"""Pool reference exports per CP subfamily (the CP-tree branch under a reference category), for the reference's sub-maps.

A density export holds one volume per category, not per CP. A system therefore contributes to a subfamily only when
all its CPs of that category belong to that one branch (in its -anchors.csv companion); mixed categories are skipped.
Each subfamily is pooled as nci_reference_regions.py pools a category: every system weighing the same (--group as
there), smoothed, saved as sub_<n>.npy. A subfamily is kept with at least --min-systems systems; summary.json lists
the kept ones with their systems and the overlap (intersection over union) of their 90 % area with their category's
in PARENT_DIR, the category map.

Branches: hydrogen bonds by donor and acceptor (N-H···O, O-H···O, ...); rings by kind (aromatic, aliphatic,
contact-closed, ...); pi interactions and other by contact class (C-H···pi, pi-stacking, H···H, cages, ...).

Usage: subfamily_regions.py [--group REGEX] [--min-systems N] PARENT_DIR OUT_DIR export1-density.csv [...]
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

from bake_reference import npy_name
from nci_reference_regions import BINS, CATEGORIES, RHO_EXTENT, S_MAX, load, mass_threshold, smooth

EDGE = 0.9


def branch(category, environment):
    """The CP's subfamily: the branch of its chemical-environment path under its reference category."""
    path = environment.split(" > ")
    if category == "Hydrogen bonds":
        # Contacts > Hydrogen bonds > N-H···O > ...
        return path[2] if len(path) > 2 else None
    # Rings > Aromatic rings > ...; Contacts > C-H···π > ...; Cage: Contact cage, 3 faces > ...
    if len(path) < 2:
        return None
    name = path[1].split(":")[0]
    return name if name != "Ring" else None


def anchors_of(density_path):
    stem = density_path[:-len("-density.csv")]
    stem = re.sub(r"-s[0-9.]+$", "", stem)
    path = stem + "-anchors.csv"
    return path if os.path.isfile(path) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", default=None)
    ap.add_argument("--min-systems", type=int, default=5)
    ap.add_argument("parent")
    ap.add_argument("out")
    ap.add_argument("exports", nargs="+")
    args = ap.parse_args()
    group = re.compile(args.group) if args.group else None
    os.makedirs(args.out, exist_ok=True)

    # {(category, branch): {group: [normalised fields]}}
    pools = defaultdict(lambda: defaultdict(list))
    for path in args.exports:
        anchors = anchors_of(path)
        if not anchors:
            print(f"{path}: no -anchors.csv companion, skipped")
            continue
        branches = defaultdict(set)
        with open(anchors, newline="") as f:
            for row in csv.DictReader(f):
                branches[row["reference_category"]].add(branch(row["reference_category"], row["chemical_environment"]))
        name = os.path.basename(path)
        match = group.search(name) if group else None
        key = match.group(1) if match else name
        for category, field in load(path).items():
            b = branches.get(category, set())
            if len(b) != 1 or None in b or field.sum() <= 0:
                continue
            pools[(category, next(iter(b)))][key].append(field / field.sum())

    parents = {}
    for category in CATEGORIES:
        p = os.path.join(args.parent, npy_name(category))
        if os.path.isfile(p):
            field = np.load(p)
            parents[category] = field >= mass_threshold(field, EDGE)

    subfamilies = []
    for (category, name), groups in sorted(pools.items(), key=lambda kv: (CATEGORIES.index(kv[0][0]), kv[0][1])):
        systems = len(groups)
        pooled = sum((sum(f) / len(f) for f in groups.values()), np.zeros((BINS, BINS)))
        field = smooth(pooled)
        mask = field >= mass_threshold(field, EDGE)
        parent = parents.get(category)
        iou = float((mask & parent).sum() / (mask | parent).sum()) if parent is not None else None
        kept = systems >= args.min_systems
        print(f"{'keep' if kept else 'skip'} {category:16} {name:24} {systems:4} systems, overlap with category "
              f"{'n/a' if iou is None else f'{100 * iou:.0f} %'}")
        if not kept:
            continue
        file = f"sub_{len(subfamilies)}.npy"
        np.save(os.path.join(args.out, file), field)
        subfamilies.append({"category": category, "name": name, "systems": systems, "file": file, "overlap": iou})
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump({"bins": BINS, "rho_extent": RHO_EXTENT, "s_max": S_MAX, "subfamilies": subfamilies}, f, indent=2,
                  ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
