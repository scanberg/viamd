#!/usr/bin/env python3
"""Pool reference exports per CP subfamily, for the reference's sub-maps.

Reads each export's -cp-density.csv (the density per CP, from viamd_nci_reference) and -anchors.csv companions, so
every CP counts in its own subfamily, whatever the other CPs of its system. Two levels:

  1. the CP-tree branch under the reference category: hydrogen bonds by donor and acceptor (N-H···O, O-H···O, ...);
     rings by kind (aromatic, aliphatic, contact-closed, ...); pi interactions and other by contact class (C-H···pi,
     pi-stacking, sigma-hole bonds, H···H, cages, ...). Written to OUT_DIR (baked with bake_reference.py).
  2. a finer split of each branch: hydrogen bonds and contacts by intra/intermolecular and element pair, rings by
     size and composition ("Aromatic 6-ring C4N2", "Closed by Hydrogen bond, 7-ring"). Written to OUT_DIR/level2.

Each subfamily is pooled as nci_reference_regions.py pools a category: every system weighing the same (--group as
there; a system's CPs of the subfamily summed), smoothed, saved as sub_<n>.npy. A subfamily is kept with at least
--min-systems systems; summary.json lists the kept ones (with the parent branch for level 2), their systems and the
overlap (intersection over union) of their 90 % area with their category's in PARENT_DIR, the category map.

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


def finer(category, row, level1):
    """The CP's level-2 subfamily within its branch, or None."""
    kind, group = row.get("kind", ""), row.get("group", "")
    if kind == "Ring":
        if group.startswith("Closed by"):
            # "Closed by Hydrogen bond, 7-ring CHNO" without the composition
            return " ".join(row.get("finer_1", "").split(" ")[:-1]) or group
        return group or None
    if kind == "Cage":
        return group.split(",")[0] or None
    ea, eb = row.get("element_a", ""), row.get("element_b", "")
    if not ea or not eb:
        return None
    ma, mb = row.get("molecule_a", ""), row.get("molecule_b", "")
    where = "" if not ma or not mb or ma == "-1" else (", intermolecular" if ma != mb else ", intramolecular")
    if category == "Hydrogen bonds":
        return level1 + where
    return f"{level1}: {'···'.join(sorted((ea, eb)))}{where}"


def companions(density_path):
    stem = density_path[:-len("-density.csv")]
    base = re.sub(r"-s[0-9.]+$", "", stem)
    anchors, cp = base + "-anchors.csv", stem + "-cp-density.csv"
    return (anchors if os.path.isfile(anchors) else None), (cp if os.path.isfile(cp) else None)


def rho_extent_of(density_path):
    """The signed-rho half-range of an export's bins (0 if none)."""
    with open(density_path, newline="") as f:
        for row in csv.DictReader(f):
            lo, hi, x = float(row["sign(lambda2)rho_min (au)"]), float(row["sign(lambda2)rho_max (au)"]), int(row["bin_x"])
            return round(-lo + x * (hi - lo), 6)
    return 0.0


def cp_fields(cp_path, bins_in):
    """{cp_index: 2D volume on the common grid}; the export's bins are an integer multiple of the common grid's."""
    factor = bins_in // BINS
    out = defaultdict(lambda: np.zeros((BINS, BINS)))
    with open(cp_path, newline="") as f:
        for row in csv.DictReader(f):
            out[int(row["cp_index"])][int(row["bin_y"]) // factor, int(row["bin_x"]) // factor] += float(row["volume (bohr^3)"])
    return out


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

    # {level: {(category, parent, name): {group: [normalised fields]}}}
    pools = {1: defaultdict(lambda: defaultdict(list)), 2: defaultdict(lambda: defaultdict(list))}
    for path in args.exports:
        anchors, cp = companions(path)
        if not anchors or not cp:
            print(f"{path}: no -anchors.csv or -cp-density.csv companion, skipped")
            continue
        extent = rho_extent_of(path)
        if abs(extent - RHO_EXTENT) > 1e-6:
            print(f"{path}: signed-rho extent {extent} differs from the map's {RHO_EXTENT}, skipped")
            continue
        keys = {}
        with open(anchors, newline="") as f:
            for row in csv.DictReader(f):
                category = row["reference_category"]
                if category not in CATEGORIES:
                    continue
                b1 = branch(category, row["chemical_environment"])
                if b1 is None:
                    continue
                b2 = finer(category, row, b1)
                keys[int(row["cp_index"])] = (category, b1, b2)
        name = os.path.basename(path)
        match = group.search(name) if group else None
        key = match.group(1) if match else name
        system = {1: defaultdict(lambda: np.zeros((BINS, BINS))), 2: defaultdict(lambda: np.zeros((BINS, BINS)))}
        for index, field in cp_fields(cp, 2 * BINS).items():
            if index not in keys:
                continue
            category, b1, b2 = keys[index]
            system[1][(category, None, b1)] += field
            if b2:
                system[2][(category, b1, b2)] += field
        for level in (1, 2):
            for sub, field in system[level].items():
                if field.sum() > 0:
                    pools[level][sub][key].append(field / field.sum())

    parents = {}
    for category in CATEGORIES:
        p = os.path.join(args.parent, npy_name(category))
        if os.path.isfile(p):
            field = np.load(p)
            parents[category] = field >= mass_threshold(field, EDGE)

    for level in (1, 2):
        out = args.out if level == 1 else os.path.join(args.out, "level2")
        os.makedirs(out, exist_ok=True)
        subfamilies = []
        order = lambda kv: (CATEGORIES.index(kv[0][0]), kv[0][1] or "", kv[0][2])
        for (category, parent, name), groups in sorted(pools[level].items(), key=order):
            systems = len(groups)
            pooled = sum((sum(f) / len(f) for f in groups.values()), np.zeros((BINS, BINS)))
            field = smooth(pooled)
            mask = field >= mass_threshold(field, EDGE)
            parent_mask = parents.get(category)
            iou = float((mask & parent_mask).sum() / (mask | parent_mask).sum()) if parent_mask is not None else None
            kept = systems >= args.min_systems
            print(f"{'keep' if kept else 'skip'} L{level} {category:16} {name:48} {systems:4} systems, overlap with category "
                  f"{'n/a' if iou is None else f'{100 * iou:.0f} %'}")
            if not kept:
                continue
            file = f"sub_{len(subfamilies)}.npy"
            np.save(os.path.join(out, file), field)
            entry = {"category": category, "name": name, "systems": systems, "file": file, "overlap": iou}
            if parent:
                entry["parent"] = parent
            subfamilies.append(entry)
        with open(os.path.join(out, "summary.json"), "w") as f:
            json.dump({"bins": BINS, "rho_extent": RHO_EXTENT, "s_max": S_MAX, "level": level, "subfamilies": subfamilies},
                      f, indent=2, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
