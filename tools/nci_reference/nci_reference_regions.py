#!/usr/bin/env python3
"""Pool VIAMD NCI -density.csv files (GUI exports or viamd_nci_reference output) into per-category expected regions (Ramachandran-like).

Each export's category reference rows (sample_scope=category, extraction_settings=category_reference_basin: each
category flooded within its CPs' sign(lambda2) basin; older category_reference_defaults exports are skipped)
are rebinned onto a common (sign(lambda2)rho, s) grid, normalised per system so each system weighs the same,
summed, smoothed lightly as the app does, and contoured at fractions of the pooled mass (the highest-density
region holding 90 % of the volume: the edge, as in the app) with its 50 % core.

With --group REGEX, exports whose file names give the same first capture of REGEX (e.g. one dimer at several
separations) share one system's weight; names that do not match are systems of their own.

Usage: nci_reference_regions.py [--group REGEX] OUT_DIR export1-density.csv [export2-density.csv ...]
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S_MAX = 2.0
BINS = 96
RHO_EXTENT = 0.07   # common signed-rho half-range, au: the largest calibrated cutoff
CATEGORIES = ["Hydrogen bonds", "Rings", "\u03c0 interactions", "Other"]
MASS_LEVELS = [0.5, 0.9]


DENSITY_COLUMN = "sign(lambda2)rho_min (au)"


def density_file(path):
    """The density export for path: path itself, or the -density companion of a samples or anchors export."""
    with open(path, newline="") as f:
        header = next(csv.reader(f), [])
    if DENSITY_COLUMN in header:
        return path
    stem = path[:-len(".csv")] if path.endswith(".csv") else path
    if stem.endswith("-anchors"):
        stem = stem[:-len("-anchors")]
    companion = stem + "-density.csv"
    if os.path.isfile(companion):
        with open(companion, newline="") as f:
            if DENSITY_COLUMN in next(csv.reader(f), []):
                return companion
    return None


def load(path):
    """{category: 2D volume on the common grid} for one export, the category reference rows only."""
    out = defaultdict(lambda: np.zeros((BINS, BINS)))
    stale = False
    xw = 2 * RHO_EXTENT / BINS
    yw = S_MAX / BINS
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["sample_scope"] != "category":
                continue
            if row["extraction_settings"] == "category_reference_defaults":
                stale = True
                continue
            if row["extraction_settings"] != "category_reference_basin":
                continue
            lo = float(row["sign(lambda2)rho_min (au)"])
            hi = float(row["sign(lambda2)rho_max (au)"])
            y = min(BINS - 1, int(float(row["s_min (dimensionless)"]) / yw + 1e-6))
            v = float(row["volume (bohr^3)"])
            # split by overlap in signed rho
            a, b = (lo + RHO_EXTENT) / xw, (hi + RHO_EXTENT) / xw
            for x in range(max(0, int(np.floor(a))), min(BINS, int(np.ceil(b)))):
                overlap = (min(b, x + 1) - max(a, x)) / (b - a)
                if overlap > 0:
                    out[row["category"]][y, x] += overlap * v
    if stale and not out:
        print(f"{path}: exported before the sign(lambda2) basin references; re-export it with the current build")
    return out


def smooth(field):
    k = np.array([1.0, 2.0, 1.0])
    k = np.outer(k, k) / 16.0
    p = np.pad(field, 1)
    return sum(k[j, i] * p[j:j + field.shape[0], i:i + field.shape[1]] for j in range(3) for i in range(3))


def mass_threshold(field, fraction):
    v = np.sort(field.ravel())[::-1]
    c = np.cumsum(v)
    if c[-1] <= 0:
        return np.inf
    return v[min(len(v) - 1, np.searchsorted(c, fraction * c[-1]))]


def main():
    args = sys.argv[1:]
    group = None
    if len(args) >= 2 and args[0] == "--group":
        group = re.compile(args[1])
        args = args[2:]
    if len(args) < 2:
        print(__doc__)
        return 1
    out_dir = args[0]
    missing = [p for p in args[1:] if not os.path.isfile(p)]
    if missing:
        print("Not found: " + ", ".join(missing))
        print("Export each system in VIAMD (Export NCI samples...) or run viamd_nci_reference / run_batch.py; both write <name>-density.csv.")
        return 1
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for path in args[1:]:
        found = density_file(path)
        if not found:
            print(f"{path}: not a density export, and no -density.csv companion was found")
            return 1
        if found != path:
            print(f"{path}: using its density companion {found}")
        if found not in paths:
            paths.append(found)
    # {category: {group: [normalised fields]}}; each group weighs one system, split evenly among its exports
    groups = {c: defaultdict(list) for c in CATEGORIES}
    for path in paths:
        name = os.path.basename(path)
        match = group.search(name) if group else None
        key = match.group(1) if match else name
        for category, field in load(path).items():
            if category not in groups or field.sum() <= 0:
                continue
            groups[category][key].append(field / field.sum())
    pooled = {c: sum((sum(f) / len(f) for f in groups[c].values()), np.zeros((BINS, BINS))) for c in CATEGORIES}
    systems = {c: len(groups[c]) for c in CATEGORIES}
    if not any(systems.values()):
        print("No category reference rows found: export with CP categories shown in the layer.")
        return 1

    extent = [-RHO_EXTENT, RHO_EXTENT, 0, S_MAX]
    xs = np.linspace(-RHO_EXTENT, RHO_EXTENT, BINS + 1)
    xc, yc = 0.5 * (xs[1:] + xs[:-1]), np.linspace(0, S_MAX, BINS + 1)
    yc = 0.5 * (yc[1:] + yc[:-1])
    summary = {}
    fig, axes = plt.subplots(1, len(CATEGORIES), figsize=(5 * len(CATEGORIES), 4.5), sharey=True)
    for ax, category in zip(axes, CATEGORIES):
        field = smooth(pooled[category])
        ax.set_title(f"{category} ({systems[category]} systems)")
        ax.set_xlabel("sign(\u03bb\u2082)\u03c1 (au)")
        ax.set_xlim(-RHO_EXTENT, RHO_EXTENT)
        ax.set_ylim(0, S_MAX)
        if field.sum() <= 0:
            continue
        density = field / field.max()
        ax.imshow(np.where(density > 0.01, density, np.nan), origin="lower", extent=extent, aspect="auto",
                  cmap="viridis")
        levels = sorted({mass_threshold(field, m) for m in MASS_LEVELS})
        ax.contour(xc, yc, field, levels=levels, colors="white", linewidths=[0.8, 1.6][:len(levels)])
        np.save(os.path.join(out_dir, f"{category.replace(' ', '_')}.npy"), field)
        summary[category] = {
            "systems": systems[category],
            "peak_signed_rho": float(xc[np.unravel_index(field.argmax(), field.shape)[1]]),
            "peak_s": float(yc[np.unravel_index(field.argmax(), field.shape)[0]]),
            "positive_share_inside_edge": float(field[:, xc >= 0][field[:, xc >= 0] >= mass_threshold(field, 0.9)].sum() / field[field >= mass_threshold(field, 0.9)].sum()),
        }
    axes[0].set_ylabel("s")
    fig.suptitle("Pooled continuous NCI density; white: 90 % volume edge, 50 % core")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "expected_regions.png"), dpi=150)
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump({"bins": BINS, "rho_extent": RHO_EXTENT, "s_max": S_MAX, "categories": summary}, f, indent=2)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
