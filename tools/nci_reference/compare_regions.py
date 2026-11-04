#!/usr/bin/env python3
"""Side by side: pooled expected regions of several runs (the .npy fields nci_reference_regions.py saves).

Usage: compare_regions.py OUT.png LABEL=REGIONS_DIR [LABEL=REGIONS_DIR ...]
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CATEGORIES = ["Hydrogen bonds", "Rings", "\u03c0 interactions", "Other"]


def mass_threshold(field, fraction):
    v = np.sort(field.ravel())[::-1]
    c = np.cumsum(v)
    return v[min(len(v) - 1, np.searchsorted(c, fraction * c[-1]))] if c[-1] > 0 else np.inf


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    runs = [a.rsplit("=", 1) for a in sys.argv[2:]]  # labels may contain "=", paths do not
    fig, axes = plt.subplots(len(runs), len(CATEGORIES), figsize=(4.2 * len(CATEGORIES), 3.6 * len(runs)),
                             sharex=True, sharey=True, squeeze=False)
    for r, (label, directory) in enumerate(runs):
        with open(os.path.join(directory, "summary.json")) as f:
            meta = json.load(f)
        ext, smax, bins = meta["rho_extent"], meta["s_max"], meta["bins"]
        xs = np.linspace(-ext, ext, bins + 1)
        xc = 0.5 * (xs[1:] + xs[:-1])
        ys = np.linspace(0, smax, bins + 1)
        yc = 0.5 * (ys[1:] + ys[:-1])
        for c, category in enumerate(CATEGORIES):
            ax = axes[r][c]
            path = os.path.join(directory, f"{category.replace(' ', '_')}.npy")
            info = meta["categories"].get(category, {})
            ax.set_title(f"{label}: {category} ({info.get('systems', 0)})", fontsize=9)
            if not os.path.isfile(path):
                continue
            field = np.load(path)
            d = field / field.max()
            ax.imshow(np.where(d > 0.01, d, np.nan), origin="lower", extent=[-ext, ext, 0, smax], aspect="auto",
                      cmap="viridis")
            levels = sorted({mass_threshold(field, m) for m in (0.5, 0.9)})
            ax.contour(xc, yc, field, levels=levels, colors="white", linewidths=[0.8, 1.6][:len(levels)])
            ax.axvline(0, color="grey", lw=0.5)
    for ax in axes[-1]:
        ax.set_xlabel("sign(\u03bb\u2082)\u03c1 (au)")
    for row in axes:
        row[0].set_ylabel("s")
    fig.suptitle("Pooled NCI category density; white: 90 % volume edge, 50 % core")
    fig.tight_layout()
    fig.savefig(sys.argv[1], dpi=130)
    return 0


if __name__ == "__main__":
    sys.exit(main())
