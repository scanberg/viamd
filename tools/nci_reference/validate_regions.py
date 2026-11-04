#!/usr/bin/env python3
"""Share of test systems' NCI volume inside a reference's expected regions, per category.

The reference is a regions directory written by nci_reference_regions.py; the test set is a list of -density.csv
exports (or their -anchors companions). Each system's category volume is placed on the same grid, and the share
inside the reference's 90 % edge (the region the app draws and counts) is reported per category: the mean over
systems, and the pooled share with every system weighed the same.

Usage: validate_regions.py REGIONS_DIR export1-density.csv [export2-density.csv ...]
"""
import os
import sys

import numpy as np

from bake_reference import npy_name
from nci_reference_regions import CATEGORIES, density_file, load, mass_threshold

EDGE = 0.9


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    regions = sys.argv[1]
    masks = {}
    for category in CATEGORIES:
        path = os.path.join(regions, npy_name(category))
        if os.path.isfile(path):
            field = np.load(path)
            masks[category] = field >= mass_threshold(field, EDGE)
    shares = {c: [] for c in masks}
    for path in sys.argv[2:]:
        found = density_file(path)
        if not found:
            continue
        for category, field in load(found).items():
            total = field.sum()
            if category in masks and total > 0:
                shares[category].append(float((field * masks[category]).sum() / total))
    print(f"{'category':<18} {'systems':>7} {'mean in':>8} {'min in':>8}")
    for category, values in shares.items():
        if values:
            print(f"{category:<18} {len(values):>7} {100 * np.mean(values):>7.1f}% {100 * np.min(values):>7.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
