#!/usr/bin/env python3
"""Render the pooled reference maps for the documentation.

From nci_reference_regions.py's REGIONS_DIR (and optionally subfamily_regions.py's SUBFAMILY_DIR) writes to OUT_DIR:

  reference_categories.png      one panel per reference category (pooled density, 50 % core and 90 % edge) and a
                                combined panel: every category in the app's colour, its 90 % edge outlined
  reference_interactions.png    one panel per kept interaction class (the level-1 subfamilies: CP-tree branches),
                                its 50 % / 90 % levels, with its category's 90 % edge dashed

The s axis stops at the extraction cutoff (s <= 1.0): the maps hold no volume above it.

Usage: plot_regions.py REGIONS_DIR OUT_DIR [--subfamilies SUBFAMILY_DIR] [--title TEXT]
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nci_reference_regions import CATEGORIES, mass_threshold

S_CUT = 1.0
# The app's interaction colours (nci_surface_model.cpp)
COLOURS = {"Hydrogen bonds": (0.20, 0.45, 0.90), "\u03c0 interactions": (0.88, 0.20, 0.20),
           "Rings": (0.20, 0.68, 0.28), "Other": (0.96, 0.55, 0.10)}
XLABEL = "sign(\u03bb\u2082)\u03c1 (au)"


def axes_grid(fields, rho_extent, s_max):
    bins = fields.shape[-1]
    xs = np.linspace(-rho_extent, rho_extent, bins + 1)
    ys = np.linspace(0, s_max, bins + 1)
    return 0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1])


def frame(ax, rho_extent, title):
    ax.set_title(title, fontsize=9)
    ax.set_xlim(-rho_extent, rho_extent)
    ax.set_ylim(0, S_CUT)
    ax.axvline(0, color="0.6", lw=0.5)
    ax.tick_params(labelsize=7)


def density(ax, field, xc, yc, rho_extent, s_max, colour):
    rgb = np.array(colour)
    d = field / field.max()
    img = np.ones(field.shape + (4,))
    img[..., :3] = 1 - (1 - rgb) * np.sqrt(d)[..., None]
    img[..., 3] = np.where(d > 0.005, 1.0, 0.0)
    ax.imshow(img, origin="lower", extent=[-rho_extent, rho_extent, 0, s_max], aspect="auto", interpolation="bilinear")
    levels = sorted({mass_threshold(field, 0.9), mass_threshold(field, 0.5)})
    ax.contour(xc, yc, field, levels=levels, colors=[tuple(0.6 * rgb)], linewidths=[0.8, 1.5][:len(levels)])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("regions")
    p.add_argument("out")
    p.add_argument("--subfamilies")
    p.add_argument("--title", default="")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    summary = json.load(open(os.path.join(a.regions, "summary.json")))
    rho_extent, s_max = summary["rho_extent"], summary["s_max"]
    fields = {}
    for c in CATEGORIES:
        path = os.path.join(a.regions, f"{c.replace(' ', '_')}.npy")
        if os.path.isfile(path):
            fields[c] = np.load(path)
    xc, yc = axes_grid(next(iter(fields.values())), rho_extent, s_max)

    fig, axes = plt.subplots(1, len(fields) + 1, figsize=(3.4 * (len(fields) + 1), 3.4), sharey=True)
    for ax, (c, field) in zip(axes, fields.items()):
        frame(ax, rho_extent, f"{c} ({summary['categories'][c]['systems']} systems)")
        density(ax, field, xc, yc, rho_extent, s_max, COLOURS[c])
        ax.set_xlabel(XLABEL, fontsize=8)
    ax = axes[-1]
    frame(ax, rho_extent, "All interactions (90 % edges)")
    total = sum(f / f.sum() for f in fields.values())
    for c, field in fields.items():
        share = (field / field.sum()) / np.maximum(total, 1e-30)
        rgb = np.array(COLOURS[c])
        img = np.zeros(field.shape + (4,))
        img[..., :3] = rgb
        img[..., 3] = 0.55 * share * (field >= mass_threshold(field, 0.9))
        ax.imshow(img, origin="lower", extent=[-rho_extent, rho_extent, 0, s_max], aspect="auto")
        ax.contour(xc, yc, field, levels=[mass_threshold(field, 0.9)], colors=[tuple(rgb)], linewidths=1.2)
    ax.set_xlabel(XLABEL, fontsize=8)
    axes[0].set_ylabel("s (cut at 1.0)", fontsize=8)
    if a.title:
        fig.suptitle(a.title, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "reference_categories.png"), dpi=130)
    plt.close(fig)

    if a.subfamilies:
        sub = json.load(open(os.path.join(a.subfamilies, "summary.json")))["subfamilies"]
        cols = 6
        rows = math.ceil(len(sub) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.5 * rows), sharex=True, sharey=True, squeeze=False)
        for ax in axes.ravel()[len(sub):]:
            ax.axis("off")
        for ax, entry in zip(axes.ravel(), sub):
            field = np.load(os.path.join(a.subfamilies, entry["file"]))
            c = entry["category"]
            frame(ax, rho_extent, f"{entry['name']} ({entry['systems']})")
            density(ax, field, xc, yc, rho_extent, s_max, COLOURS[c])
            if c in fields:
                ax.contour(xc, yc, fields[c], levels=[mass_threshold(fields[c], 0.9)], colors="0.35",
                           linewidths=0.7, linestyles="dashed")
        for ax in axes[-1]:
            ax.set_xlabel(XLABEL, fontsize=8)
        for ax in axes[:, 0]:
            ax.set_ylabel("s", fontsize=8)
        fig.suptitle("Per interaction class (number of systems); dashed: its category's 90 % edge"
                     + (f"\n{a.title}" if a.title else ""), fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "reference_interactions.png"), dpi=110)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
