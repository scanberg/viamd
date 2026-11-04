#!/usr/bin/env python3
"""Render publication figures from the pooled NCI reference maps.

The inputs are the ``regions`` and ``subfamilies`` directories produced by
the archived ``v3-extended/pool_extended.py`` pipeline.
No source CSVs are required for plotting.  The script writes:

* ``reference_categories.{png,svg}``: four category maps and their combined
  90%-edge overlay;
* ``reference_combined.{png,svg}``: the category 90%-edge and 50%-core lines;
* ``guides_<category>.{png,svg}``: every broad and fine map for that category,
  with common axes and the number of supporting pooled systems.

The densities are shown as relative colour within each map; line levels are
highest-density thresholds containing 90% (solid) and 50% (dashed) of the
map's pooled mass.  Figures use the map's actual rho and s extents.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap


CATEGORIES = ("Hydrogen bonds", "Rings", "π interactions", "Other")
COLORS = {
    "Hydrogen bonds": "#3473d1",
    "Rings": "#35ad54",
    "π interactions": "#d83e3e",
    "Other": "#e68a17",
}
LABEL_X = "sign(λ₂)ρ (a.u.)"
LABEL_Y = "Reduced density gradient, s"


def load_fields(regions: Path, subfamilies: Path):
    region_summary = json.loads((regions / "summary.json").read_text())
    family_summary = json.loads((subfamilies / "summary.json").read_text())
    fine_summary_path = subfamilies / "level2" / "summary.json"
    fine_summary = json.loads(fine_summary_path.read_text()) if fine_summary_path.exists() else {"subfamilies": []}
    bins = int(region_summary["bins"])
    rho_extent = float(region_summary["rho_extent"])
    s_max = float(region_summary["s_max"])
    if int(family_summary["bins"]) != bins:
        raise ValueError("category and subfamily map bin counts differ")
    if any(int(summary.get("bins", bins)) != bins for summary in (family_summary, fine_summary)):
        raise ValueError("subfamily map bin counts differ")
    for summary in (family_summary, fine_summary):
        if (float(summary.get("rho_extent", rho_extent)) != rho_extent
                or float(summary.get("s_max", s_max)) != s_max):
            raise ValueError("category and subfamily map axes differ")

    fields = {}
    for category in CATEGORIES:
        filename = category.replace(" ", "_") + ".npy"
        field = np.load(regions / filename)
        validate_field(filename, field, bins)
        fields[category] = field

    guides = {category: [] for category in CATEGORIES}
    for directory, entries in (
        (subfamilies, family_summary["subfamilies"]),
        (subfamilies / "level2", fine_summary["subfamilies"]),
    ):
        for entry in entries:
            category = entry["category"]
            if category not in guides:
                raise ValueError(f"unknown category in subfamily summary: {category}")
            field = np.load(directory / entry["file"])
            validate_field(entry["file"], field, bins)
            guides[category].append((entry["name"], int(entry["systems"]), field))
    for category in guides:
        guides[category].sort(key=lambda item: (item[0].count(":"), item[0].casefold()))
    return region_summary, fields, guides, rho_extent, s_max


def validate_field(name: str, field: np.ndarray, bins: int) -> None:
    if field.shape != (bins, bins) or not np.isfinite(field).all() or (field < 0).any():
        raise ValueError(f"{name}: expected a finite non-negative {bins}×{bins} map")
    if field.sum() <= 0:
        raise ValueError(f"{name}: map has no pooled mass")


def mass_threshold(field: np.ndarray, fraction: float) -> float:
    descending = np.sort(field.ravel())[::-1]
    return float(descending[min(len(descending) - 1,
                                np.searchsorted(np.cumsum(descending), fraction * descending.sum()))])


def grid_centres(bins: int, rho_extent: float, s_max: float):
    return (np.linspace(-rho_extent, rho_extent, bins + 1)[:-1]
            + rho_extent / bins,
            np.linspace(0.0, s_max, bins + 1)[:-1] + s_max / (2 * bins))


def draw_density(ax, field: np.ndarray, xc: np.ndarray, yc: np.ndarray, extent, color: str,
                 *, contours: bool = True):
    normalized = field / field.max()
    rgba = np.zeros((*field.shape, 4), dtype=float)
    rgb = matplotlib.colors.to_rgb(color)
    rgba[..., :3] = rgb
    rgba[..., 3] = np.where(normalized > 0.003, 0.78 * np.sqrt(normalized), 0.0)
    ax.imshow(rgba, origin="lower", extent=extent, aspect="auto", interpolation="nearest",
              rasterized=True, zorder=1)
    if contours:
        ax.contour(xc, yc, field, levels=[mass_threshold(field, 0.9)], colors=[color],
                   linewidths=1.05, zorder=3)
        ax.contour(xc, yc, field, levels=[mass_threshold(field, 0.5)], colors=[color],
                   linewidths=0.8, linestyles="dashed", zorder=4)


def style_axes(ax, rho_extent: float, s_max: float):
    ax.set_xlim(-rho_extent, rho_extent)
    ax.set_ylim(0, s_max)
    ax.set_xticks(np.linspace(-rho_extent, rho_extent, 5))
    ax.axvline(0, color="#666666", lw=0.45, alpha=0.6, zorder=0)
    ax.tick_params(labelsize=7, direction="out", length=2)
    ax.set_facecolor("white")


def save_pair(fig, out: Path, stem: str, dpi: int = 300):
    fig.savefig(out / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_category_summary(out: Path, summary, fields, xc, yc, rho_extent, s_max):
    fig, axes = plt.subplots(1, 5, figsize=(16.2, 3.65), sharex=True, sharey=True)
    extent = [-rho_extent, rho_extent, 0, s_max]
    for ax, category in zip(axes[:4], CATEGORIES):
        draw_density(ax, fields[category], xc, yc, extent, COLORS[category])
        count = summary["categories"][category]["systems"]
        ax.set_title(f"{category}\n{count} pooled systems", fontsize=9, pad=5)
        style_axes(ax, rho_extent, s_max)
    combined = axes[4]
    for category in CATEGORIES:
        field = fields[category]
        color = COLORS[category]
        combined.contour(xc, yc, field, levels=[mass_threshold(field, 0.9)],
                         colors=[color], linewidths=1.15, zorder=3)
        combined.contour(xc, yc, field, levels=[mass_threshold(field, 0.5)],
                         colors=[color], linewidths=0.85, linestyles="dashed", zorder=4)
    combined.set_title("Combined 90% edges / 50% cores", fontsize=9, pad=5)
    style_axes(combined, rho_extent, s_max)
    handles = [
        plt.Line2D([0], [0], color=COLORS[c], lw=1.3,
                   label=f"{c} (n={summary['categories'][c]['systems']})")
        for c in CATEGORIES
    ]
    combined.legend(handles=handles, loc="upper right", fontsize=6.5, frameon=True,
                    framealpha=0.92, borderpad=0.35, handlelength=1.2)
    fig.supxlabel(LABEL_X, fontsize=9)
    fig.supylabel(LABEL_Y, fontsize=9, x=0.015)
    fig.suptitle("Pooled NCI reference densities and highest-density mass contours",
                 fontsize=11, y=1.02)
    fig.tight_layout(pad=0.7, rect=(0.045, 0.02, 1.0, 0.97))
    save_pair(fig, out, "reference_categories")


def plot_combined(out: Path, summary, fields, xc, yc, rho_extent, s_max):
    fig, ax = plt.subplots(figsize=(8.1, 6.1))
    for category in CATEGORIES:
        field = fields[category]
        color = COLORS[category]
        ax.contour(xc, yc, field, levels=[mass_threshold(field, 0.9)],
                   colors=[color], linewidths=1.7)
        ax.contour(xc, yc, field, levels=[mass_threshold(field, 0.5)],
                   colors=[color], linewidths=1.1, linestyles="dashed")
    style_axes(ax, rho_extent, s_max)
    handles = [
        plt.Line2D([0], [0], color=COLORS[c], lw=1.7,
                   label=f"{c}: 90% edge / 50% core; n={summary['categories'][c]['systems']}")
        for c in CATEGORIES
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True, framealpha=0.94)
    ax.set_xlabel(LABEL_X)
    ax.set_ylabel(LABEL_Y)
    ax.set_title("Combined category map (color identifies category)", fontsize=11)
    fig.tight_layout()
    save_pair(fig, out, "reference_combined")


def plot_guide_sheet(out: Path, category: str, guides, xc, yc, rho_extent, s_max,
                     columns: int = 5):
    entries = guides[category]
    columns = min(columns, max(1, len(entries)))
    rows = math.ceil(len(entries) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(columns * 2.75, rows * 2.33),
                             sharex=True, sharey=True, squeeze=False)
    extent = [-rho_extent, rho_extent, 0, s_max]
    color = COLORS[category]
    for ax in axes.ravel()[len(entries):]:
        ax.axis("off")
    for ax, (name, systems, field) in zip(axes.ravel(), entries):
        draw_density(ax, field, xc, yc, extent, color)
        wrapped = "\n".join(textwrap.wrap(name, width=34, break_long_words=False,
                                          break_on_hyphens=False))
        ax.set_title(f"{wrapped}\nn = {systems} systems", fontsize=7.1, pad=3)
        style_axes(ax, rho_extent, s_max)
    for ax in axes[-1]:
        if ax.axison:
            ax.set_xlabel(LABEL_X, fontsize=8)
    for ax in axes[:, 0]:
        if ax.axison:
            ax.set_ylabel(LABEL_Y, fontsize=8)
    fig.suptitle(f"{category}: all {len(entries)} supported broad and fine guides",
                 fontsize=11, y=0.998)
    fig.tight_layout(rect=(0.015, 0.015, 0.995, 0.985), h_pad=0.85, w_pad=0.5)
    save_pair(fig, out, "guides_" + category.lower().replace(" ", "_").replace("π", "pi"),
              dpi=300)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", required=True, type=Path,
                        help="pooled category-map directory")
    parser.add_argument("--subfamilies", required=True, type=Path,
                        help="pooled broad/fine guide directory")
    parser.add_argument("--out", required=True, type=Path, help="figure output directory")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary, fields, guides, rho_extent, s_max = load_fields(args.regions, args.subfamilies)
    bins = int(summary["bins"])
    xc, yc = grid_centres(bins, rho_extent, s_max)
    plot_category_summary(args.out, summary, fields, xc, yc, rho_extent, s_max)
    plot_combined(args.out, summary, fields, xc, yc, rho_extent, s_max)
    for category in CATEGORIES:
        plot_guide_sheet(args.out, category, guides, xc, yc, rho_extent, s_max)
    n_guides = sum(len(entries) for entries in guides.values())
    if n_guides != 106:
        raise ValueError(f"expected all 106 supported guides (18 broad + 88 fine); found {n_guides}")
    print(json.dumps({
        "category_systems": {c: int(summary["categories"][c]["systems"]) for c in CATEGORIES},
        "guide_counts": {c: len(guides[c]) for c in CATEGORIES},
        "maps": n_guides,
        "grid": {"bins": bins, "signed_rho_extent_au": rho_extent, "s_max": s_max},
        "figures": sorted(path.name for path in args.out.glob("*")),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
