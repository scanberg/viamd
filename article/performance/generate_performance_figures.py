#!/usr/bin/env python3
# pyright: reportMissingModuleSource=false
"""Generate performance figures for the VIAMD--VeloxChem article.

Inputs:
    article/performance/performance-timings-clean.csv
Outputs:
    article/figure/performance-timings.png
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = ROOT / "article" / "performance" / "performance-timings-clean.csv"
FIG_DIR = ROOT / "article" / "figure"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TIMINGS_FIG = FIG_DIR / "performance-timings.png"


def read_rows() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="") as fh:
        return list(csv.DictReader(fh))


def generate_timing_figure(rows: list[dict[str, str]]) -> None:
    # Keep operations ordered by article relevance.
    operation_order = [
        "HDF5 parse",
        "Basis extraction",
        "MO volume",
        "Transition density",
    ]
    datasets = []
    for row in rows:
        if row["dataset"] not in datasets:
            datasets.append(row["dataset"])

    data = {(row["dataset"], row["operation"]): float(row["median_ms"]) for row in rows}

    colors = {
        "HDF5 parse": "#4C78A8",
        "Basis extraction": "#72B7B2",
        "MO volume": "#F58518",
        "Transition density": "#E45756",
    }

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    x_base = list(range(len(datasets)))
    width = 0.18
    offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

    for op, off in zip(operation_order, offsets):
        values = [data.get((dataset, op), 0.0) for dataset in datasets]
        labels = ["" if v == 0.0 else f"{v:.1f}" for v in values]
        bars = ax.bar([x + off for x in x_base], values, width=width, label=op, color=colors[op])
        for bar, label in zip(bars, labels):
            if label:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    rotation=90 if bar.get_height() > 12 else 0,
                )

    ax.set_xticks(x_base)
    ax.set_xticklabels(datasets)
    ax.set_ylabel("Median wall time (ms)")
    ax.set_title("First-draft CPU benchmark of VeloxChem data operations in VIAMD")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.legend(ncol=2, frameon=False)
    ax.set_ylim(0, max(float(row["median_ms"]) for row in rows) * 1.25)

    # Add compact metadata as an in-figure note.
    meta = []
    seen = set()
    for row in rows:
        dataset = row["dataset"]
        if dataset in seen:
            continue
        seen.add(dataset)
        meta.append(f"{dataset}: {row['atoms']} atoms, {row['aos']} AOs, {row['states']} states")
    note = "64³ grid, 9 repeats, Release CPU build\n" + "\n".join(meta)
    ax.text(
        0.01,
        0.98,
        note,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#CCCCCC", "alpha": 0.9},
    )

    fig.tight_layout()
    fig.savefig(TIMINGS_FIG, dpi=300)
    plt.close(fig)


def main() -> None:
    rows = read_rows()
    generate_timing_figure(rows)
    print(TIMINGS_FIG)


if __name__ == "__main__":
    main()
