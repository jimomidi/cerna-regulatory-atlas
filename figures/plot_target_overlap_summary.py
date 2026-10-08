"""
Bar charts summarizing what fraction of a candidate target-gene set meets
combined correlation / FDR / fold-change criteria.

These functions plot *precomputed* percentages (typically produced by
``analysis/target_gene_overlap.py``) rather than recomputing them, so they
can be reused for any gene set / threshold combination.

Usage
-----
    python -m figures.plot_target_overlap_summary --demo results/demo_overlap.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np

from config import FIGURE_DPI
from utils import apply_default_style


def plot_criteria_percentages(
    percentages: Mapping[str, float],
    output_path: Path,
    title: str = "",
) -> None:
    """Simple bar chart of percentages, one bar per criterion label."""
    apply_default_style(font_size=18)
    labels = list(percentages.keys())
    values = list(percentages.values())

    fig, ax = plt.subplots(figsize=(7, 2.8))
    bars = ax.bar(range(len(labels)), values, width=0.5, color="black", edgecolor="black")
    if title:
        fig.suptitle(title, fontsize=18, weight="bold", y=0.98)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 1.5, f"{height:.2f}%", ha="center", va="bottom", fontsize=18)
    ax.set_ylim(0, max(values) + 10)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)


def plot_logfc_threshold_summary(
    labels: Sequence[str],
    values: Sequence[float],
    group_legend: Mapping[str, str],
    output_path: Path,
    xlabel: str = "logFC",
    colors: Sequence[str] | None = None,
) -> None:
    """Grouped bar chart of percentages passing a series of logFC thresholds,
    under two correlation/FDR regimes (or any two conditions supplied via
    ``group_legend``, e.g. {"r < -0.3 & FDR < 0.05": color1, ...}).
    """
    from matplotlib.patches import Patch

    apply_default_style(font_size=18)
    colors = colors or ["#006DAB", "#A30055", "#1B5E20", "#B3DDF2", "#F4B6D2", "#A5D6A7"]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8.5, 3.5))
    bars = ax.bar(x, values, color=colors[: len(labels)], width=0.6, edgecolor="black")

    ax.set_title(xlabel, pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel(xlabel, labelpad=10)
    ax.set_ylim(0, max(values) + 10)
    ax.set_yticks(np.arange(0, max(values) + 11, 10))

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 0.8, f"{height:.2f}%", ha="center", va="bottom", fontsize=14)

    legend_elements = [Patch(facecolor=color, edgecolor="black", label=label) for label, color in group_legend.items()]
    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=14, labelspacing=0.4)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)


def _demo(output_path: Path) -> None:
    """Render the summary chart with arbitrary, clearly synthetic values so
    the plotting code can be verified without a full analysis run. These
    numbers are not derived from any analysis or input file -- do not
    interpret them as a result; use ``plot_criteria_percentages`` with your
    own computed percentages instead."""
    plot_criteria_percentages(
        {
            "r < -0.3\nFDR < 0.05": 25.0,
            "r < -0.2\nFDR < 0.05": 40.0,
            "r < -0.3": 15.0,
            "r < -0.2": 30.0,
        },
        output_path,
        title="Example (synthetic data): % of Target Gene Set Meeting Criteria",
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", type=Path, help="Render a demo figure with placeholder values to this path.")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.demo:
        _demo(args.demo)
        print(f"Saved demo figure to {args.demo}")
    else:
        print("Import this module and call plot_criteria_percentages / "
              "plot_logfc_threshold_summary with your own computed percentages, "
              "or pass --demo <path> to render an illustrative example.")


if __name__ == "__main__":
    main()
