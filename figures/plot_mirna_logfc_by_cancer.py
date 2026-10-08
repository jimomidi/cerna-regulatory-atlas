"""
Bar plot of a target miRNA's log2 fold-change across multiple cancer types.

Expected input
--------------
csv (or DataFrame) with columns ``miRNA``, ``cancer_type`` and ``log2FC``.

Usage
-----
    python -m figures.plot_mirna_logfc_by_cancer \
        --input data/mirna_logfc_by_cancer.csv \
        --mirna-id target_mirna \
        --output results/mirna_logfc_by_cancer.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from config import FIGURE_DPI, TARGET_MIRNA_ID
from utils import apply_default_style

UP_COLOR = "#67001f"
DOWN_COLOR = "#053061"


def plot_mirna_logfc_by_cancer(df: pd.DataFrame, mirna_id: str, output_path: Path) -> None:
    subset = df[df["miRNA"] == mirna_id].copy()
    if subset.empty:
        raise ValueError(f"No rows found for miRNA '{mirna_id}'.")
    subset["color"] = subset["log2FC"].apply(lambda x: UP_COLOR if x > 0 else DOWN_COLOR)

    apply_default_style(font_size=18)
    x_pos = [i * 0.6 for i in range(len(subset))]

    fig, ax = plt.subplots(figsize=(5, 3))
    bars = ax.bar(x_pos, subset["log2FC"], color=subset["color"], width=0.35, edgecolor="black", linewidth=2.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(subset["cancer_type"], fontsize=18)
    ax.axhline(0, color="black", linewidth=2.5)
    ax.set_ylim(-2, 3)
    ax.set_ylabel("~ log2 (FC)", fontsize=18, fontweight="bold")
    ax.tick_params(axis="y", labelsize=18)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--mirna-id", default=TARGET_MIRNA_ID)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    df = pd.read_csv(args.input)
    plot_mirna_logfc_by_cancer(df, args.mirna_id, args.output)
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
