"""
Time-dependent ROC AUC plot with a bootstrap confidence band.

This module only visualizes precomputed time-dependent AUC values; it does
not compute them. Time-dependent AUC with bootstrap confidence intervals is
typically produced by a dedicated survival-ROC method (e.g. R's ``timeROC``
package, or ``scikit-survival`` in Python) -- run that analysis first and
save its output as a CSV with the columns below.

Expected input
--------------
csv with columns ``time`` (follow-up time, e.g. months), ``auc``,
``ci_low`` and ``ci_high``.

Usage
-----
    python -m figures.plot_time_dependent_auc \
        --input data/time_dependent_auc.csv \
        --output results/time_dependent_auc.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from config import FIGURE_DPI
from utils import apply_default_style


def plot_time_dependent_auc(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "Time-dependent ROC AUC (Bootstrap 95% CI)",
    band_color: str = "#5AB5F6",
) -> None:
    apply_default_style()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.fill_between(df["time"], df["ci_low"], df["ci_high"], color=band_color, alpha=0.2)
    ax.plot(df["time"], df["auc"], color="black", linewidth=2, marker="o", markersize=10, label="AUC")
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1)

    ax.set_xlabel("Time (Months)", fontsize=18, weight="bold", labelpad=10)
    ax.set_ylabel("Time-dependent AUC", fontsize=18, weight="bold", labelpad=10)
    ax.set_title(title, fontsize=16, weight="bold", pad=12)
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xticks(df["time"])
    ax.set_ylim(0.4, 1.02)
    ax.set_xlim(df["time"].min(), df["time"].max())

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI)
    plt.close(fig)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--title", default="Time-dependent ROC AUC (Bootstrap 95% CI)")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    df = pd.read_csv(args.input)
    plot_time_dependent_auc(df, args.output, title=args.title)
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
