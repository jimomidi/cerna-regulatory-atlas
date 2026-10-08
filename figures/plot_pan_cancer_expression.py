"""
Pan-cancer expression boxplot: target miRNA expression across TCGA tumor
types, split by Normal vs. Tumor sample group.

Expected input
--------------
csv with columns ``project`` (e.g. "TCGA-XXXX"), ``group`` ("normal"/"tumor",
case-insensitive) and ``expr`` (log2(CPM) expression value).

Usage
-----
    python -m figures.plot_pan_cancer_expression \
        --input data/pan_cancer_expression.csv \
        --output results/pan_cancer_expression.png \
        --exclude-project PROJECT_X
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from config import FIGURE_DPI
from utils import apply_default_style

GROUP_PALETTE = {"Normal": "#87CEEB", "Tumor": "#E78AC3"}


def mark_outliers(df: pd.DataFrame, value_col: str = "expr") -> pd.DataFrame:
    """Flag Tukey-fence outliers within each (project, group) subset.

    Implemented with ``groupby().transform`` (rather than ``groupby().apply``)
    so the original columns are always preserved regardless of pandas version.
    """
    df = df.copy()
    grouped = df.groupby(["project", "group"])[value_col]
    q1 = grouped.transform(lambda s: s.quantile(0.25))
    q3 = grouped.transform(lambda s: s.quantile(0.75))
    iqr = q3 - q1
    lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    df["is_outlier"] = ~df[value_col].between(lower, upper)
    return df


def plot_pan_cancer_expression(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "Pan-Cancer Expression",
    exclude_projects: tuple[str, ...] = (),
    reference_line: float | None = None,
) -> None:
    df = df.copy()
    df["project"] = df["project"].str.replace("TCGA-", "", regex=False)
    df["group"] = df["group"].str.capitalize()
    df = df[~df["project"].isin(exclude_projects)]
    df = mark_outliers(df)

    apply_default_style(font_size=20)
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(16, 6))

    sns.boxplot(
        data=df,
        x="project",
        y="expr",
        hue="group",
        linewidth=1.5,
        fliersize=5,
        palette=GROUP_PALETTE,
        medianprops=dict(color="black"),
        boxprops=dict(edgecolor="black"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        flierprops=dict(marker="o", markerfacecolor="black", markeredgecolor="black", markersize=8),
        ax=ax,
    )

    ax.set_xlabel("")
    ax.grid(False)
    ax.spines["bottom"].set_color("black")
    ax.set_title(title, fontsize=24, weight="bold")
    ax.set_ylabel("log2(CPM)", fontsize=24, weight="bold", labelpad=10)
    ax.tick_params(axis="x", labelrotation=90, labelsize=22)
    ax.tick_params(axis="y", labelsize=22)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin - 0.5, xmax + 0.5)
    ax.legend(fontsize=20, loc="upper left", bbox_to_anchor=(0, 1), frameon=False)

    if reference_line is not None:
        ax.axhline(y=reference_line, color="black", linewidth=1)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.close(fig)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--title", default="Pan-Cancer Expression")
    parser.add_argument("--exclude-project", action="append", default=[])
    parser.add_argument("--reference-line", type=float, default=None)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    df = pd.read_csv(args.input)
    plot_pan_cancer_expression(
        df,
        args.output,
        title=args.title,
        exclude_projects=tuple(args.exclude_project),
        reference_line=args.reference_line,
    )
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
