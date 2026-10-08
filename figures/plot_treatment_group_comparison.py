"""
Figures comparing target miRNA expression between treatment groups: a
box + strip plot, and a raincloud plot (distribution + boxplot combined).

Expected input
--------------
Output of ``analysis.treatment_group_comparison.prepare_treatment_comparison_table``:
a DataFrame with ``miRNA``, ``treatment_group``, ``expression`` and
``log_expression`` columns.

Design note
-----------
This module only plots data-driven comparisons from a caller-supplied
table. It does not embed sample identifiers, and it does not draw
per-sample pairings unless real pairing information is provided; the
box/strip and raincloud views below make no pairing assumption.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from config import FIGURE_DPI
from utils import apply_default_style

TREATMENT_PALETTE = ["#8b002a", "#053061"]


def plot_treatment_boxplot(df: pd.DataFrame, output_path: Path) -> None:
    apply_default_style()
    fig = plt.figure(figsize=(4, 5))

    sns.boxplot(data=df, x="miRNA", y="log_expression", hue="treatment_group", palette=TREATMENT_PALETTE, showfliers=False)
    sns.stripplot(data=df, x="miRNA", y="log_expression", hue="treatment_group", dodge=True, jitter=True, color="black", alpha=0.7, size=8, legend=False)

    plt.ylabel("log10(RPM)", fontsize=22, weight="bold", labelpad=10)
    plt.xlabel("")
    plt.legend(fontsize=15, frameon=False)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close(fig)


def plot_treatment_raincloud(df: pd.DataFrame, output_path: Path) -> None:
    """Requires the optional ``ptitprince`` package (see requirements.txt)."""
    import ptitprince as pt

    apply_default_style()
    fig = plt.figure(figsize=(5, 6))
    ax = pt.RainCloud(
        x="treatment_group", y="log_expression", data=df,
        palette=["#1B9E77", "#D95F9A"], width_viol=0.6, orient="v", alpha=0.7, move=0.2,
        rain_alpha=0, box_showfliers=False,
    )

    for artist in ax.artists:
        artist.set_edgecolor("black")
        artist.set_linewidth(2)
    for line in ax.lines:
        line.set_color("black")
        line.set_linewidth(2)

    plt.ylabel("log10(CPM)", fontsize=20, weight="bold")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=FIGURE_DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)
