"""
Figures for the drug-gene interaction summary in
``analysis/drug_gene_interaction.py``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from config import FIGURE_DPI
from utils import apply_default_style

APPROVED_COLOR = "#1B9E77"
NOT_APPROVED_COLOR = "#D95F9A"


def plot_approval_summary(percentages: pd.Series, output_path: Path) -> None:
    """Grouped bar chart of Cancer / Inhibitor / Cancer+Inhibitor interactions
    by regulatory approval status.

    ``percentages`` is expected to have the index produced by
    ``drug_gene_interaction.summarize_percentages``.
    """
    apply_default_style(font_size=24)

    values = [
        percentages["cancer_approved"], percentages["cancer_not_approved"],
        percentages["inhibitor_approved"], percentages["inhibitor_not_approved"],
        percentages["cancer_inhibitor_approved"], percentages["cancer_inhibitor_not_approved"],
    ]
    x = np.array([0, 0.25, 0.55, 0.8, 1.1, 1.35])
    bar_width = 0.18
    colors = [APPROVED_COLOR, NOT_APPROVED_COLOR] * 3

    fig = plt.figure(figsize=(6.3, 3.5))
    bars = plt.bar(x, values, width=bar_width, color=colors)
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height + 0.3, f" {height:.2f}%", ha="center", fontsize=20, rotation=90)

    xtick_positions = [(x[0] + x[1]) / 2, (x[2] + x[3]) / 2, (x[4] + x[5]) / 2]
    plt.xticks(xtick_positions, ["Cancer", "Inhibitor", "Cancer + Inhibitor"], fontsize=18)
    plt.yticks(fontsize=24)
    plt.ylabel("Percentage (%)", fontsize=26, weight="bold", labelpad=10)
    plt.ylim(0, max(values) + 5.0)

    legend_elements = [Patch(facecolor=APPROVED_COLOR, label="Approved"), Patch(facecolor=NOT_APPROVED_COLOR, label="Not Approved")]
    plt.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0, 1.04), fontsize=20, frameon=False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.close(fig)


def plot_interaction_type_counts(counts: pd.Series, output_path: Path) -> None:
    """Bar chart of interaction counts by mechanism type."""
    apply_default_style(font_size=18)
    fig = plt.figure(figsize=(4, 5.5))
    x = np.arange(len(counts))
    bars = plt.bar(x, counts.values, color="black", width=0.6)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height + 3, f"{int(height)}", ha="center", fontsize=16)

    plt.xticks(x, counts.index, fontsize=18, rotation=20)
    plt.yticks(fontsize=18)
    plt.ylabel("Count", fontsize=22, weight="bold", labelpad=10)
    plt.ylim(0, counts.max() + 25)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.close(fig)
