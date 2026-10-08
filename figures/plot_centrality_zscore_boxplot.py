"""
Boxplot comparing normalized (z-scored) network betweenness centrality
between two network conditions (e.g. Tumor vs. Normal ceRNA networks).

Expected input
--------------
The output of ``analysis.network_centrality.betweenness_zscore_summary``:
a DataFrame with columns ``network_type`` and ``zscore``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from config import FIGURE_DPI
from utils import apply_default_style

GROUP_PALETTE = {"Tumor": "#E78AC3", "Normal": "#87CEEB"}


def plot_centrality_zscore_boxplot(df: pd.DataFrame, output_path: Path) -> None:
    apply_default_style(font_size=24)

    fig = plt.figure(figsize=(6.5, 3))
    sns.boxplot(
        data=df, y="network_type", x="zscore", showfliers=False, width=0.6,
        boxprops=dict(facecolor="black", edgecolor="black"),
        medianprops=dict(color="white", linewidth=1),
    )
    sns.stripplot(
        data=df, y="network_type", x="zscore", jitter=0.25, size=15,
        palette=GROUP_PALETTE, alpha=0.85, edgecolor="black", linewidth=0.8,
    )

    plt.xlabel("Z-score (Betweenness Centrality)", labelpad=12)
    plt.ylabel("")
    ax = plt.gca()
    ax.tick_params(axis="y", labelsize=25, width=1.5)
    for label in ax.get_yticklabels():
        label.set_fontweight("bold")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)
