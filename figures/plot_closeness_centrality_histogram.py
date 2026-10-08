"""
Histogram (with kernel density estimate) comparing closeness centrality
distributions between two ceRNA networks (e.g. Tumor vs. Normal).

Expected input
--------------
Two array-likes of closeness centrality values, e.g. from
``analysis.network_centrality.closeness_centrality``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import seaborn as sns

from config import FIGURE_DPI
from utils import apply_default_style


def plot_closeness_centrality_histogram(
    tumor_values: Sequence[float],
    normal_values: Sequence[float],
    output_path: Path,
    bins: int = 40,
) -> None:
    apply_default_style(font_size=24)

    fig = plt.figure(figsize=(7, 5))
    sns.histplot(tumor_values, bins=bins, kde=True, stat="frequency", color="#E78AC3", label="Tumor", alpha=0.8)
    sns.histplot(normal_values, bins=bins, kde=True, stat="frequency", color="#87CEEB", label="Normal", alpha=0.8)

    plt.xscale("log")
    plt.xlabel("Closeness Centrality (log scale)", labelpad=15)
    plt.ylabel("Frequency", labelpad=15)
    plt.legend(frameon=False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)
