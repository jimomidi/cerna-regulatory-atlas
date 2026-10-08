"""
Log-log degree distribution scatter plot with a smoothing spline, comparing
two ceRNA networks (e.g. Tumor vs. Normal).

Expected input
--------------
Two DataFrames as returned by ``analysis.network_centrality.degree_distribution``
(columns ``degree``, ``node_count``).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

from config import FIGURE_DPI
from utils import apply_default_style


def _fit_spline(degree_df: pd.DataFrame, smoothing: float = 0.8, n_points: int = 400):
    log_x = np.log10(degree_df["degree"])
    log_y = np.log10(degree_df["node_count"])
    spline = UnivariateSpline(log_x, log_y, s=smoothing)
    x_fit = np.logspace(log_x.min(), log_x.max(), n_points)
    y_fit = 10 ** spline(np.log10(x_fit))
    return x_fit, y_fit


def plot_degree_distribution(
    normal_degrees: pd.DataFrame,
    tumor_degrees: pd.DataFrame,
    output_path: Path,
) -> None:
    apply_default_style(font_size=24)

    x_fit_n, y_fit_n = _fit_spline(normal_degrees)
    x_fit_t, y_fit_t = _fit_spline(tumor_degrees)

    fig = plt.figure(figsize=(6, 4))
    plt.plot(x_fit_n, y_fit_n, color="#4682B4", linewidth=3, zorder=1)
    plt.plot(x_fit_t, y_fit_t, color="#C71585", linewidth=3, zorder=1)
    plt.scatter(normal_degrees["degree"], normal_degrees["node_count"], s=200, c="#87CEEB", edgecolors="black", label="Normal", zorder=2)
    plt.scatter(tumor_degrees["degree"], tumor_degrees["node_count"], s=200, c="#E78AC3", edgecolors="black", label="Tumor", zorder=2)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Degree (log scale)", labelpad=15)
    plt.ylabel("Number of Nodes\n(log scale)", labelpad=15)
    plt.legend(frameon=False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)
