"""
Figures for the survival analyses in ``analysis/survival_analysis.py``:
a Kaplan-Meier curve, a Cox model hazard-ratio forest plot, and a bootstrap
C-index validation histogram.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.stats import norm

from config import FIGURE_DPI, TARGET_MIRNA_ID
from utils import apply_default_style

GROUP_COLORS = {"High": "#67001f", "Low": "#053061"}


def plot_kaplan_meier(km_result: dict, output_path: Path, mirna_label: str = TARGET_MIRNA_ID) -> None:
    """Plot Kaplan-Meier curves from ``survival_analysis.fit_kaplan_meier`` output."""
    apply_default_style()
    fig, ax = plt.subplots(figsize=(5, 5))

    for label, kmf in km_result["fitters"].items():
        kmf.plot(ax=ax, ci_show=True, color=GROUP_COLORS[label], linewidth=2, linestyle="-", legend=False)

    ax.set_title(f"{mirna_label} Survival", fontsize=25, fontweight="bold")
    ax.set_xlabel("Months", fontsize=25, fontweight="bold", labelpad=10)
    ax.set_ylabel("Percent Survival", fontsize=25, fontweight="bold", labelpad=10)
    ax.set_xlim(0, 160)
    ax.set_ylim(0, 1)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.tick_params(axis="both", labelsize=25)

    def _fmt_median(value: float) -> str:
        return f"{value:.1f}m" if not np.isinf(value) else "inf m"

    stats_text = (
        f"Logrank p = {km_result['logrank_p']:.1e}\n"
        f"HR(high) = {km_result['hazard_ratio']:.2f}\n"
        f"Median High = {_fmt_median(km_result['median_high'])}\n"
        f"Median Low  = {_fmt_median(km_result['median_low'])}"
    )
    fig.subplots_adjust(right=0.85)
    fig.text(1.01, 0.55, stats_text, fontsize=18, va="top", fontweight="normal")

    custom_lines = [Line2D([0], [0], color=GROUP_COLORS["Low"], lw=2), Line2D([0], [0], color=GROUP_COLORS["High"], lw=2)]
    fig.legend(custom_lines, ["Low Group", "High Group"], loc="upper left", bbox_to_anchor=(1.01, 0.73), frameon=False, fontsize=18)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=FIGURE_DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)


def plot_cox_hazard_ratios(
    cph: CoxPHFitter,
    output_path: Path,
    covariate_labels: dict[str, str] | None = None,
) -> None:
    """Forest plot of hazard ratios from a fitted Cox model."""
    apply_default_style()
    fig = plt.figure(figsize=(5, 3.5))
    ax = cph.plot(hazard_ratios=True)

    x_min, x_max = ax.get_xlim()
    ax.set_xticks(np.arange(np.floor(x_min), np.ceil(x_max) + 1, 1))
    plt.title("Multivariate Cox Regression", fontsize=20, weight="bold")
    ax.set_xlabel("HR (95% CI)", fontsize=15, weight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)

    if covariate_labels:
        current_labels = [covariate_labels.get(t.get_text(), t.get_text()) for t in ax.get_yticklabels()]
        ax.set_yticklabels(current_labels, fontsize=15, weight="bold")

    plt.subplots_adjust(left=0.35)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.close(fig)


def plot_bootstrap_cindex(c_indices: list[float], output_path: Path) -> None:
    """Histogram of bootstrap C-index values with a fitted normal curve overlay."""
    if not c_indices:
        raise ValueError("No bootstrap C-index values were produced.")

    apply_default_style()
    mean_c = float(np.mean(c_indices))
    std_c = float(np.std(c_indices))

    fig = plt.figure(figsize=(5, 3.5))
    sns.histplot(c_indices, bins=30, color="#C39BD3", edgecolor="black", stat="density")

    x_vals = np.linspace(min(c_indices), max(c_indices), 300)
    plt.plot(x_vals, norm.pdf(x_vals, mean_c, std_c), color="#6A0DAD", linewidth=2.5, label="Normal Curve")
    plt.axvline(mean_c, color="#FFD700", linestyle="--", linewidth=2, label="Mean C-index")

    plt.title("Bootstrap C-index (Internal Validation)", fontsize=18, weight="bold")
    plt.xlabel("C-index", fontsize=18, weight="bold")
    plt.ylabel("Density", fontsize=18, weight="bold")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI)
    plt.close(fig)
