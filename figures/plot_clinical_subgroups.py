"""
Box + strip plots of target miRNA expression across clinical subgroups
(age, gender, pathological stage) and a Normal-vs-Tumor comparison.

Design note
-----------
The shared boxplot style (``_BOX_KWS``) and the y-axis upper limit
(``ymax``, derived from the data with 10% headroom) are defined explicitly
in this module. They affect plot appearance only; no statistical
calculation depends on them.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from analysis.clinical_association import AGE_ORDER, STAGE_ORDER
from config import FIGURE_DPI
from utils import apply_default_style

AGE_COLOR = "#4682B4"
GENDER_COLOR = "#B23A6F"
STAGE_COLOR = "#228B22"
NORMAL_COLOR = "#4682B4"
TUMOR_COLOR = "#B23A6F"

# Fixed boxplot style shared by the three clinical-subgroup panels below.
_BOX_KWS = dict(color="black", width=0.5, fliersize=0, medianprops={"color": "white", "linewidth": 1.4})


def _panel(df: pd.DataFrame, x: str, order: list[str], color: str, title: str, ymax: float, ax) -> None:
    sns.boxplot(data=df, x=x, y="RPM", order=order, ax=ax, **_BOX_KWS)
    sns.stripplot(data=df, x=x, y="RPM", order=order, color=color, size=8, jitter=True, alpha=0.6, ax=ax)
    ax.set_title(title, fontweight="bold", fontsize=20)
    ax.set_xlabel("")
    ax.set_ylabel("RPM", fontsize=22, fontweight="bold", labelpad=12)
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax + 1, 10))
    ax.grid(False)


def plot_clinical_subgroups(df: pd.DataFrame, output_dir: Path, file_prefix: str = "expression") -> None:
    """Save three panels (age, gender, stage) as separate PNG files."""
    apply_default_style(font_size=20)
    output_dir.mkdir(parents=True, exist_ok=True)
    ymax = float(df["RPM"].max()) * 1.1

    fig1, ax1 = plt.subplots(figsize=(4, 4))
    _panel(df, "age_group", AGE_ORDER, AGE_COLOR, "Age Group", ymax, ax1)
    ax1.set_xticklabels(["21-40", "41-60", "61-80"], fontsize=22)
    fig1.tight_layout()
    fig1.savefig(output_dir / f"{file_prefix}_age.png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(4, 4))
    _panel(df, "gender", ["female", "male"], GENDER_COLOR, "Gender Group", ymax, ax2)
    ax2.set_xticklabels(["Female", "Male"], fontsize=22)
    fig2.tight_layout()
    fig2.savefig(output_dir / f"{file_prefix}_gender.png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig2)

    df_stage = df[df["stage_group"].isin(STAGE_ORDER)]
    fig3, ax3 = plt.subplots(figsize=(4, 4))
    _panel(df_stage, "stage_group", STAGE_ORDER, STAGE_COLOR, "Stage Group", ymax, ax3)
    ax3.set_xticklabels([f"S #{i+1}" for i in range(len(STAGE_ORDER))], fontsize=22)
    fig3.tight_layout()
    fig3.savefig(output_dir / f"{file_prefix}_stage.png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig3)


def plot_group_comparison(df: pd.DataFrame, output_path: Path, value_col: str = "log_RPM_adj") -> None:
    """Normal vs Tumor box + strip plot."""
    apply_default_style(font_size=20)
    group_order = ["Normal", "Tumor"]

    fig = plt.figure(figsize=(4, 4))
    ax = plt.gca()
    sns.boxplot(data=df, x="Group", y=value_col, order=group_order, ax=ax, **_BOX_KWS)
    sns.stripplot(data=df[df["Group"] == "Normal"], x="Group", y=value_col, order=group_order, color=NORMAL_COLOR, size=8, jitter=True, alpha=0.6, ax=ax)
    sns.stripplot(data=df[df["Group"] == "Tumor"], x="Group", y=value_col, order=group_order, color=TUMOR_COLOR, size=8, jitter=True, alpha=0.6, ax=ax)

    ax.set_title("Normal vs Tumor", fontweight="bold", fontsize=20)
    ax.set_xlabel("")
    ax.set_ylabel("log2(RPM+0.01)", fontsize=22, fontweight="bold", labelpad=12)
    ax.set_xticklabels(["Normal", "Tumor"], fontsize=22)
    ax.grid(False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)
