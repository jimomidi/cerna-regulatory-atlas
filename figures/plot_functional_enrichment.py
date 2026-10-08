"""
Bubble plots for functional enrichment results (e.g. GO or Reactome
over-representation analysis of ceRNA network target genes).

Expected input
--------------
csv with columns ``Term``, ``Count``, ``List Total`` and ``FDR``; for the
GO variant, also a ``Category`` column identifying each term as Biological
Process (BP), Molecular Function (MF) or Cellular Component (CC).

Usage
-----
    python -m figures.plot_functional_enrichment \
        --input data/reactome_enrichment.csv --kind reactome \
        --output results/reactome_enrichment.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import ListedColormap

from config import FIGURE_DPI
from utils import apply_default_style

GO_CATEGORY_BACKGROUND = {"BP": "#D0F5E8", "CC": "#FBE4EC", "MF": "#D6ECFF"}


def _prepare(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.str.strip()
    df["GeneRatio"] = df["Count"] / df["List Total"]
    df["Term"] = df["Term"].str.strip()
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    return df


def _fdr_colormap() -> ListedColormap:
    base = cm.get_cmap("hsv", 256)
    return ListedColormap(base(np.linspace(0.3, 1, 256)))


def plot_reactome_bubble(df: pd.DataFrame, output_path: Path, title: str = "Reactome Enrichment") -> None:
    apply_default_style()
    df = _prepare(df)

    fig = plt.figure(figsize=(4, 6.5))
    bubble = sns.scatterplot(
        data=df, x="GeneRatio", y="Term", size="Count", hue="FDR",
        sizes=(100, 600), palette=_fdr_colormap(), edgecolor="black", alpha=0.8, legend="brief",
    )
    bubble.set_xlabel("GeneRatio", fontsize=20, weight="bold", labelpad=10)
    bubble.set_ylabel("")
    plt.title(title, fontsize=20, weight="bold")
    bubble.tick_params(axis="both", labelsize=18)
    bubble.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
    plt.legend(title="FDR", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.close(fig)


def plot_go_bubble(df: pd.DataFrame, output_path: Path, title: str = "GO Enrichment") -> None:
    apply_default_style()
    df = _prepare(df)

    def _go_type(category: str) -> str:
        for key in ("BP", "MF", "CC"):
            if key in category:
                return key
        return "Other"

    df["GO_Type"] = df["Category"].apply(_go_type)

    fig, ax = plt.subplots(figsize=(9.5, 6))
    for i, term in enumerate(df["Term"]):
        go_type = df.loc[df["Term"] == term, "GO_Type"].values[0]
        ax.axhspan(i - 0.5, i + 0.5, color=GO_CATEGORY_BACKGROUND.get(go_type, "#FFFFFF"), alpha=0.6, zorder=0)

    sns.scatterplot(
        data=df, x="GeneRatio", y="Term", size="Count", hue="FDR",
        sizes=(120, 650), palette=_fdr_colormap(), edgecolor="black", alpha=0.85, legend="brief", zorder=2, ax=ax,
    )

    ax.set_xlabel("GeneRatio", fontsize=20, weight="bold", labelpad=10)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=20, weight="bold")
    ax.tick_params(axis="both", labelsize=18)
    ax.set_ylim(-0.5, len(df["Term"].unique()) - 0.5)
    ax.legend(title="FDR", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--kind", choices=["reactome", "go"], required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    df = pd.read_csv(args.input)
    if args.kind == "reactome":
        plot_reactome_bubble(df, args.output)
    else:
        plot_go_bubble(df, args.output)
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
