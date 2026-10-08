"""
Summary bar charts for a curated list of target genes:

1. ``plot_reported_vs_computed_logfc`` -- compares a reported logFC (e.g.
   from a published differential-expression table) against logFC
   recomputed directly from tumor/normal expression matrices, with the
   discrepancy shown as an error bar.
2. ``plot_median_expression`` -- median log2(TPM + 1) expression (with
   standard deviation error bars) for the same gene list in a single
   sample group (call once for tumor, once for normal).

Expected inputs
----------------
tumor_expression : csv
    Long-format table with ``gene_name``, ``case_id`` and
    ``tpm_unstranded`` columns.
normal_expression : csv
    Wide-format table with a gene-name column (e.g. GTEx-style ``Name`` /
    ``Description`` export) and one column per sample.
reported_logfc : csv
    Table with ``gene`` and ``logFC`` columns (e.g. concatenated
    up/down-regulated differential expression results).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import FIGURE_DPI
from utils import apply_default_style


def compute_logfc_comparison(
    reported_logfc: pd.DataFrame,
    tumor_expr: pd.DataFrame,
    normal_expr: pd.DataFrame,
    gene_list: Sequence[str],
    normal_gene_col: str = "Description",
    normal_drop_cols: Sequence[str] = ("Name",),
) -> pd.DataFrame:
    """Combine reported and recomputed logFC for a gene list, returning a
    DataFrame with columns ``logFC`` (reported), ``Computed_logFC`` and
    ``Deviation``.
    """
    reported = reported_logfc.copy()
    reported.columns = reported.columns.str.strip()
    reported["gene"] = reported["gene"].str.upper()
    reported = reported[reported["gene"].isin(gene_list)][["gene", "logFC"]]
    reported = reported.groupby("gene").median().reindex(gene_list)

    tumor_pivot = tumor_expr.pivot_table(index="gene_name", columns="case_id", values="tpm_unstranded", aggfunc="sum")
    tumor_pivot.index = tumor_pivot.index.str.upper()
    tumor_log = np.log2(tumor_pivot + 1)

    normal = normal_expr.copy()
    normal.columns = normal.columns.str.strip()
    normal[normal_gene_col] = normal[normal_gene_col].str.upper()
    normal = normal[normal[normal_gene_col].isin(gene_list)].set_index(normal_gene_col)
    normal = normal.drop(columns=[c for c in normal_drop_cols if c in normal.columns])
    normal_log = np.log2(normal.T.apply(pd.to_numeric, errors="coerce") + 1)

    computed = {}
    for gene in gene_list:
        if gene in tumor_log.index and gene in normal_log.columns:
            computed[gene] = tumor_log.loc[gene].median(skipna=True) - normal_log[gene].median(skipna=True)
        else:
            computed[gene] = np.nan

    reported["Computed_logFC"] = pd.Series(computed)
    reported["Deviation"] = reported["logFC"] - reported["Computed_logFC"]
    return reported


def plot_reported_vs_computed_logfc(summary: pd.DataFrame, output_path: Path, title: str = "logFC") -> None:
    apply_default_style()
    lower_errors = np.where(summary["logFC"] < 0, summary["Deviation"].abs(), 0)
    upper_errors = np.where(summary["logFC"] > 0, summary["Deviation"].abs(), 0)

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.bar(
        x=summary.index,
        height=summary["logFC"],
        yerr=[lower_errors, upper_errors],
        capsize=5,
        color="black",
        edgecolor="black",
        linewidth=2.5,
    )
    ax.axhline(0, color="black", linewidth=1.5, linestyle="--")
    ax.set_title(title, fontsize=30, weight="bold", pad=10)
    ax.set_xticks(np.arange(len(summary.index)))
    ax.set_xticklabels(summary.index, rotation=90, ha="center", fontsize=30)
    ax.tick_params(axis="y", labelsize=22)
    plt.margins(x=0.01)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_median_expression(
    expr_log: pd.DataFrame,
    gene_list: Sequence[str],
    output_path: Path,
    color: str,
    title: str = "Median log2(TPM + 1)",
    low_expr_threshold: float = 0.5,
    min_sample_fraction: float = 0.1,
) -> None:
    """``expr_log`` is a samples-by-genes matrix (NOT yet log-transformed);
    low-expression values are masked before log2 transformation, using the
    ``low_expr_threshold`` and ``min_sample_fraction`` filtering rule.
    """
    apply_default_style()
    expr = expr_log[list(gene_list)].copy()
    min_samples = int(expr.shape[0] * min_sample_fraction)
    for gene in expr.columns:
        low_expr = expr[gene] <= low_expr_threshold
        if low_expr.sum() > (expr.shape[0] - min_samples):
            expr.loc[low_expr, gene] = np.nan
    expr_log2 = np.log2(expr + 1)

    medians = expr_log2.median(skipna=True)
    stds = expr_log2.std(skipna=True)

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.bar(
        x=medians.index,
        height=medians.values,
        yerr=[np.zeros_like(stds.values), stds.values],
        capsize=5,
        color=color,
        edgecolor="black",
        linewidth=1.8,
    )
    ax.set_title(title, fontsize=30, weight="bold", pad=15)
    ax.set_xticks(np.arange(len(medians.index)))
    ax.set_xticklabels(medians.index, rotation=90, ha="center", fontsize=30)
    ax.tick_params(axis="y", labelsize=22)
    ax.set_ylim(bottom=0)
    plt.margins(x=0.01)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
