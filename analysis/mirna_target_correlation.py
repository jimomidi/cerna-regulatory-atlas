"""
miRNA-mRNA expression correlation with multiple-testing correction.

For a miRNA of interest and a candidate list of target genes, this script
computes the Pearson correlation between log2-transformed miRNA and mRNA
expression across paired tumor samples, applies Benjamini-Hochberg FDR
correction across all tested pairs, and reports pairs passing a negative
correlation / FDR threshold (the expected signature of miRNA-mediated
target repression).

The candidate gene list is generic: it can be predicted targets from a
target-prediction database, a differential-expression gene set, or any
other gene list of interest -- simply point ``--gene-list`` at the
appropriate CSV file.

Expected inputs
----------------
mirna_expression : csv
    Long-format table with at least the columns ``miRNA_ID``, ``case_id``
    and ``reads_per_million_miRNA_mapped``.
mrna_expression : csv
    Long-format table with at least the columns ``gene_name``, ``case_id``
    and ``fpkm_unstranded``.
gene_list : csv
    Table with a ``gene`` column listing candidate target genes.

Output
------
A CSV with one row per tested (miRNA, gene) pair, including the Pearson
correlation coefficient, p-value, FDR-adjusted p-value and sample size.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

from config import (
    FDR_THRESHOLD,
    MIN_EXPRESSED_FRACTION,
    MIN_PAIRED_SAMPLES,
    TARGET_MIRNA_ID,
)


def _filter_low_expression(
    df: pd.DataFrame,
    id_col: str,
    value_col: str,
    min_expressed_fraction: float = MIN_EXPRESSED_FRACTION,
) -> pd.DataFrame:
    """Drop features expressed above 1 unit in fewer than the given fraction of samples."""
    wide = df.pivot_table(index=id_col, columns="case_id", values=value_col)
    keep = wide.index[wide.gt(1).sum(axis=1) >= min_expressed_fraction * wide.shape[1]]
    return df[df[id_col].isin(keep)]


def compute_mirna_target_correlations(
    mirna_expr: pd.DataFrame,
    mrna_expr: pd.DataFrame,
    gene_list: Iterable[str],
    mirna_ids: Iterable[str] = (TARGET_MIRNA_ID,),
    min_paired_samples: int = MIN_PAIRED_SAMPLES,
) -> pd.DataFrame:
    """Compute Pearson correlations between miRNA and candidate target-gene expression.

    Parameters
    ----------
    mirna_expr, mrna_expr
        Long-format expression tables as described in the module docstring.
    gene_list
        Candidate target gene symbols to test.
    mirna_ids
        miRNA identifiers to test (usually a single target miRNA).
    min_paired_samples
        Minimum number of samples with both miRNA and mRNA measurements
        required to attempt a correlation.

    Returns
    -------
    DataFrame with columns ``miRNA``, ``gene``, ``r``, ``p``, ``n`` and,
    after FDR correction, ``FDR``.
    """
    mirna_expr = mirna_expr.copy()
    mrna_expr = mrna_expr.copy()

    mirna_expr = _filter_low_expression(mirna_expr, "miRNA_ID", "reads_per_million_miRNA_mapped")
    mrna_expr = _filter_low_expression(mrna_expr, "gene_name", "fpkm_unstranded")

    mirna_expr["log2_rpm"] = np.log2(mirna_expr["reads_per_million_miRNA_mapped"] + 0.01)
    mrna_expr["log2_fpkm"] = np.log2(mrna_expr["fpkm_unstranded"] + 0.01)

    gene_list = list(dict.fromkeys(gene_list))  # de-duplicate, preserve order
    results: List[dict] = []

    for mirna in mirna_ids:
        mir_df = (
            mirna_expr.loc[mirna_expr["miRNA_ID"] == mirna, ["case_id", "log2_rpm"]]
            .drop_duplicates()
        )
        for gene in gene_list:
            gene_df = (
                mrna_expr.loc[mrna_expr["gene_name"] == gene, ["case_id", "log2_fpkm"]]
                .drop_duplicates()
            )
            merged = mir_df.merge(gene_df, on="case_id")
            if len(merged) < min_paired_samples:
                continue

            r, p = pearsonr(merged["log2_rpm"], merged["log2_fpkm"])
            results.append({"miRNA": mirna, "gene": gene, "r": r, "p": p, "n": len(merged)})

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df["FDR"] = multipletests(result_df["p"], method="fdr_bh")[1]
    return result_df


def significant_targets(
    result_df: pd.DataFrame,
    correlation_threshold: float,
    fdr_threshold: float = FDR_THRESHOLD,
) -> pd.DataFrame:
    """Filter correlation results to pairs meeting the significance criteria."""
    if result_df.empty:
        return result_df
    return result_df[(result_df["r"] < correlation_threshold) & (result_df["FDR"] < fdr_threshold)]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mirna-expression", required=True, type=Path)
    parser.add_argument("--mrna-expression", required=True, type=Path)
    parser.add_argument("--gene-list", required=True, type=Path, help="CSV with a 'gene' column")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mirna-id", default=TARGET_MIRNA_ID)
    parser.add_argument("--correlation-threshold", type=float, default=-0.3)
    parser.add_argument("--fdr-threshold", type=float, default=FDR_THRESHOLD)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    mirna_expr = pd.read_csv(args.mirna_expression)
    mrna_expr = pd.read_csv(args.mrna_expression)
    gene_list = pd.read_csv(args.gene_list)["gene"].dropna().astype(str).unique().tolist()

    result_df = compute_mirna_target_correlations(
        mirna_expr, mrna_expr, gene_list, mirna_ids=[args.mirna_id]
    )

    if result_df.empty:
        print("No correlation results were produced. Check that the input files "
              "contain overlapping sample IDs and that the gene list is non-empty.")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(args.output, index=False)
    print(f"Saved {len(result_df)} correlation results to {args.output}")

    sig = significant_targets(result_df, args.correlation_threshold, args.fdr_threshold)
    print(
        f"{len(sig)} gene(s) pass r < {args.correlation_threshold} and "
        f"FDR < {args.fdr_threshold}:"
    )
    if not sig.empty:
        print(", ".join(sig["gene"].tolist()))


if __name__ == "__main__":
    main()
