"""
Overlap between correlation-derived target genes and differential-expression
results.

Two related comparisons are provided:

1. ``genes_common_to_lists`` -- set overlap between one or more gene lists
   (e.g. correlation-significant targets vs. an independent target-prediction
   list), useful for sanity-checking agreement between methods.
2. ``count_by_logfc_threshold`` -- among genes with a correlation coefficient
   below a chosen cutoff, how many also show an up-regulation logFC above a
   series of thresholds (1, 2, 3), both as raw counts and as a percentage of
   a reference gene-set size.

Expected inputs
----------------
correlation_results : csv
    Output of ``mirna_target_correlation.py`` (must contain ``gene`` and
    ``r`` columns).
differential_expression : csv
    Differential expression results with ``gene`` and ``logFC`` columns.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd


def genes_common_to_lists(*gene_lists: Iterable[str]) -> List[str]:
    """Return genes present in every supplied list (case-insensitive)."""
    normalized = [set(str(g).strip().upper() for g in genes) for genes in gene_lists]
    common = set.intersection(*normalized) if normalized else set()
    return sorted(common)


def _merge_on_gene(correlation_df: pd.DataFrame, deg_df: pd.DataFrame) -> pd.DataFrame:
    correlation_df = correlation_df.dropna(subset=["r", "gene"]).copy()
    deg_df = deg_df.dropna(subset=["logFC", "gene"]).copy()
    correlation_df["gene_key"] = correlation_df["gene"].str.lower()
    deg_df["gene_key"] = deg_df["gene"].str.lower()
    return pd.merge(correlation_df, deg_df, on="gene_key", how="inner", suffixes=("_corr", "_deg"))


def count_by_logfc_threshold(
    correlation_df: pd.DataFrame,
    deg_df: pd.DataFrame,
    correlation_thresholds: Iterable[float] = (-0.3, -0.2),
    logfc_thresholds: Iterable[float] = (1, 2, 3),
    reference_total: int | None = None,
) -> pd.DataFrame:
    """Count (and optionally express as a percentage of ``reference_total``) genes
    that pass both a correlation cutoff and a series of logFC cutoffs.

    Returns a tidy DataFrame with one row per (correlation_threshold, logfc_threshold).
    """
    merged = _merge_on_gene(correlation_df, deg_df)

    rows = []
    for r_thresh in correlation_thresholds:
        subset = merged[merged["r"] < r_thresh]
        for fc_thresh in logfc_thresholds:
            count = int((subset["logFC"] > fc_thresh).sum())
            row = {
                "correlation_threshold": r_thresh,
                "logfc_threshold": fc_thresh,
                "count": count,
            }
            if reference_total:
                row["percent"] = round(count * 100 / reference_total, 2)
            rows.append(row)
    return pd.DataFrame(rows)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--correlation-results", required=True, type=Path)
    parser.add_argument("--differential-expression", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--reference-total",
        type=int,
        default=None,
        help="Total candidate gene count to express results as percentages.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    correlation_df = pd.read_csv(args.correlation_results)
    deg_df = pd.read_csv(args.differential_expression)

    summary = count_by_logfc_threshold(
        correlation_df, deg_df, reference_total=args.reference_total
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, index=False)
    print(summary.to_string(index=False))
    print(f"\nSaved summary to {args.output}")


if __name__ == "__main__":
    main()
