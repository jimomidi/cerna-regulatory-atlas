"""
Summarize drug-gene interaction / druggability annotations for a candidate
target gene set (e.g. from a resource such as DGIdb).

Computes the fraction of annotated interactions that are cancer-related,
regulatory-approval status, and whether the interacting drug acts as an
inhibitor/antagonist/blocker, then reports candidate druggable targets.

Expected input
--------------
interactions : csv
    Table with (at least) columns ``gene``, ``drug``, ``indication``,
    ``regulatory approval`` and ``interaction score``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CANCER_PATTERN = r"cancer|neoplasm|antineoplastic|tumor|oncolog"
INHIBITOR_DRUG_PATTERN = r"inhibitor|antagonist|blocker|anti-|mab$|umab$|zumab$|ib$"
INHIBITOR_INDICATION_PATTERN = r"inhibitor|antagonist|blocker|suppress|suppressor"


def load_and_clean(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    for col in ("drug", "indication", "regulatory approval"):
        df[col] = df[col].astype(str).str.strip().str.lower()
    return df


def annotate_categories(df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean columns flagging cancer relevance, approval status and inhibitor mechanism."""
    df = df.copy()
    df["is_cancer"] = df["indication"].str.contains(CANCER_PATTERN, na=False)
    df["is_approved"] = df["regulatory approval"] == "approved"
    df["is_not_approved"] = df["regulatory approval"] == "not approved"
    df["is_inhibitor"] = df["drug"].str.contains(INHIBITOR_DRUG_PATTERN, na=False) | df[
        "indication"
    ].str.contains(INHIBITOR_INDICATION_PATTERN, na=False)
    df["is_cancer_and_inhibitor"] = df["is_cancer"] & df["is_inhibitor"]
    return df


def summarize_percentages(df: pd.DataFrame) -> pd.Series:
    """Percentage of all interactions falling into each category combination."""
    total = len(df)
    categories = {
        "cancer_approved": df["is_cancer"] & df["is_approved"],
        "cancer_not_approved": df["is_cancer"] & df["is_not_approved"],
        "inhibitor_approved": df["is_inhibitor"] & df["is_approved"],
        "inhibitor_not_approved": df["is_inhibitor"] & df["is_not_approved"],
        "cancer_inhibitor_approved": df["is_cancer_and_inhibitor"] & df["is_approved"],
        "cancer_inhibitor_not_approved": df["is_cancer_and_inhibitor"] & df["is_not_approved"],
    }
    return pd.Series({name: round(mask.sum() / total * 100, 2) for name, mask in categories.items()})


def cancer_inhibitor_targets(df: pd.DataFrame, min_score: float | None = None) -> pd.DataFrame:
    """List candidate cancer-relevant inhibitor interactions, optionally filtered by score."""
    subset = df[df["is_cancer_and_inhibitor"]].copy()
    if min_score is not None:
        subset = subset[subset["interaction score"] > min_score]
    columns = [c for c in ("gene", "drug", "interaction score", "regulatory approval") if c in subset.columns]
    result = subset[columns].copy()
    result["regulatory approval"] = result["regulatory approval"].str.title()
    return result


def interaction_type_counts(df: pd.DataFrame, type_col: str = "interaction_type") -> pd.Series:
    """Count interactions by mechanism type, if the input table provides one."""
    return df[type_col].value_counts()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interactions", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-score", type=float, default=None)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = annotate_categories(load_and_clean(args.interactions))

    percentages = summarize_percentages(df)
    print(percentages.to_string())
    percentages.to_csv(args.output_dir / "interaction_summary_percentages.csv")

    targets = cancer_inhibitor_targets(df, min_score=args.min_score)
    suffix = f"_score_gt_{args.min_score}" if args.min_score is not None else ""
    out_path = args.output_dir / f"cancer_inhibitor_targets{suffix}.csv"
    targets.to_csv(out_path, index=False)
    print(f"\nSaved {len(targets)} candidate targets to {out_path}")


if __name__ == "__main__":
    main()
