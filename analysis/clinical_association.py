"""
Association between target miRNA expression and clinical variables.

Two comparisons are provided:

1. ``prepare_subgroup_table`` -- merges miRNA expression with clinical
   metadata and buckets samples into age, gender and pathological-stage
   subgroups, for visual inspection of expression differences across
   clinical categories.
2. ``prepare_group_table`` -- a simple Normal-vs-Tumor comparison from a
   two-group expression table.

Plotting for both lives in ``figures/plot_clinical_subgroups.py``.

Expected inputs
----------------
clinical : csv
    Table with ``case_id``, ``age_at_index``, ``gender`` and
    ``ajcc_pathologic_stage`` (or ``ajcc_clinical_stage``) columns.
mirna_expression : csv
    Long-format table with ``miRNA_ID``, ``case_id`` and
    ``reads_per_million_miRNA_mapped`` columns.

Note
----
The shared plotting style and axis limits used for these comparisons are
defined in ``figures/plot_clinical_subgroups.py``; see that module's
docstring.
"""

from __future__ import annotations

import re

import pandas as pd

from config import TARGET_MIRNA_ID

AGE_ORDER = ["21-40 Yrs", "41-60 Yrs", "61-80 Yrs"]
STAGE_ORDER = ["Stage 1", "Stage 2", "Stage 3", "Stage 4"]


def _categorize_age(age: float) -> str:
    if 21 <= age <= 40:
        return "21-40 Yrs"
    if 41 <= age <= 60:
        return "41-60 Yrs"
    if 61 <= age <= 80:
        return "61-80 Yrs"
    return "Other"


def _normalize_stage(stage: object) -> str | None:
    if stage is None:
        return None
    raw = str(stage).strip()
    if raw in ("--", "", "nan", "None"):
        return None
    lowered = raw.lower().replace("stage", "").strip()
    match = re.search(r"(iv|iii|ii|i|[1-4])", lowered)
    if not match:
        return None
    token = match.group(1)
    roman_to_stage = {"i": "Stage 1", "ii": "Stage 2", "iii": "Stage 3", "iv": "Stage 4"}
    if token in roman_to_stage:
        return roman_to_stage[token]
    if token.isdigit():
        return f"Stage {token}"
    return None


def prepare_subgroup_table(
    clinical: pd.DataFrame,
    mirna_expr: pd.DataFrame,
    mirna_id: str = TARGET_MIRNA_ID,
) -> pd.DataFrame:
    """Merge expression with clinical metadata and derive age/stage subgroup labels."""
    clinical = clinical.copy()
    clinical.columns = clinical.columns.str.strip()

    stage_col = "ajcc_pathologic_stage" if "ajcc_pathologic_stage" in clinical.columns else "ajcc_clinical_stage"
    clinical = clinical[["case_id", "age_at_index", "gender", stage_col]].rename(columns={stage_col: "stage"})

    mirna_expr = mirna_expr.copy()
    mirna_expr.columns = mirna_expr.columns.str.strip()
    mirna_expr = mirna_expr[mirna_expr["miRNA_ID"].astype(str).str.strip().str.lower() == mirna_id.lower()]
    if mirna_expr.empty:
        raise ValueError(f"miRNA '{mirna_id}' not found in the expression table.")

    mirna_expr = (
        mirna_expr[["case_id", "reads_per_million_miRNA_mapped"]]
        .rename(columns={"reads_per_million_miRNA_mapped": "RPM"})
        .groupby("case_id", as_index=False)["RPM"]
        .mean()
    )

    df = pd.merge(mirna_expr, clinical, on="case_id", how="inner")
    df["age_at_index"] = pd.to_numeric(df["age_at_index"], errors="coerce")
    df["RPM"] = pd.to_numeric(df["RPM"], errors="coerce")
    df["gender"] = df["gender"].astype(str).str.strip().str.lower()
    df["stage"] = df["stage"].astype(str).str.strip()
    df = df.dropna(subset=["age_at_index", "RPM", "gender"])
    df = df[df["gender"].isin(["male", "female"])]

    df["age_group"] = df["age_at_index"].apply(_categorize_age)
    df["stage_group"] = df["stage"].apply(_normalize_stage)
    return df


def prepare_group_table(expr_df: pd.DataFrame, group_col: str = "Group", value_col: str = "log_RPM_adj") -> pd.DataFrame:
    """Clean a two-group (e.g. Normal vs Tumor) expression table for plotting."""
    df = expr_df.copy()
    df.columns = df.columns.str.strip()
    df = df.dropna(subset=[group_col, value_col])
    df[group_col] = df[group_col].astype(str).str.strip()
    return df[df[group_col].isin(["Normal", "Tumor"])]
