"""
Survival analysis for a target miRNA's expression level.

Three related analyses are provided:

1. ``find_optimal_cutoff`` -- scans a range of expression percentiles to
   find the High/Low expression split that minimizes the log-rank p-value
   (a data-driven alternative to a fixed median split), subject to a
   minimum group size.
2. ``fit_kaplan_meier`` -- fits Kaplan-Meier survival curves for the
   resulting High/Low groups and returns the fitted models plus summary
   statistics (log-rank p-value, hazard ratio, median survival times).
3. ``fit_cox_model`` / ``bootstrap_concordance`` -- fits a multivariate Cox
   proportional-hazards model (expression + clinical covariates) and
   validates it internally via bootstrap resampling of the concordance
   index (C-index).

Expected inputs
----------------
mirna_expression : csv
    Long-format table with ``miRNA_ID``, ``case_id`` and
    ``reads_per_million_miRNA_mapped`` columns.
clinical : csv
    Table with ``case_id``, ``days_to_death``, ``days_to_last_follow_up``,
    ``vital_status``, and (for the Cox model) ``age_at_diagnosis``,
    ``gender`` and ``ajcc_pathologic_stage``.

Overall survival time is defined as min(days to death, days to last
follow-up) / 30, in months; the event indicator is 1 for deceased patients
and 0 for censored (alive) patients.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

from config import MIN_EXPRESSED_FRACTION, N_BOOTSTRAP, TARGET_MIRNA_ID

STAGE_TO_NUMERIC = {"I": 1, "II": 2, "III": 3, "IV": 4}


def prepare_survival_table(
    mirna_expr: pd.DataFrame,
    clinical: pd.DataFrame,
    mirna_id: str = TARGET_MIRNA_ID,
    min_expressed_fraction: float = MIN_EXPRESSED_FRACTION,
) -> pd.DataFrame:
    """Merge miRNA expression with clinical outcome data into one analysis-ready table.

    Returns a DataFrame with columns ``case_id``, ``log2_rpm``, ``OS_time``
    (months) and ``OS_event`` (1 = death, 0 = censored).
    """
    expr_wide = mirna_expr.pivot(
        index="miRNA_ID", columns="case_id", values="reads_per_million_miRNA_mapped"
    )
    low_expression_mask = (expr_wide < 1).sum(axis=1) <= (1 - min_expressed_fraction) * expr_wide.shape[1]
    keep_mirnas = expr_wide.index[low_expression_mask]
    mirna_expr = mirna_expr[mirna_expr["miRNA_ID"].isin(keep_mirnas)]

    expr = mirna_expr.loc[mirna_expr["miRNA_ID"] == mirna_id, ["case_id", "reads_per_million_miRNA_mapped"]].copy()
    expr["log2_rpm"] = np.log2(expr["reads_per_million_miRNA_mapped"] + 0.01)

    clin = clinical[["case_id", "days_to_death", "days_to_last_follow_up", "vital_status"]].drop_duplicates("case_id").copy()
    clin["days_to_death"] = pd.to_numeric(clin["days_to_death"], errors="coerce")
    clin["days_to_last_follow_up"] = pd.to_numeric(clin["days_to_last_follow_up"], errors="coerce")
    clin["OS_time"] = clin[["days_to_death", "days_to_last_follow_up"]].min(axis=1) / 30.0
    clin["OS_event"] = clin["vital_status"].apply(lambda x: 1 if str(x).strip().lower() == "dead" else 0)

    merged = pd.merge(expr, clin, on="case_id").dropna(subset=["OS_time", "OS_event", "log2_rpm"])
    return merged


def find_optimal_cutoff(
    df: pd.DataFrame,
    value_col: str = "log2_rpm",
    percentile_range: Tuple[int, int] = (20, 80),
    min_group_size: int = 25,
) -> float:
    """Search for the expression percentile cutoff that minimizes the
    log-rank p-value between the resulting High/Low survival groups.
    """
    best_p, best_cutoff = np.inf, None
    for percentile in np.arange(*percentile_range, 1):
        cutoff = np.percentile(df[value_col], percentile)
        group = np.where(df[value_col] >= cutoff, "High", "Low")
        if pd.Series(group).value_counts().min() < min_group_size:
            continue

        result = logrank_test(
            df.loc[group == "Low", "OS_time"],
            df.loc[group == "High", "OS_time"],
            event_observed_A=df.loc[group == "Low", "OS_event"],
            event_observed_B=df.loc[group == "High", "OS_event"],
        )
        if result.p_value < best_p:
            best_p, best_cutoff = result.p_value, cutoff

    if best_cutoff is None:
        raise ValueError("No cutoff produced two groups meeting the minimum group size.")
    return best_cutoff


def fit_kaplan_meier(df: pd.DataFrame, cutoff: float, value_col: str = "log2_rpm") -> dict:
    """Fit Kaplan-Meier curves for the High/Low expression groups defined by ``cutoff``.

    Returns a dict with the fitted ``KaplanMeierFitter`` objects, log-rank
    test result, hazard ratio and group median survival times.
    """
    df = df.copy()
    df["group"] = np.where(df[value_col] >= cutoff, "High", "Low")

    fitters = {}
    for label in ("Low", "High"):
        kmf = KaplanMeierFitter()
        subset = df[df["group"] == label]
        kmf.fit(durations=subset["OS_time"], event_observed=subset["OS_event"], label=label)
        fitters[label] = kmf

    logrank_result = logrank_test(
        df.loc[df["group"] == "Low", "OS_time"],
        df.loc[df["group"] == "High", "OS_time"],
        event_observed_A=df.loc[df["group"] == "Low", "OS_event"],
        event_observed_B=df.loc[df["group"] == "High", "OS_event"],
    )

    n_high = (df["group"] == "High").sum()
    n_low = (df["group"] == "Low").sum()
    hazard_ratio = (df.loc[df["group"] == "High", "OS_event"].sum() / n_high) / (
        df.loc[df["group"] == "Low", "OS_event"].sum() / n_low
    )

    return {
        "fitters": fitters,
        "logrank_p": logrank_result.p_value,
        "hazard_ratio": hazard_ratio,
        "median_high": df.loc[df["group"] == "High", "OS_time"].median(),
        "median_low": df.loc[df["group"] == "Low", "OS_time"].median(),
        "n_high": n_high,
        "n_low": n_low,
    }


def prepare_cox_table(
    mirna_expr: pd.DataFrame,
    clinical: pd.DataFrame,
    mirna_id: str = TARGET_MIRNA_ID,
) -> pd.DataFrame:
    """Merge expression with clinical covariates for a multivariate Cox model."""
    clinical = clinical.replace("--", np.nan)

    expr = mirna_expr.loc[mirna_expr["miRNA_ID"] == mirna_id, ["case_id", "reads_per_million_miRNA_mapped"]].copy()
    expr["log2_rpm"] = np.log2(expr["reads_per_million_miRNA_mapped"] + 0.01)

    clin = clinical[
        [
            "case_id",
            "days_to_death",
            "days_to_last_follow_up",
            "vital_status",
            "age_at_diagnosis",
            "gender",
            "ajcc_pathologic_stage",
        ]
    ].drop_duplicates("case_id").copy()

    clin["days_to_death"] = pd.to_numeric(clin["days_to_death"], errors="coerce")
    clin["days_to_last_follow_up"] = pd.to_numeric(clin["days_to_last_follow_up"], errors="coerce")
    clin["age_at_diagnosis"] = pd.to_numeric(clin["age_at_diagnosis"], errors="coerce")
    clin["OS_time"] = clin[["days_to_death", "days_to_last_follow_up"]].min(axis=1) / 30.0
    clin["vital_status"] = clin["vital_status"].str.lower().str.strip()
    clin["OS_event"] = clin["vital_status"].map({"dead": 1, "alive": 0})

    df = pd.merge(clin, expr[["case_id", "log2_rpm"]], on="case_id", how="inner")
    df = df.dropna(
        subset=["OS_time", "OS_event", "log2_rpm", "age_at_diagnosis", "gender", "ajcc_pathologic_stage"]
    )

    df["gender"] = df["gender"].str.lower().map({"male": 1, "female": 0})
    df["stage_num"] = df["ajcc_pathologic_stage"].str.extract(r"Stage\s*([IVX]+)", expand=False).map(STAGE_TO_NUMERIC)
    df = df.dropna(subset=["stage_num"])

    return df[["OS_time", "OS_event", "log2_rpm", "age_at_diagnosis", "gender", "stage_num"]]


def fit_cox_model(cox_table: pd.DataFrame) -> CoxPHFitter:
    """Fit a multivariate Cox proportional-hazards model."""
    cph = CoxPHFitter()
    cph.fit(cox_table, duration_col="OS_time", event_col="OS_event")
    return cph


def bootstrap_concordance(
    cox_table: pd.DataFrame,
    n_bootstrap: int = N_BOOTSTRAP,
    min_events_per_class: int = 5,
) -> list[float]:
    """Internally validate a Cox model via bootstrap resampling of the C-index."""
    c_indices = []
    for i in range(n_bootstrap):
        sample = cox_table.sample(frac=1, replace=True, random_state=i)
        event_counts = sample["OS_event"].value_counts()
        if 0 not in event_counts or 1 not in event_counts:
            continue
        if event_counts[0] < min_events_per_class or event_counts[1] < min_events_per_class:
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(sample, duration_col="OS_time", event_col="OS_event")
            predicted_risk = cph.predict_partial_hazard(sample)
            c_indices.append(concordance_index(sample["OS_time"], -predicted_risk, sample["OS_event"]))
        except Exception:
            continue
    return c_indices
