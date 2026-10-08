"""
Compare target miRNA expression between two treatment groups (e.g. patients
who received chemotherapy vs. radiotherapy), using group membership
supplied by the caller rather than hard-coded sample identifiers.

Expected inputs
----------------
mirna_expression : csv
    Long-format table with ``miRNA_ID``, ``case_id`` and
    ``reads_per_million_miRNA_mapped`` columns.
treatment_groups : csv
    Table with ``case_id`` and ``treatment_group`` columns, mapping each
    sample to its treatment category (e.g. "Pharmaceutical", "Radiation").
    Build this from your own clinical/treatment records; sample
    identifiers are not shipped with this repository.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from config import TARGET_MIRNA_ID


def prepare_treatment_comparison_table(
    mirna_expr: pd.DataFrame,
    treatment_groups: pd.DataFrame,
    mirna_ids: Sequence[str] = (TARGET_MIRNA_ID,),
) -> pd.DataFrame:
    """Build a long-format table of (miRNA, treatment_group, expression, log_expression)
    for the requested miRNA(s), restricted to samples with a known treatment group
    and non-zero expression.
    """
    pivot = mirna_expr.pivot_table(index="case_id", columns="miRNA_ID", values="reads_per_million_miRNA_mapped")

    records = []
    for mirna in mirna_ids:
        if mirna not in pivot.columns:
            continue
        for _, row in treatment_groups.iterrows():
            case_id, group = row["case_id"], row["treatment_group"]
            if case_id not in pivot.index:
                continue
            expr = pivot.loc[case_id, mirna]
            if pd.notna(expr) and expr > 0:
                records.append({"miRNA": mirna, "treatment_group": group, "expression": expr})

    df = pd.DataFrame(records)
    if not df.empty:
        df["log_expression"] = np.log10(df["expression"])
    return df
