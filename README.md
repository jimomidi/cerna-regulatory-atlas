# ceRNA Regulatory Atlas

A collection of scripts (Python, plus one R script) for analyzing a
candidate miRNA's role in a competing endogenous RNA (ceRNA) regulatory
network: expression correlation with target genes, network centrality
analysis, survival association, time-dependent predictive performance,
clinical subgroup comparison, drug-target druggability summaries, and the
corresponding publication-style figures.

The toolkit is data-agnostic: it does not ship with, or assume, any
particular disease, dataset, or gene/miRNA identity. Point it at your own
expression matrices, clinical tables, and interaction/annotation files and
configure the target miRNA and gene list in `config.py`.

## Workflow

The pipeline mirrors a typical miRNA-ceRNA characterization study:

1. **Correlation analysis** (`analysis/mirna_target_correlation.py`) --
   test whether a candidate miRNA is negatively correlated with predicted
   or differentially-expressed target genes across tumor samples.
2. **Target overlap** (`analysis/target_gene_overlap.py`) -- quantify
   agreement between correlation-derived targets and independent gene
   lists (e.g. differential expression, alternative prediction tools).
3. **Network centrality** (`analysis/network_centrality.py`) -- build a
   ceRNA interaction network from a curated interaction table filtered by
   expression correlation, and rank nodes by betweenness/closeness
   centrality and degree.
4. **Survival analysis** (`analysis/survival_analysis.py`) -- Kaplan-Meier
   analysis with a data-driven optimal expression cutoff, and a
   multivariate Cox proportional-hazards model with bootstrap internal
   validation.
5. **Time-dependent predictive performance** (`analysis/time_dependent_auc.R`)
   -- a spline-adjusted Cox risk score evaluated as a time-dependent ROC
   AUC at several follow-up horizons, with bootstrap 95% confidence bands.
   Written in R (using `timeROC`); its CSV output matches what
   `figures/plot_time_dependent_auc.py` expects, so the two can be chained
   across languages.
6. **Clinical association** (`analysis/clinical_association.py`) --
   compare expression across age, gender and pathological-stage subgroups,
   and between normal and tumor samples.
7. **Treatment group comparison** (`analysis/treatment_group_comparison.py`)
   -- compare expression between clinically defined treatment groups.
8. **Drug-gene interaction summary** (`analysis/drug_gene_interaction.py`)
   -- summarize druggability and regulatory-approval status of network
   target genes.

Each `analysis/*.py` module returns plain DataFrames/objects; the matching
script in `figures/` turns those results into a saved plot. The one R
script (`analysis/time_dependent_auc.R`) instead writes its results to a
CSV file, which `figures/plot_time_dependent_auc.py` can then read and
plot. Keeping analysis and plotting separate means you can rerun
statistics without regenerating figures, or restyle figures without
recomputing statistics.

## Requirements

* Python 3.10+
* See `requirements.txt` for package versions. Install with:

```bash
pip install -r requirements.txt
```

`ptitprince` (used only by the raincloud plot in
`figures/plot_treatment_group_comparison.py`) can be tricky to install on
some platforms; if you don't need that specific figure, remove it from
`requirements.txt`.

`analysis/time_dependent_auc.R` additionally requires R (4.0+) with the
`survival`, `timeROC`, `ggplot2`, `dplyr` and `splines` packages. The
script installs any of these that are missing on first run; to install
them manually instead:

```r
install.packages(c("survival", "timeROC", "ggplot2", "dplyr", "splines"))
```

## Project layout

```
cerna-regulatory-atlas/
├── config.py           # Central paths and analysis parameters
├── utils.py             # Shared plotting/IO helpers
├── analysis/             # Statistics and data processing (no plotting)
│   ├── mirna_target_correlation.py
│   ├── target_gene_overlap.py
│   ├── network_centrality.py
│   ├── survival_analysis.py
│   ├── time_dependent_auc.R
│   ├── clinical_association.py
│   ├── treatment_group_comparison.py
│   └── drug_gene_interaction.py
├── figures/               # Figure-generation scripts (take analysis output, save plots)
│   ├── plot_pan_cancer_expression.py
│   ├── plot_mirna_logfc_by_cancer.py
│   ├── plot_target_gene_expression_summary.py
│   ├── plot_target_overlap_summary.py
│   ├── plot_network_graph.py
│   ├── plot_centrality_zscore_boxplot.py
│   ├── plot_closeness_centrality_histogram.py
│   ├── plot_degree_distribution.py
│   ├── plot_survival_results.py
│   ├── plot_clinical_subgroups.py
│   ├── plot_treatment_group_comparison.py
│   ├── plot_drug_interaction_summary.py
│   └── plot_functional_enrichment.py
├── requirements.txt
├── LICENSE
└── .gitignore
```

`data/` and `results/` are not included in the repository (see
`.gitignore`); create them locally, or point the scripts at wherever your
own files live.

## Expected input format

None of the scripts ship with example data. Each module's docstring
documents its exact expected columns; in summary:

| File | Key columns |
|---|---|
| miRNA expression | `miRNA_ID`, `case_id`, `reads_per_million_miRNA_mapped` |
| mRNA expression (long format) | `gene_name`, `case_id`, `fpkm_unstranded` / `tpm_unstranded` |
| mRNA expression (wide format, e.g. GTEx-style) | a gene-name column, one column per sample |
| Clinical | `case_id`, `days_to_death`, `days_to_last_follow_up`, `vital_status`, `age_at_diagnosis`/`age_at_index`, `gender`, `ajcc_pathologic_stage` |
| Gene list | `gene` |
| Interaction table | columns renameable to `miRNA_ID`, `gene_id`, `weight` (`Var1`/`Var2`/`Freq` handled automatically) |
| Drug interaction table | `gene`, `drug`, `indication`, `regulatory approval`, `interaction score` |
| Treatment groups | `case_id`, `treatment_group` |
| Functional enrichment | `Term`, `Count`, `List Total`, `FDR` (+ `Category` for GO) |

## Usage

Run modules from the repository root using `python -m` so relative
imports (`config`, `utils`, `analysis.*`) resolve correctly. For example:

```bash
python -m analysis.mirna_target_correlation \
    --mirna-expression data/tumor_mirna.csv \
    --mrna-expression data/tumor_mrna.csv \
    --gene-list data/candidate_targets.csv \
    --output results/correlation_results.csv \
    --mirna-id hsa-miR-XXX-3p

python -m figures.plot_pan_cancer_expression \
    --input data/pan_cancer_expression.csv \
    --output results/pan_cancer_expression.png \
    --exclude-project PROJECT_X
```

`analysis/time_dependent_auc.R` is a standalone R script, run separately
from the Python modules. It takes four positional arguments (miRNA
expression CSV, clinical CSV, target miRNA ID, output path prefix) and
writes both a CSV and a PNG; the CSV can then be restyled with the
Python figure script if desired:

```bash
Rscript analysis/time_dependent_auc.R \
    data/tumor_mirna.csv data/clinical.csv target_mirna \
    results/time_dependent_auc

python -m figures.plot_time_dependent_auc \
    --input results/time_dependent_auc.csv \
    --output results/time_dependent_auc.png
```

Scripts that don't expose a command-line interface (most of the
`network_centrality`, `survival_analysis` and `plot_*` figure modules,
which take DataFrames/graphs rather than single file paths) are meant to
be imported into a short driver script or notebook of your own, e.g.:

```python
import pandas as pd
from analysis.survival_analysis import prepare_survival_table, find_optimal_cutoff, fit_kaplan_meier
from figures.plot_survival_results import plot_kaplan_meier

mirna_expr = pd.read_csv("data/tumor_mirna.csv")
clinical = pd.read_csv("data/clinical.csv")

df = prepare_survival_table(mirna_expr, clinical, mirna_id="hsa-miR-XXX-3p")
cutoff = find_optimal_cutoff(df)
result = fit_kaplan_meier(df, cutoff)
plot_kaplan_meier(result, "results/survival_curve.png")
```

## Outputs

Each script writes CSV summary tables and/or PNG/SVG figures to whatever
path you pass via `--output` / `output_path`. No output paths are
hard-coded; nothing is written outside the location you specify.

## Reproducibility

* All file paths are supplied by the caller (CLI arguments or function
  parameters) -- nothing is hard-coded to a particular machine.
* Statistical parameters (correlation/FDR thresholds, random seed,
  bootstrap iterations, minimum sample fractions, etc.) are centralized in
  `config.py`.
* The random seed used for permutation-based centrality z-scores
  (`config.RANDOM_SEED`) is fixed for reproducibility; the Cox bootstrap
  validation uses a deterministic per-iteration seed (`random_state=i`).
* `analysis/time_dependent_auc.R` is a standalone R script (not tied to
  `config.py`); its random seed, bootstrap count and evaluation horizons
  are set in its own "Configuration" section at the top of the file.

## Customization

Edit `config.py` to set:

* `TARGET_MIRNA_ID` -- the miRNA identifier used across the correlation,
  survival, clinical-association and treatment-comparison modules.
* `TARGET_GENE_LIST` -- default candidate target genes.
* Statistical thresholds (`CORRELATION_THRESHOLD`, `FDR_THRESHOLD`,
  `MIN_INTERACTION_WEIGHT`, etc.).
* Plot styling (`FONT_FAMILY`, `FIGURE_DPI`).

## Design notes

* **No manual edits to results.** The interaction network is built only
  from the interaction table and the sign/significance of the paired
  expression correlation. No node or edge is added, removed, or relabeled
  by hand, and hubs are ranked purely from the computed centrality values.
  `figures/plot_network_graph.py` always labels a node with its own
  identifier and has no parameter for substituting another node's label.
* **No hard-coded results.** Every plotted value is either computed by
  the analysis modules or read from a caller-supplied file. The `--demo`
  option of `figures/plot_target_overlap_summary.py` uses arbitrary
  synthetic numbers purely to smoke-test the plotting code.
* **Plot styling is separate from statistics.** Shared box-plot styles and
  axis limits (e.g. in `figures/plot_clinical_subgroups.py`) affect
  appearance only.
* **Batch-effect correction is out of scope.** No script here performs
  batch correction; any such step belongs upstream, on the files these
  scripts read. `analysis/clinical_association.prepare_group_table` and
  `figures/plot_clinical_subgroups.plot_group_comparison` default to a
  value column named `log_RPM_adj`, so they can consume an already
  adjusted expression column if your own pipeline produces one.
