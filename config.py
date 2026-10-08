"""
Central configuration for the ceRNA Regulatory Atlas toolkit.

Every script in this repository imports its default paths and analysis
parameters from here. Edit the values below (or override them from the
command line, where supported) to point the pipeline at your own data
without touching the analysis code itself.

All paths are relative to the repository root so the project remains
portable across machines and operating systems.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# ---------------------------------------------------------------------------
# Study identifiers (replace with your own miRNA / gene of interest)
# ---------------------------------------------------------------------------
# The pipeline focuses on a single candidate miRNA. Set this to the mature
# miRNA identifier used in your own expression matrices
# (e.g. "hsa-miR-XXX-3p").
TARGET_MIRNA_ID = "target_mirna"

# Candidate target genes to evaluate (e.g. predicted targets from a
# target-prediction database, or a differential expression gene list).
# Replace with your own gene symbols.
TARGET_GENE_LIST = ["GENE_1", "GENE_2", "GENE_3"]

# ---------------------------------------------------------------------------
# Statistical parameters
# ---------------------------------------------------------------------------
RANDOM_SEED = 42

# Minimum fraction of samples in which a feature must exceed a low-expression
# threshold to be retained during filtering.
MIN_EXPRESSED_FRACTION = 0.2

# Pearson correlation coefficient and FDR (Benjamini-Hochberg) cutoffs used
# to call a miRNA-gene pair a significant negative correlation.
CORRELATION_THRESHOLD = -0.3
FDR_THRESHOLD = 0.05

# Minimum number of paired samples required to compute a correlation.
MIN_PAIRED_SAMPLES = 5

# Minimum interaction-database support (e.g. number of prediction tools or
# evidence count) required to keep a miRNA-gene edge in the ceRNA network.
MIN_INTERACTION_WEIGHT = 5

# Number of randomized graphs used to build the null distribution for
# centrality z-scores.
N_RANDOM_GRAPHS = 100

# Number of bootstrap resamples used for Cox model internal validation.
N_BOOTSTRAP = 1000

# ---------------------------------------------------------------------------
# Plot styling
# ---------------------------------------------------------------------------
FONT_FAMILY = "Times New Roman"
FIGURE_DPI = 300
