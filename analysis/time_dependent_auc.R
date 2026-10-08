#!/usr/bin/env Rscript
#
# Time-dependent ROC AUC for a target miRNA's expression, with a
# spline-adjusted Cox risk score and bootstrap 95% confidence bands.
#
# For a series of follow-up horizons, this script fits a Cox proportional
# hazards model (miRNA expression via a natural spline, plus age and
# gender), derives a linear-predictor risk score, and computes the
# time-dependent AUC of that risk score against the observed outcome using
# timeROC::timeROC(). Confidence bands come from bootstrap resampling of
# the full pipeline (refit + re-score + re-estimate AUC per resample),
# which remains valid even when the analytic (iid-based) standard error is
# unavailable.
#
# This script only performs the statistical computation and exports its
# results; the companion Python script `figures/plot_time_dependent_auc.py`
# renders the resulting CSV as a styled figure consistent with the rest of
# this repository. A local plot is also produced here for convenience.
#
# Expected inputs
# ----------------
# mirna_expression : csv
#   Long-format table with `miRNA_ID`, `case_id` and
#   `reads_per_million_miRNA_mapped` columns.
# clinical : csv
#   Table with `case_id`, `days_to_death`, `days_to_last_follow_up`,
#   `vital_status`, `age_at_diagnosis` and `gender` columns.
#
# Output
# ------
# A CSV with columns `time`, `auc`, `ci_low`, `ci_high` (matching the input
# format expected by `figures/plot_time_dependent_auc.py`), and a PNG plot.
#
# Usage
# -----
#   Rscript analysis/time_dependent_auc.R \
#     data/tumor_mirna.csv data/clinical.csv target_mirna results/time_dependent_auc

suppressWarnings(suppressMessages({

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
required_packages <- c("survival", "timeROC", "ggplot2", "dplyr", "splines")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

}))

# ---------------------------------------------------------------------------
# Configuration (override via command-line arguments; see Usage above)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

mirna_expression_path <- if (length(args) >= 1) args[[1]] else "data/tumor_mirna.csv"
clinical_path         <- if (length(args) >= 2) args[[2]] else "data/clinical.csv"
target_mirna_id       <- if (length(args) >= 3) args[[3]] else "target_mirna"
output_prefix         <- if (length(args) >= 4) args[[4]] else "results/time_dependent_auc"

evaluation_times <- seq(12, 108, by = 12)  # follow-up horizons, in months
spline_df <- 5                              # degrees of freedom for the expression spline
n_bootstrap <- 500                          # lower (e.g. 200) for speed, or raise (e.g. 1000) for tighter CIs
random_seed <- 42

# ---------------------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------------------
expr_all <- read.csv(mirna_expression_path)

expr <- expr_all |>
  dplyr::filter(miRNA_ID == target_mirna_id) |>
  dplyr::select(case_id, reads_per_million_miRNA_mapped) |>
  dplyr::mutate(log2_rpm = log2(reads_per_million_miRNA_mapped + 0.01))

clinical <- read.csv(clinical_path)
clinical[clinical == "--"] <- NA

clinical <- clinical |>
  dplyr::mutate(
    days_to_death          = suppressWarnings(as.numeric(days_to_death)),
    days_to_last_follow_up = suppressWarnings(as.numeric(days_to_last_follow_up)),
    OS_time = pmin(days_to_death, days_to_last_follow_up, na.rm = TRUE) / 30,
    vital_status = tolower(trimws(vital_status)),
    OS_event = ifelse(vital_status == "dead", 1, 0),
    age = suppressWarnings(as.numeric(age_at_diagnosis)),
    gender = as.factor(trimws(gender))
  )

survival_data <- dplyr::inner_join(clinical, expr, by = "case_id") |>
  dplyr::select(OS_time, OS_event, log2_rpm, age, gender) |>
  dplyr::filter(
    is.finite(OS_time),
    is.finite(OS_event),
    is.finite(log2_rpm),
    OS_time >= 1
  ) |>
  dplyr::filter(!is.na(age), !is.na(gender))

# Remove extreme expression outliers (top 1%)
upper_threshold <- as.numeric(stats::quantile(survival_data$log2_rpm, 0.99, na.rm = TRUE))
survival_data <- survival_data |> dplyr::filter(log2_rpm <= upper_threshold)

# ---------------------------------------------------------------------------
# Core computation: fit Cox model -> risk score -> time-dependent AUC
# ---------------------------------------------------------------------------
compute_auc_vector <- function(data, times, spline_df = 5) {
  fit <- survival::coxph(
    survival::Surv(OS_time, OS_event) ~ splines::ns(log2_rpm, df = spline_df) + age + gender,
    data = data,
    x = TRUE
  )
  data$risk_score <- stats::predict(fit, type = "lp")

  roc <- timeROC::timeROC(
    T = data$OS_time,
    delta = data$OS_event,
    marker = data$risk_score,
    cause = 1,
    times = times,
    iid = FALSE
  )
  as.numeric(roc$AUC)
}

# Point estimate
auc_point <- compute_auc_vector(survival_data, evaluation_times, spline_df = spline_df)

# Bootstrap confidence interval (full pipeline re-run per resample)
set.seed(random_seed)
n_samples <- nrow(survival_data)
bootstrap_auc <- matrix(NA_real_, nrow = n_bootstrap, ncol = length(evaluation_times))

for (b in seq_len(n_bootstrap)) {
  resample_idx <- sample.int(n_samples, size = n_samples, replace = TRUE)
  resampled_data <- survival_data[resample_idx, , drop = FALSE]

  # A resample can occasionally collapse to a single event class or gender
  # level; treat such failures as missing rather than aborting the run.
  bootstrap_auc[b, ] <- tryCatch(
    compute_auc_vector(resampled_data, evaluation_times, spline_df = spline_df),
    error = function(e) rep(NA_real_, length(evaluation_times))
  )

  if (b %% 50 == 0) message("Bootstrap: ", b, "/", n_bootstrap)
}

successful_runs <- apply(bootstrap_auc, 1, function(x) all(is.finite(x)))
bootstrap_auc_ok <- bootstrap_auc[successful_runs, , drop = FALSE]

if (nrow(bootstrap_auc_ok) < 30) {
  warning("Fewer than 30 successful bootstrap runs. Consider relaxing filters or reducing spline_df.")
}

ci_low  <- apply(bootstrap_auc_ok, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
ci_high <- apply(bootstrap_auc_ok, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
results <- data.frame(
  time    = evaluation_times,
  auc     = auc_point,
  ci_low  = pmax(0, ci_low),
  ci_high = pmin(1, ci_high)
)

output_dir <- dirname(output_prefix)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

csv_path <- paste0(output_prefix, ".csv")
write.csv(results, csv_path, row.names = FALSE)

plot <- ggplot2::ggplot(results, ggplot2::aes(x = time, y = auc)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                        fill = "gray80", alpha = 0.4) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
  ggplot2::coord_cartesian(ylim = c(0.4, 1.0)) +
  ggplot2::labs(
    title = paste0("Time-dependent ROC AUC (Spline Cox Model, Bootstrap 95% CI; B=", n_bootstrap, ")"),
    x = "Time (months)",
    y = "Time-dependent AUC"
  ) +
  ggplot2::theme_classic(base_family = "serif")

png_path <- paste0(output_prefix, ".png")
ggplot2::ggsave(png_path, plot = plot, width = 6, height = 4, dpi = 300)

print(round(results, 3))
cat("\nSaved results table to:", csv_path, "\n")
cat("Saved figure to:", png_path, "\n")
cat("Bootstrap successful runs:", nrow(bootstrap_auc_ok), "out of", n_bootstrap, "\n")
