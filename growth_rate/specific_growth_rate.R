#!/usr/bin/env Rscript
# ==============================================================================
# Bacterial Growth Rate Analysis
# ------------------------------------------------------------------------------
# Calculates the maximum specific growth rate (mu) and doubling time for each
# biological replicate by identifying the optimal exponential growth interval
# via linear regression on log-transformed OD values.
#
# Usage:
#   Rscript growth_rate_analysis.R <excel_path> <sheet> [od_min] [od_max] [min_points] [min_fold_change] [output_dir]
#
# Required arguments:
#   excel_path      Path to the input Excel file
#   sheet           Sheet name to read (e.g. "1B")
#
# Optional arguments (positional, provide in order):
#   od_min          Minimum OD to include in exponential window  (default: 0.15)
#   od_max          Maximum OD to include in exponential window  (default: 1.0)
#   min_points      Minimum number of time points per interval   (default: 5)
#   min_fold_change Minimum OD fold-change across interval       (default: 1.5)
#   output_dir      Directory for the output Excel file          (default: Desktop)
#
# Example:
#   Rscript growth_rate_analysis.R data/experiment.xlsx "1B" 0.1 0.9 4 1.5 results/
# ==============================================================================

# ── 1. Parse command-line arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript growth_rate_analysis.R <excel_path> <sheet>",
      "[od_min] [od_max] [min_points] [min_fold_change] [output_dir]\n")
  quit(status = 1)
}

excel_path     <- args[1]
sheet          <- args[2]
od_min         <- if (length(args) >= 3) as.numeric(args[3]) else 0.15
od_max         <- if (length(args) >= 4) as.numeric(args[4]) else 1.0
min_points     <- if (length(args) >= 5) as.integer(args[5]) else 5L
min_fold_change <- if (length(args) >= 6) as.numeric(args[6]) else 1.5
output_dir     <- if (length(args) >= 7) args[7] else file.path(
  Sys.getenv("USERPROFILE", unset = "~"), "Desktop"
)

cat("─────────────────────────────────────────────────────\n")
cat("Input file  :", excel_path, "\n")
cat("Sheet       :", sheet,      "\n")
cat("OD window   :", od_min, "–", od_max, "\n")
cat("Min points  :", min_points, "\n")
cat("Min FC      :", min_fold_change, "\n")
cat("Output dir  :", output_dir, "\n")
cat("─────────────────────────────────────────────────────\n")

# ── 2. Load libraries ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  library(writexl)
  library(stringr)
})

# ── 3. Import data ────────────────────────────────────────────────────────────
if (!file.exists(excel_path)) {
  stop("Excel file not found: ", excel_path)
}

raw_data <- as.data.frame(read_excel(excel_path, sheet = sheet))
cat("Rows read from sheet:", nrow(raw_data), "\n")

# ── 4. Reshape to long format and derive helper columns ───────────────────────
# Each column (except Time) represents one replicate well; pivot to tidy format.
data_long <- raw_data %>%
  pivot_longer(
    cols      = -Time,
    names_to  = "Sample",
    values_to = "OD"
  ) %>%
  filter(!is.na(OD), OD > 0) %>%
  mutate(
    Strain = str_remove(Sample, "_[0-9]+$"),   # strip trailing _<number>
    Rep    = str_extract(Sample, "[0-9]+$"),    # keep only the replicate number
    lnOD   = log(OD)                            # natural log for linear regression
  )

cat("Samples detected:", n_distinct(data_long$Sample), "\n")
cat("Strains detected:", n_distinct(data_long$Strain), "\n")

# ── 5. Function: find the optimal exponential-growth interval ─────────────────
# Strategy: exhaustive search over all contiguous sub-windows within the
# OD range [od_min, od_max]. Each window is scored by adj-R², then by number
# of points (longer windows preferred at equal fit), then by mu. The top-
# ranked window is returned for downstream summarisation.
find_interval <- function(df,
                          min_pts      = min_points,
                          od_lo        = od_min,
                          od_hi        = od_max,
                          min_fc       = min_fold_change) {

  df <- df %>%
    arrange(Time) %>%
    filter(OD >= od_lo, OD <= od_hi)

  n <- nrow(df)
  if (n < min_pts) return(NULL)   # not enough points in the valid OD range

  # Iterate over all start (i) / end (j) index pairs
  intervals <- map_dfr(seq_len(n - min_pts + 1L), function(i) {
    map_dfr(seq(i + min_pts - 1L, n), function(j) {

      sub   <- df[i:j, ]
      model <- lm(lnOD ~ Time, data = sub)

      mu  <- coef(model)[["Time"]]
      r2  <- summary(model)$r.squared
      ar2 <- summary(model)$adj.r.squared
      fc  <- max(sub$OD) / min(sub$OD)

      tibble(
        t_start           = min(sub$Time),
        t_end             = max(sub$Time),
        n_points          = nrow(sub),
        od_start          = first(sub$OD),
        od_end            = last(sub$OD),
        fold_change       = fc,
        mu                = mu,            # h⁻¹  (specific growth rate)
        doubling_time_h   = log(2) / mu,   # hours
        doubling_time_min = log(2) / mu * 60,
        r2                = r2,
        adj_r2            = ar2
      )
    })
  }) %>%
    filter(mu > 0, fold_change >= min_fc) %>%
    arrange(desc(adj_r2), desc(n_points), desc(mu))

  intervals %>% slice(1)   # return the single best-fitting window
}

# ── 6. Apply to every biological replicate ────────────────────────────────────
results_replicates <- data_long %>%
  group_by(Strain, Rep) %>%
  group_modify(~ find_interval(.x)) %>%
  ungroup()

cat("Replicates processed:", nrow(results_replicates), "\n")

# ── 7. Export results ─────────────────────────────────────────────────────────
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

output_file <- file.path(output_dir,
                         paste0("growth_results_", sheet, ".xlsx"))

write_xlsx(results_replicates, output_file)
cat("Results saved to:", output_file, "\n")
