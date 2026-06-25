# Growth Rate Analysis

R script for calculating the maximum specific growth rate (μ) and doubling time of bacterial cultures from OD measurements in microplate format. For each biological replicate, the script identifies the optimal exponential growth window via exhaustive linear regression on log-transformed OD values and exports the results to Excel.

---

## Key script

| Script | Description |
|--------|-------------|
| `growth_rate_analysis.R` | Reads a wide-format OD table from Excel, reshapes it, identifies the best-fit exponential interval per replicate, and writes results to `.xlsx` |

---

## Workflow overview

```
Excel file (wide format: Time | Sample_1 | Sample_2 | ...)
       │
       ▼
Reshape to long format
       │  pivot_longer → Strain + Rep + lnOD
       ▼
Exponential window search   (per replicate)
       │  filter OD ∈ [od_min, od_max]
       │  exhaustive start/end index pairs
       │  lm(lnOD ~ Time) for each window
       │  rank by adj-R², n_points, μ
       ▼
Best window selected
       │  μ (h⁻¹), doubling time (h and min), R², fold-change
       ▼
growth_results_<sheet>.xlsx
```

---

## Usage

```bash
Rscript growth_rate_analysis.R <excel_path> <sheet> [od_min] [od_max] [min_points] [min_fold_change] [output_dir]
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `excel_path` | ✓ | — | Path to the input Excel file |
| `sheet` | ✓ | — | Sheet name to read (e.g. `"1B"`) |
| `od_min` | — | `0.15` | Lower OD bound for the exponential window |
| `od_max` | — | `1.0` | Upper OD bound for the exponential window |
| `min_points` | — | `5` | Minimum number of time points per window |
| `min_fold_change` | — | `1.5` | Minimum OD fold-change across the window |
| `output_dir` | — | Desktop | Directory for the output Excel file |

### Examples

```bash
# Minimal — use all defaults
Rscript growth_rate_analysis.R data/experiment.xlsx "1B"

# Custom OD window and output folder
Rscript growth_rate_analysis.R data/experiment.xlsx "1B" 0.1 0.9 4 1.5 results/
```

---

## Input format

The Excel sheet must have a `Time` column (in hours) followed by one column per replicate well. Column names should follow the pattern `<Strain>_<Rep>` (e.g. `WT_1`, `WT_2`, `mutA_1`).

```
Time | WT_1  | WT_2  | mutA_1 | mutA_2
0    | 0.05  | 0.06  | 0.05   | 0.05
1    | 0.07  | 0.08  | 0.06   | 0.07
...
```

---

## Output

Single Excel file `growth_results_<sheet>.xlsx` with one row per replicate:

| Column | Description |
|--------|-------------|
| `Strain` | Strain name (prefix before `_<number>`) |
| `Rep` | Replicate number |
| `t_start` / `t_end` | Time bounds of the exponential window (h) |
| `n_points` | Number of time points in the window |
| `od_start` / `od_end` | OD at the window boundaries |
| `fold_change` | OD fold-change across the window |
| `mu` | Specific growth rate μ (h⁻¹) |
| `doubling_time_h` | Doubling time in hours |
| `doubling_time_min` | Doubling time in minutes |
| `r2` / `adj_r2` | R² and adjusted R² of the linear fit |

---

## Dependencies

```r
install.packages(c("readxl", "dplyr", "tidyr", "purrr", "broom", "writexl", "stringr"))
```
