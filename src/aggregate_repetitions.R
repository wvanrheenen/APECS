# src/aggregate_repetitions.R
# Aggregate prediction_model.csv AND combined_metrics.csv (ONLY "overall" row) across repetitions
source("src/libraries_simPed.R")

# Get input files
pred_files <- snakemake@input[["prediction_models"]]
metrics_files <- snakemake@input[["metrics_files"]]
set_name <- snakemake@params[["set"]]

cat("Aggregating for set:", set_name, "\n")
cat("- Prediction files:", length(pred_files), "\n")
cat("- Metrics files:", length(metrics_files), "\n\n")

# Safe read function
safe_read_csv <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    cat("Skipping empty/missing file:", basename(file), "\n")
    return(NULL)
  }
  tryCatch({
    df <- read.csv(file, stringsAsFactors = FALSE)
    if (nrow(df) == 0) {
      cat("Skipping empty data file:", basename(file), "\n")
      return(NULL)
    }
    cat("Read", nrow(df), "rows from:", basename(file), "\n")
    df
  }, error = function(e) {
    cat("Failed to read", basename(file), ":", e$message, "\n")
    NULL
  })
}

# === AGGREGATE PREDICTION MODELS (unchanged) ===
models_list <- lapply(pred_files, safe_read_csv)
models_list <- models_list[!sapply(models_list, is.null)]
if (length(models_list) > 0) {
  all_models <- rbindlist(models_list, idcol = "rep", fill = TRUE)
  aggregated_predictions <- all_models[, .(
    mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
    sd_Sensitivity = sd(Sensitivity, na.rm = TRUE),
    mean_Specificity = mean(Specificity, na.rm = TRUE),
    sd_Specificity = sd(Specificity, na.rm = TRUE),
    mean_PPV = mean(PPV, na.rm = TRUE),
    sd_PPV = sd(PPV, na.rm = TRUE),
    mean_NPV = mean(NPV, na.rm = TRUE),
    sd_NPV = sd(NPV, na.rm = TRUE),
    mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    sd_Accuracy = sd(Accuracy, na.rm = TRUE),
    mean_Odds_Ratio = mean(Odds_Ratio, na.rm = TRUE),
    sd_Odds_Ratio = sd(Odds_Ratio, na.rm = TRUE),
    n_reps = .N
  ), by = Scenario]
  write.csv(aggregated_predictions, snakemake@output[["prediction_model"]], row.names = FALSE)
  cat("✓ Wrote aggregated predictions:", nrow(aggregated_predictions), "scenarios\n")
} else {
  cat("⚠ No valid prediction files found\n")
}

# === AGGREGATE METRICS (ONLY "overall" row, preserve column order) ===
metrics_list <- lapply(metrics_files, safe_read_csv)
metrics_list <- metrics_list[!sapply(metrics_list, is.null)]

if (length(metrics_list) == 0) {
  stop("No valid combined_metrics files found")
}

# Filter to ONLY "overall" rows
overall_metrics <- lapply(metrics_list, function(df) {
  if ("yob_interval" %in% names(df)) {
    overall_row <- df[df$yob_interval == "overall", ]
  } else {
    overall_row <- df[1, ]  # fallback
  }
  if (nrow(overall_row) == 0) return(NULL)
  overall_row
})

overall_metrics <- overall_metrics[!sapply(overall_metrics, is.null)]
if (length(overall_metrics) == 0) {
  stop("No 'overall' rows found in metrics files")
}

# Combine all overall rows (use FIRST file's column order as reference)
all_overall_metrics <- rbindlist(overall_metrics, idcol = "rep", fill = TRUE)
cat("Combined", nrow(all_overall_metrics), "'overall' rows from", length(unique(all_overall_metrics$rep)), "repetitions\n")

# Get column order from FIRST metrics file (preserve original order)
reference_cols <- names(metrics_list[[1]])
numeric_cols <- reference_cols[sapply(reference_cols, function(col) is.numeric(all_overall_metrics[[col]]))]

cat("Preserving column order from reference:", paste(head(numeric_cols, 3), collapse = ", "), "...\n")

# Create LONG FORMAT output: one row per original column
final_metrics <- data.table(
  metric = character(),
  mean_value = numeric(),
  sd_value = numeric(),
  min_value = numeric(),
  max_value = numeric(),
  n_reps = integer()
)

for (col in numeric_cols) {
  col_values <- all_overall_metrics[[col]]
  mean_val <- round(mean(col_values, na.rm = TRUE), 5)
  sd_val <- round(sd(col_values, na.rm = TRUE), 5)
  min_val <- round(min(col_values, na.rm = TRUE), 5)
  max_val <- round(max(col_values, na.rm = TRUE), 5)
  
  final_metrics <- rbind(final_metrics, data.table(
    metric = col,
    mean_value = mean_val,
    sd_value = sd_val,
    min_value = min_val,
    max_value = max_val,
    n_reps = nrow(all_overall_metrics)
  ))
}

# Add yob_interval row first
yob_row <- data.table(
  metric = "yob_interval",
  mean_value = NA_real_,
  sd_value = NA_real_,
  min_value = NA_real_,
  max_value = NA_real_,
  n_reps = nrow(all_overall_metrics)
)
final_metrics <- rbind(yob_row, final_metrics, fill = TRUE)

# Insert n_reps row second (after yob_interval)
nreps_row <- data.table(
  metric = "n_reps",
  mean_value = nrow(all_overall_metrics),
  sd_value = NA_real_,
  min_value = NA_real_,
  max_value = NA_real_,
  n_reps = nrow(all_overall_metrics)
)
final_metrics <- rbind(final_metrics[1, ], nreps_row, final_metrics[-1, ], fill = TRUE)

write.csv(final_metrics, snakemake@output[["aggregate_metrics"]], row.names = FALSE)
cat("✓ Wrote aggregate_metrics.csv:", nrow(final_metrics), "rows,", ncol(final_metrics), "columns\n")
cat("Format: metric,mean_value,sd_value,min_value,max_value,n_reps (preserving original column order)\n")
