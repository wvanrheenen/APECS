# ../src/aggregate_repetitions.R
# Aggregate prediction_model.csv across 5 repetitions for specific scenario
source("../src/libraries_simPed.R")

# Get input files (all 5 prediction_model.csv files for this {set})
pred_files <- snakemake@input[["prediction_models"]]
set_name <- snakemake@params[["set"]]

cat("Aggregating", length(pred_files), "repetition files for set:", set_name, "\n")

# Read all prediction model files safely
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

models_list <- lapply(pred_files, safe_read_csv)
models_list <- models_list[!sapply(models_list, is.null)]
if (length(models_list) == 0) stop("No valid prediction_model files found")

# Combine all models with rep identifier
all_models <- rbindlist(models_list, idcol = "rep", fill = TRUE)

# Filter to target scenario: "4.2) ≥2 1st/2nd ALS/FTD"
target_scenario <- "4.2) ≥2 1st/2nd ALS/FTD"
scenario_data <- all_models[Scenario == target_scenario, ]

if (nrow(scenario_data) == 0) {
  stop("Target scenario '", target_scenario, "' not found in any repetition files")
}

cat("Found scenario data for", nrow(scenario_data), "repetitions\n")

# Compute mean and SD for key metrics
summary_stats <- scenario_data[, .(
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
), ]

cat("Summary statistics computed for", summary_stats$n_reps, "repetitions\n")

# Create final output with Scenario column preserved
final_result <- data.table(
  Scenario = target_scenario,
  summary_stats[, c("mean_Sensitivity", "sd_Sensitivity", "mean_Specificity", "sd_Specificity", 
                   "mean_PPV", "sd_PPV", "mean_NPV", "sd_NPV", "mean_Accuracy", "sd_Accuracy",
                   "mean_Odds_Ratio", "sd_Odds_Ratio", "n_reps"), with = FALSE]
)

# Write final aggregated result
write.csv(final_result, snakemake@output[["prediction_model"]], row.names = FALSE)

cat("Wrote aggregated results to:", snakemake@output[["prediction_model"]], "\n")
cat("Key metrics for", target_scenario, ":\n")
cat("  Sensitivity: ", round(summary_stats$mean_Sensitivity, 3), " ±", round(summary_stats$sd_Sensitivity, 3), "\n")
cat("  PPV:        ", round(summary_stats$mean_PPV, 3), " ±", round(summary_stats$sd_PPV, 3), "\n")
cat("  Specificity : ", round(summary_stats$mean_Specificity, 3), " ±", round(summary_stats$sd_Specificity, 3), "\n")
cat("  NPV: ", round(summary_stats$mean_NPV, 3), " ±", round(summary_stats$sd_NPV, 3), "\n")
