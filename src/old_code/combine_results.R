# Combine and summarize metrics
source("src/libraries_simPed.R")

# Combine simulation results
result_files <- snakemake@input[["results"]]
combined_results <- do.call(rbind, lapply(result_files, read.csv))
write.csv(combined_results, snakemake@output[["combined_results"]], row.names = FALSE)

# Combine metrics
metric_files <- snakemake@input[["metrics"]]
combined_metrics <- do.call(rbind, lapply(metric_files, read.csv))

# Sum all columns per yob_interval
summed_metrics <- combined_metrics %>%
  group_by(yob_interval) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

overall_row <- summed_metrics %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(yob_interval = "overall") %>%
  dplyr::select(yob_interval, everything())

add_metrics <- function(df) {
  df %>%
    mutate(
      lifetime_risk_ALS = ifelse(n_individuals > 0, n_ALS / n_individuals, NA),
      lifetime_risk_FTD = ifelse(n_individuals > 0, n_FTD / n_individuals, NA),
      lifetime_risk_dementia = ifelse(n_individuals > 0, n_dementia / n_individuals, NA),
      penetrance_c9_ALS = ifelse(n_c9_alleles > 0, n_mendelian_ALS / n_c9_alleles, NA),
      penetrance_c9_FTD = ifelse(n_c9_alleles > 0, n_c9_FTD / n_c9_alleles, NA),
      penetrance_c9_dementia = ifelse(n_c9_alleles > 0, n_c9_dementia / n_c9_alleles, NA),
      proportion_monogenic_ALS = ifelse(n_ALS > 0, n_mendelian_ALS / n_ALS, NA),
      proportion_monogenic_FTD = ifelse(n_FTD > 0, n_mendelian_FTD / n_FTD, NA),
      proportion_c9_FTD = ifelse(n_mendelian_FTD > 0, n_c9_FTD / n_mendelian_FTD, NA)
    )
}
summed_metrics <- add_metrics(summed_metrics)
overall_row <- add_metrics(overall_row)

# Write summarized metrics to CSV
final_metrics <- bind_rows(summed_metrics, overall_row)
write.csv(final_metrics, snakemake@output[["combined_metrics"]], row.names = FALSE)

## Predictor model
# Read the combined simulations data
results_df <- read.csv(snakemake@output[["combined_results"]])

# Source the prediction model functions
source("src/prediction_model.R")

# Run the prediction model analyses
scenarios_1_8 <- define_scenarios_1_8(results_df)
filtered_results_df <- results_df[results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als <= 1, ]
scenarios_9_10 <- define_scenarios_9_10(filtered_results_df)
all_scenarios <- c(scenarios_1_8, scenarios_9_10)

# Calculate metrics for all scenarios
metrics_list_1_8 <- lapply(scenarios_1_8, calculate_metrics_vector, mendelian = results_df$mendel_ALS_Y == 1)
metrics_list_9_10 <- lapply(scenarios_9_10, calculate_metrics_vector, mendelian = filtered_results_df$mendel_ALS_Y == 1)
metrics_list <- c(metrics_list_1_8, metrics_list_9_10)

# Convert the list to a dataframe
metrics_df <- do.call(rbind, metrics_list)

# Add row names as a column
metrics_df <- cbind(Scenario = rownames(metrics_df), as.data.frame(metrics_df))
rownames(metrics_df) <- NULL

# Round numeric columns to 3 decimal places
metrics_df[, 2:7] <- round(metrics_df[, 2:7], 3)

# Write the prediction model results
write.csv(metrics_df, file = snakemake@output[["prediction_model"]], row.names = FALSE)