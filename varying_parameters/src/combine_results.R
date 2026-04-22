# Combine and summarize metrics
source("../src/libraries_simPed.R")

safe_read_csv <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    cat("Skipping empty/missing file:", basename(file), "\n")
    return(NULL)
  }
  tryCatch({
    df <- read.csv(file)
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

# Combine simulation results
result_files <- snakemake@input[["results"]]
result_list <- lapply(result_files, safe_read_csv)
result_list <- result_list[!sapply(result_list, is.null)]
if (length(result_list) == 0) stop("No valid simulation files found")
combined_results <- do.call(rbind, result_list)
write.csv(combined_results, snakemake@output[["combined_results"]], row.names = FALSE)


# Calculate phenocopy metrics from combined_results
monogenic_index <- combined_results[combined_results$mendel_ALS_Y == 1, ]
n_monogenic_index <- nrow(monogenic_index)
n_phenocopies <- sum(monogenic_index$phenocopy, na.rm = TRUE)
n_phenocopies_1st <- sum(monogenic_index$phenocopy_1st, na.rm = TRUE)
n_phenocopies_2nd <- sum(monogenic_index$phenocopy_2nd, na.rm = TRUE)
n_phenocopies_3rd <- sum(monogenic_index$phenocopy_3rd, na.rm = TRUE)
phenocopy_rate <- ifelse(n_monogenic_index > 0, n_phenocopies / n_monogenic_index, NA)
phenocopy_rate_1st <- ifelse(n_monogenic_index > 0, n_phenocopies_1st / n_monogenic_index, NA)
phenocopy_rate_2nd <- ifelse(n_monogenic_index > 0, n_phenocopies_2nd / n_monogenic_index, NA)
phenocopy_rate_3rd <- ifelse(n_monogenic_index > 0, n_phenocopies_3rd / n_monogenic_index, NA)

cat("\nDEBUG: Total monogenic index cases:", n_monogenic_index, "\n")
cat("DEBUG: Pedigrees with phenocopies among them:", n_phenocopies, "\n")
cat("DEBUG: Overall phenocopy rate:", round(phenocopy_rate, 4), "\n")
cat("DEBUG: Phenocopy rate in 1st degree:", round(phenocopy_rate_1st, 4), "\n")
cat("DEBUG: Phenocopy rate in 2nd degree:", round(phenocopy_rate_2nd, 4), "\n")
cat("DEBUG: Phenocopy rate in 3rd degree:", round(phenocopy_rate_3rd, 4), "\n\n")

n_affected_relatives_1st_dgr <- sum(monogenic_index$relatives_1st_als)
n_affected_relatives_1st_dgr_poly <- sum(monogenic_index$relatives_1st_als_polygenic)
n_affected_relatives_2nd_dgr <- sum(monogenic_index$relatives_2nd_als)
n_affected_relatives_2nd_dgr_poly <- sum(monogenic_index$relatives_2nd_als_polygenic)
n_affected_relatives_3rd_dgr <- sum(monogenic_index$relatives_3rd_als)
n_affected_relatives_3rd_dgr_poly <- sum(monogenic_index$relatives_3rd_als_polygenic)

phenocopy_rate_relatives_1st <- ifelse(n_affected_relatives_1st_dgr > 0, n_affected_relatives_1st_dgr_poly / n_affected_relatives_1st_dgr, NA)
phenocopy_rate_relatives_2nd <- ifelse(n_affected_relatives_2nd_dgr > 0, n_affected_relatives_2nd_dgr_poly / n_affected_relatives_2nd_dgr, NA)
phenocopy_rate_relatives_3rd <- ifelse(n_affected_relatives_3rd_dgr > 0, n_affected_relatives_3rd_dgr_poly / n_affected_relatives_3rd_dgr, NA)

cat("\n\nDEBUG: Phenocopy individuals among affected in 1st degree:", round(phenocopy_rate_relatives_1st, 4), "\n")
cat("DEBUG: Phenocopy individuals among affected in 2nd degree:", round(phenocopy_rate_relatives_2nd, 4), "\n")
cat("DEBUG: Phenocopy individuals among affected in 3rd degree:", round(phenocopy_rate_relatives_3rd, 4), "\n\n")

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

phenocopy_metrics <- tibble(
  n_phenocopies = n_phenocopies,
  n_phenocopies_1st = n_phenocopies_1st,
  n_phenocopies_2nd = n_phenocopies_2nd,
  n_phenocopies_3rd = n_phenocopies_3rd,
  phenocopy_rate = phenocopy_rate,
  phenocopy_rate_1st = phenocopy_rate_1st,
  phenocopy_rate_2nd = phenocopy_rate_2nd,
  phenocopy_rate_3rd = phenocopy_rate_3rd
)

add_metrics <- function(df) {
  df %>%
    mutate(
      lifetime_risk_ALS = ifelse(n_individuals > 0, n_ALS / n_individuals, NA),
      lifetime_risk_FTD = ifelse(n_individuals > 0, n_FTD / n_individuals, NA),
      lifetime_risk_dementia = ifelse(n_individuals > 0, n_dementia / n_individuals, NA),
      penetrance_c9_ALS = ifelse(n_c9_alleles > 0, n_mendelian_c9_ALS / n_c9_alleles, NA),
      penetrance_c9_FTD = ifelse(n_c9_alleles > 0, n_c9_FTD / n_c9_alleles, NA),
      penetrance_c9_dementia = ifelse(n_c9_alleles > 0, n_c9_dementia / n_c9_alleles, NA),
      proportion_monogenic_ALS = ifelse(n_ALS > 0, n_mendelian_ALS / n_ALS, NA),
      proportion_monogenic_FTD = ifelse(n_FTD > 0, n_mendelian_FTD / n_FTD, NA),
      proportion_C9SOD1FUS_FTD = ifelse(n_mendelian_FTD > 0, (n_c9_FTD + n_SOD1FUS_FTD) / n_mendelian_FTD, NA)        
    )
}

# Then call it with the scalar values:
summed_metrics <- summed_metrics %>% 
  add_metrics() %>% 
  bind_cols(phenocopy_metrics)

overall_row <- overall_row %>% 
  add_metrics() %>% 
  bind_cols(phenocopy_metrics)
  
# Write summarized metrics to CSV
final_metrics <- bind_rows(summed_metrics, overall_row)
write.csv(final_metrics, snakemake@output[["combined_metrics"]], row.names = FALSE)

## Predictor model
# Read the combined simulations data
results_df <- read.csv(snakemake@output[["combined_results"]])

# Source the prediction model functions
source("../src/prediction_model.R")

# Run the prediction model analyses
scenarios_1_9 <- define_scenarios_1_9(results_df)

# Calculate metrics for all scenarios
metrics_list_1_9 <- lapply(scenarios_1_9, calculate_metrics_vector, mendelian = results_df$mendel_ALS_Y == 1)
metrics_list <- c(metrics_list_1_9)

# Convert the list to a dataframe
metrics_df <- do.call(rbind, metrics_list)

# Add row names as a column
metrics_df <- cbind(Scenario = rownames(metrics_df), as.data.frame(metrics_df))
rownames(metrics_df) <- NULL

# Round numeric columns to 3 decimal places
metrics_df <- metrics_df %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Write the prediction model results
write.csv(metrics_df, file = snakemake@output[["prediction_model"]], row.names = FALSE)