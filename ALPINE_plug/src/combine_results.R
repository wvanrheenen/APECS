# Combine and summarize metrics
source("src/libraries_simPed.R")

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

# Add ALPINE columns for 4th/5th degree relatives
combined_results <- combined_results %>%
  mutate(
    SPORADIC = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0) * 1,
    # No closer relatives affected AND at least one 4th degree ALS/FTD case
    ALPINE_3rd = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0 &
                  (relatives_3rd_als > 0 | relatives_3rd_ftd > 0)) * 1,
    # No closer relatives affected AND at least one 4th degree ALS/FTD case
    ALPINE_4th = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0 &
                  (relatives_4th_als > 0 | relatives_4th_ftd > 0)) * 1,
    
    # No closer relatives affected AND at least one 5th degree ALS/FTD case
    ALPINE_5th = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0 &
                  (relatives_5th_als > 0 | relatives_5th_ftd > 0)) * 1,

    # Any distant relative, without overlap 
    ALPINE_4th_5th = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0 &
                  (relatives_4th_als > 0 | relatives_4th_ftd > 0 |
                   relatives_5th_als > 0 | relatives_5th_ftd > 0)) * 1,    

    # Any distant relative, without overlap 
    ALPINE_any = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0 &
                  (relatives_3rd_als > 0 | relatives_3rd_ftd > 0 |
                   relatives_4th_als > 0 | relatives_4th_ftd > 0 |
                   relatives_5th_als > 0 | relatives_5th_ftd > 0)) * 1    
  )

# Overwrite the combined results file with new columns
write.csv(combined_results, snakemake@output[["combined_results_ALPINE"]], row.names = FALSE)

## Calculate rate of ALPINE cases among monogenic ALS; sporadic within 3rd degree, but more distant affected relatives
monogenic_index <- combined_results[combined_results$mendel_ALS_Y == 1, ]

n_monogenic_sporadic <- sum(monogenic_index$SPORADIC)
n_alpine_3rd <- sum(monogenic_index$ALPINE_3rd)
n_alpine_4th <- sum(monogenic_index$ALPINE_4th)
n_alpine_5th <- sum(monogenic_index$ALPINE_5th) 
n_alpine_4th_5th <- sum(monogenic_index$ALPINE_4th_5th) 
n_alpine_any <- sum(monogenic_index$ALPINE_any)

cat("Number of monogenic cases, sporadic within 3rd degree:", n_monogenic_sporadic, "\n")
cat("Number of sporadic monogenic cases, with affected 3rd degree relative:", n_alpine_3rd, "\n")
cat("Number of sporadic monogenic cases, with affected 4th degree relative:", n_alpine_4th, "\n")
cat("Number of sporadic monogenic cases, with affected 5th degree relative:", n_alpine_5th, "\n")
cat("Number of sporadic monogenic cases, with affected 4th or 5th degree relative:", n_alpine_4th_5th, "\n")
cat("Number of sporadic monogenic cases, with any affected relative >3rd degree:", n_alpine_any, "\n")
cat("ALPINE_3rd rate among sporadic monogenic:", round(n_alpine_3rd/n_monogenic_sporadic, 4), "\n")
cat("ALPINE_4th rate among sporadic monogenic:", round(n_alpine_4th/n_monogenic_sporadic, 4), "\n")
cat("ALPINE_4th_5th rate among sporadic monogenic:", round(n_alpine_4th_5th/n_monogenic_sporadic, 4), "\n")
cat("ALPINE_5th rate among sporadic monogenic:", round(n_alpine_5th/n_monogenic_sporadic, 4), "\n")
cat("ALPINE_any rate among sporadic monogenic:", round(n_alpine_any/n_monogenic_sporadic, 4), "\n\n")


# Calculate phenocopy metrics from combined_results
monogenic_index <- combined_results[combined_results$mendel_ALS_Y == 1, ]
n_monogenic_index <- nrow(monogenic_index)
n_phenocopies <- sum(monogenic_index$phenocopy, na.rm = TRUE)
phenocopy_rate <- ifelse(n_monogenic_index > 0, n_phenocopies / n_monogenic_index, NA)

cat("DEBUG: Total monogenic index cases:", n_monogenic_index, "\n")
cat("DEBUG: Phenocopies among them:", n_phenocopies, "\n")
cat("DEBUG: Overall phenocopy rate:", round(phenocopy_rate, 4), "\n\n")


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


add_metrics <- function(df, n_phenocopies, phenocopy_rate) {
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
      proportion_C9SOD1FUS_FTD = ifelse(n_mendelian_FTD > 0, (n_c9_FTD + n_SOD1FUS_FTD) / n_mendelian_FTD, NA),
      n_phenocopies = n_phenocopies,
      phenocopy_rate = phenocopy_rate      
    )
}

# Then call it with the scalar values:
summed_metrics <- add_metrics(summed_metrics, n_phenocopies, phenocopy_rate)
overall_row <- add_metrics(overall_row, n_phenocopies, phenocopy_rate)


# Write summarized metrics to CSV
final_metrics <- bind_rows(summed_metrics, overall_row)
write.csv(final_metrics, snakemake@output[["combined_metrics"]], row.names = FALSE)

## Predictor model
# Read the combined simulations data
results_df <- read.csv(snakemake@output[["combined_results"]])

# Source the prediction model functions
source("src/prediction_model.R")

# Run the prediction model analyses
scenarios_1_6 <- define_scenarios_1_6(results_df)

# Calculate metrics for all scenarios
metrics_list_1_6 <- lapply(scenarios_1_6, calculate_metrics_vector, mendelian = results_df$mendel_ALS_Y == 1)
metrics_list <- c(metrics_list_1_6)

# Convert the list to a dataframe
metrics_df <- do.call(rbind, metrics_list)

# Add row names as a column
metrics_df <- cbind(Scenario = rownames(metrics_df), as.data.frame(metrics_df))
rownames(metrics_df) <- NULL

# Round numeric columns to 3 decimal places
metrics_df[, 2:7] <- round(metrics_df[, 2:7], 3)

# Write the prediction model results
write.csv(metrics_df, file = snakemake@output[["prediction_model"]], row.names = FALSE)