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

# Add SPORADIC column
combined_results <- combined_results %>%
  mutate(
    SPORADIC_close = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_close_3rd_als == 0 & relatives_close_3rd_ftd == 0) * 1,
    SPORADIC_all = (relatives_1st_als == 0 & relatives_1st_ftd == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd == 0 &
                  relatives_3rd_als == 0 & relatives_3rd_ftd == 0) * 1,
    SPORADIC_unique = (relatives_1st_als == 0 & relatives_1st_ftd_unique == 0 &
                  relatives_2nd_als == 0 & relatives_2nd_ftd_unique == 0 &
                  relatives_3rd_als == 0 & relatives_3rd_ftd_unique == 0) * 1,                                        
  )


# TOTAL COUNTS
n_total_index <- nrow(combined_results)

# 1) Monogenic/polygenic proportions (among ALL index patients)
n_monogenic <- sum(combined_results$mendel_ALS_Y == 1, na.rm = TRUE)
n_polygenic <- sum(combined_results$polygenicY_ALS == 1, na.rm = TRUE)
prop_monogenic_total <- n_monogenic / n_total_index
prop_polygenic_total <- n_polygenic / n_total_index

# 2) Sporadic/non-sporadic proportions (among ALL index patients)  
n_sporadic_close <- sum(combined_results$SPORADIC_close == 1)
n_familial_close <- sum(combined_results$SPORADIC_close == 0)
prop_sporadic_close <- n_sporadic_close / n_total_index
prop_familial_close <- n_familial_close / n_total_index

# 3) Sporadic proportions WITHIN monogenic/polygenic
monogenic_index <- combined_results[combined_results$mendel_ALS_Y == 1, ]
polygenic_index <- combined_results[combined_results$polygenicY_ALS == 1, ]

prop_sporadic_close_mono <- sum(monogenic_index$SPORADIC_close) / nrow(monogenic_index)
prop_sporadic_close_poly <- sum(polygenic_index$SPORADIC_close) / nrow(polygenic_index)

# 4) Monogenic/polygenic proportions WITHIN sporadic/familial
sporadic_index <- combined_results[combined_results$SPORADIC_close == 1, ]
familial_index <- combined_results[combined_results$SPORADIC_close == 0, ]

prop_monogenic_sporadic <- sum(sporadic_index$mendel_ALS_Y == 1) / nrow(sporadic_index)
prop_polygenic_sporadic <- sum(sporadic_index$polygenicY_ALS == 1) / nrow(sporadic_index)
prop_monogenic_familial <- sum(familial_index$mendel_ALS_Y == 1) / nrow(familial_index)
prop_polygenic_familial <- sum(familial_index$polygenicY_ALS == 1) / nrow(familial_index)

# Count sums of affected relatives
total_relatives           <- sum(combined_results$relatives_1st) + sum(combined_results$relatives_2nd) + sum(combined_results$relatives_3rd)
total_relatives_als       <- sum(combined_results$relatives_1st_als) + sum(combined_results$relatives_2nd_als) + sum(combined_results$relatives_3rd_als)
total_relatives_ftd       <- sum(combined_results$relatives_1st_ftd) + sum(combined_results$relatives_2nd_ftd) + sum(combined_results$relatives_3rd_ftd)
total_relatives_dementia  <- sum(combined_results$relatives_1st_dementia) + sum(combined_results$relatives_2nd_dementia) + sum(combined_results$relatives_3rd_dementia)

# PRINT RESULTS
cat("=== SIMULATION SUMMARY ===\n")
cat(sprintf("Total index patients: %d\n", n_total_index))

cat(sprintf("\nTotal number of relatives: %d\n", total_relatives))
cat(sprintf("\nTotal number of relatives with ALS: %d\n", total_relatives_als))
cat(sprintf("\nTotal number of relatives with FTD: %d\n", total_relatives_ftd))
cat(sprintf("\nTotal number of relatives with any dementia: %d\n", total_relatives_dementia))


cat(sprintf("\n1) PROPORTION BY ETIOLOGY:\n"))
cat(sprintf("   Monogenic: %.1f%% (%d/%d)\n", prop_monogenic_total*100, n_monogenic, n_total_index))
cat(sprintf("   Polygenic: %.1f%% (%d/%d)\n", prop_polygenic_total*100, n_polygenic, n_total_index))

cat(sprintf("\n2) PROPORTION BY FAMILY HISTORY:\n"))
cat(sprintf("   Sporadic (close 3rd): %.1f%% (%d/%d)\n", prop_sporadic_close*100, n_sporadic_close, n_total_index))
cat(sprintf("   Familial: %.1f%% (%d/%d)\n", prop_familial_close*100, n_familial_close, n_total_index))

cat(sprintf("\n3) SPORADIC PATIENTS AMONG ALL MONOGENIC/POLYGENIC INDEX PATIENTS:\n"))
cat(sprintf("   Monogenic sporadic: %.1f%% (%d/%d)\n", prop_sporadic_close_mono*100, 
            sum(monogenic_index$SPORADIC_close), nrow(monogenic_index)))
cat(sprintf("   Polygenic sporadic: %.1f%% (%d/%d)\n", prop_sporadic_close_poly*100, 
            sum(polygenic_index$SPORADIC_close), nrow(polygenic_index)))

cat(sprintf("\n4) ETIOLOGY AMONG ALL SPORADIC/FAMILIAL INDEX PATIENTS:\n"))
cat(sprintf("   Sporadic → Monogenic: %.1f%% (%d/%d)\n", prop_monogenic_sporadic*100, 
            sum(sporadic_index$mendel_ALS_Y == 1), nrow(sporadic_index)))
cat(sprintf("   Sporadic → Polygenic: %.1f%% (%d/%d)\n", prop_polygenic_sporadic*100, 
            sum(sporadic_index$polygenicY_ALS == 1), nrow(sporadic_index)))
cat(sprintf("   Familial → Monogenic: %.1f%% (%d/%d)\n", prop_monogenic_familial*100, 
            sum(familial_index$mendel_ALS_Y == 1), nrow(familial_index)))
cat(sprintf("   Familial → Polygenic: %.1f%% (%d/%d)\n", prop_polygenic_familial*100, 
            sum(familial_index$polygenicY_ALS == 1), nrow(familial_index)))
cat("==========================\n")

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

cat("\nPhenocopy rate as number of pedigrees with any phenocopy among total number of pedigrees:\n")
cat("\nTotal monogenic index cases:", n_monogenic_index, "\n")
cat("Pedigrees with phenocopies among them:", n_phenocopies, "\n")
cat(sprintf(
  "Overall phenocopy rate: %.4f (%d/%d)\n",
  phenocopy_rate, n_phenocopies, n_monogenic_index
))
cat(sprintf(
  "Phenocopy rate in 1st degree: %.4f (%d/%d)\n",
  phenocopy_rate_1st, n_phenocopies_1st, n_monogenic_index
))
cat(sprintf(
  "Phenocopy rate in 2nd degree: %.4f (%d/%d)\n",
  phenocopy_rate_2nd, n_phenocopies_2nd, n_monogenic_index
))
cat(sprintf(
  "Phenocopy rate in 3rd degree: %.4f (%d/%d)\n",
  phenocopy_rate_3rd, n_phenocopies_3rd, n_monogenic_index
))

# Total affected relatives per degree
n_affected_relatives_1st_dgr <- sum(monogenic_index$relatives_1st_als)
n_affected_relatives_2nd_dgr <- sum(monogenic_index$relatives_2nd_als)
n_affected_relatives_3rd_dgr <- sum(monogenic_index$relatives_3rd_als)

# Phenocopy components per degree
n_affected_relatives_1st_dgr_polygenic <- sum(monogenic_index$relatives_1st_als_polygenic)
n_affected_relatives_1st_dgr_diff_ancestor <- sum(monogenic_index$relatives_1st_different_monogenic_ancestor)

n_affected_relatives_2nd_dgr_polygenic <- sum(monogenic_index$relatives_2nd_als_polygenic)
n_affected_relatives_2nd_dgr_diff_ancestor <- sum(monogenic_index$relatives_2nd_different_monogenic_ancestor)

n_affected_relatives_3rd_dgr_polygenic <- sum(monogenic_index$relatives_3rd_als_polygenic)
n_affected_relatives_3rd_dgr_diff_ancestor <- sum(monogenic_index$relatives_3rd_different_monogenic_ancestor)

# Total phenocopy counts per degree (polygenic + different-ancestor monogenic)
n_affected_relatives_1st_dgr_phenocopy <- n_affected_relatives_1st_dgr_polygenic +
  n_affected_relatives_1st_dgr_diff_ancestor

n_affected_relatives_2nd_dgr_phenocopy <- n_affected_relatives_2nd_dgr_polygenic +
  n_affected_relatives_2nd_dgr_diff_ancestor

n_affected_relatives_3rd_dgr_phenocopy <- n_affected_relatives_3rd_dgr_polygenic +
  n_affected_relatives_3rd_dgr_diff_ancestor

# Phencopy rate as number of phencopies among number of affected relatives
phenocopy_rate_relatives_1st <- ifelse(n_affected_relatives_1st_dgr > 0, n_affected_relatives_1st_dgr_phenocopy / n_affected_relatives_1st_dgr, NA)
phenocopy_rate_relatives_2nd <- ifelse(n_affected_relatives_2nd_dgr > 0, n_affected_relatives_2nd_dgr_phenocopy / n_affected_relatives_2nd_dgr, NA)
phenocopy_rate_relatives_3rd <- ifelse(n_affected_relatives_3rd_dgr > 0, n_affected_relatives_3rd_dgr_phenocopy / n_affected_relatives_3rd_dgr, NA)

cat("\nPhenocopy rate as number of phencopies among number of affected relatives:\n")
cat("\nTotal number of affected 1st degree relatives among monogenic index patients:",
    n_affected_relatives_1st_dgr, "\n")
cat("Total number of affected 2nd degree relatives among monogenic index patients:",
    n_affected_relatives_2nd_dgr, "\n")
cat("Total number of affected 3rd degree relatives among monogenic index patients:",
    n_affected_relatives_3rd_dgr, "\n")

cat(sprintf(
  "\nPhenocopy individuals among affected in 1st degree: %.4f (%d/%d) [polygenic: %d, diff-ancestor monogenic: %d]\n",
  phenocopy_rate_relatives_1st,
  n_affected_relatives_1st_dgr_phenocopy,
  n_affected_relatives_1st_dgr,
  n_affected_relatives_1st_dgr_polygenic,
  n_affected_relatives_1st_dgr_diff_ancestor
))

cat(sprintf(
  "Phenocopy individuals among affected in 2nd degree: %.4f (%d/%d) [polygenic: %d, diff-ancestor monogenic: %d]\n",
  phenocopy_rate_relatives_2nd,
  n_affected_relatives_2nd_dgr_phenocopy,
  n_affected_relatives_2nd_dgr,
  n_affected_relatives_2nd_dgr_polygenic,
  n_affected_relatives_2nd_dgr_diff_ancestor
))

cat(sprintf(
  "Phenocopy individuals among affected in 3rd degree: %.4f (%d/%d) [polygenic: %d, diff-ancestor monogenic: %d]\n\n",
  phenocopy_rate_relatives_3rd,
  n_affected_relatives_3rd_dgr_phenocopy,
  n_affected_relatives_3rd_dgr,
  n_affected_relatives_3rd_dgr_polygenic,
  n_affected_relatives_3rd_dgr_diff_ancestor
))

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
      lifetime_risk_ALS = round(ifelse(n_individuals > 0, n_ALS / n_individuals, NA), 4),
      lifetime_risk_FTD = round(ifelse(n_individuals > 0, n_FTD / n_individuals, NA), 4),
      lifetime_risk_dementia = round(ifelse(n_individuals > 0, n_dementia / n_individuals, NA), 4),
      penetrance_c9_ALS = round(ifelse(n_c9_alleles > 0, n_mendelian_c9_ALS / n_c9_alleles, NA), 4),
      penetrance_c9_FTD = round(ifelse(n_c9_alleles > 0, n_c9_FTD / n_c9_alleles, NA), 4),
      penetrance_c9_dementia = round(ifelse(n_c9_alleles > 0, n_c9_dementia / n_c9_alleles, NA), 4),
      proportion_monogenic_ALS = round(ifelse(n_ALS > 0, n_mendelian_ALS / n_ALS, NA), 4),
      proportion_monogenic_FTD = round(ifelse(n_FTD > 0, n_mendelian_FTD / n_FTD, NA), 4),
      proportion_C9SOD1FUS_FTD = round(ifelse(n_mendelian_FTD > 0, (n_c9_FTD + n_SOD1FUS_FTD) / n_mendelian_FTD, NA), 4)        
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
source("src/prediction_model.R")

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