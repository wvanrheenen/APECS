# Combine and summarize metrics
source("src/libraries_simPed.R")

# Combine simulation results
result_files <- snakemake@input[["results"]]
combined_results <- do.call(rbind, lapply(result_files, read.csv))
write.csv(combined_results, snakemake@output[["combined_results"]], row.names = FALSE)


# Combine metrics
metric_files <- snakemake@input[["metrics"]]
combined_metrics <- do.call(rbind, lapply(metric_files, read.csv))

# Define columns to be summed and averaged
sum_columns <- c("Total_Simulations", "Total_Individuals", "Total_ALS_Cases", "Total_Mendelian_ALS_Cases",
                 "Total_Mendelian_ALS_Cases_common", "Total_Mendelian_ALS_Cases_patho", "Total_Polygenic_ALS_Cases",
                 "Total_FTD_Cases", "Total_Mendelian_FTD_Cases", "Total_Mendelian_FTD_Cases_common",
                 "Total_Mendelian_FTD_Cases_patho", "Total_Mendelian_FTD_Cases_ftd", "Total_Polygenic_FTD_Cases", "Total_Carriers",
                 "Total_ALS_FTD_Cases", "Total_Dementia_Cases", "Total_other_dementia_Cases",
                 "LastGen_Mendelian_ALS_Cases", "LastGen_Polygenic_ALS_Cases", "LastGen_Mendelian_FTD_Cases",
                 "LastGen_Polygenic_FTD_Cases", "LastGen_Dementia_Cases", "LastGen_Carriers", "LastGen_Sporadic", "LastGen_Sporadic_Mendel")

avg_columns <- c("Simulated_ALS_LifetimeRisk", "Mendelian_ALS_Percentage", "Polygenic_ALS_Percentage",
                 "Simulated_FTD_LifetimeRisk", "Mendelian_FTD_Percentage", "Polygenic_FTD_Percentage",
                 "Total_Penetrance_ALS", "Total_Penetrance_FTD", "Simulated_other_dementia_LifetimeRisk", "LastGen_Sporadic_Mendel_perc")

# Summarize metrics
summarized_metrics <- combined_metrics %>%
  summarise(
    across(all_of(sum_columns), sum),
    across(all_of(avg_columns), mean)
  )

# Calculate overall percentages and risks
total_als_cases <- summarized_metrics$Total_ALS_Cases
total_ftd_cases <- summarized_metrics$Total_FTD_Cases
total_individuals <- summarized_metrics$Total_Individuals

summarized_metrics <- combined_metrics %>%
  summarise(
    across(all_of(sum_columns), sum),
    across(all_of(avg_columns), ~mean(., na.rm = TRUE))
  ) %>%
  mutate(
    Overall_Penetrance_ALS = Total_Mendelian_ALS_Cases / Total_Carriers,
    Overall_Penetrance_FTD = Total_Mendelian_FTD_Cases / Total_Carriers,
    Overall_ALS_LifetimeRisk = Total_ALS_Cases / Total_Individuals,
    Overall_FTD_LifetimeRisk = Total_FTD_Cases / Total_Individuals,
    Overall_Mendelian_ALS_Percentage = Total_Mendelian_ALS_Cases / Total_ALS_Cases,
    Overall_Polygenic_ALS_Percentage = Total_Polygenic_ALS_Cases / Total_ALS_Cases,
    Overall_Mendelian_FTD_Percentage = Total_Mendelian_FTD_Cases / Total_FTD_Cases,
    Overall_Polygenic_FTD_Percentage = Total_Polygenic_FTD_Cases / Total_FTD_Cases
  )

# Write summarized metrics to CSV
write.csv(summarized_metrics, snakemake@output[["combined_metrics"]], row.names = FALSE)


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