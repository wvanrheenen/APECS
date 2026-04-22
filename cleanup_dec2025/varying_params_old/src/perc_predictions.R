
# Combine and summarize metrics
input_file <- snakemake@input[["filtered"]]
output_file <- snakemake@output[["model"]]

# Read the filtered simulation file
results_df <- read.csv(input_file)

source("src/libraries_simPed.R")
source("src/prediction_model.R")

print(input_file)
print(file.exists(input_file))
print(file.info(input_file)$size)

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

# Write metrics to output file
write.csv(metrics_df, output_file, row.names=FALSE)

# Read the filtered simulation file
results_df <- read.csv(input_file)

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

# Write metrics to output file
write.csv(metrics_df, output_file, row.names=FALSE)