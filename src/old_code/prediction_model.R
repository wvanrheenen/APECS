# Get the script's path
options( warn = 2 )
source("src/libraries_simPed.R")

## ANALYSES

# Function to calculate metrics
calculate_metrics_vector <- function(condition, mendelian) {
  TP <- sum(condition & mendelian)
  FP <- sum(condition & !mendelian)
  TN <- sum(!condition & !mendelian)
  FN <- sum(!condition & mendelian)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  odds_ratio <- (TP * TN) / (FP * FN)
  
  return(c(Sensitivity = sensitivity, Specificity = specificity, 
           PPV = PPV, NPV = NPV, Accuracy = accuracy, Odds_Ratio = odds_ratio))
}

# Function to define scenarios
define_scenarios_1_8 <- function(results_df) {
  scenarios_1_8 <- list(
    "1.1) ≥1 1st ALS" = results_df$relatives_1st_als >= 1,
    "1.2) ≥1 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als >= 1,
    "1.3) ≥1 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 1,
    "2.1) ≥2 1st ALS" = results_df$relatives_1st_als >= 2,
    "2.2) ≥2 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als >= 2,
    "2.3) ≥2 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 2,
    "3.1) ≥1 1st ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_1st_ftd_unique >= 1,
    "3.2) ≥1 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 1,
    "3.3) ≥1 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 1,
    "4.1) ≥2 1st ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_1st_ftd_unique >= 2,
    "4.2) ≥2 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 2,
    "4.3) ≥2 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 2,
    "5.1) ≥1 1st ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_1st_dementia_unique>= 1,
    "5.2) ≥1 1st/2nd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 1,
    "5.3) ≥1 1st/2nd/3rd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                         results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 1,
    "6.1) ≥2 1st ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_1st_dementia_unique>= 2,
    "6.2) ≥2 1st/2nd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 2,
    "6.3) ≥2 1st/2nd/3rd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                         results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 2,
    "7.0) 1 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1, 
    "7.1) 1 1st/2nd/3rd ALS, ≥1 1st FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & results_df$relatives_1st_ftd_unique >= 1,
    "7.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 1),
    "7.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 1),
    "8.1) 1 1st/2nd/3rd ALS, ≥1 1st dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & results_df$relatives_1st_dementia_unique>= 1,
    "8.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 1),
    "8.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 1),
    "Paper 1) 1 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als == 1,
    "Paper 2) 2 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als == 2,
    "Paper 3) 1 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique == 1,
    "Paper 4) 2 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique == 2,
    "Questionnaire 1) 1 1st degree ALS" = results_df$relatives_1st_als == 1,
    "Questionnaire 2) 2 1st degree ALS" = results_df$relatives_1st_als == 2,
    "Questionnaire 1) 1 1st degree ALS, 1 2nd degree FTD" = results_df$relatives_1st_als == 1 & results_df$relatives_2nd_ftd_unique == 1,
    "Questionnaire 1) 1 1st degree ALS, 1 2nd degree dementia" = results_df$relatives_1st_als == 1 & results_df$relatives_2nd_dementia_unique == 1
  )
  
  return(scenarios_1_8)
}

# Define scenarios 1.0 to 6.3
scenarios_1_6 <- list(
    "1.1) ≥1 1st ALS" = results_df$relatives_1st_als >= 1,
    "1.2) ≥1 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als >= 1,
    "1.3) ≥1 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 1,
    "2.1) ≥2 1st ALS" = results_df$relatives_1st_als >= 2,
    "2.2) ≥2 1st/2nd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als >= 2,
    "2.3) ≥2 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 2,
    "3.1) ≥1 1st ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_1st_ftd_unique >= 1,
    "3.2) ≥1 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 1,
    "3.3) ≥1 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 1,
    "4.1) ≥2 1st ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_1st_ftd_unique >= 2,
    "4.2) ≥2 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 2,
    "4.3) ≥2 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 2,
    "5.1) ≥1 1st ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_1st_dementia_unique>= 1,
    "5.2) ≥1 1st/2nd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 1,
    "5.3) ≥1 1st/2nd/3rd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                         results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 1,
    "6.1) ≥2 1st ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_1st_dementia_unique>= 2,
    "6.2) ≥2 1st/2nd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 2,
    "6.3) ≥2 1st/2nd/3rd ALS/dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                         results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 2
)

# Define scenarios 7.0 to 8.3
scenarios_7_8 <- list(
    "7.0) 1 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1, 
    "7.1) 1 1st/2nd/3rd ALS, ≥1 1st FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & results_df$relatives_1st_ftd_unique >= 1,
    "7.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 1),
    "7.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 1),
    "8.1) 1 1st/2nd/3rd ALS, ≥1 1st dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & results_df$relatives_1st_dementia_unique>= 1,
    "8.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique>= 1),
    "8.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd dementia" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als == 1 & (results_df$relatives_1st_dementia_unique+ results_df$relatives_2nd_dementia_unique+ results_df$relatives_3rd_dementia_unique>= 1)
)

# Function to create and print crosstab
print_crosstab <- function(scenario, name) {
  crosstab <- table(Scenario = scenario, Mendelian = results_df$mendel_ALS_Y == 1)
  cat(paste("Crosstab for Scenario", name, ":\n"))
  print(crosstab)
  cat("\n\n")
}

# Print crosstabs for scenarios 1.0 to 6.3
for (scenario_name in names(scenarios_1_6)) {
  print_crosstab(scenarios_1_6[[scenario_name]], scenario_name)
}

# Print crosstabs for scenarios 7.0 to 8.3
for (scenario_name in names(scenarios_7_8)) {
  print_crosstab(scenarios_7_8[[scenario_name]], scenario_name)
}

# Define scenarios_paper
scenarios_paper <- list(
  "Paper 1) 1 1st/2nd ALS" = (results_df$relatives_1st_als + results_df$relatives_2nd_als) == 1,
  "Paper 2) 2 1st/2nd ALS" = (results_df$relatives_1st_als + results_df$relatives_2nd_als) == 2,
  "Paper 3) 1 1st/2nd ALS/FTD" = (results_df$relatives_1st_als + results_df$relatives_2nd_als +
                                  results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique) == 1,
  "Paper 4) 2 1st/2nd ALS/FTD" = (results_df$relatives_1st_als + results_df$relatives_2nd_als +
                                  results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique) == 2
)

# Define scenarios_poster
scenarios_poster <- list(
  "Questionnaire 1) 1 1st degree ALS" = results_df$relatives_1st_als == 1,
  "Questionnaire 2) 2 1st degree ALS" = results_df$relatives_1st_als == 2,
  "Questionnaire 1) 1 1st degree ALS, 1 2nd degree FTD" = results_df$relatives_1st_als == 1 & results_df$relatives_2nd_ftd_unique == 1,
  "Questionnaire 1) 1 1st degree ALS, 1 2nd degree dementia" = results_df$relatives_1st_als == 1 & results_df$relatives_2nd_dementia_unique == 1
)

# Function to create and print crosstab
print_crosstab <- function(scenario, name) {
  crosstab <- table(Scenario = scenario, Mendelian = results_df$mendel_ALS_Y == 1)
  cat(paste("Crosstab for Scenario", name, ":\n"))
  print(crosstab)
  cat("\n\n")
}

# Print crosstabs for scenarios in paper
for (scenario_name in names(scenarios_paper)) {
  print_crosstab(scenarios_paper[[scenario_name]], scenario_name)
}

# Filter results_df to keep only pedigrees with 0 or 1 ALS relative
# Checks for multiple ALS relatives
checks <- list(
  "More than one 1st degree ALS relative" = sum(results_df$relatives_1st_als > 1),
  "More than one 2nd degree ALS relative" = sum(results_df$relatives_2nd_als > 1),
  "More than one 3rd degree ALS relative" = sum(results_df$relatives_3rd_als > 1),
  "More than one 1st/2nd/3rd degree ALS relative" = sum(results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als > 1)
)

# Calculate counts and percentages of mendelian ALS for each category
mendelian_counts <- list(
  "More than one 1st degree ALS relative" = sum(results_df$relatives_1st_als > 1 & results_df$mendel_ALS_Y == 1),
  "More than one 2nd degree ALS relative" = sum(results_df$relatives_2nd_als > 1 & results_df$mendel_ALS_Y == 1),
  "More than one 3rd degree ALS relative" = sum(results_df$relatives_3rd_als > 1 & results_df$mendel_ALS_Y == 1),
  "More than one 1st/2nd/3rd degree ALS relative" = sum((results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als > 1) & results_df$mendel_ALS_Y == 1)
)

# Print results
cat("Counts of individuals with multiple ALS relatives and mendelian ALS:\n")
for (category in names(checks)) {
  total_count <- checks[[category]]
  mendelian_count <- mendelian_counts[[category]]
  percentage <- ifelse(total_count > 0, (mendelian_count / total_count) * 100, 0)
  
  cat(paste0(category, ":\n"))
  cat(paste0("  Total count: ", total_count, "\n"))
  cat(paste0("  Count with mendelian ALS: ", mendelian_count, "\n"))
  cat(paste0("  Percentage with mendelian ALS: ", round(percentage, 2), "%\n\n"))
}

# Checks for multiple FTD relatives
ftd_checks <- list(
  "More than one 1st degree FTD relative" = sum(results_df$relatives_1st_ftd > 1),
  "More than one 2nd degree FTD relative" = sum(results_df$relatives_2nd_ftd > 1),
  "More than one 3rd degree FTD relative" = sum(results_df$relatives_3rd_ftd > 1),
  "More than one 1st/2nd/3rd degree FTD relative" = sum(results_df$relatives_1st_ftd + results_df$relatives_2nd_ftd + results_df$relatives_3rd_ftd > 1)
)

# Calculate counts and percentages of mendelian FTD for each category
mendelian_ftd_counts <- list(
  "More than one 1st degree FTD relative" = sum(results_df$relatives_1st_ftd > 1 & results_df$mendel_FTD_Y == 1),
  "More than one 2nd degree FTD relative" = sum(results_df$relatives_2nd_ftd > 1 & results_df$mendel_FTD_Y == 1),
  "More than one 3rd degree FTD relative" = sum(results_df$relatives_3rd_ftd > 1 & results_df$mendel_FTD_Y == 1),
  "More than one 1st/2nd/3rd degree FTD relative" = sum((results_df$relatives_1st_ftd + results_df$relatives_2nd_ftd + results_df$relatives_3rd_ftd > 1) & results_df$mendel_FTD_Y == 1)
)

# Print results
cat("Counts of individuals with multiple FTD relatives and mendelian FTD:\n")
for (category in names(ftd_checks)) {
  total_count <- ftd_checks[[category]]
  mendelian_count <- mendelian_ftd_counts[[category]]
  percentage <- ifelse(total_count > 0, (mendelian_count / total_count) * 100, 0)
  
  cat(paste0(category, ":\n"))
  cat(paste0("  Total count: ", total_count, "\n"))
  cat(paste0("  Count with mendelian FTD: ", mendelian_count, "\n"))
  cat(paste0("  Percentage with mendelian FTD: ", round(percentage, 2), "%\n\n"))
}


filtered_results_df <- results_df[results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als <= 1, ]



define_scenarios_9_10 <- function(filtered_results_df) {
  scenarios_9_10 <- list(
    "9.0) 1 1st/2nd/3rd ALS" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1, 
    "9.1) 1 1st/2nd/3rd ALS, ≥1 1st FTD" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & filtered_results_df$relatives_1st_ftd_unique >= 1,
    "9.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd FTD" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & (filtered_results_df$relatives_1st_ftd_unique + filtered_results_df$relatives_2nd_ftd_unique >= 1),
    "9.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd FTD" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & (filtered_results_df$relatives_1st_ftd_unique + filtered_results_df$relatives_2nd_ftd_unique + filtered_results_df$relatives_3rd_ftd_unique >= 1),
    "10.1) 1 1st/2nd/3rd ALS, ≥1 1st dementia" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & filtered_results_df$relatives_1st_dementia_unique>= 1,
    "10.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd dementia" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & (filtered_results_df$relatives_1st_dementia_unique+ filtered_results_df$relatives_2nd_dementia_unique>= 1),
    "10.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd dementia" = filtered_results_df$relatives_1st_als + filtered_results_df$relatives_2nd_als + filtered_results_df$relatives_3rd_als == 1 & (filtered_results_df$relatives_1st_dementia_unique+ filtered_results_df$relatives_2nd_dementia_unique+ filtered_results_df$relatives_3rd_dementia_unique>= 1)
  )
  return(scenarios_9_10)
}

# Crosstab for scenario 9.0
# Define scenarios 9.0 to 10.3
scenarios_9_10 <- define_scenarios_9_10(filtered_results_df)

# Function to create and print crosstab
print_crosstab <- function(scenario, name) {
  crosstab <- table(Scenario = scenario, Mendelian = filtered_results_df$mendel_ALS_Y == 1)
  cat(paste("Crosstab for Scenario", name, ":\n"))
  print(crosstab)
  cat("\n\n")
}

# Print crosstabs for scenarios 9.0 to 10.3
scenario_names <- c("9.0) 1 1st/2nd/3rd ALS", 
                    "9.1) 1 1st/2nd/3rd ALS, ≥1 1st FTD",
                    "9.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd FTD",
                    "9.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd FTD",
                    "10.1) 1 1st/2nd/3rd ALS, ≥1 1st dementia",
                    "10.2) 1 1st/2nd/3rd ALS, ≥1 1st/2nd dementia",
                    "10.3) 1 1st/2nd/3rd ALS, ≥1 1st/2nd/3rd dementia")

for (scenario_name in scenario_names) {
  print_crosstab(scenarios_9_10[[scenario_name]], scenario_name)
}
