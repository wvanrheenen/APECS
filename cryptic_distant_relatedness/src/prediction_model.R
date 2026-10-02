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
define_scenarios_1_6 <- function(results_df) {
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
  
  return(scenarios_1_6)
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
