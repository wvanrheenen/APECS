# Get the script's path
options( warn = 2 )
source("../../src/libraries_simPed.R")

## ANALYSES

# Add this Wilson CI function after source("src/libraries_simPed.R")
wilson_ci <- function(x, n, conf_level = 0.95) {
  if (n == 0 || x < 0 || x > n) return(c(NA, NA, NA))
  
  z <- qnorm(1 - (1 - conf_level) / 2)
  z2 <- z^2
  p_hat <- x / n
  center <- (x + z2 / 2) / (n + z2)
  margin <- z * sqrt((x * (n - x) + z2 / 4) / (n + z2)^2) / sqrt(n)
  
  ci_lower <- center - margin
  ci_upper <- center + margin
  
  # Clamp to [0,1]
  ci_lower <- pmax(0, ci_lower)
  ci_upper <- pmin(1, ci_upper)
  
  return(c(estimate = p_hat, ci_lower = ci_lower, ci_upper = ci_upper))
}

# Function to calculate metrics
calculate_metrics_vector <- function(condition, mendelian) {
  TP <- sum(condition & mendelian)
  FP <- sum(condition & !mendelian)
  TN <- sum(!condition & !mendelian)
  FN <- sum(!condition & mendelian)
  total <- TP + FP + TN + FN
  
  # Point estimates
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  accuracy <- (TP + TN) / total
  odds_ratio <- (TP * TN) / (FP * FN)
  
  # Wilson 95% CIs
  sens_ci <- wilson_ci(TP, TP + FN)
  spec_ci <- wilson_ci(TN, TN + FP)
  ppv_ci <- wilson_ci(TP, TP + FP)
  npv_ci <- wilson_ci(TN, TN + FN)
  acc_ci <- wilson_ci(TP + TN, total)
  
  return(c(
    Sensitivity = sensitivity, Sensitivity_CI_low = sens_ci[2], Sensitivity_CI_high = sens_ci[3],
    Specificity = specificity, Specificity_CI_low = spec_ci[2], Specificity_CI_high = spec_ci[3],
    PPV = PPV, PPV_CI_low = ppv_ci[2], PPV_CI_high = ppv_ci[3],
    NPV = NPV, NPV_CI_low = npv_ci[2], NPV_CI_high = npv_ci[3],
    Accuracy = accuracy, Accuracy_CI_low = acc_ci[2], Accuracy_CI_high = acc_ci[3],
    Odds_Ratio = odds_ratio
  ))
}

# Function to define scenarios
define_scenarios_1_9 <- function(results_df) {
  scenarios_1_9 <- list(
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
    "7.1) ≥3 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 3,
    "7.2) ≥3 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 3,
    "7.3) ≥3 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 3,                                    
    "8.1) ≥4 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 4,
    "8.2) ≥4 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 4,
    "8.3) ≥4 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 4,
    "9.1) ≥5 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 5,
    "9.2) ≥5 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 5,
    "9.3) ≥5 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 5
  )
  
  return(scenarios_1_9)
}

# Define scenarios 1.0 to 6.3
scenarios_1_9 <- list(
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
    "7.1) ≥3 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 3,
    "7.2) ≥3 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 3,
    "7.3) ≥3 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 3,                                    
    "8.1) ≥4 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 4,
    "8.2) ≥4 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 4,
    "8.3) ≥4 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 4,
    "9.1) ≥5 1st/2nd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique >= 5,
    "9.2) ≥5 1st/2nd/3rd ALS" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als >= 5,
    "9.3) ≥5 1st/2nd/3rd ALS/FTD" = results_df$relatives_1st_als + results_df$relatives_2nd_als + results_df$relatives_3rd_als + 
                                    results_df$relatives_1st_ftd_unique + results_df$relatives_2nd_ftd_unique + results_df$relatives_3rd_ftd_unique >= 5
)

# Function to create and print crosstab
print_crosstab <- function(scenario, name) {
  crosstab <- table(Scenario = scenario, Mendelian = results_df$mendel_ALS_Y == 1)
  cat(paste("Crosstab for Scenario", name, ":\n"))
  print(crosstab)
  cat("\n\n")
}

# Print crosstabs for scenarios 1.0 to 6.3
for (scenario_name in names(scenarios_1_9)) {
  print_crosstab(scenarios_1_9[[scenario_name]], scenario_name)
}
