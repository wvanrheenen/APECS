# Load necessary libraries
library(tidyverse)
library(ggplot2) 
library(scales)

# Read ALS incidence data
als <- read_csv("lifetime_ALS_risk_processed.csv") %>%
  rename(als_incidence = incidence_per_100k_PY)

# Read mortality data
mortality <- read_tsv("danish_population_deaths.txt") %>%
  mutate(age = as.numeric(gsub("[^0-9]", "", age))) %>%
  rename(mortality_rate = weighted_mortality_rate)

# Create a full age range (from min to max in both datasets)
full_ages <- tibble(age = min(c(als$age, mortality$age)):max(c(als$age, mortality$age)))

# Join and fill missing values with 0
df <- full_ages %>%
  left_join(als, by = "age") %>%
  left_join(mortality, by = "age") %>%
  mutate(
    als_incidence = replace_na(als_incidence, 0),
    mortality_rate = replace_na(mortality_rate, 0)
  )

# Calculate survival and disease-free probabilities
df <- df %>%
  mutate(
    # Probability of surviving to this age (cumulative product up to previous age)
    surv_prob = cumprod(1 - lag(mortality_rate, default = 0)),
    # Probability of being ALS-free up to this age (cumulative product up to previous age)
    alsfree_prob = cumprod(1 - lag(als_incidence / 100000, default = 0)),
    # Age-specific risk contribution
    age_contribution = (als_incidence / 100000) * surv_prob * alsfree_prob,
    cumulative_risk = cumsum(age_contribution),
    
    # New: Uncorrected calculations (without mortality adjustment)
    uncorrected_age_contribution = (als_incidence / 100000) * alsfree_prob,
    uncorrected_cumulative_risk = cumsum(uncorrected_age_contribution)  
  )

# Create combined plot
combined_plot <- ggplot(df) +
  geom_line(aes(x = age, y = cumulative_risk, color = "Competing Risk Adjusted"), 
            linewidth = 1) +
  geom_line(aes(x = age, y = uncorrected_cumulative_risk, 
                color = "No Mortality Adjustment"), 
            linewidth = 1, linetype = "dashed") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_color_manual(
    name = "Method",
    values = c("Competing Risk Adjusted" = "#1f78b4", 
               "No Mortality Adjustment" = "#e31a1c")
  ) +
  labs(
    title = "Cumulative ALS Risk Comparison",
    subtitle = "Effect of mortality adjustment on lifetime risk estimates",
    x = "Age (years)",
    y = "Cumulative Probability",
    caption = "Data source: Danish population statistics"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Display plot
print(combined_plot)

# Sum contributions to get lifetime risk
lifetime_risk <- sum(df$age_contribution)

# Print both risk estimates
cat(sprintf("Competing risk-adjusted lifetime risk: %.2f%%\n", lifetime_risk * 100))
cat(sprintf("Uncorrected lifetime risk: %.2f%%\n", 
            max(df$uncorrected_cumulative_risk) * 100))


# Define scaling factor based on known real-world risk vs model risk
scaling_factor <- (1/350) / sum(df$age_contribution)

# Generate age-specific risk data with adjusted cumulative risk
age_risk_data <- df %>%
  mutate(
    cum_inc = uncorrected_cumulative_risk, 
    cum_risk = cumulative_risk,
    cum_risk_adjusted = cumulative_risk * scaling_factor
  ) %>%
  select(age, cum_inc, cum_risk, cum_risk_adjusted) %>%
  arrange(age)

# Save as tab-separated file
write_tsv(age_risk_data, "lifetime_ALS_risk.txt")
write_tsv(age_risk_data, "/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_risk.txt")  