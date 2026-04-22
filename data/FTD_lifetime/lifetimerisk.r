library(tidyverse)
library(ggplot2)
library(scales)

# Load FTD incidence data (5-year bins)
ftd <- read_csv("lifetime_ftd_risk_formatted.csv") %>%
  separate(`Age Bin`, into = c("age_start", "age_end"), convert = TRUE) %>%
  mutate(
    age_start = as.numeric(age_start),
    age_end = as.numeric(age_end),
    incidence_rate = `Incidence Rate per 100k PY` / 100000
  ) %>%
  select(age_start, age_end, incidence_rate)

# Load mortality data (5-year bins)
mortality <- read_csv("death-rate-by-age-group-in-england-and-wales.csv") %>%
  select(-Entity, -Code, -Year) %>%
  select(1:(ncol(.) - 18)) %>%   # Drop the last 18 columns
  pivot_longer(
    cols = everything(),
    names_to = "age_group",
    values_to = "death_rate_per_1000"
  ) %>%
  mutate(
    death_rate = death_rate_per_1000 / 1000,  # Convert to per-person rate
    age_group = case_when(
      age_group == "<1 year old" ~ "0-0",
      age_group == "1-4 years old" ~ "1-4",
      age_group == "5-9 years old" ~ "5-9",
      TRUE ~ str_replace(age_group, " years old|\\+", "")
    )
  ) %>%
  separate(age_group, into = c("age_start", "age_end"), convert = TRUE) %>%
  mutate(
    age_end = ifelse(is.na(age_end), age_start, age_end),
    age_end = case_when(
      age_start == 0 & age_end == 0 ~ 0,
      age_start == 1 & age_end == 4 ~ 4,
      age_start == 5 & age_end == 9 ~ 9,
      TRUE ~ age_start + 4
    )
  )

# Create full age grid (0-4,5-9,...,85-89)
age_grid <- tibble(
  age_start = seq(0, 85, 5),
  age_end = age_start + 4
) %>% 
  add_row(age_start = 90, age_end = 99)  # Adjust for mortality's 80+ group

# Merge data with full age grid
df <- age_grid %>%
  left_join(ftd, by = c("age_start", "age_end")) %>%
  left_join(mortality, by = c("age_start", "age_end")) %>%
  mutate(
    incidence_rate = replace_na(incidence_rate, 0),
    death_rate = approx(
      x = age_start + 2.5,  # Use mid-bin for interpolation
      y = death_rate,
      xout = age_start + 2.5,
      rule = 2
    )$y
  )

# Calculate cumulative probabilities with and without mortality adjustment
df <- df %>%
  mutate(
    # Probability of surviving this 5-year bin
    survival_this_bin = (1 - death_rate)^5,
    # Cumulative survival to start of bin
    cumulative_survival = cumprod(lag(survival_this_bin, default = 1)),
    # Probability of remaining FTD-free to start of bin
    disease_free = cumprod(1 - 5 * lag(incidence_rate, default = 0)),
    # Risk contribution from this bin (adjusted for mortality)
    risk_contribution = 5 * incidence_rate * cumulative_survival * disease_free,
    cumulative_risk = cumsum(risk_contribution),
    
    # Uncorrected risk contribution (no mortality adjustment)
    uncorrected_risk_contribution = 5 * incidence_rate * disease_free,
    uncorrected_cumulative_risk = cumsum(uncorrected_risk_contribution)
  )

# Plot cumulative risks: adjusted vs unadjusted
cumulative_plot <- ggplot(df) +
  geom_line(aes(x = age_start, y = cumulative_risk, color = "Adjusted for Mortality"), size = 1) +
  geom_line(aes(x = age_start, y = uncorrected_cumulative_risk, color = "Unadjusted"), size = 1, linetype = "dashed") +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
  scale_color_manual(
    name = "Risk Type",
    values = c("Adjusted for Mortality" = "#1f78b4", "Unadjusted" = "#e31a1c")
  ) +
  labs(
    title = "Cumulative Lifetime Risk of FTD by Age",
    subtitle = "Comparison of mortality-adjusted and unadjusted estimates",
    x = "Age (years)",
    y = "Cumulative Probability (%)",
    caption = "Data source: UK mortality and FTD incidence data"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(cumulative_plot)

# Print lifetime risk estimates
cat(sprintf("Mortality-adjusted lifetime FTD risk: %.3f%%\n", max(df$cumulative_risk) * 100))
cat(sprintf("Unadjusted lifetime FTD risk: %.3f%%\n", max(df$uncorrected_cumulative_risk) * 100))

# Define scaling factor based on known real-world risk vs model risk
scaling_factor <- (1/750) / 0.002708

# Generate age-specific risk data with adjusted cumulative risk
age_risk_data <- df %>%
  rowwise() %>%
  mutate(age = list(seq(age_start, age_end))) %>%
  unnest(age) %>%
  mutate(cum_risk_adjusted = cumulative_risk * scaling_factor) %>%
  select(age, cum_risk = cumulative_risk, cum_risk_adjusted) %>%
  arrange(age)

# Save as tab-separated file
write_tsv(age_risk_data, "ftd_lifetime_risk_by_age.txt")
write_tsv(age_risk_data, "/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_FTD_risk.txt")