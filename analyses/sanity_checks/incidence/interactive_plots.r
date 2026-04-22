library(tidyverse)
library(ggplot2)

setwd("/hpc/hers_en/pbeele/simPed/analyses/sanity_checks/incidence")

# Load your simulated lifetime risk data
sim_data <- read.csv("results/population_metrics.csv") %>%
  dplyr::select(age, prob_als_at_age) %>%
    arrange(age) %>%
  mutate(
    cum_risk_sim = cumsum(prob_als_at_age)  # cumulative sum over age
  )

# Load the reference cumulative incidence data
ref_data <- read.table("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt", header = TRUE) %>%
  dplyr::select(age, cum_incidence_c9_ALS)

ref_data <- read.table("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt", header = TRUE) %>%
  dplyr::select(age, cum_incidence_patho_ALS)

ref_data <- read.table("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_risk.txt", header = TRUE) %>%
  dplyr::select(age, cum_risk_adjusted)

# merge sets
plot_data <- inner_join(sim_data, ref_data, by = "age") %>%
  rename(
    Simulated = cum_risk_sim,
    Reference = cum_risk_adjusted
  )

# plot
ggplot(plot_data, aes(x = age)) +
  geom_line(aes(y = 100 * Simulated, color = "Simulated Lifetime Risk"), size = 1) +
  geom_line(aes(y = 100 * Reference, color = "Reference Lifetime Risk"), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("Simulated Lifetime Risk" = "blue", "Reference Lifetime Risk" = "red")) +
  labs(
    title = "Cumulative Lifetime Risk of ALS by Age",
    x = "Age",
    y = "Cumulative Lifetime Risk (%)",
    color = "Legend"
  ) +
  theme_minimal()
