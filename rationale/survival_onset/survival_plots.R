#!/usr/bin/env Rscript

options(warn = 1)

library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)
library(RColorBrewer)


# Read command line args: input_file output_dir
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Time grid
time_grid <- seq(0, 20, by = 0.1)

# Four-class ColorBrewer RdBu palette
rdbu4 <- RColorBrewer::brewer.pal(4, "RdBu")
cols_emp_sim <- c(
  "Empirical" = rdbu4[1],
  "Simulated" = rdbu4[4]
)

########################################
# (A) ALS
########################################
als_emp <- data.frame(
  time = time_grid,
  density = dlnorm(time_grid, meanlog = log(2.5), sdlog = 0.68),
  type = "Empirical"
)

als_sim <- data %>%
  filter(!is.na(post_surv_ALS)) %>%
  mutate(time_bin = round(post_surv_ALS)) %>%
  count(time_bin, name = "n") %>%
  mutate(density = n / sum(n), type = "Simulated") %>%
  rename(time = time_bin)

als_df <- bind_rows(als_emp, als_sim)

p_als <- ggplot(als_df, aes(x = time, y = density, color = type)) +
  geom_line(linewidth = 1) +
  labs(title = "(A) Empirical vs. Simulated ALS Survival",
       x = "Time Since Onset (years)", y = "Density") +
  scale_color_manual(values = cols_emp_sim, name = "") +
  coord_cartesian(xlim = c(0, 20)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt"))

ggsave(file.path(output_dir, "als_survival_3x3.pdf"), p_als, 
       width = 4.5, height = 3, units = "in", dpi = 300)

########################################
# (B) FTD 
########################################
ftd_emp <- data.frame(
  time = time_grid,
  density = dlnorm(time_grid, meanlog = log(8.7), sdlog = 0.37),
  type = "Empirical"
)

ftd_sim <- data %>%
  filter(!is.na(post_surv_FTD)) %>%
  mutate(time_bin = round(post_surv_FTD)) %>%
  count(time_bin, name = "n") %>%
  mutate(density = n / sum(n), type = "Simulated") %>%
  rename(time = time_bin)

ftd_df <- bind_rows(ftd_emp, ftd_sim)

p_ftd <- ggplot(ftd_df, aes(x = time, y = density, color = type)) +
  geom_line(linewidth = 1) +
  labs(title = "(B) Empirical vs. Simulated FTD Survival",
       x = "Time Since Onset (years)", y = "Density") +
  scale_color_manual(values = cols_emp_sim, name = "") +
  coord_cartesian(xlim = c(0, 20)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt"))

ggsave(file.path(output_dir, "ftd_survival_3x3.pdf"), p_ftd, 
       width = 4.5, height = 3, units = "in", dpi = 300)

########################################
# (C) Dementia
########################################

# Eleven-class ColorBrewer RdBu palette
rdbu11 <- RColorBrewer::brewer.pal(11, "RdBu")

data_dem <- data %>%
  filter(!is.na(post_surv_dem), !is.na(age_dementia)) %>%
  mutate(age_group = case_when(
    age_dementia < 70 ~ "<70",
    age_dementia < 80 ~ "70–79", 
    age_dementia < 90 ~ "80–89",
    TRUE ~ "90+"
  ))

dem_emp <- bind_rows(
  data.frame(time = time_grid, density = dlnorm(time_grid, log(10.7), 0.42), age_group = "<70"),
  data.frame(time = time_grid, density = dlnorm(time_grid, log(5.4), 0.44), age_group = "70–79"),
  data.frame(time = time_grid, density = dlnorm(time_grid, log(4.3), 0.48), age_group = "80–89"),
  data.frame(time = time_grid, density = dlnorm(time_grid, log(3.8), 0.40), age_group = "90+")
) %>% mutate(type = "Empirical")

dem_sim <- data_dem %>%
  mutate(time_bin = round(post_surv_dem)) %>%
  count(age_group, time_bin, name = "n") %>%
  group_by(age_group) %>%
  mutate(density = n / sum(n)) %>%
  ungroup() %>%
  rename(time = time_bin) %>%
  mutate(type = "Simulated")

dem_df <- bind_rows(dem_emp, dem_sim) %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c("<70", "70–79", "80–89", "90+")
    ),
    type = factor(
      type,
      levels = c("Empirical", "Simulated")
    ),
    group_type = factor(
      paste(age_group, type, sep = " - "),
      levels = c(
        "<70 - Empirical", "<70 - Simulated",
        "70–79 - Empirical", "70–79 - Simulated",
        "80–89 - Empirical", "80–89 - Simulated",
        "90+ - Empirical", "90+ - Simulated"
      )
    )
  )


# Dementia: empirical = red shades; simulated = blue shades.
# Colours darken with later onset age.
cols_dem <- c(
  "<70 - Empirical"   = rdbu11[5],
  "<70 - Simulated"   = rdbu11[8],
  "70–79 - Empirical" = rdbu11[4],
  "70–79 - Simulated" = rdbu11[9],
  "80–89 - Empirical" = rdbu11[3],
  "80–89 - Simulated" = rdbu11[10],
  "90+ - Empirical"   = rdbu11[2],
  "90+ - Simulated"   = rdbu11[11]
)


p_dem <- ggplot(dem_df, aes(x = time, y = density, color = group_type, linetype = type)) +
  geom_line(linewidth = 1) +
  labs(title = "(C) Empirical vs. Simulated Dementia Survival",
       subtitle = "Survival after onset dependent on age at onset",  # Added subtitle
       x = "Time Since Onset (years)", y = "Density") +
  scale_color_manual(
    values = cols_dem,
    name = "Age at onset"
  ) +
  scale_linetype_manual(
    values = c("Empirical" = "solid", "Simulated" = "22"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 20)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 10)),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt"))

ggsave(file.path(output_dir, "dementia_survival_6x3.pdf"), p_dem, 
       width = 9, height = 3, units = "in", dpi = 300)

cat("All survival plots saved to:", output_dir, "\n")