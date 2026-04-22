#!/usr/bin/env Rscript
# yob_aao_analysis.R: Analyze year of birth, onset year, and age-at-onset distributions

# Load libraries
library(ggplot2)
library(dplyr)

# Read data (ONCE)
df <- read.csv("current_year/combined_simulations.csv", stringsAsFactors = FALSE)

# Counts (ONCE)
cat("Total ALS cases (Y_ALS == 1):\n"); print(sum(df$Y_ALS == 1, na.rm = TRUE))
cat("\nMonogenic ALS (mendel_ALS_Y == 1):\n"); print(sum(df$mendel_ALS_Y == 1, na.rm = TRUE))
cat("\nPolygenic ALS (polygenicY_ALS == 1):\n"); print(sum(df$polygenicY_ALS == 1, na.rm = TRUE))
cat("\nContingency table:\n"); print(table(Monogenic = df$mendel_ALS_Y, Polygenic = df$polygenicY_ALS, useNA = "ifany"))

## YEAR OF BIRTH ANALYSIS (1925-2005) ##
monogenic_yob <- df[df$mendel_ALS_Y == 1, ]
polygenic_yob <- df[df$polygenicY_ALS == 1, ]
n_mono_yob <- nrow(monogenic_yob); n_poly_yob <- nrow(polygenic_yob)

p_yob <- ggplot(df, aes(x = yob)) + geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Year of Birth - Affected Index Patients", x = "Year of Birth", y = "Count") +
  xlim(1925, 2005) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave("yob_histo.pdf", plot = p_yob, width = 10, height = 6, dpi = 300)

p_mono_yob <- ggplot(monogenic_yob, aes(x = yob)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "red", color = "black", alpha = 0.8) +
  labs(title = "Monogenic ALS - Year of Birth", subtitle = paste("n =", n_mono_yob), x = "Year of Birth", y = "Percentage (%)") + 
  xlim(1925, 2005) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p_poly_yob <- ggplot(polygenic_yob, aes(x = yob)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Polygenic ALS - Year of Birth", subtitle = paste("n =", n_poly_yob), x = "Year of Birth", y = "Percentage (%)") + 
  xlim(1925, 2005) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

df_yob_combined <- rbind(
  data.frame(value = monogenic_yob$yob, group = "Monogenic"),
  data.frame(value = polygenic_yob$yob, group = "Polygenic")
)
p_yob_density <- ggplot(df_yob_combined, aes(x = value, color = group, fill = group)) +
  geom_density(aes(y = ..density.. * 100), size = 1.2, alpha = 0.3) +
#   geom_histogram(aes(y = ..density.. * 100), binwidth = 1, alpha = 0.4, position = "identity") +
  scale_color_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  scale_fill_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  labs(title = "Year of Birth: Monogenic vs Polygenic Density Overlay", x = "Year of Birth", y = "Percentage (%)") +
  xlim(1925, 2005) + theme_minimal() + theme(legend.position = "top")

ggsave("yob_densityplot.pdf", p_yob_density, width = 10, height = 6, dpi = 300)
ggsave("yob_mono_percent.pdf", p_mono_yob, width = 10, height = 6, dpi = 300)
ggsave("yob_poly_percent.pdf", p_poly_yob, width = 10, height = 6, dpi = 300)


cat("\n=== YEAR OF BIRTH ===\n")
cat(sprintf("Monogenic (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_mono_yob, median(monogenic_yob$yob, na.rm=TRUE),
            quantile(monogenic_yob$yob, 0.25, na.rm=TRUE), quantile(monogenic_yob$yob, 0.75, na.rm=TRUE)))
cat(sprintf("Polygenic  (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_poly_yob, median(polygenic_yob$yob, na.rm=TRUE),
            quantile(polygenic_yob$yob, 0.25, na.rm=TRUE), quantile(polygenic_yob$yob, 0.75, na.rm=TRUE)))

## YEAR OF ONSET ANALYSIS (1940-2025) ##
monogenic_aao <- df[df$mendel_ALS_Y == 1 & !is.na(df$year_onset_ALS), ]
polygenic_aao <- df[df$polygenicY_ALS == 1 & !is.na(df$year_onset_ALS), ]
n_mono_aao <- nrow(monogenic_aao); n_poly_aao <- nrow(polygenic_aao)

p_aao <- ggplot(df[!is.na(df$year_onset_ALS), ], aes(x = year_onset_ALS)) + geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Year of ALS Onset - Affected Index Patients", x = "Year of ALS Onset", y = "Count") +
  xlim(1940, 2025) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave("onset_year_histo.pdf", plot = p_aao, width = 10, height = 6, dpi = 300)

p_mono_aao <- ggplot(monogenic_aao, aes(x = year_onset_ALS)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "red", color = "black", alpha = 0.8) +
  labs(title = "Monogenic ALS - Onset Year", subtitle = paste("n =", n_mono_aao), x = "Year of ALS Onset", y = "Percentage (%)") + 
  xlim(1940, 2025) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p_poly_aao <- ggplot(polygenic_aao, aes(x = year_onset_ALS)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Polygenic ALS - Onset Year", subtitle = paste("n =", n_poly_aao), x = "Year of ALS Onset", y = "Percentage (%)") + 
  xlim(1940, 2025) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

df_aao_combined <- rbind(
  data.frame(value = monogenic_aao$year_onset_ALS, group = "Monogenic"),
  data.frame(value = polygenic_aao$year_onset_ALS, group = "Polygenic")
)
p_aao_density <- ggplot(df_aao_combined, aes(x = value, color = group, fill = group)) +
  geom_density(aes(y = ..density.. * 100), size = 1.2, alpha = 0.3) +
#   geom_histogram(aes(y = ..density.. * 100), binwidth = 1, alpha = 0.4, position = "identity") +
  scale_color_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  scale_fill_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  labs(title = "Year of onset: Monogenic vs Polygenic Density Overlay", x = "Year of ALS Onset", y = "Percentage (%)") +
  xlim(1940, 2025) + theme_minimal() + theme(legend.position = "top")

ggsave("onset_year_densityplot.pdf", p_aao_density, width = 10, height = 6, dpi = 300)
ggsave("onset_year_mono_percent.pdf", p_mono_aao, width = 10, height = 6, dpi = 300)
ggsave("onset_year_poly_percent.pdf", p_poly_aao, width = 10, height = 6, dpi = 300)

cat("\n=== YEAR OF ALS ONSET ===\n")
cat(sprintf("Monogenic (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_mono_aao, median(monogenic_aao$year_onset_ALS, na.rm=TRUE),
            quantile(monogenic_aao$year_onset_ALS, 0.25, na.rm=TRUE), quantile(monogenic_aao$year_onset_ALS, 0.75, na.rm=TRUE)))
cat(sprintf("Polygenic  (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_poly_aao, median(polygenic_aao$year_onset_ALS, na.rm=TRUE),
            quantile(polygenic_aao$year_onset_ALS, 0.25, na.rm=TRUE), quantile(polygenic_aao$year_onset_ALS, 0.75, na.rm=TRUE)))

## AGE AT ALS ONSET ANALYSIS (20-90) ##
monogenic_age <- df[df$mendel_ALS_Y == 1 & !is.na(df$age_ALS), ]
polygenic_age <- df[df$polygenicY_ALS == 1 & !is.na(df$age_ALS), ]
n_mono_age <- nrow(monogenic_age); n_poly_age <- nrow(polygenic_age)

p_age <- ggplot(df[!is.na(df$age_ALS), ], aes(x = age_ALS)) + geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Age at ALS Onset - Affected Index Patients", x = "Age at ALS Onset (years)", y = "Count") +
  xlim(20, 90) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave("age_ALS_histo.pdf", plot = p_age, width = 10, height = 6, dpi = 300)

p_mono_age <- ggplot(monogenic_age, aes(x = age_ALS)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "red", color = "black", alpha = 0.8) +
  labs(title = "Monogenic ALS - Age at Onset", subtitle = paste("n =", n_mono_age), x = "Age at ALS Onset (years)", y = "Percentage (%)") + 
  xlim(20, 90) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p_poly_age <- ggplot(polygenic_age, aes(x = age_ALS)) + geom_histogram(aes(y = ..density.. * 100), binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Polygenic ALS - Age at Onset", subtitle = paste("n =", n_poly_age), x = "Age at ALS Onset (years)", y = "Percentage (%)") + 
  xlim(20, 90) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

df_age_combined <- rbind(
  data.frame(value = monogenic_age$age_ALS, group = "Monogenic"),
  data.frame(value = polygenic_age$age_ALS, group = "Polygenic")
)
p_age_density <- ggplot(df_age_combined, aes(x = value, color = group, fill = group)) +
  geom_density(aes(y = ..density.. * 100), size = 1.2, alpha = 0.3) +
#   geom_histogram(aes(y = ..density.. * 100), binwidth = 1, alpha = 0.4, position = "identity") +
  scale_color_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  scale_fill_manual(values = c("Monogenic" = "red", "Polygenic" = "steelblue")) +
  labs(title = "Age at onset: Monogenic vs Polygenic Density Overlay", x = "Age at ALS Onset (years)", y = "Percentage (%)") +
  xlim(20, 90) + theme_minimal() + theme(legend.position = "top")

ggsave("aao_densityplot.pdf", p_age_density, width = 10, height = 6, dpi = 300)
ggsave("aao_percent.pdf", p_mono_age, width = 10, height = 6, dpi = 300)
ggsave("aao_poly_percent.pdf", p_poly_age, width = 10, height = 6, dpi = 300)

cat("\n=== AGE AT ALS ONSET ===\n")
cat(sprintf("Monogenic (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_mono_age, median(monogenic_age$age_ALS, na.rm=TRUE),
            quantile(monogenic_age$age_ALS, 0.25, na.rm=TRUE), quantile(monogenic_age$age_ALS, 0.75, na.rm=TRUE)))
cat(sprintf("Polygenic  (n=%d): median = %.0f, IQR = %.0f-%.0f\n", n_poly_age, median(polygenic_age$age_ALS, na.rm=TRUE),
            quantile(polygenic_age$age_ALS, 0.25, na.rm=TRUE), quantile(polygenic_age$age_ALS, 0.75, na.rm=TRUE)))
