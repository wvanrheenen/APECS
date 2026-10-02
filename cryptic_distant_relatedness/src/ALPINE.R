#!/usr/bin/env Rscript
# analyze_alpine.R - Updated to save filtered sporadic monogenic patients

# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Get command line argument (file path)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript analyze_alpine.R results/alpine/combined_simulations_ALPINE.csv")
}
input_file <- args[1]

# Read the data
cat("Reading:", input_file, "\n")
df <- read_csv(input_file, show_col_types = FALSE)

# Create SPORADIC_Monogenic column
df <- df %>%
  mutate(
    SPORADIC_Monogenic = (mendel_ALS_Y == 1 & SPORADIC == 1) * 1,
    SPORADIC_Polygenic = (polygenicY_ALS == 1 & SPORADIC == 1) * 1
  )

cat("Loaded", nrow(df), "rows and", ncol(df), "columns\n")
cat("Created SPORADIC_Monogenic column\n")

## 1. Summary Statistics (UPDATED)
cat("\n=== ALPINE CASE SUMMARY ===\n")
cat("Total cases:", nrow(df), "\n")
cat("Monogenic ALS cases (mendel_ALS_Y == 1):", sum(df$mendel_ALS_Y == 1, na.rm = TRUE), "\n")
cat("Polygenic ALS cases (polygenicY_ALS == 1):", sum(df$polygenicY_ALS == 1, na.rm = TRUE), "\n")
cat("Sporadic within 3rd degree:", sum(df$SPORADIC, na.rm = TRUE), "\n")
cat("Sporadic Monogenic cases:", sum(df$SPORADIC_Monogenic, na.rm = TRUE), "\n")
cat("Sporadic Polygenic cases:", sum(df$SPORADIC_Polygenic, na.rm = TRUE), "\n")
cat("ALPINE 3rd degree:", sum(df$ALPINE_3rd, na.rm = TRUE), "\n")
cat("ALPINE 4th degree:", sum(df$ALPINE_4th, na.rm = TRUE), "\n")
cat("ALPINE 5th degree:", sum(df$ALPINE_5th, na.rm = TRUE), "\n")
cat("ALPINE 4th/5th degree:", sum(df$ALPINE_4th_5th, na.rm = TRUE), "\n")
cat("ALPINE any >3rd degree:", sum(df$ALPINE_any, na.rm = TRUE), "\n")

## 2. ALPINE Rates Among Sporadic Monogenic
sporadic_mono <- df[df$SPORADIC_Monogenic == 1, ]
if (nrow(sporadic_mono) > 0) {
  n_sporadic_mono <- nrow(sporadic_mono)
  cat("\n=== RATES AMONG SPORADIC MONOGENIC (", n_sporadic_mono, " cases (", round((n_sporadic_mono/sum(df$mendel_ALS_Y == 1, na.rm = TRUE))*100, 2), "% of all monogenic)) ===\n", sep="")  
  cat("ALPINE_3rd rate:", round(mean(df$ALPINE_3rd[df$SPORADIC_Monogenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_4th rate:", round(mean(df$ALPINE_4th[df$SPORADIC_Monogenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_5th rate:", round(mean(df$ALPINE_5th[df$SPORADIC_Monogenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_4th_5th rate:", round(mean(df$ALPINE_4th_5th[df$SPORADIC_Monogenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_any rate:", round(mean(df$ALPINE_any[df$SPORADIC_Monogenic == 1], na.rm=TRUE), 4), "\n")
  
  # **NEW: Save FILTERED sporadic monogenic patients**
  write_csv(sporadic_mono, "sporadic_monogenic_counts.csv")
  cat("✓ Filtered sporadic monogenic patients saved to: sporadic_monogenic_counts.csv (", n_sporadic_mono, " rows)\n", sep="")
} else {
  cat("No sporadic monogenic cases found - skipping filtered export\n")
}

## 3. ALPINE Rates Among Sporadic Polygenic
sporadic_poly <- df[df$SPORADIC_Polygenic == 1, ]
if (nrow(sporadic_poly) > 0) {
  n_sporadic_poly <- nrow(sporadic_poly)
  cat("\n=== RATES AMONG SPORADIC POLYGENIC (", n_sporadic_poly, " cases (", round((n_sporadic_poly/sum(df$polygenicY_ALS == 1, na.rm = TRUE))*100, 2), "% of all polygenic)) ===\n", sep="")  
  cat("ALPINE_3rd rate:", round(mean(df$ALPINE_3rd[df$SPORADIC_Polygenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_4th rate:", round(mean(df$ALPINE_4th[df$SPORADIC_Polygenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_5th rate:", round(mean(df$ALPINE_5th[df$SPORADIC_Polygenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_4th_5th rate:", round(mean(df$ALPINE_4th_5th[df$SPORADIC_Polygenic == 1], na.rm=TRUE), 4), "\n")
  cat("ALPINE_any rate:", round(mean(df$ALPINE_any[df$SPORADIC_Polygenic == 1], na.rm=TRUE), 4), "\n")
  
  # **NEW: Save FILTERED sporadic polygenic patients**
  write_csv(sporadic_poly, "sporadic_polygenic_counts.csv")
  cat("✓ Filtered sporadic polygenic patients saved to: sporadic_polygenic_counts.csv (", n_sporadic_poly, " rows)\n", sep="")
} else {
  cat("No sporadic polygenic cases found - skipping filtered export\n")
}

## 4. Save detailed summary table (for full dataset)
summary_table <- df %>%
  summarise(
    total_cases = n(),
    monogenic_ALS = sum(mendel_ALS_Y == 1, na.rm = TRUE),
    polygenic_ALS = sum(polygenicY_ALS == 1, na.rm = TRUE),
    sporadic_3rd = sum(SPORADIC, na.rm = TRUE),
    Sporadic_Monogenic = sum(SPORADIC_Monogenic, na.rm = TRUE),
    Sporadic_Polygenic = sum(SPORADIC_Polygenic, na.rm = TRUE),
    ALPINE_3rd = sum(ALPINE_3rd, na.rm = TRUE),
    ALPINE_4th = sum(ALPINE_4th, na.rm = TRUE),
    ALPINE_5th = sum(ALPINE_5th, na.rm = TRUE),
    ALPINE_4th_5th = sum(ALPINE_4th_5th, na.rm = TRUE),
    ALPINE_any = sum(ALPINE_any, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "count")

output_file <- gsub(".csv$", "_analysis_summary.csv", input_file)
write_csv(summary_table, output_file)
cat("Detailed summary saved to:", output_file, "\n")

