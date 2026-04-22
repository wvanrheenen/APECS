#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(scales)
library(wesanderson)  # Wes Anderson color palettes

# Access Snakemake objects
sets <- snakemake@params$sets
input_files <- snakemake@input
output_file <- snakemake@output[[1]]

# Create output directory
output_dir <- dirname(output_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read and combine all prediction model files across sets
metrics <- map_dfr(seq_along(input_files), function(i) {
  set_name <- sets[i]
  file_path <- input_files[[i]]
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Filter for the specific fALS scenario
  df_scen <- df %>% 
    filter(Scenario == "4.2) ≥2 1st/2nd ALS/FTD")
  
  if (nrow(df_scen) == 0) return(NULL)
  
  # Add set identifier
  df_scen %>%
    mutate(set = set_name)
})

# Check if we found data
if (nrow(metrics) == 0) {
  stop("No data found for scenario '4.2) ≥2 1st/2nd ALS/FTD' across all sets")
}

# Prepare data for plotting (pivot to long format)
metrics_long <- metrics %>%
  select(set, Sensitivity, Specificity, PPV, NPV) %>%
  pivot_longer(
    cols = c(Sensitivity, Specificity, PPV, NPV),
    names_to = "Metric", 
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, 
                    levels = c("Sensitivity", "Specificity", "PPV", "NPV")),
    # Convert sets to numeric for smoother lines (assuming chronological order)
    set_num = as.numeric(factor(set, levels = sets))
  )

# Identify key timepoints for labeling (1st, 4th, 7th positions for 2005, 2035, 2065)
key_positions <- c(1, 4, 7)  # Indices for 2005, 2035, 2065 in your SETS

# Darjeeling1 color palette (beautiful muted tones)
darjeeling_colors <- wes_palette("Darjeeling1")[1:4]

# Create the line plot with Wes Anderson colors
p <- ggplot(metrics_long, aes(x = set_num, y = Value, color = Metric)) +
  geom_line(linewidth = 1.4, alpha = 0.9) +
  geom_point(size = 3.5, alpha = 1, stroke = 0.3) +
  # Add labels only at key timepoints (2005, 2035, 2065)
  geom_text(
    data = . %>% filter(set_num %in% key_positions),
    aes(label = scales::percent(Value, accuracy = 0.1)),  # 1 decimal place (89.2% → 89.2%)
    vjust = -1.3,        # HIGHER: -0.8 → -1.2 (more negative = higher up)
    hjust = 0.5, 
    size = 3.2, 
    fontface = "bold",
    color = "black", 
    stroke = 0.3, 
    family = "sans"
  ) +
  scale_x_continuous(
    breaks = 1:length(sets),
    labels = sets,
    name = "Year"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1), 
    limits = c(0.20, 1),  # 20% to 100%
    expand = expansion(mult = c(0, 0.02))  # Small padding at top, none at bottom
  ) +
  scale_color_manual(
    values = darjeeling_colors,
    name = "Metric"
  ) +
  labs(
    title = "fALS Prediction Metrics: ≥2 ALS/FTD relatives within two degrees",
    subtitle = "Metrics simulated for index patients aged 25-85 years old",
    x = "Year",
    y = "Metric Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.key.width = unit(1.5, "cm")
  )

# Save plot
ggsave(output_file, p, width = 11, height = 7.5, dpi = 300)  # Slightly taller for labels
cat("Metrics line plot saved to:", output_file, "\n")
cat("Using Wes Anderson 'Darjeeling1' color palette\n")
cat("Labels added at 2005, 2035, and 2065 timepoints\n")
cat("Data summary:\n")
print(metrics_long %>% 
      group_by(set, Metric) %>% 
      summarise(mean_value = mean(Value), .groups = "drop"))
