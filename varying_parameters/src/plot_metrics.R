#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(scales)
library(wesanderson)

# Access Snakemake objects
sets <- snakemake@params$sets
input_files <- snakemake@input
output_file <- snakemake@output[[1]]

# Get current wd (i.e. variable name)
current_dir <- basename(getwd())

# Dynamic title and x_label mapping
labels <- list(
  "age_onset_monogenic" = list(title = "(C) Varying onset of monogenic ALS", x_label = "Shift in onset of monogenic ALS (years)"),
  "age_onset_polygenic" = list(title = "(D) Varying onset of polygenic ALS", x_label = "Shift in onset of polygenic ALS (years)"),
  "c9_daf" = list(title = "(E) Varying ALS-FTD-common disease allele frequency", x_label = "ALS-FTD-common disease allele frequency"),
  "c9_penetrance" = list(title = "(F) Varying ALS-FTD-common disease allele penetrance", x_label = "ALS-FTD-common disease allele penetrance"),
  "fussod_daf" = list(title = "(G) Varying ALS-FTD-rare disease allele frequency", x_label = "ALS-FTD-rare disease allele frequency"),
  "fussod_penetrance" = list(title = "(H) Varying ALS-FTD-rare disease allele penetrance", x_label = "ALS-FTD-rare disease allele penetrance"),
  "fert_rate" = list(title = "(A) Varying fertility rate", x_label = "Fertility rate per generation"),
  "genetic_correlation" = list(title = "(K) Varying genetic correlation between ALS~FTD", x_label = "Genetic correlation ALS~FTD"),
  "heritability_ALS" = list(title = "(J) Varying ALS heritability", x_label = "ALS heritability (h2)"),
  "life_exp" = list(title = "(B) Varying mean life expectancy", x_label = "Mean life expectancy per individual"),
  "lifetime_risk_ALS" = list(title = "(I) Varying ALS lifetime risk", x_label = "ALS lifetime risk (K)"),
  "simulation_year" = list(title = "(L) Varying censoring years when running simulation", x_label = "Simulation censoring year")
)

# Get labels for current directory (with fallback)
label_info <- labels[[current_dir]]
if (is.null(label_info)) {
  title <- paste("Varying", current_dir)
  x_label <- paste(current_dir, "Value")
} else {
  title <- label_info$title
  x_label <- label_info$x_label
}

# Create output directory
output_dir <- dirname(output_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read and process each file (keep all mean/sd columns)
plot_data <- map2_dfr(input_files, sets, function(file_path, set_name) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df %>% mutate(set = set_name)
})

# Check if we found data
if (nrow(plot_data) == 0) {
  stop("No data found across all sets")
}

cat("Plotting data for", nrow(plot_data), "sets\n")
print(plot_data)

# Create plotting data for main metrics
main_metrics <- plot_data %>%
  select(set, mean_Sensitivity, mean_PPV, mean_Specificity, mean_NPV, mean_Accuracy) %>%
  pivot_longer(-set, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = case_when(
      grepl("Sensitivity", Metric) ~ "Sensitivity",
      grepl("PPV", Metric) ~ "PPV", 
      TRUE ~ NA_character_
    ),
    Metric = factor(Metric, levels = c("Sensitivity", "PPV")),
    set_num = as.numeric(factor(set, levels = sets))
  ) %>%
  filter(!is.na(Metric))

# Colors (only for Sensitivity and PPV since that's what error bars show)
darjeeling_colors <- wes_palette("Darjeeling1")[c(1, 2)]

# Create the plot
p <- ggplot(main_metrics, aes(x = set_num, y = Value, color = Metric)) +
  # Error bars for Sensitivity using plot_data directly
  geom_errorbar(
    data = plot_data,
    aes(x = match(set, sets), 
        ymin = pmax(mean_Sensitivity - sd_Sensitivity, 0), 
        ymax = pmin(mean_Sensitivity + sd_Sensitivity, 1),
        color = "Sensitivity"),
    width = 0.2, linewidth = 1, alpha = 0.7,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Error bars for PPV using plot_data directly
  geom_errorbar(
    data = plot_data,
    aes(x = match(set, sets), 
        ymin = pmax(mean_PPV - sd_PPV, 0), 
        ymax = pmin(mean_PPV + sd_PPV, 1),
        color = "PPV"),
    width = 0.2, linewidth = 1, alpha = 0.7,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Main points and lines (only Sensitivity and PPV will have error bars)
  geom_line(linewidth = 1) +
  geom_point(size = 2, stroke = 0.5) +
  # Value labels
  geom_text(
    aes(label = scales::percent(Value, accuracy = 0.1)),
    nudge_y = 0.04,           
    vjust = -1.5,                
    hjust = 0.5,
    size = 2.0, color = "black",
    family = "sans", check_overlap = TRUE
  ) + 
  scale_x_continuous(
    breaks = 1:length(sets),
    labels = sets,
    name = x_label
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1), 
    limits = c(0.0, 1.0),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_color_manual(
    values = darjeeling_colors, 
    breaks = c("Sensitivity", "PPV")
  ) +
  labs(
    title = title,
    x = x_label,
    y = "Predictive Metric Value",
    color = "Metric"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0),
    plot.margin = margin(5, 5, 5, 5, "pt")    
  )

# Save plot
ggsave(output_file, p, width = 3, height = 3, dpi = 300)
cat("Metrics plot saved to:", output_file, "\n")