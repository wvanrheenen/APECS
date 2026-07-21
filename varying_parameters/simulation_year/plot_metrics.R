#!/usr/bin/env Rscript

# =========================================================
# Standalone interactive plotting script
# Styled to match poster figures
# =========================================================

library(tidyverse)
library(ggplot2)
library(scales)

# ---------------------------------------------------------
# DEFINE SETS
# ---------------------------------------------------------

sets <- c(
  "1985",
  "1995",
  "2005",
  "2015",
  "2025",
  "2035",
  "2045",
  "2055",
  "2065"
)

# ---------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------

input_files <- file.path(
  "results",
  sets,
  "aggregate_prediction_model.csv"
)

output_file <- "results/plots/fALS_metrics_poster_style.pdf"

# ---------------------------------------------------------
# LABELS
# ---------------------------------------------------------

title <- "Simulation censoring year sensitivity analysis"
x_label <- "Simulation year"

# ---------------------------------------------------------
# POSTER COLORS
# ---------------------------------------------------------

metric_colors <- c(
  "PPV" = "#7CC8A3",
  "Sensitivity" = "#f39100"
)

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

plot_data <- map2_dfr(
  input_files,
  sets,
  function(file_path, set_name) {

    if (!file.exists(file_path)) {
      stop(paste("Missing file:", file_path))
    }

    cat("Reading:", file_path, "\n")

    read_csv(
      file_path,
      show_col_types = FALSE
    ) %>%
      mutate(set = set_name)
  }
)

# ---------------------------------------------------------
# PREPARE PLOTTING DATA
# ---------------------------------------------------------

main_metrics <- plot_data %>%
  select(
    set,
    mean_Sensitivity,
    mean_PPV
  ) %>%
  pivot_longer(
    -set,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      grepl("Sensitivity", Metric) ~ "Sensitivity",
      grepl("PPV", Metric) ~ "PPV"
    ),

    Metric = factor(
      Metric,
      levels = c("Sensitivity", "PPV")
    ),

    set_num = as.numeric(
      factor(set, levels = sets)
    )
  )

# ---------------------------------------------------------
# CREATE PLOT
# ---------------------------------------------------------

p <- ggplot(
  main_metrics,
  aes(
    x = set_num,
    y = Value,
    color = Metric,
    group = Metric
  )
) +

  # -------------------------------------------------------
  # ERROR BARS
  # -------------------------------------------------------

  geom_errorbar(
    data = plot_data,
    aes(
      x = match(set, sets),
      ymin = pmax(mean_Sensitivity - sd_Sensitivity, 0),
      ymax = pmin(mean_Sensitivity + sd_Sensitivity, 1),
      color = "Sensitivity"
    ),
    width = 0.15,
    linewidth = 0.6,
    alpha = 0.7,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +

  geom_errorbar(
    data = plot_data,
    aes(
      x = match(set, sets),
      ymin = pmax(mean_PPV - sd_PPV, 0),
      ymax = pmin(mean_PPV + sd_PPV, 1),
      color = "PPV"
    ),
    width = 0.15,
    linewidth = 0.6,
    alpha = 0.7,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +

  # -------------------------------------------------------
  # LINES
  # -------------------------------------------------------

  geom_line(
    linewidth = 0.6
  ) +

  # -------------------------------------------------------
  # POINTS
  # -------------------------------------------------------

  geom_point(
    size = 3,
    stroke = 0.4
  ) +

  # -------------------------------------------------------
  # DIRECT LINE LABELS
  # -------------------------------------------------------

  geom_text(
    data = main_metrics %>%
      group_by(Metric) %>%
      slice_max(set_num, n = 1),

    aes(label = Metric),

    size = 2.5,
    fontface = "bold",

    show.legend = FALSE
  ) +

  # -------------------------------------------------------
  # X AXIS
  # -------------------------------------------------------

  scale_x_continuous(
    breaks = 1:length(sets),
    labels = sets,
    name = x_label
  ) +

  # -------------------------------------------------------
  # Y AXIS
  # -------------------------------------------------------

  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.2, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.08))
  ) +

  # -------------------------------------------------------
  # COLORS
  # -------------------------------------------------------

  scale_color_manual(
    values = metric_colors
  ) +

  # -------------------------------------------------------
  # LABELS
  # -------------------------------------------------------

  labs(
    x = NULL,
    y = "fALS accuracy",
    color = NULL
  ) +

  # -------------------------------------------------------
  # THEME
  # -------------------------------------------------------

  theme_bw(base_size = 10) +

  theme(

    # Remove title
    plot.title = element_blank(),

    # Axis text
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    # Axis titles
    axis.title.x = element_blank(),

    axis.title.y = element_text(
      color = "black",
      size = 11
    ),

    # Remove legend
    legend.position = "none",

    # Grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    # Background
    panel.background = element_blank(),
    panel.border = element_blank(),

    # Axis ticks
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

# ---------------------------------------------------------
# DISPLAY
# ---------------------------------------------------------

print(p)

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  output_file,
  p,
  width = 5,
  height = 2,
  dpi = 300
)

cat("\nSaved plot to:\n")
cat(output_file, "\n")