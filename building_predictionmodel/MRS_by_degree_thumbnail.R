suppressPackageStartupMessages(suppressWarnings({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(patchwork); library(tibble); library(purrr); library(wesanderson)
}))

# Input / output setup
df <- read_csv("../results/phenocopies/combined_simulations.csv")
if (!dir.exists("mrs_relatives")) dir.create("mrs_relatives")

# Group classification
df_grouped <- df %>%
  mutate(
    group = case_when(
      polygenicY_ALS == 1 & mendel_ALS_Y != 1 ~ "Polygenic ALS",
      mendel_ALS_Y_common == 1 ~ "Common Monogenic ALS",
      mendel_ALS_Y_patho == 1 ~ "Rare Monogenic ALS"
    )
  ) %>%
  filter(!is.na(group))

# Case counts for y-axis labels
group_counts <- df_grouped %>%
  count(group) %>%
  deframe()

group_names <- c(
  "Polygenic ALS" = "Polygenic ALS",
  "Common Monogenic ALS" = "Moderate penetrance\nMonogenic ALS",
  "Rare Monogenic ALS" = "High penetrance\nMonogenic ALS"
)

group_labels <- setNames(
  paste0(
    group_names[names(group_counts)],
    "\n(n=", scales::comma(group_counts), ")"
  ),
  names(group_counts)
)

# Disease definitions: 1st + 2nd degree relatives only
disease_scenarios <- list(
  ALS_FTD = list(
    name = "ALS or FTD-affected relatives",
    cols = c(
      "relatives_1st_als", "relatives_1st_ftd_unique",
      "relatives_2nd_als", "relatives_2nd_ftd_unique"
    )
  )
)

# Data preparation
create_mrs_data <- function(df, cols) {
  df %>%
    rowwise() %>%
    mutate(total_affected = sum(c_across(all_of(cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      total_affected_cat = case_when(
        total_affected >= 5 ~ 5,
        TRUE ~ total_affected
      ),
      total_affected_cat_stack = factor(
        total_affected_cat,
        levels = 5:0,
        labels = c(">4", "4", "3", "2", "1", "0")
      ),
      total_affected_cat_legend = factor(
        total_affected_cat,
        levels = 0:5,
        labels = c("0", "1", "2", "3", "4", ">4")
      )
    ) %>%
    count(
      group,
      total_affected_cat,
      total_affected_cat_stack,
      total_affected_cat_legend
    ) %>%
    group_by(group) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
}

# Plotting function
create_horizontal_mrs_plot <- function(data, title) {
  ggplot(
    data,
    aes(y = group, x = proportion, fill = total_affected_cat_stack)
  ) +
    geom_col(width = 0.8, color = "black", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(proportion >= 0.03, percent(proportion, 0.1), "")),
      position = position_fill(vjust = 0.5),
      size = 4,
      color = "black"
    ) +
    scale_fill_manual(
      values = setNames(
        RColorBrewer::brewer.pal(n = 6, name = "RdBu"),
        c(">4", "4", "3", "2", "1", "0")
      ),
      breaks = c("0", "1", "2", "3", "4", ">4"),
      labels = c("0", "1", "2", "3", "4", ">4"),
      name = "No. affected relatives",
      guide = guide_legend(ncol = 6, title.position = "left")
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    scale_y_discrete(labels = group_labels) +
    labs(
      title = title,
      x = "Percent of index patients",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 14,
      hjust = 0.5,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    legend.title = element_text(
      color = "black"
    ),
    legend.text = element_text(
      color = "black"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )
}

# Create and save one 1st + 2nd-degree plot per disease definition
for (disease_key in names(disease_scenarios)) {
  disease_info <- disease_scenarios[[disease_key]]

  mrs_data <- create_mrs_data(
    df = df_grouped,
    cols = disease_info$cols
  ) %>%
    mutate(
      group = factor(
        group,
        levels = c(
          "Rare Monogenic ALS",
          "Common Monogenic ALS",
          "Polygenic ALS"
        )
      )
    )

  plot_title <- paste0(
    "Distribution of ", disease_info$name,
    " in ALS index patients\n(1st and 2nd degree relatives)"
  )

  p <- create_horizontal_mrs_plot(mrs_data, plot_title)

  filename <- file.path(
    "mrs_relatives",
    "als_ftd_1st_2nd_degree_thumbnail_plot.pdf"
  )

  ggsave(
    filename,
    p,
    width = 8.5,
    height = 5.5
  )

  cat("✔", disease_info$name, "plot saved to", filename, "\n")
}