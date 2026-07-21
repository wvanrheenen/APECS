suppressPackageStartupMessages(suppressWarnings({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(tibble)
}))

# Input / output setup
df <- read_csv("../results/testset/combined_simulations.csv")
if (!dir.exists("mrs_relatives")) dir.create("mrs_relatives")

# Group classification and labels
df_grouped <- df %>%
  mutate(
    group = case_when(
      (polygenicY_ALS == 1 & mendel_ALS_Y != 1) ~ "Polygenic ALS",
      mendel_ALS_Y_common == 1 ~ "Common Monogenic ALS",
      mendel_ALS_Y_patho == 1 ~ "Rare Monogenic ALS"
    )
  ) %>%
  filter(!is.na(group))

group_counts <- df_grouped %>% count(group) %>% deframe()

group_names <- c(
  "Polygenic ALS" = "Polygenic ALS",
  "Common Monogenic ALS" = "Common\nMonogenic ALS",
  "Rare Monogenic ALS" = "Rare, pathogenic\nMonogenic ALS"
)

group_labels <- setNames(
  paste0(group_names[names(group_counts)], "\n(n=", group_counts, ")"),
  names(group_counts)
)

# Only 3rd-degree relatives
disease_scenarios <- list(
  ALS_FTD = list(
    name = "ALS or FTD-affected relatives up to 2nd degree",
    cols = c("relatives_1st_als", "relatives_1st_ftd_unique", "relatives_2nd_als", "relatives_2nd_ftd_unique")
  )
)

create_mrs_data <- function(df, cols, scenario) {
  df %>%
    rowwise() %>%
    mutate(total_affected = sum(c_across(all_of(cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      total_affected_cat = case_when(total_affected >= 5 ~ 5, TRUE ~ total_affected),
      scenario = scenario,
      total_affected_cat_stack = factor(total_affected_cat, levels = 5:0, labels = c(">4", "4", "3", "2", "1", "0")),
      total_affected_cat_legend = factor(total_affected_cat, levels = 0:5, labels = c("0", "1", "2", "3", "4", ">4"))
    ) %>%
    count(group, total_affected_cat, total_affected_cat_stack, total_affected_cat_legend, scenario) %>%
    group_by(group, scenario) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
}

# Updated 6-color theme matching your poster
mrs_colors <- c(
  "0"  = "#EAF5F4",  # Light blue (your first color)
  "1"  = "#ADD8D7",  # Light cyan transition
  "2"  = "#489F9D",  # Your second color (teal) 
  "3"  = "#38A7EA",  # Your third color (blue)
  "4"  = "#016FB7",  # Your fourth color (dark blue)
  ">4" = "#1F4E79"   # Darker blue completion
)

create_horizontal_mrs_plot <- function(data, scenario_name, show_legend = TRUE) {
  p <- ggplot(data, aes(y = group, x = proportion, fill = total_affected_cat_stack)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(proportion >= 0.07, percent(proportion, 0.1), "")),
      position = position_fill(vjust = 0.5),
      size = 4.5,
      color = "black"
    ) +
    scale_fill_manual(
      values = mrs_colors,
      breaks = levels(data$total_affected_cat_legend),
      labels = levels(data$total_affected_cat_legend),
      name = "No. affected relatives",
      guide = guide_legend(ncol = 6, title.position = "left")
    ) +
    scale_x_continuous(labels = percent_format(1L), limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_discrete(labels = group_labels) +
    labs(
      title    = scenario_name,
      x        = "Percent of index patients",
      y        = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title   = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
      plot.subtitle = element_text(size = 16, color = "black", hjust = 0.5),  # 16 pt subtitle
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(size = 14, margin = margin(r = -20), color = "black"),
      axis.title.x = element_text(color = "black"),
      legend.text  = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background   = element_blank(),
      panel.border       = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.ticks.y       = element_blank()
    )

  if (!show_legend) {
    p <- p + theme(legend.position = "top")
  }

  p
}

# Generate one plot per disease
for (disease_key in names(disease_scenarios)) {
  disease_info <- disease_scenarios[[disease_key]]

  plot_data <- create_mrs_data(df_grouped, disease_info$cols, "2nd degree relatives") %>%
    mutate(group = factor(group, levels = c("Rare Monogenic ALS", "Common Monogenic ALS", "Polygenic ALS")))

  p <- create_horizontal_mrs_plot(
    plot_data,
    paste0(disease_info$name),
    show_legend = FALSE
  )

  # Updated dimensions: 75.1 cm wide, 34.3 cm high
  filename <- paste0("mrs_relatives/", tolower(disease_key), "_2nd_degree_ESHG_plot.pdf")
  ggsave(filename, p, width = 53.5/2.54, height = 34.5/2.54, units = "cm", device = "pdf")
  cat("✔", disease_info$name, "plot saved to", filename, "\n")
}