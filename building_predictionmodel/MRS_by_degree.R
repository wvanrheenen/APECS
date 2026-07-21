# COMPLETE MRS_by_degree.R script with case counts in titles 
suppressPackageStartupMessages(suppressWarnings({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(patchwork); library(tibble); library(purrr); library(wesanderson)
}))


# Input / output setup
df <- read_csv("../results/phenocopies/combined_simulations.csv")
if (!dir.exists("mrs_relatives")) dir.create("mrs_relatives")


# Group classification and labels (unchanged)
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
  paste0(group_names[names(group_counts)], "\n(n=", scales::comma(group_counts), ")"),
  names(group_counts)
)


# Define 4 scenario sets
disease_scenarios <- list(
  ALS = list(
    name = "ALS-affected relatives",
    scenarios = list(
      "1st degree relatives" = c("relatives_1st_als"),
      "1st and 2nd degree relatives" = c("relatives_1st_als", "relatives_2nd_als"), 
      "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_als", "relatives_2nd_als", "relatives_3rd_als")
    )
  ),
  FTD = list(
    name = "FTD-affected relatives",
    scenarios = list(
      "1st degree relatives" = c("relatives_1st_ftd"),
      "1st and 2nd degree relatives" = c("relatives_1st_ftd", "relatives_2nd_ftd"), 
      "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_ftd", "relatives_2nd_ftd", "relatives_3rd_ftd")
    )
  ),
  dementia = list(
    name = "Dementia-affected relatives",
    scenarios = list(
      "1st degree relatives" = c("relatives_1st_dementia"),
      "1st and 2nd degree relatives" = c("relatives_1st_dementia", "relatives_2nd_dementia"), 
      "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_dementia", "relatives_2nd_dementia", "relatives_3rd_dementia")
    )
  ),
  ALS_FTD = list(
    name = "ALS or FTD-affected relatives",
    scenarios = list(
      "1st degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique"),
      "1st and 2nd degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique", "relatives_2nd_als", "relatives_2nd_ftd_unique"), 
      "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique", "relatives_2nd_als", "relatives_2nd_ftd_unique", "relatives_3rd_als", "relatives_3rd_ftd_unique")
    )
  )
)


# Data prep and plotting functions (modified)
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

# Plotting function
create_horizontal_mrs_plot <- function(data, scenario_name, show_legend = TRUE) {
  p <- ggplot(data, aes(y = group, x = proportion, fill = total_affected_cat_stack)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(proportion >= 0.03, percent(proportion, 0.1), "")),
      position = position_fill(vjust = 0.5),
      size = 3,
      color = "black"
    ) +
    # Legend: keep categorical labels, but ensure any numbers use commas
    scale_fill_manual(
      values = c(
        ">4" = "#DC143C", "4" = "#FF6347", "3" = "#FFA500",
        "2" = "#FFD700", "1" = "#90EE90", "0" = "#2E8B57"
      ),
      breaks = levels(data$total_affected_cat_legend),
      labels = levels(data$total_affected_cat_legend),
      name = "No. affected relatives",
      guide = guide_legend(ncol = 6, title.position = "left")
    ) +
    # X-axis: no expansion, exactly 0–1 (0–100%)
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    scale_y_discrete(labels = group_labels) +
    labs(
      title = scenario_name,
      x = "Percent of index patients",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background   = element_blank(),
      panel.border       = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.ticks.y       = element_blank()
    )

  if (show_legend) {
    p <- p + theme(legend.position = "bottom")
  } else {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}


# Generate 4 plots
for (disease_key in names(disease_scenarios)) {
  disease_info <- disease_scenarios[[disease_key]]
  scenario_names <- names(disease_info$scenarios)
  
  # Prepare data for this disease
  all_scenarios_data <- bind_rows(
    imap_dfr(disease_info$scenarios, ~create_mrs_data(df_grouped, .x, .y))
  ) %>%
    mutate(group = factor(group, levels = c("Rare Monogenic ALS", "Common Monogenic ALS", "Polygenic ALS")))
  
  # Create individual scenario plots - only bottom plot gets legend
  plots <- lapply(seq_along(scenario_names), function(i) {
    show_legend <- (scenario_names[i] == "1st, 2nd, and 3rd degree relatives")
    create_horizontal_mrs_plot(filter(all_scenarios_data, scenario == scenario_names[i]), 
                               scenario_names[i], 
                               show_legend = show_legend)
  })
  
  # Combine into single figure
  combined_plot <- wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = paste0("Figure 4: Distribution of ", disease_info$name, " in index patients"),
      subtitle = "Proportion of affected relatives across increasing kinship distance",
      theme = theme(
        plot.title = element_text(size = 13, hjust = 0.6, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.6)
      )
    )
  
  # Save
  filename <- paste0("mrs_relatives/", tolower(disease_key), "_degree_mrs_plots.pdf")
  ggsave(filename, combined_plot, width = 12, height = 11)
  cat("✔", disease_info$name, "plots saved to", filename, "\n")
}