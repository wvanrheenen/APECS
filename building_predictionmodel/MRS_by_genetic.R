# COMPLETE MRS_by_genetic.R script with case counts in titles 
suppressPackageStartupMessages(suppressWarnings({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(patchwork); library(tibble); library(purrr); library(wesanderson)
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

# Case counts for titles (exactly like MRS_by_degree.R)
group_counts <- df_grouped %>% count(group) %>% deframe()
group_labels <- setNames(
  paste0(names(group_counts), "\n(n=", group_counts, ")"),
  names(group_counts)
)

# Define 4 scenario sets
disease_scenarios <- list(
  ALS = list(name = "ALS-affected relatives", scenarios = list(
    "1st degree relatives" = c("relatives_1st_als"),
    "1st and 2nd degree relatives" = c("relatives_1st_als", "relatives_2nd_als"), 
    "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_als", "relatives_2nd_als", "relatives_3rd_als")
  )),
  FTD = list(name = "FTD-affected relatives", scenarios = list(
    "1st degree relatives" = c("relatives_1st_ftd"),
    "1st and 2nd degree relatives" = c("relatives_1st_ftd", "relatives_2nd_ftd"), 
    "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_ftd", "relatives_2nd_ftd", "relatives_3rd_ftd")
  )),
  dementia = list(name = "Dementia-affected relatives", scenarios = list(
    "1st degree relatives" = c("relatives_1st_dementia"),
    "1st and 2nd degree relatives" = c("relatives_1st_dementia", "relatives_2nd_dementia"), 
    "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_dementia", "relatives_2nd_dementia", "relatives_3rd_dementia")
  )),
  ALS_FTD = list(name = "ALS or FTD-affected relatives", scenarios = list(
    "1st degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique"),
    "1st and 2nd degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique", "relatives_2nd_als", "relatives_2nd_ftd_unique"), 
    "1st, 2nd, and 3rd degree relatives" = c("relatives_1st_als", "relatives_1st_ftd_unique", "relatives_2nd_als", "relatives_2nd_ftd_unique", "relatives_3rd_als", "relatives_3rd_ftd_unique")
  ))
)

# Data prep function
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

# YOUR EXACT PLOT FUNCTION - now uses group_labels with case counts
create_architecture_mrs_plot <- function(data, group_name) {
  # Get case count for this group
  n_cases <- group_counts[[group_name]]
  title_with_count <- paste0(group_name, " (n=", n_cases, ")")
  
  ggplot(data, aes(y = scenario, x = proportion, fill = total_affected_cat_stack)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.3) +
    geom_text(aes(label = ifelse(proportion >= 0.03, percent(proportion, 0.1), "")), 
              position = position_fill(vjust = 0.5), size = 3, color = "black") +
    geom_vline(xintercept = seq(0.2, 0.8, 0.2), color = "grey85", alpha=0.5, linewidth = 0.3) +
    scale_fill_manual(
      values = c(">4" = "#DC143C", "4" = "#FF6347", "3" = "#FFA500", "2" = "#FFD700", "1" = "#90EE90", "0" = "#2E8B57"),
      breaks = levels(data$total_affected_cat_legend),
      labels = levels(data$total_affected_cat_legend),
      name = "No. affected relatives",
      guide = guide_legend(ncol = 6, title.position = "left")
    ) +
    scale_x_continuous(labels = percent_format(1L), limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = title_with_count, x = "Percent of index patients", y = NULL) +  # FIXED: case count in TITLE
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(size = 11, hjust = 0.5, face = "bold"))
}

# Generate 4 plots - GENETIC ARCHITECTURE panels
for (disease_key in names(disease_scenarios)) {
  disease_info <- disease_scenarios[[disease_key]]
  
  # Prepare data - 1st degree TOP → 3rd degree BOTTOM
  all_scenarios_data <- bind_rows(
    imap_dfr(disease_info$scenarios, ~create_mrs_data(df_grouped, .x, .y))
  ) %>%
    mutate(
      group = factor(group, levels = c("Polygenic ALS", "C9orf72-like Monogenic ALS", "FUS/SOD1-like Monogenic ALS")),
      scenario = factor(scenario, levels = rev(names(disease_info$scenarios)))
    )
  
  # Create 3 plots using YOUR function - Polygenic → C9 → FUS (with case counts)
  polygenic_plot <- create_architecture_mrs_plot(
    filter(all_scenarios_data, group == "Polygenic ALS"), 
    "Polygenic ALS"
  )
  
  c9_plot <- create_architecture_mrs_plot(
    filter(all_scenarios_data, group == "C9orf72-like Monogenic ALS"), 
    "C9orf72-like Monogenic ALS"
  )
  
  fus_plot <- create_architecture_mrs_plot(
    filter(all_scenarios_data, group == "FUS/SOD1-like Monogenic ALS"), 
    "FUS/SOD1-like Monogenic ALS"
  )
  
  # Combine VERTICALLY exactly like your original script
  combined_plot <- wrap_plots(list(polygenic_plot, c9_plot, fus_plot), ncol = 1) +
    plot_annotation(
      title = paste0("Distribution of ", disease_info$name, " in index patients"),
      subtitle = "Genetic architecture across increasing kinship distance",
      theme = theme(
        plot.title = element_text(size = 13, hjust = 0.6, face = "bold"),
        plot.subtitle = element_text(size = 11, hjust = 0.6)
      )
    )
  
  filename <- paste0("plot_c9FUS/", tolower(disease_key), "_architecture_mrs_plots.pdf")
  ggsave(filename, combined_plot, width = 12, height = 11)
  cat("✔", disease_info$name, "NEW architecture plots saved to", filename, "\n")
}
