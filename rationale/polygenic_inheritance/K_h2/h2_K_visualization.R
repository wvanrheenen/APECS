library(readr)
library(dplyr)
library(ggplot2)
library(wesanderson)

outcome_dirs <- c("h2_als", "K_als")
varying_params <- c("varying_h2_als", "varying_K_als", "varying_h2_ftd", 
                    "varying_K_ftd", "varying_rg", "varying_lambda", "varying_gen")

read_combined_files <- function(outcome_dir, varying_param) {
  file_path <- file.path("results", outcome_dir, paste0(varying_param, ".all"))
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  read_tsv(file_path, show_col_types = FALSE) %>%
    mutate(outcome = outcome_dir, varying_param = varying_param)
}

df_list <- lapply(varying_params, function(param) {
  lapply(outcome_dirs, function(outcome) read_combined_files(outcome, param)) %>% bind_rows()
})
names(df_list) <- varying_params

param_to_column <- list(
  varying_h2_als = "h2_als",
  varying_K_als = "K_als",
  varying_h2_ftd = "h2_ftd",
  varying_K_ftd = "K_ftd",
  varying_rg = "rg",
  varying_lambda = "lambda",
  varying_gen = "gen"
)

# Custom labeller to prepend "K: " for facets of K_als (adjust as needed)
k_als_labeller <- function(x) paste0("K: ", x)

for (param_name in names(df_list)) {
  df <- df_list[[param_name]] %>%
    filter(!is.na(get(param_to_column[[param_name]]))) %>%
    mutate(
      h2_als = factor(h2_als),
      K_als = factor(K_als),
      rg = factor(rg),
      lambda = factor(lambda)
    )
  
  cat(sprintf("Plotting results for %s\n", param_name))
  
  varying_col <- param_to_column[[param_name]]
  df <- df %>% mutate(varying_value = factor(df[[varying_col]]))
  
  # Prepare expected segments for dashed line (replace expected points)
  expected_segments_K <- df %>%
    arrange(varying_value, K_als) %>%
    group_by(varying_value, K_als) %>%
    summarize(
      y = unique(K_als), # numeric for expected K_als
      .groups = "drop"
    ) %>%
    mutate(
      x = as.numeric(K_als) - 0.5,
      xend = as.numeric(K_als) + 0.5,
      yend = y
    )
  
  p1 <- ggplot(df, aes(x = K_als, y = K_empirical, fill = varying_value)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_segment(data = expected_segments_K,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linetype = "dashed",
                 inherit.aes = FALSE) +
    facet_wrap(~ varying_value, scales = "free_x") +
    scale_fill_manual(values = wes_palette("Darjeeling1", n = nlevels(df$varying_value)),
                      name = varying_col) +
    labs(
      title = paste("Simulated vs Expected Prevalence for", param_name),
      fill = varying_col,
      x = "K ALS",
      y = "Simulated Prevalence"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  png(paste0("K_", param_name, ".png"), width = 7, height = 5, units = "in", res = 300)
  print(p1)
  dev.off()
  
  # Prepare expected segments for h2_als plot
  expected_segments_h2 <- df %>%
    arrange(varying_value, h2_als) %>%
    group_by(varying_value, h2_als) %>%
    summarize(
      y = unique(h2_empirical), 
      .groups = "drop"
    ) %>%
    mutate(
      x = as.numeric(h2_als) - 0.5,
      xend = as.numeric(h2_als) + 0.5,
      yend = y
    )
  
  p2 <- ggplot(df, aes(x = h2_als, y = h2_empirical, fill = varying_value)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_segment(data = expected_segments_h2,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linetype = "dashed",
                 inherit.aes = FALSE) +
    facet_wrap(~ varying_value, scales = "free_x") +
    scale_fill_manual(values = wes_palette("Darjeeling1", n = nlevels(df$varying_value)),
                      name = varying_col) +
    labs(
      title = paste("Simulated vs Expected Heritability for", param_name),
      fill = varying_col,
      x = "Heritability (h2 ALS)",
      y = "Simulated Heritability"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  png(paste0("h2_", param_name, ".png"), width = 7, height = 5, units = "in", res = 300)
  print(p2)
  dev.off()
}
