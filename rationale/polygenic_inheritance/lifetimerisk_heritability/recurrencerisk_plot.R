#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(wesanderson)
  library(RColorBrewer)
})

# ---- argument parsing ----
parser <- ArgumentParser(description = "Recurrence risk λ₁ validation plots")
parser$add_argument("-i", "--input", type = "character", default = "polygenic_inheritance_all.tsv")
parser$add_argument("-p", "--plotdir", type = "character", default = "plots")
args <- parser$parse_args()

# ---- data ----
df <- read_tsv(args$input, show_col_types = FALSE) %>%
  mutate(
    lambda1_emp = as.numeric(lambda1_emp),
    lambda1_theor = as.numeric(lambda1_theor),
    lambda1_theor = round(lambda1_theor, 2)
  ) %>%
  filter(!is.na(lambda1_emp) & !is.na(lambda1_theor))

dir.create(file.path(args$plotdir, "lambda1"), recursive = TRUE, showWarnings = FALSE)

# ---- label mapping for fill variables ----
label_map <- c(
  "K_ALS_input"          = "Lifetime risk ALS",
  "rg_ALSFTD_input"      = "Correlation ALS~FTD",
  "h2_ALS_input"         = "h² ALS",
  "lambda"               = "Fertility rate"
)
# fill_options will be names of the ones you actually want to loop over
fill_options <- c("lambda", "K_FTD_input", "h2_FTD_input", "rg_ALSFTD_input")


# Custom labeller function using label_map
nice_facet_label <- function(variable, value) {
  var_label <- if (variable %in% names(label_map)) label_map[variable] else variable
  paste0(var_label, " = ", value)
}

rdbu_cols <- brewer.pal(4, "RdBu")[c(1, 3, 4)]

for(fill_var in fill_options) {
  # K_ALS row (1 row × 3 cols)
  df_kals <- df %>%
    mutate(x_label = paste0("RR=", lambda1_theor, "\n(h2 ALS=", h2_ALS_input, ")"))
  
  p_kals <- ggplot(df_kals, aes(x = factor(lambda1_theor), y = lambda1_emp, fill = factor(.data[[fill_var]]))) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 1.0,
                 outlier.shape = 16, outlier.size = 1) +
    facet_wrap(~ K_ALS_input, nrow = 1, ncol = 3, scales = "free",
               labeller = nice_facet_label) +
    scale_fill_manual(values = rdbu_cols) +
    scale_x_discrete(labels = function(x) unique(df_kals$x_label[as.character(df_kals$lambda1_theor) == x])) +
    labs(title = "(B) Polygenic Recurrence Risk Validation",
         x = "Theoretical Recurrence Risk", y = "Observed RR",
         fill = label_map[fill_var]) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.margin = margin(b = 2),
      axis.text.x = element_text(color = "black", size = 7),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      plot.margin = unit(c(5,5,5,5), "pt")
    )
  
  ggsave(sprintf("%s/lambda1/lambda1_theor_by_K_ALS_%s.pdf", args$plotdir, fill_var), 
         p_kals, width = 9, height = 3, units = "in", dpi = 300)
  
  # h2_ALS row (1 row × 3 cols)
  df_h2als <- df %>%
    mutate(x_label = paste0("RR=", lambda1_theor, "\n(K ALS=", K_ALS_input, ")"))
  
  p_h2als <- ggplot(df_h2als, aes(x = factor(lambda1_theor), y = lambda1_emp, fill = factor(.data[[fill_var]]))) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 1.0,
                 outlier.shape = 16, outlier.size = 1) +
    facet_wrap(~ h2_ALS_input, nrow = 1, ncol = 3, scales = "free",
               labeller = nice_facet_label) +
    scale_fill_manual(values = rdbu_cols) +
    scale_x_discrete(labels = function(x) unique(df_h2als$x_label[as.character(df_h2als$lambda1_theor) == x])) +
    labs(title = "(B) Polygenic Recurrence Risk Validation",
         x = "Theoretical Recurrence Risk", y = "Observed RR",
         fill = label_map[fill_var]) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 5)),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.margin = margin(b = 2),
      axis.text.x = element_text(color = "black", size = 7),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      plot.margin = unit(c(5,5,5,5), "pt")
    )
  
  ggsave(sprintf("%s/lambda1/lambda1_theor_by_h2_ALS_%s.pdf", args$plotdir, fill_var), 
         p_h2als, width = 9, height = 3, units = "in", dpi = 300)
  
  cat(sprintf("✓ lambda1_theor_by_[K_ALS|h2_ALS]_%s.pdf\n", fill_var))
}


cat("\nAll recurrence risk λ₁ plots saved as 9x3 PDFs\n")
cat("RdBu palette used throughout (colours 1, 3, and 4 of RdBu-4)\n")
