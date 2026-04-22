#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(wesanderson)
})


# ---- argument parsing ----
parser <- ArgumentParser(description = "rg_ALSFTD validation plots")
parser$add_argument("-i", "--input", type = "character", default = "polygenic_inheritance_all.tsv")
parser$add_argument("-f", "--facets", type = "character", nargs = '+', 
                    default = c("K_ALS_input", "h2_FTD_input", "K_FTD_input", "lambda"))
parser$add_argument("-p", "--plotdir", type = "character", default = "plots")
args <- parser$parse_args()


# ---- data ----
df <- read_tsv(args$input, show_col_types = FALSE) %>%
  mutate(rg_ALS_FTD_obs = as.numeric(rg_ALS_FTD_obs)) %>%
  filter(!is.na(rg_ALS_FTD_obs))

dir.create(file.path(args$plotdir, "rg_ALSFTD"), recursive = TRUE, showWarnings = FALSE)


# ---- label mapping ----
label_map <- c(
  "K_ALS_input"          = "Lifetime risk ALS",
  "h2_ALS_input"         = "h² ALS",
  "h2_FTD_input"         = "h² FTD",
  "K_FTD_input"          = "Lifetime risk FTD",
  "lambda"               = "Fertility rate",
  "rg_ALSFTD_input"      = "Genetic correlation ALS~FTD"
)

# vector‑safe version of nice_label
nice_label <- function(var) {
  var <- as.character(var)
  out <- ifelse(var %in% names(label_map), label_map[var], var)
  structure(out, names = NULL)
}

# Custom labeller that uses nice_label
facet_labeller <- function(variable, value) {
  vlabel <- nice_label(variable)
  paste0(vlabel, " = ", value)
}


# ---- plotting ----
facet_vars <- args$facets

for (facet_var in facet_vars) {
  df_plot <- df
  y_max <- max(df_plot$rg_ALS_FTD_obs, na.rm = TRUE) * 1.1
  n_levels <- length(unique(df_plot[[facet_var]]))
  
  output_file <- file.path(args$plotdir, "rg_ALSFTD", paste0(facet_var, "_validation.pdf"))
  fill_levels <- sort(unique(df_plot[[facet_var]]))
  darjeeling_cols <- wes_palette("Darjeeling1", length(fill_levels))
  
  p <- ggplot(df_plot, aes(x = factor(rg_ALSFTD_input), y = rg_ALS_FTD_obs, 
                          fill = factor(.data[[facet_var]]))) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7,
                 outlier.shape = 16, outlier.size = 1) +
    facet_wrap(as.formula(paste("~ h2_ALS_input")), ncol = 4,
               labeller = facet_labeller) +    # not labeller(…=…); just pass the function
    scale_y_continuous(limits = c(0, y_max), breaks = pretty(c(0, y_max), n = 6)) +
    scale_fill_manual(values = darjeeling_cols) +
    labs(title = "(E) Genetic Correlation ALS~FTD Validation",
         x = "Input genetic correlation",
         y = "Observed correlation",
         fill = nice_label(facet_var)) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.margin = margin(b = 2),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      plot.margin = unit(c(5,5,5,5), "pt")
    )
  
  ggsave(output_file, p, width = 9, height = 3, units = "in", dpi = 300)
  cat(sprintf("✓ %s\n", basename(output_file)))
}


cat("\nAll genetic correlation plots saved as 9x3 PDFs\n")
cat("Darjeeling1 palette used throughout\n")