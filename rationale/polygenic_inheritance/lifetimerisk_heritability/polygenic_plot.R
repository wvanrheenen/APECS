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
parser <- ArgumentParser(description = "Polygenic simulation validation plots")

parser$add_argument(
  "-t", "--type", type = "character", default = "h2_ALS",
  choices = c("h2_ALS", "K_ALS"),
  help = "Plot type: 'h2_ALS' or 'K_ALS' [default %(default)s]"
)

parser$add_argument(
  "-i", "--input", type = "character", 
  default = "polygenic_inheritance_all.tsv",
  help = "Input TSV file [default %(default)s]"
)

parser$add_argument(
  "-f", "--facets", type = "character", nargs = '+',
  default = c("h2_FTD_input"),
  help = "Variables to facet by (repeatable) [default %(default)s]"
)

parser$add_argument(
  "-p", "--plotdir", type = "character", default = "plots",
  help = "Plot output directory [default %(default)s]"
)

args <- parser$parse_args()


# ---- data ----
df <- read_tsv(args$input, show_col_types = FALSE) %>%
  mutate(h2_ALS_obs = as.numeric(h2_ALS_obs), K_ALS_obs = as.numeric(K_ALS_obs))


# ---- label mapping for variables ----
label_map <- c(
  "K_ALS"                = "Lifetime risk ALS",
  "K_ALS_input"          = "Lifetime risk ALS",
  "h2_ALS_input"         = "h² ALS",
  "h2_FTD_input"         = "h² FTD",
  "lambda"               = "Fertility rate",
  "rg_ALSFTD_input"      = "Genetic correlation ALS~FTD"
)

# default label fallback if not in map
nice_label <- function(var) {
  if (var %in% names(label_map)) label_map[var] else var
}


# ---- plotting ----
facet_vars <- args$facets
type_dir <- if (args$type == "h2_ALS") "h2_ALS" else "K_ALS"


# Custom labeller that uses nice_label
facet_labeller <- function(variable, value) {
  vlabel <- nice_label(variable)
  paste0(vlabel, " = ", value)
}


for (facet_var in facet_vars) {
  df_plot <- df %>% filter(!is.na(.data[[paste0(args$type, "_obs")]]))
  
  y_var <- if (args$type == "h2_ALS") "h2_ALS_obs" else "K_ALS_obs"
  y_max <- max(df_plot[[y_var]], na.rm = TRUE) * 1.1
  n_levels <- length(unique(df_plot[[facet_var]]))
  ncol_val <- n_levels
  
  output_file <- file.path(args$plotdir, type_dir, paste0(facet_var, "_validation.pdf"))
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  if (args$type == "h2_ALS") {
    fill_var <- "K_ALS_input"
    fill_levels <- sort(unique(df_plot[[fill_var]]))
    darjeeling_cols <- wes_palette("Darjeeling1", length(fill_levels))
    
    p <- ggplot(df_plot, aes(x = factor(h2_ALS_input), y = h2_ALS_obs, fill = factor(.data[[fill_var]]))) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7,
                   outlier.shape = 16, outlier.size = 1) +
      facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol_val,
                 labeller = facet_labeller) +
      scale_y_continuous(limits = c(0, y_max), breaks = pretty(c(0, y_max), n = 6)) +
      scale_fill_manual(values = darjeeling_cols) +
      labs(title = "(C) Polygenic h² ALS Validation",
           x = "Input h² ALS",
           y = "Observed h² ALS",
           fill = nice_label(fill_var)) +
      theme_bw() +
      theme(
        legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
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
      
  } else {  # K_ALS
    fill_var <- "h2_ALS_input"
    fill_levels <- sort(unique(df_plot[[fill_var]]))
    darjeeling_cols <- wes_palette("Darjeeling1", length(fill_levels))
    
    p <- ggplot(df_plot, aes(x = factor(K_ALS_input), y = K_ALS_obs, fill = factor(.data[[fill_var]]))) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7,
                   outlier.shape = 16, outlier.size = 1) +
      facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol_val,
                 labeller = facet_labeller) +
      scale_y_continuous(limits = c(0, y_max), breaks = pretty(c(0, y_max), n = 6)) +
      scale_fill_manual(values = darjeeling_cols) +
      labs(title = "(D) Polygenic Lifetime Risk ALS Validation",
           x = "Input K ALS",
           y = "Observed K ALS",
           fill = nice_label(fill_var)) +
      theme_bw() +
      theme(
        legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
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
  }
  
  ggsave(output_file, p, width = 9, height = 3, units = "in", dpi = 300)
  cat(sprintf("✓ %s\n", basename(output_file)))
}


cat("\nAll polygenic validation plots saved as 9x3 PDFs\n")
cat("Darjeeling1 palette used throughout\n")