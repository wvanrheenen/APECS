#!/usr/bin/env Rscript
# options(warn = 2)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(wesanderson)  # ADD THIS
})

# ---- argument parsing ----
parser <- ArgumentParser(
  description = "Polygenic simulation validation plots"
)

parser$add_argument(
  "-t", "--type",
  type = "character",
  default = "h2_ALS",
  choices = c("h2_ALS", "K_ALS"),
  help = "Plot type: 'h2_ALS' or 'K_ALS' [default %(default)s]"
)

parser$add_argument(
  "-i", "--input",
  type = "character",
  default = "polygenic_inheritance_all.tsv",
  help = "Input TSV file [default %(default)s]"
)

parser$add_argument(
  "-f", "--facets", 
  type = "character", 
  nargs = '+',
  default = c("h2_FTD_input"),
  help = "Variables to facet by (repeatable) [default %(default)s]"
)

parser$add_argument(
  "-p", "--plotdir", 
  type = "character",
  default = "plots",
  help = "Plot output directory [default %(default)s]"
)

args <- parser$parse_args()

# ---- data ----
df <- read_tsv(args$input, show_col_types = FALSE) %>%
  mutate(
    h2_ALS_obs = as.numeric(h2_ALS_obs),
    K_ALS_obs  = as.numeric(K_ALS_obs)
  )

#!/usr/bin/env Rscript
# options(warn = 2)


suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})


# ---- argument parsing ----
parser <- ArgumentParser(
  description = "Polygenic simulation validation plots"
)


parser$add_argument(
  "-t", "--type",
  type = "character",
  default = "h2_ALS",
  choices = c("h2_ALS", "K_ALS"),
  help = "Plot type: 'h2_ALS' or 'K_ALS' [default %(default)s]"
)


parser$add_argument(
  "-i", "--input",
  type = "character",
  default = "polygenic_inheritance_all.tsv",
  help = "Input TSV file [default %(default)s]"
)


parser$add_argument(
  "-f", "--facets", 
  type = "character", 
  nargs = '+',
  default = c("h2_FTD_input"),
  help = "Variables to facet by (repeatable) [default %(default)s]"
)


parser$add_argument(
  "-p", "--plotdir", 
  type = "character",
  default = "plots",
  help = "Plot output directory [default %(default)s]"
)


args <- parser$parse_args()


# ---- data ----
df <- read_tsv(args$input, show_col_types = FALSE) %>%
  mutate(
    h2_ALS_obs = as.numeric(h2_ALS_obs),
    K_ALS_obs  = as.numeric(K_ALS_obs)
  )

# ---- plotting ----
facet_vars <- args$facets
type_dir <- if (args$type == "h2_ALS") "h2_ALS" else "K_ALS"

for (facet_var in facet_vars) {
  df_plot <- df %>% filter(!is.na(.data[[paste0(args$type, "_obs")]]))
  
  # DYNAMIC Y-AXIS
  y_var <- if (args$type == "h2_ALS") "h2_ALS_obs" else "K_ALS_obs"
  y_max <- max(df_plot[[y_var]], na.rm = TRUE) * 1.1
  
  n_levels <- length(unique(df_plot[[facet_var]]))
  ncol_val <- n_levels
  
  output_file <- file.path(args$plotdir, type_dir, paste0(facet_var, "_validation.png"))
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  if (args$type == "h2_ALS") {
    fill_levels <- sort(unique(df_plot$K_ALS_input))
    darjeeling_cols <- wes_palette("Darjeeling1", n = length(fill_levels), type = "discrete")
    
    p <- ggplot(df_plot, aes(x = factor(h2_ALS_input), y = h2_ALS_obs, fill = factor(K_ALS_input))) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7,
                   outlier.shape = 16, outlier.size = 1) +
      facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol_val,
                 labeller = labeller(.default = label_both)) +
      scale_y_continuous(limits = c(0, y_max), breaks = pretty(c(0, y_max), n = 6)) +
      # AUTO-MATCH Darjeeling1 colors to K_ALS_input levels
      scale_fill_manual(values = setNames(darjeeling_cols, fill_levels)) +
      labs(x = "Input h²_ALS", y = "Observed h²_ALS (Falconer)", 
           fill = "Prevalence\nK_ALS", 
           title = sprintf("Polygenic h²_ALS validation (max=%.2f)", round(y_max, 2))) +
      theme_bw(base_size = 12) + theme(legend.position = "top")
      
  } else {  # K_ALS
    fill_levels <- sort(unique(df_plot$h2_ALS_input))
    darjeeling_cols <- wes_palette("Darjeeling1", n = length(fill_levels), type = "discrete")
    
    p <- ggplot(df_plot, aes(x = factor(K_ALS_input), y = K_ALS_obs, fill = factor(h2_ALS_input))) +
      geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7,
                   outlier.shape = 16, outlier.size = 1) +
      facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol_val,
                 labeller = labeller(.default = label_both)) +
      scale_y_continuous(limits = c(0, y_max), breaks = pretty(c(0, y_max), n = 6)) +
      # AUTO-MATCH Darjeeling1 colors to h2_ALS_input levels
      scale_fill_manual(values = setNames(darjeeling_cols, fill_levels)) +
      labs(x = "Input K_ALS", y = "Observed K_ALS", 
           fill = "Heritability\nh²_ALS", 
           title = sprintf("Polygenic K_ALS validation (max=%.3f)", round(y_max, 3))) +
      theme_bw(base_size = 12) + theme(legend.position = "top")
  }
  
  ggsave(output_file, p, height = 5, width = max(10, ncol_val * 1.5), dpi = 300, bg = "white")
  cat(sprintf("Plot saved: %s (ncol=%d, y_max=%.3f, %d fill levels)\n", 
              output_file, ncol_val, y_max, length(fill_levels)))
}
