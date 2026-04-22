library(tidyverse)
library(ggplot2)

# Use snakemake@params for arguments
output_dir <- snakemake@params$output_dir
sweep_type <- snakemake@params$sweep_type
metrics_out <- file.path(output_dir, sweep_type, paste0("metrics_plots_", sweep_type, ".pdf"))

# List all prediction model files for this sweep
files <- list.files(
  path = file.path(output_dir, sweep_type),
  pattern = "prediction_model_perc.csv",
  recursive = TRUE,
  full.names = TRUE
)


# Read and combine all files, extract scenario "4.2) ≥2 1st/2nd ALS/FTD"
metrics <- map_dfr(files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  path_parts <- strsplit(f, "/")[[1]]
  if (length(path_parts) < 4) return(NULL)
  param_dir <- path_parts[3]

  # Extract parameters from param_dir
  frmod <- as.numeric(gsub(".*frmod_([0-9.]+).*", "\\1", param_dir))
  lemod <- as.numeric(gsub(".*lemod_([0-9.]+).*", "\\1", param_dir))
  penalsmod <- as.numeric(gsub(".*penalsmod_([0-9.]+).*", "\\1", param_dir))
  penftdmod <- as.numeric(gsub(".*penftdmod_([0-9.]+).*", "\\1", param_dir))

  # Extract perc from path_parts
  perc_pos <- grep("^[0-9]+perc$", path_parts)
  if (length(perc_pos) == 0) return(NULL)
  perc <- as.numeric(gsub("^([0-9]+)perc$", "\\1", path_parts[perc_pos[1]]))
  if (length(perc) != 1 || is.na(perc)) return(NULL)

  # Filter for scenario and check for NA in parameters
  if (any(is.na(c(frmod, lemod, penalsmod, penftdmod)))) return(NULL)
  df_scen <- df %>% filter(Scenario == "4.2) ≥2 1st/2nd ALS/FTD")
  if (nrow(df_scen) == 0) return(NULL)

  df_scen %>%
    mutate(
      frmod = frmod,
      lemod = lemod,
      penalsmod = penalsmod,
      penftdmod = penftdmod,
      perc = perc
    )
})

# Determine which parameter is varying
if (sweep_type == "varying_fertility") {
  x_var <- "frmod"
  x_lab <- "Fertility Rate"
} else if (sweep_type == "varying_life_expectancy") {
  x_var <- "lemod"
  x_lab <- "Life Expectancy at birth (yrs)"
} else if (sweep_type == "varying_penetrance_als") {
  x_var <- "penalsmod"
  x_lab <- "Cumulative ALS Penetrance"
} else if (sweep_type == "varying_penetrance_ftd") {
  x_var <- "penftdmod"
  x_lab <- "Cumulative FTD Penetrance"
} else if (sweep_type == "varying_percentage") {
  x_var <- "perc"
  x_lab <- "Mendelian Proportion (%)"
}

# Plot each metric
metrics_long <- metrics %>%
  pivot_longer(cols = c(Sensitivity, Specificity, PPV, NPV),
               names_to = "Metric", values_to = "Value")

print(metrics_long)

p <- ggplot(metrics_long, aes_string(x = x_var, y = "Value", color = "Metric")) +
  geom_line() +
  geom_point() +
  labs(
    title = "Predictive metrics for fALS criteria:\n≥2 ALS/FTD relatives within two degrees",
    x = gsub(" Modifier", "", x_lab),  # Remove "Modifier" from x-axis label
    y = "Percentage (%)"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # Show y-axis as percentages
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title[1][2][6]

pdf(metrics_out, width = 8, height = 6)
print(p)
dev.off()
