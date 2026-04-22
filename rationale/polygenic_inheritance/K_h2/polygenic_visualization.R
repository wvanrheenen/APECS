library(readr)
library(dplyr)
library(ggplot2)

# Define outcome directories and variation parameters you want to analyze
outcome_dirs <- c("h2_als", "K_als")
varying_params <- c("varying_h2_als", "varying_K_als", "varying_h2_ftd", 
                    "varying_K_ftd", "varying_rg", "varying_lambda", "varying_gen")

# Function to read combined .all files and add identifying columns
read_combined_files <- function(outcome_dir, varying_param) {
  file_path <- file.path("results", outcome_dir, paste0(varying_param, ".all"))
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  df <- read_tsv(file_path, show_col_types = FALSE) %>%
    mutate(outcome = outcome_dir,
           varying_param = varying_param)
  return(df)
}

# Read all combined files into one data frame
df_list <- lapply(varying_params, function(param) {
  lapply(outcome_dirs, function(outcome) read_combined_files(outcome, param)) %>%
    bind_rows()
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



# Here is the plotting loop with corrected axis mappings
for (param_name in names(df_list)) {
  df <- df_list[[param_name]] %>%
        mutate(
          h2_als = factor(h2_als),
          K_als = factor(K_als),
          rg = factor(rg),
          lambda = factor(lambda)
        )
  cat(sprintf("Plotting results for %s\n", param_name))
  varying_col <- param_to_column[[param_name]]
  df <- df %>% mutate(varying_value = factor(df[[varying_col]]))
  
    p1 <- ggplot(df, aes_string(x = "K_als", y = "K_empirical", fill = "varying_value")) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
        geom_point(aes(y = as.numeric(as.character(K_als))),
                    color = "black", shape = 18, size = 2,
                    position = position_dodge(width = 0.8)) +
        theme_minimal() +
        ggtitle(paste("Simulated vs Expected Prevalence for", param_name)) +
        labs(fill = varying_col)
    print(p1)

    p2 <- ggplot(df, aes_string(x = "h2_als", y = "h2_empirical", fill = "varying_value")) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
        geom_point(aes(y = as.numeric(as.character(h2_als))),
                    color = "black", shape = 18, size = 2,
                    position = position_dodge(width = 0.8)) +
        theme_minimal() +
        ggtitle(paste("Simulated vs Expected Heritability for", param_name)) +
        labs(fill = varying_col)
    print(p2)    
}