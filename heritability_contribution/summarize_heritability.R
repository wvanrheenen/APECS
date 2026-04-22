library(dplyr)

# Use snakemake object - tsvs is ALREADY a list of file paths
input_files <- snakemake@input[["tsvs"]]  # List of 15 files
output_file <- snakemake@output[["tsv"]]

print(paste("Reading", length(input_files), "files:"))
print(input_files)

# Read all results (input_files is already a vector of paths)
all_results <- do.call(rbind, lapply(input_files, function(f) {
  read.table(f, header=TRUE, sep="\t")
}))

print(paste("Read", nrow(all_results), "rows total"))
print(table(all_results$set))

# Summary across replicates
summary <- all_results %>%
  group_by(set) %>%
  summarise(
    h2_obs_mean = mean(h2_ALS_obs),
    h2_obs_range = paste0(round(min(h2_ALS_obs),3), "-", round(max(h2_ALS_obs),3)),
    K_obs_mean = mean(K_ALS_obs),
    n_reps = n(),
    total_peds = sum(n_peds),
    .groups="drop"
  ) %>%
  mutate(
    h2_contribution = case_when(
      set == "monogenic" ~ h2_obs_mean,
      set == "polygenic" ~ h2_obs_mean,
      set == "combined" ~ h2_obs_mean,
      TRUE ~ NA_real_
    )
  )

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.table(summary, output_file, sep="\t", quote=FALSE, row.names=FALSE)
print(summary)
