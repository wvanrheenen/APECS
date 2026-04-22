#!/usr/bin/env Rscript

## Load packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

## Input path ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  in_file <- args[1]
} else {
  in_file <- "combined_simulations.csv"
}

message("Reading: ", in_file)

# Load & filter to monogenic ALS cases
df <- read.csv(in_file, header = TRUE)
df_als <- df %>% filter(mendel_ALS_Y == 1)

# Define column names as CHARACTER VECTOR (this is key!)
target_cols <- c(
  "mendel_ALS_Y", "polygenicY_ALS", 
  "relatives_1st", "relatives_1st_als", "relatives_1st_als_monogenic", "relatives_1st_als_polygenic",
  "relatives_2nd", "relatives_2nd_als", "relatives_2nd_als_monogenic", "relatives_2nd_als_polygenic",
  "relatives_3rd", "relatives_3rd_als", "relatives_3rd_als_monogenic", "relatives_3rd_als_polygenic"
)

## Filter to phenocopies + select columns
df_phenocopy <- df_als %>% 
  filter(phenocopy == 1) %>%
  select(id, all_of(target_cols))  # Now works! target_cols is character vector

# Output stats
message("Found ", nrow(df_phenocopy), " phenocopy cases")

# Write CSV
out_file <- "phenocopy_relatives.csv"
write_csv(df_phenocopy, out_file)
message("Wrote to: ", out_file)

# Preview
print(head(df_phenocopy))
