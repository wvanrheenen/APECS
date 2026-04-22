# src/filter_mendelian_proportion.R

args <- commandArgs(trailingOnly=TRUE)
input_file <- snakemake@input[["combined_results"]]
output_file <- snakemake@output[["filtered"]]
perc <- as.numeric(gsub(".*_(\\d+)perc_mend.csv", "\\1", output_file)) # Extract percentage from filename

# Read the data
df <- read.csv(input_file)

# Identify Mendelian and polygenic cases
mendel_rows <- which(df$mendel_ALS_Y == 1)
polygenic_rows <- which(df$mendel_ALS_Y == 0)
n_mendel <- length(mendel_rows)
n_polygenic <- length(polygenic_rows)
n_total <- n_mendel + n_polygenic

set.seed(42)

# Start by prioritizing Mendelian cases
x <- min(n_mendel, round((perc / 100) * n_total)) # Mendelian cases
y <- round(x * (100 - perc) / perc) # Polygenic cases

# Adjust if there are not enough polygenic cases to meet the target
if (y > n_polygenic) {
    y <- n_polygenic
    x <- min(n_mendel, round(y * perc / (100 - perc)))
}

# Adjust if there are still not enough mendelian cases
if (x > n_mendel) {
   x <- n_mendel
   y <- min(n_polygenic, round(x * (100 - perc) / perc))
}

selected_mendel <- sample(mendel_rows, x)
selected_polygenic <- sample(polygenic_rows, y)

filtered_df <- df[c(selected_mendel, selected_polygenic), ]
filtered_df <- filtered_df[sample(nrow(filtered_df)), ]

write.csv(filtered_df, output_file, row.names=FALSE)
