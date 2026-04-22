#!/usr/bin/env Rscript

options(warn = 2)

# filter_AF_sum.R - Filters gnomAD CSV to Pathogenic variants and sums allele frequencies
# Usage: Rscript filter_AF_sum.R GENE1_DATE.csv GENE2_DATE.csv ...

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript filter_AF_sum.R {GENE_DATE.csv} [more GENE_DATE.csv files]")
}

for (input_file in args) {

  cat("Loading", input_file, "\n")

  # Extract gene name from filename (GENE_DATE.csv -> GENE)
  gene_name <- sub("_\\d+.*\\.csv$", "", basename(input_file))

  data <- read.csv(input_file, stringsAsFactors = FALSE)

  # Filter to Pathogenic variants (column 15: ClinVar Germline Classification)
  clinvar_col <- 15
  pathogenic_terms <- c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic")
  filtered_data <- data[data[[clinvar_col]] %in% pathogenic_terms, ]

  # Allele Frequency column
  af_raw <- filtered_data[["Allele.Frequency"]]
  af_values <- as.numeric(af_raw)
  total_af <- sum(af_values, na.rm = TRUE)

  # AC/AN surrogate SNP method
  ac_col <- which(colnames(filtered_data) == "Allele.Count")
  an_col <- which(colnames(filtered_data) == "Allele.Number")
  if (length(ac_col) > 0 && length(an_col) > 0) {
    ac_values <- as.numeric(filtered_data[[ac_col]])
    an_values <- as.numeric(filtered_data[[an_col]])
    total_ac <- sum(ac_values, na.rm = TRUE)
    total_an <- max(an_values, na.rm = TRUE)
    calculated_af <- total_ac / total_an
  } else {
    calculated_af <- NA
  }

  cat("Number of pathogenic variants:", nrow(filtered_data), "\n")
  cat("Number of pathogenic alleles:", total_ac, "\n")
  cat("Highest number of total alleles:", total_an, "\n")

  # Write per-gene text file
  output_txt <- paste0(gene_name, "_allelefreq.txt")
  cat("mean_DAF =", sprintf("%.8f", total_af), "\n", "AC_AN_DAF =", sprintf("%.8f", calculated_af), "\n", file = output_txt)
  cat("Written:", output_txt, "\n")
}
