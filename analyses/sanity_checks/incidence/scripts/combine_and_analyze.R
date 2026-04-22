# scripts/combine_and_analyze.R

suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(tidyverse))
suppressMessages(library(MASS))
suppressMessages(library(mvnfast))
suppressMessages(library(kinship2))
suppressMessages(library(pedtools))
suppressMessages(library(ribd))
suppressMessages(library(RColorBrewer))
suppressMessages(library(this.path))
suppressMessages(library(wesanderson))
suppressMessages(library(broom))
suppressMessages(library(fmsb))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


# Load all expanded_years chunks
chunk_files <- snakemake@input
expanded_years <- purrr::map_df(chunk_files, readRDS)
start_year <- as.integer(snakemake@params$start_year)
end_year <- as.integer(snakemake@params$end_year)
years_of_interest <- as.character(seq(start_year, end_year))

dir.create("results/debugging", showWarnings = FALSE, recursive = TRUE)


# functions
create_age_year_population <- function(indiv, current_year) {
  # Expand individual's timeline
  df = expand_individual_years(indiv, current_year)
  # Create population count matrix
  df %>%
    group_by(age, year) %>%
    summarise(population = sum(alive), .groups = "drop") %>%
    pivot_wider(
      names_from = year,
      values_from = population,
      values_fill = 0
    ) %>%
    complete(age = full_seq(age, 1), fill = list(population = 0)) %>%
    arrange(age)
}

reorder_year_columns <- function(df) {
  # Extract all columns except 'age'
  year_cols <- setdiff(names(df), "age")
  # Convert year column names to numeric for proper ordering
  ordered_cols <- year_cols[order(as.numeric(year_cols))]
  # Select 'age' first, then ordered year columns
  dplyr::select(df, age, dplyr::all_of(ordered_cols))
}

calculate_lifetime_risk_period <- function(age, mortality_rate, incidence_rate, total_population, total_deaths, total_cases) {
  survival_to_age <- cumprod(1 - dplyr::lag(mortality_rate, default = 0))
  prob_als_at_age <- survival_to_age * incidence_rate
  
  lifetime_risk <- sum(prob_als_at_age, na.rm = TRUE)
  
  tibble(
    age = age,
    total_population = total_population,
    total_deaths = total_deaths,
    total_cases = total_cases,          # Add total ALS cases here
    mortality_rate = mortality_rate,
    incidence_rate = incidence_rate,
    survival_to_age = survival_to_age,
    prob_als_at_age = prob_als_at_age
  ) -> risk_table
  
  cat(sprintf("Estimated lifetime risk of ALS over period: %.2f%%\n", 100 * lifetime_risk))
  return(list(lifetime_risk = lifetime_risk, risk_table = risk_table))
}

calculate_lifetime_risk_bins <- function(age_year_population, age_year_death, age_year_als_onset, start_year, end_year, bin_size = 10) {
  bins <- seq(start_year, end_year, by = bin_size)
  bin_labels <- paste0(bins, "-", bins + bin_size - 1)
  
  binned_results <- purrr::map_df(seq_along(bins), function(i) {
    bin_start <- bins[i]
    bin_end <- min(bins[i] + bin_size - 1, end_year)
    years_in_bin <- as.character(seq(bin_start, bin_end))
    
    # Aggregate over bin years
    total_deaths <- age_year_death %>%
      mutate(total_deaths = rowSums(across(all_of(years_in_bin)))) %>%
      dplyr::select(age, total_deaths)
    
    total_als_onsets <- age_year_als_onset %>%
      mutate(total_als_onsets = rowSums(across(all_of(years_in_bin)))) %>%
      dplyr::select(age, total_als_onsets)
    
    total_population <- age_year_population %>%
      mutate(total_population = rowSums(across(all_of(years_in_bin)))) %>%
      dplyr::select(age, total_population)
    
    age_specific_rates <- total_deaths %>%
      left_join(total_als_onsets, by = "age") %>%
      left_join(total_population, by = "age") %>%
      mutate(
        mortality_rate = ifelse(total_population > 0, total_deaths / total_population, 0),
        als_incidence_rate = ifelse(total_population > 0, total_als_onsets / total_population, 0)
      )
    
    # Calculate lifetime risk for this bin
    risk_result <- calculate_lifetime_risk_period(
      age = age_specific_rates$age,
      mortality_rate = age_specific_rates$mortality_rate,
      incidence_rate = age_specific_rates$als_incidence_rate,
      total_population = age_specific_rates$total_population,
      total_deaths = age_specific_rates$total_deaths,
      total_cases = age_specific_rates$total_als_onsets
    )
    
    tibble(
      bin = bin_labels[i],
      start_year = bin_start,
      end_year = bin_end,
      lifetime_risk = risk_result$lifetime_risk,
      lifetime_risk_percent = 100 * risk_result$lifetime_risk
    )
  })
  return(binned_results)
}

# Build population matrices 
age_year_population_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(population = sum(alive), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = population, values_fill = 0) %>%
  reorder_year_columns()

age_year_death_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(deaths = sum(death), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = deaths, values_fill = 0) %>%
  reorder_year_columns()

age_year_als_onset_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(als_onsets = sum(als_onset), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = als_onsets, values_fill = 0) %>%
  reorder_year_columns()

## aggregate deaths, ALS onsets, and populations over period of interest
total_deaths <- age_year_death_ordered %>%
  mutate(total_deaths = rowSums(across(all_of(years_of_interest)))) %>%
  dplyr::select(age, total_deaths)

total_als_onsets <- age_year_als_onset_ordered %>%
  mutate(total_als_onsets = rowSums(across(all_of(years_of_interest)))) %>%
  dplyr::select(age, total_als_onsets)

total_population <- age_year_population_ordered %>%
  mutate(total_population = rowSums(across(all_of(years_of_interest)))) %>%
  dplyr::select(age, total_population)

age_specific_rates <- total_deaths %>%
  left_join(total_als_onsets, by = "age") %>%
  left_join(total_population, by = "age") %>%
  mutate(
    mortality_rate = ifelse(total_population > 0, total_deaths / total_population, 0),
    als_incidence_rate = ifelse(total_population > 0, total_als_onsets / total_population, 0)
  )  

# Save for debugging
write_csv(age_year_population_ordered, "results/debugging/age_year_population.csv")
write_csv(age_year_death_ordered, "results/debugging/age_year_death.csv")
write_csv(age_year_als_onset_ordered, "results/debugging/age_year_als_onset.csv")

write_csv(total_deaths, "results/debugging/total_deaths.csv")
write_csv(total_als_onsets, "results/debugging/total_als_onsets.csv")
write_csv(total_population, "results/debugging/total_population.csv")

write_csv(age_specific_rates, "results/debugging/age_specific_rates.csv")
  
# Lifetime risk calculation for period
result <- calculate_lifetime_risk_period(
  age = age_specific_rates$age,
  mortality_rate = age_specific_rates$mortality_rate,
  incidence_rate = age_specific_rates$als_incidence_rate,
  total_population = age_specific_rates$total_population,
  total_deaths = age_specific_rates$total_deaths,
  total_cases = age_specific_rates$total_als_onsets   # Pass total ALS cases here
)

binned_risks <- calculate_lifetime_risk_bins(
  age_year_population_ordered,
  age_year_death_ordered,
  age_year_als_onset_ordered,
  start_year = start_year,
  end_year = end_year,
  bin_size = 5
)

cat(sprintf("Start period: %s\n", start_year))
cat(sprintf("End period: %s\n", end_year))
cat(sprintf("Estimated lifetime risk over period: %.4f (%.2f%%)\n", result$lifetime_risk, 100 * result$lifetime_risk))

# Get output file path from Snakemake
saveRDS(result, snakemake@output$lifetime_risk)
rds_outfile <- snakemake@output$lifetime_risk[[1]]

# Create CSV filename by replacing .rds with .csv
csv_outfile <- sub("\\.rds$", ".csv", rds_outfile)

# Save CSV summary
write.csv(result$risk_table, file = csv_outfile, row.names = FALSE)
write.csv(binned_risks, file = snakemake@output$binned_risk, row.names = FALSE)

cat("Saved 5-year binned lifetime risks to", snakemake@output$binned_risk, "\n")

