## Rscript to run various simulation in parallel

args <- commandArgs(trailingOnly = TRUE)
set <- args[1]
rep <- args[2]
chunk_start <- as.numeric(args[3])
simulations_per_job <- as.numeric(args[4])

source("../src/libraries_simPed.R")
source("../src/functions_simPed.R")

# Parameters for simulation
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table(file.path(script_dir, "../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025
disease_onset = read.csv(file.path(script_dir, "../../data/age_of_onset/age_of_onset.csv"), header = T)

# Sensitivity analysis: shift ALL onset curves right by 5 years
shift_years <- as.integer(set)  # e.g., "5" 
if (!is.na(shift_years) && shift_years != 0) {
  ages <- disease_onset$age
  min_age <- min(ages)
  max_age <- max(ages)
  n <- nrow(disease_onset)
  
  # Choose which curves to shift (monogenic only here)
  col_names <- c("C9_ALS")
  for (col in col_names) {
    if (!(col %in% names(disease_onset))) next
    old_vals <- disease_onset[[col]]
    new_vals <- numeric(n)
    for (i in seq_len(n)) {
      a <- ages[i]
      src_age <- a - shift_years  # where this age should take its incidence from
      if (src_age < min_age) {
        # earlier than table range: set to 0 (no incidence before first age)
        new_vals[i] <- 0
      } else if (src_age > max_age) {
        # beyond last age: keep plateau at last value (or 0 if you prefer)
        new_vals[i] <- old_vals[ages == max_age]
      } else {
        new_vals[i] <- old_vals[ages == src_age]
      }
    }
    disease_onset[[col]] <- new_vals
  }
}
# write.csv(disease_onset, "disease_onset_shifted.csv", row.names = FALSE)

DAF_common = 0.00075 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.00012 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0.00015 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
h2_ALS = 0.40 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.00267 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
penetrance_FTD_common = 0.1  # derived from Gao 2025
penetrance_FTD_patho = 0.1 # copied from Gao 2025
penetrance_FTD_nonALS = 0.9 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
penetrance_dem_common = 0.5 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0.1 # guestimate, no data
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024
plot = FALSE # Do not plot pedigree

yob_index_min = 1940
yob_index_max = 2000

yob_min = yob_index_min - 70
yob_max = yob_index_max + 10
interval = 10

# Run simulations
simulation_output = pedigree_simulations(n_simulations = simulations_per_job, k=k, yob_index_min = yob_index_min, yob_index_max = yob_index_max, lambda=lambda, mean_gen_yr=mean_gen_yr, 
                                  fert_rate=fert_rate, disease_onset=disease_onset, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd,
                                  penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                                  penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                                  K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, h2_dementia=h2_dementia, life_expectancy=life_expectancy, 
                                  current_year=current_year, rg_ALSFTD=rg_ALSFTD, rg_ALSdem=rg_ALSdem, rg_FTDdem=rg_FTDdem, yob_min = yob_min, yob_max = yob_max, interval = interval, plot=plot)

results_df <- simulation_output$results_df
metrics_df <- simulation_output$interval_metrics_combined

# Write results and metrics to separate CSV files
sim_dir <- file.path("results", set, rep, "simulations")
metrics_dir <- file.path("results", set, rep, "metrics")

if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
if (!dir.exists(metrics_dir)) dir.create(metrics_dir, recursive = TRUE)

write.csv(results_df, file.path(sim_dir, sprintf("simulation_%d.csv", chunk_start)), row.names = FALSE)
write.csv(metrics_df, file.path(metrics_dir, sprintf("metrics_%d.csv", chunk_start)), row.names = FALSE)
