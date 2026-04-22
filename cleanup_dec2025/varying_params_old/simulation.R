## Rscript to run various simulation in parallel

args <- commandArgs(trailingOnly = TRUE)
chunk_start <- as.numeric(args[1])
simulations_per_job <- as.numeric(args[2])
fertility_mod <- as.numeric(args[3])
lifeexp_mod <- as.numeric(args[4])
penetrance_als_mod <- as.numeric(args[5])
penetrance_ftd_mod <- as.numeric(args[6])
sweep_type <- args[7]
results_out <- args[8]
metrics_out <- args[9]

# Ensure output directories exist
dir.create(dirname(results_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_out), recursive = TRUE, showWarnings = FALSE)

source("src/libraries_simPed.R")
source("src/functions_simPed.R")

# Parameters for simulation
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "/data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table(file.path(script_dir, "/data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2100
DAF_common = 0.0011 # Disease allele frequency, lower than reported in Van Wijk 2024
DAF_patho = 0.0
DAF_ftd = 0.0001
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = 0.45 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = 0.0033 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
penetrance_FTD_common = 0.1 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 
penetrance_FTD_patho = 0.1 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_FTD_nonALS = 0.9 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_FTD = 0.375 # additive polygenic heritability
K_FTD = 0.0010 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_FTD_risk.txt" ), header=T) # derived from Coyle-Gilchrist 2016
penetrance_dem_common = 0.45
penetrance_dem_patho = 0.2
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.4 # heritabiliy is 0.5 for AD according to Gatz 2006; lowered to compensate for prevalence of e.g. vascular dementia
rg = 0.8
plot = FALSE # Do not plot pedigree


# Apply modifiers
if (sweep_type == "varying_fertility") {
  fert_rate$mean_fertility <- fertility_mod
  # Do NOT modify life_expectancy or penetrances
} 
if (sweep_type == "varying_life_expectancy") {
  life_expectancy$life_expectancy <- lifeexp_mod
} 
if (sweep_type == "varying_year_of_birth") {
  life_expectancy$life_expectancy <- lifeexp_mod
} 
if (sweep_type == "varying_penetrance_als") {
  penetrance_ALS_common$cum_incidence_c9_ALS <- penetrance_ALS_common$cum_incidence_c9_ALS * (penetrance_als_mod / max(penetrance_ALS_common$cum_incidence_c9_ALS))
  penetrance_ALS_patho$cum_incidence_patho_ALS <- penetrance_ALS_patho$cum_incidence_patho_ALS * (penetrance_als_mod / max(penetrance_ALS_patho$cum_incidence_patho_ALS))
}
if (sweep_type == "varying_penetrance_ftd") {
  penetrance_FTD_common$cum_incidence_c9_FTD <- penetrance_FTD_common$cum_incidence_c9_FTD * (penetrance_ftd_mod / max(penetrance_FTD_common$cum_incidence_c9_FTD))
  penetrance_FTD_common$cum_incidence_patho_FTD <- penetrance_FTD_common$cum_incidence_patho_FTD * (penetrance_ftd_mod / max(penetrance_FTD_common$cum_incidence_patho_FTD))
}


# Run simulations
simulation_output = pedigree_simulations(n_simulations = simulations_per_job, k=k, lambda=lambda, mean_gen_yr=mean_gen_yr, 
                                  fert_rate=fert_rate, disease_onset=disease_onset, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd,
                                  penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                                  penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                                  K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, h2_dementia=h2_dementia, life_expectancy=life_expectancy, 
                                  current_year=current_year, rg=rg, yob_min = yob_min, yob_max = yob_max, interval = interval, plot=plot)

results_df <- simulation_output$results_df
metrics_df <- simulation_output$metrics

# Write results and metrics to specified output files
write.csv(results_df, results_out, row.names = FALSE)
write.csv(metrics_df, metrics_out, row.names = FALSE)