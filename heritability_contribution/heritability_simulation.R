#!/usr/bin/env Rscript
# heritability_simulation.R - Monogenic/Polygenic heritability decomposition

options(warn = 2)
source("../src/libraries_simPed.R")
source("adapted_functions_simPed.R")  # Contains pedigree_simulations(), add_pheno(), etc.

## Main parameters (from your main analysis)
k <- 2
lambda <- NA
fert_rate <- read.table("../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr <- 30
life_expectancy <- read.table("../data/life_expectancy_at_fifteen/OWID_life_expectancy_NL_100.txt", header=T)
current_year <- 2200
disease_onset <- read.csv("../data/age_of_onset/age_of_onset.csv", header=T)

yob_index_min <- 1940
yob_index_max <- 2000
yob_min <- yob_index_min - 70
yob_max <- yob_index_max + 10

# Argument parser
args <- commandArgs(trailingOnly = TRUE)
set <- args[1]  # "monogenic", "polygenic", "combined"
n_peds <- as.integer(args[2])
rep_id <- as.integer(args[3])
h2_ALS <- as.numeric(args[4])
output <- args[5]

cat(sprintf("Running %s: %d peds, rep %d → %s\n", set, n_peds, rep_id, output))

# Scenario-specific parameters
if (set == "monogenic") {
  DAF_common <- 0.0008; DAF_patho <- 0.00012; DAF_ftd <- 0.00015
  penetrance_ALS_common <- 0.21; penetrance_ALS_patho <- 0.35
  h2_ALS <- h2_ALS; K_ALS <- 0.0027 * 0  # No polygenic
  penetrance_FTD_common <- 0.1; penetrance_FTD_patho <- 0.1; penetrance_FTD_nonALS <- 0.9
  h2_FTD <- 0.45; K_FTD <- 0.0013 * 0
  penetrance_dem_common <- 0.5; penetrance_dem_patho <- 0.1
  h2_dementia <- 0.55; K_dementia <- 0.418 * 0;
  rg_ALSFTD <- 0.6; rg_ALSdem <- 0.25; rg_FTDdem <- 0.35
} else if (set == "polygenic") {
  DAF_common <- 0; DAF_patho <- 0; DAF_ftd <- 0
  penetrance_ALS_common <- 0; penetrance_ALS_patho <- 0
  h2_ALS <- h2_ALS; K_ALS <- 0.0027 * 0.9
  penetrance_FTD_common <- 0; penetrance_FTD_patho <- 0; penetrance_FTD_nonALS <- 0
  h2_FTD <- 0.45; K_FTD <- 0.0013 * 0.75
  penetrance_dem_common <- 0; penetrance_dem_patho <- 0
  h2_dementia <- 0.55; K_dementia <- 0.418;
  rg_ALSFTD <- 0.6; rg_ALSdem <- 0.25; rg_FTDdem <- 0.35
} else {  # combined (default)
  DAF_common <- 0.0008; DAF_patho <- 0.00012; DAF_ftd <- 0.00015
  penetrance_ALS_common <- 0.21; penetrance_ALS_patho <- 0.35
  h2_ALS <- h2_ALS; K_ALS <- 0.0027 * 0.9
  penetrance_FTD_common <- 0.1; penetrance_FTD_patho <- 0.1; penetrance_FTD_nonALS <- 0.9
  h2_FTD <- 0.45; K_FTD <- 0.0013 * 0.75
  penetrance_dem_common <- 0.5; penetrance_dem_patho <- 0.1
  h2_dementia <- 0.55; K_dementia <- 0.418;
  rg_ALSFTD <- 0.6; rg_ALSdem <- 0.25; rg_FTDdem <- 0.35
}

# Run pedigree simulations for this scenario
# Helper function to simulate heritability for one scenario/replicate
sim_heritability <- function(set, k, n_peds, rep_id) {
  
  # AGGREGATE counters across ALL pedigrees (like sim_polygenic)
  total_population_N <- 0
  total_affected_N <- 0
  total_offspring_N <- 0  
  total_offspring_Y1 <- 0
  
  cat(sprintf("Aggregating %d pedigrees for %s...\n", n_peds, set))
  
  for (sim in 1:n_peds) {
    if(sim %% 50 == 0) cat(sprintf("Sim %d/%d\n", sim, n_peds))
    
    yob_index <- sample(yob_index_min:yob_index_max, 1)
    
    # Build pedigree pipeline (EXACTLY like sim_polygenic)
    core_ped <- init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, 
                        k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, 
                        life_expectancy=life_expectancy, current_year=current_year)
    
    core_ped <- add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, 
                       DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, 
                       mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, 
                       current_year=current_year)
    
    core_ped <- add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, 
                          DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, 
                          life_expectancy=life_expectancy, current_year=current_year) 
    
    core_ped <- add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common,
                          DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr,
                          fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)

    core_ped <- add_pheno(core_ped, disease_onset=disease_onset, 
                         penetrance_ALS_common=penetrance_ALS_common, 
                         penetrance_ALS_patho=penetrance_ALS_patho, h2_ALS=h2_ALS, K_ALS=K_ALS, 
                         penetrance_FTD_common=penetrance_FTD_common, 
                         penetrance_FTD_patho=penetrance_FTD_patho, 
                         penetrance_FTD_nonALS=penetrance_FTD_nonALS, h2_FTD=h2_FTD, K_FTD=K_FTD, 
                         penetrance_dem_common=penetrance_dem_common, 
                         penetrance_dem_patho=penetrance_dem_patho, 
                         h2_dementia=h2_dementia, K_dementia=K_dementia, 
                         rg_ALSFTD=rg_ALSFTD, rg_ALSdem=rg_ALSdem, rg_FTDdem=rg_FTDdem)

    # Aggregate population totals (EXACTLY like sim_polygenic)
    total_population_N <- total_population_N + nrow(core_ped)
    total_affected_N <- total_affected_N + sum(core_ped$Y_ALS == 1, na.rm = TRUE)
    
    # Parent-offspring transmission (EXACTLY like sim_polygenic)
    affected_parents <- core_ped$id[core_ped$Y_ALS == 1]
    for(affected_parent in affected_parents) {
      children <- core_ped[(core_ped$pid == affected_parent | core_ped$mid == affected_parent) & 
                          !is.na(core_ped$pid) & !is.na(core_ped$mid), ]
      total_offspring_N <- total_offspring_N + nrow(children)
      total_offspring_Y1 <- total_offspring_Y1 + sum(children$Y_ALS == 1, na.rm = TRUE)
    }
  }
  
  # FINAL aggregate estimates (EXACTLY like sim_polygenic)
  obs_K <- total_affected_N / max(1, total_population_N)
  
  # Published population prevalence for h2 calculation
  K_ALS_published <- 0.0028 # 1/360
  LT <- -qnorm(obs_K, 0, 1)
  Z <- dnorm(LT)
  i <- Z / obs_K
  
  Kr <- total_offspring_Y1 / max(1, total_offspring_N)
  obs_h2 <- if(Kr > 0 && Kr < 1) {
    LTr <- -qnorm(Kr, 0, 1)
    Zr <- dnorm(LTr)
    aR <- 0.5
    inside <- 1 - (1 - LT/i) * (LT^2 - LTr^2)
    if(is.na(inside) || inside < 0) inside <- 0
    (LT - LTr * sqrt(inside)) / (aR * (i + (i - LT) * LTr^2))
  } else NA
  
  lambda1_emp <- Kr / obs_K  
  
  cat(sprintf("FINAL: h2=%.3f, K=%.4f, λ1=%.2f\n", obs_h2, obs_K, lambda1_emp))
  
  # Output matching sim_polygenic structure (no lambda/rg)
  data.frame(
    set = set,
    k = k,
    n_peds = n_peds,
    total_individuals = total_population_N,
    total_affected = total_affected_N,
    h2_ALS_obs = obs_h2,
    K_ALS_obs = obs_K,
    offspring_N = total_offspring_N,
    offspring_affected = total_offspring_Y1,
    lambda1_emp = lambda1_emp,
    run = rep_id
  )
}

# Run simulation (single row output per rep)
cat(sprintf("Running rep %d: %d peds for %s\n", rep_id, n_peds, set))
results <- sim_heritability(set, k, n_peds, rep_id)

# Write TSV for Snakemake
write.table(results, output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("Results written to: %s\n", output))