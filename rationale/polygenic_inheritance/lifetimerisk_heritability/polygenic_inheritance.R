## This is a script to provide a rationale for the way the simulations work under a polygenic disease assumption. 

### Helper functions ### 

#' make sure warnings are treated as errors
options( warn = 2 )

source("../../../src/libraries_simPed.R")
source("adapted_functions_simPed.R")

check_duplicates <- function(df, stage) {
  n_total <- nrow(df)
  n_unique <- length(unique(df$id))
  n_dupes <- sum(duplicated(df$id))
  
  cat(sprintf("%s: %d rows, %d unique IDs, %d duplicates\n", stage, n_total, n_unique, n_dupes))
  
  if (n_dupes > 0) {
    dupes <- names(sort(table(df$id), decreasing = TRUE))[1:10]  # Top 10 offenders
    cat("DUPLICATE IDs:", paste(dupes, table(df$id)[dupes]), "\n")
    print(head(df[df$id %in% dupes[1:3], c("gen", "id", "pid", "mid")], 10))  # Show first 3 dupes
    stop(sprintf("DUPLICATES FOUND at %s! Fix generation code.", stage))
  }
}


## Parameters for simulation
## For this analysis, we're only interested in simulating polygenic disease. 

# script_path <- this.path()
# script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table("../../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("../../../data/life_expectancy_at_fifteen/OWID_life_expectancy_NL_100.txt", header=T)
current_year = 2200
disease_onset = read.csv("../../../data/age_of_onset/age_of_onset.csv", header = T)

DAF_common = 0
DAF_patho = 0
DAF_ftd = 0
penetrance_ALS_common = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019
K_ALS = (0.0027 * 0.85) # derived from Levison 2025, 1/360 = 0.0028. Corrected for 10% monogenic. 
penetrance_FTD_common = 0
penetrance_FTD_patho = 0
penetrance_FTD_nonALS = 0
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
penetrance_dem_common = 0
penetrance_dem_patho = 0
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024
plot = FALSE # Do not plot pedigree

yob_index_min = 1960
yob_index_max = 2000

yob_min = yob_index_min - 70
yob_max = yob_index_max + 10

### Actual analysis ### 

# Argument parser ----
args <- commandArgs(trailingOnly = TRUE)
h2_ALS <- as.numeric(args[1])  
K_ALS <- as.numeric(args[2])                 
h2_FTD <- as.numeric(args[3])                      
K_FTD <- as.numeric(args[4])                 
rg_ALSFTD <- as.numeric(args[5])    
lambda <- as.numeric(args[6])
k <- as.integer(args[7])                      
n_peds <- as.integer(args[8])                 
single_rep <- as.integer(args[9])    # Always 1, ignore
rep_id <- as.integer(args[10])        # Rep number 1-5
output <- args[11]                    # Output path

if (length(args) != 11) stop("Expected 11 arguments")
cat(sprintf("Starting: h2_ALS=%.2f, K_ALS=%.4f, h2_FTD=%.2f, K_FTD=%.4f, rg_ALSFTD=%.2f, lambda=%.2f, k=%d, n_peds=%d, rep=%d, output=%s\n", 
            h2_ALS, K_ALS, h2_FTD, K_FTD, rg_ALSFTD, lambda, k, n_peds, rep_id, output))
            
# Helper function to simulate theoretical recurrence risk ----            
lambda1_theory <- function(K, h2) {
  T  <- -qnorm(K)
  z  <- dnorm(T)
  i  <- z / K
  numer  <- T - h2 * i / 2
  denom  <- sqrt(1 - h2^2 / 4)
  K1     <- 1 - pnorm(numer / denom)
  return(K1 / K)
}

# Helper function to simulate one parameter set ----
sim_polygenic <- function(lambda, h2_ALS, K_ALS, h2_FTD, K_FTD, rg_ALSFTD, k, n_peds) {
  
  # AGGREGATE counters across ALL pedigrees
  total_population_N <- 0
  total_affected_N <- 0
  total_offspring_N <- 0  
  total_offspring_Y1 <- 0
  rg_ALS_FTD_obs <- numeric(n_peds)

  cat(sprintf("Aggregating %d pedigrees...\n", n_peds))
  
  for (sim in 1:n_peds) {
    if(sim %% 50 == 0) cat(sprintf("Sim %d/%d\n", sim, n_peds))
    
    yob_index = sample(yob_index_min:yob_index_max, 1)
    
    # Build pedigree pipeline (unchanged)
    core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, 
                        k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, 
                        life_expectancy=life_expectancy, current_year=current_year)
    
    core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, 
                       DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, 
                       mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, 
                       current_year=current_year)
    
    core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, 
                          DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, 
                          life_expectancy=life_expectancy, current_year=current_year) 
    
    core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common,
                          DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr,
                          fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)

    core_ped = add_pheno(core_ped, disease_onset=disease_onset, 
                         penetrance_ALS_common=penetrance_ALS_common, 
                         penetrance_ALS_patho=penetrance_ALS_patho, h2_ALS=h2_ALS, K_ALS=K_ALS, 
                         penetrance_FTD_common=penetrance_FTD_common, 
                         penetrance_FTD_patho=penetrance_FTD_patho, 
                         penetrance_FTD_nonALS=penetrance_FTD_nonALS, h2_FTD=h2_FTD, K_FTD=K_FTD, 
                         penetrance_dem_common=penetrance_dem_common, 
                         penetrance_dem_patho=penetrance_dem_patho, 
                         h2_dementia=h2_dementia, K_dementia=K_dementia, 
                         rg_ALSFTD=rg_ALSFTD, rg_ALSdem=rg_ALSdem, rg_FTDdem=rg_FTDdem)

    # Aggregate population totals
    total_population_N <- total_population_N + nrow(core_ped)
    total_affected_N <- total_affected_N + sum(core_ped$Y_ALS == 1, na.rm = TRUE)
    
    # Parent-offspring h2
    affected_parents <- core_ped$id[core_ped$Y_ALS == 1]
    for(affected_parent in affected_parents) {
      children <- core_ped[(core_ped$pid == affected_parent | core_ped$mid == affected_parent) & 
                          !is.na(core_ped$pid) & !is.na(core_ped$mid), ]
      total_offspring_N <- total_offspring_N + nrow(children)
      total_offspring_Y1 <- total_offspring_Y1 + sum(children$Y_ALS == 1, na.rm = TRUE)
    
    }
    rg_ALS_FTD_obs[sim] <- cor(core_ped$G_ALS, core_ped$G_FTD, use = "complete.obs")

  }
  
  # FINAL aggregate estimates
  obs_K <- total_affected_N / max(1, total_population_N)
  
  # Existing h2 calculation (unchanged)
  LT <- -qnorm(K_ALS, 0, 1)
  Z <- dnorm(LT)
  i <- Z / K_ALS
  
  Kr <- total_offspring_Y1 / max(1, total_offspring_N)
  obs_h2 <- if(Kr > 0 && Kr < 1) {
    LTr <- -qnorm(Kr, 0, 1)
    Zr <- dnorm(LTr)
    aR <- 0.5
    inside <- 1 - (1 - LT/i) * (LT^2 - LTr^2)
    if(inside < 0) inside <- 0
    (LT - LTr * sqrt(inside)) / (aR * (i + (i - LT) * LTr^2))
  } else NA
  
  lambda1_emp <- Kr / obs_K  
  lambda1_theor <- lambda1_theory(K_ALS, h2_ALS)

  mean_rg_ALS_FTD_obs <- mean(rg_ALS_FTD_obs, na.rm = TRUE)

  cat(sprintf("FINAL: h2=%.3f, K=%.4f, λ1=%.2f(theor:%.2f)\n", 
              obs_h2, obs_K, lambda1_emp, lambda1_theor))
  cat(sprintf("RG ALS-FTD: input=%.3f vs obs=%.3f (diff=%.3f)\n", 
              rg_ALSFTD, mean_rg_ALS_FTD_obs, mean_rg_ALS_FTD_obs - rg_ALSFTD))  # ADD THIS

  
  # EXPANDED output with lambda1
  data.frame(
    lambda = lambda,
    h2_ALS_input = h2_ALS,
    K_ALS_input = K_ALS,
    h2_FTD_input = h2_FTD,
    K_FTD_input = K_FTD,
    rg_ALSFTD_input = rg_ALSFTD,
    k = k,
    n_peds = n_peds,
    total_individuals = total_population_N,
    total_affected = total_affected_N,
    h2_ALS_obs = obs_h2,
    K_ALS_obs = obs_K,
    offspring_N = total_offspring_N,
    offspring_affected = total_offspring_Y1,
    lambda1_emp = lambda1_emp,
    lambda1_theor = lambda1_theor,
    rg_ALS_FTD_obs = mean_rg_ALS_FTD_obs,
    run = rep_id
  )
}

# Main execution - single parameter set
cat(sprintf("Running rep %d: %d peds\n", rep_id, n_peds))
results <- sim_polygenic(lambda = lambda, h2_ALS = h2_ALS, K_ALS = K_ALS, 
                        h2_FTD = h2_FTD, K_FTD = K_FTD, rg_ALSFTD = rg_ALSFTD,
                        k = k, n_peds = n_peds)
results <- results %>% mutate(run = rep_id)

# Write TSV for Snakemake
write.table(results, output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("Results written to: %s\n", output))