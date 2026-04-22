## This is a script to provide a rationale for the way the simulations work under a monogenic disease assumption. 
## Use this is as a rationale and evidence that inheritance patterns follow actual mendelian inheritance 

### Helper functions ### 

#' make sure warnings are treated as errors
options( warn = 2 )

source("../../src/libraries_simPed.R")
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
## For this analysis, we're only interested in simulating monogenic disease. For analysis purposes, we give the founder always 1 disease allele, and 1 normal allele.  
## Since the c9-like and fus/sod1-like monogenic disease are modelled the same way, we just simulate c9-like disease
## We also can't have people dying due to FTD/dementia so we set penetrance and polygenic lifetime risk to zero. 
## Also, for the analysis to work, people need to become old enough to 'live up' to the age disease becomes penetrant, thus we set life-expectancy and the current year to full lifespans

# script_path <- this.path()
# script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table("../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("../../data/life_expectancy_at_fifteen/OWID_life_expectancy_NL_100.txt", header=T)
current_year = 2200
disease_onset = read.csv("../../data/age_of_onset/age_of_onset.csv", header = T)


DAF_common = 0 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = 0 # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
penetrance_FTD_common = 0  # derived from Gao 2025
penetrance_FTD_patho = 0 # copied from Gao 2025
penetrance_FTD_nonALS = 0 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = 0 # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
penetrance_dem_common = 0 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0 # guestimate, no data
K_dementia = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
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
penetrance_ALS_common <- as.numeric(args[1])  
lambda <- as.numeric(args[2])                 
k <- as.integer(args[3])                      
n_peds <- as.integer(args[4])                 
single_rep <- as.integer(args[5])    # Always 1, ignore
rep_id <- as.integer(args[6])        # Rep number 1-5
output <- args[7]                    # Output path

if (length(args) != 7) stop("Expected 7 arguments")
cat(sprintf("Starting: pen=%.2f, lambda=%.1f, k=%d, n_peds=%d, rep=%d\n", 
            penetrance_ALS_common, lambda, k, n_peds, rep_id))

# Helper function to simulate one parameter set ----
sim_monogenic <- function(lambda, p, k, n_peds) {
  penetrance_ALS_common = p # Local override
  affected_counts <- numeric(k)
  total_counts    <- numeric(k)
  
  for (sim in 1:n_peds) {
    cat(sprintf("=== Sim %d/%d ===\n", sim, n_peds))
    yob_index = sample(yob_index_min:yob_index_max, 1)
    cat(sprintf("Year of birth: %d\n", yob_index))

    core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, 
                        k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, 
                        life_expectancy=life_expectancy, current_year=current_year)
    # check_duplicates(core_ped, "after init_ped")

    core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, 
                       DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, 
                       mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, 
                       current_year=current_year)
    # check_duplicates(core_ped, "after add_gen")  

    core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, 
                          DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, 
                          life_expectancy=life_expectancy, current_year=current_year) 
    # check_duplicates(core_ped, "after add_inlaws")  

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

    for (gen in 1:k) {
      inds <- core_ped[core_ped$gen == gen & grepl("^C", core_ped$id), ]
      affected_counts[gen] <- affected_counts[gen] + sum(inds$mendel_ALS_Y == 1, na.rm = TRUE)
      total_counts[gen]    <- total_counts[gen] + nrow(inds)
    }
  }

  empirical <- affected_counts / total_counts
  theoretical <- p * (0.5)^(1:k)
  data.frame(
    penetrance = rep(p, k),
    lambda = rep(lambda, k),
    peds = rep(n_peds, k),
    generation = 1:k,
    expected_proportion = theoretical,
    simulated_proportion = empirical
  )
}


# Run all combinations with multiple runs:
# Main execution - single parameter set
cat(sprintf("Running single rep %d: %d peds\n", rep_id, n_peds))
results <- sim_monogenic(lambda = lambda, p = penetrance_ALS_common, k = k, n_peds = n_peds)
results <- results %>% mutate(run = rep_id)

# Write to stdout (Snakemake redirects to TSV)
write.table(results, output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("Results written to output: %s\n", output))