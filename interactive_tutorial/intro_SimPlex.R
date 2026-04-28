#' ---
#' title: "SimPlex: a simulation scheme for pedigrees with Mendelian and complex traits"
#' output: github_document
#' author: Wouter van Rheenen, Paul Beele
#' date: "`r format(Sys.time(), '%d %B %Y')`"
#' ---
#'

options( warn = 2 )

#' ## R code and functions.
source("../src/libraries_simPed.R")
source("../src/functions_simPed.R")

#' ### Function to initiate pedigree with one founder with a mutation
print(init_ped)

#' ### Function to simulate a next generation
print(add_gen)

#' ### Function to simulate the "inlaws"
#' These are the ancestors for those who married into this pedigree
print(add_inlaws)

#' ### Function to simulate all external branches of the pedigree 
#' These are the branches with individuals unlinked to the core pedigree.
print(add_ext_branches)

#' ### Function to add genetic values for Mendelian and polygenic inheritance
print(add_pheno)

#' functions to count number of affected relatives
print(add_1st_relatives)
print(add_2nd_and_3rd_relatives)

#' wrapper function to simulate full pedigree
print(sim_ped)

#' function to simulate a desired number of pedigrees with/without monogenic inheritance
print(pedigree_simulations)

#' ## step-wise simulation with plots 
# Parameters for simulation
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table(file.path(script_dir, "../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025
disease_onset = read.csv(file.path(script_dir, "../data/age_of_onset/age_of_onset.csv"), header = T)

# Mendelian parameters
DAF_common = 0.00075 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.00012 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0.00015 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
penetrance_FTD_common = 0.1  # derived from Gao 2025
penetrance_FTD_patho = 0.1 # copied from Gao 2025
penetrance_FTD_nonALS = 0.9 
penetrance_dem_common = 0.5 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0.1 # guestimate, no data

# polygenic parameters:
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.0027 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024

#' Step 1. simulate the core pedigree (with monogenic disease)
#+ step1_core_pedigree
core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=2010, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
# run the loop to add generations of offspring
core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
ped_plt1 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
color = "white"

png("step1_core_pedigree.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt1,
  title="Step 1 - Core Pedigree",
  cex=1.1,
  carrier = carriers,
  fill=color
)
dev.off()

#' step 2. Simulate the inlaws, ancestors of the spouses married into this pedigree
#+ step2_inlaws
core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
ped_plt2 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
color = ifelse(ped_plt2$ID %in% ped_plt1$ID, "white", "orange")

plot(ped_plt2, title="Step 2 - simulated the in-laws", cex=1.5, carrier = carriers, fill=color)

png("step2_inlaws.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt2,
  title="Step 2 - In-laws",
  cex=1.1,
  carrier = carriers,
  fill=color
)
dev.off()

#' step 3. Simulate external branches, unlinked to founder
#+ step3_unlinked
core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, 
                            life_expectancy=life_expectancy, current_year=current_year)
ped_plt3 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id

color = ifelse(ped_plt3$ID %in% ped_plt1$ID, "white", ifelse(ped_plt3$ID %in% ped_plt2$ID, "orange", "red"))

png("step3_external_branches.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt3,
  title="Step 3 - External branches",
  cex=1.0,
  carrier = carriers,
  fill=color
)
dev.off()

#' step 4. Simulate phenotypes
#+ step4_pheno
core_ped = add_pheno(core_ped, disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, 
                     penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem)

carriers_ALS_FTD = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
affected_ALS = filter(core_ped, Y_ALS == 1)$id
affected_FTD = filter(core_ped, Y_FTD == 1)$id
affected_ALS_FTD = filter(core_ped, Y_ALS == 1 | Y_FTD == 1)$id
affected_dementia = filter(core_ped, Y_dementia_other == 1)$id
affected_either = filter(core_ped, Y_ALS == 1 | Y_dementia == 1)$id
status = filter(core_ped, status == "dead")$id
# make plots
png("plot_ALS.png", width=800, height=600)  
plot(ped_plt3, title="Binary ALS phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_ALS, cex=1.0)
dev.off()
png("plot_FTD.png", width=800, height=600)  
plot(ped_plt3, title="Binary ALS+FTD phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_ALS_FTD, cex=1.0)
dev.off()
png("plot_dem.png", width=800, height=600) 
plot(ped_plt3, title="Binary ALS+dementia phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_either, cex=1.0)
dev.off()

#' step 5. Count relatives
core_ped = add_1st_relatives(core_ped, k)
core_ped = add_2nd_and_3rd_relatives(core_ped, k)

write.csv(core_ped, "core_ped.csv", row.names = FALSE)

#' step 6. Wrapper function
simulated_ped = sim_ped(i=100, k=2, plot=TRUE,
                        DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, disease_onset = disease_onset, mean_gen_yr=mean_gen_yr,
                        penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                        penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                        K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, 
                        h2_dementia=h2_dementia, life_expectancy=life_expectancy, current_year=current_year, rg_ALSFTD=rg_ALSFTD, rg_ALSdem=rg_ALSdem, rg_FTDdem=rg_FTDdem
                        )

#' step 7. Simulate 100 pedigrees 
#' Variables can be altered 

# Demographic parameters
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table(file.path(script_dir, "../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025
disease_onset = read.csv(file.path(script_dir, "../data/age_of_onset/age_of_onset.csv"), header = T)
plot = FALSE

# Mendelian parameters
DAF_common = 0.5 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.00012 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0.00015 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
penetrance_FTD_common = 0.1  # derived from Gao 2025
penetrance_FTD_patho = 0.1 # copied from Gao 2025
penetrance_FTD_nonALS = 0.9 
penetrance_dem_common = 0.5 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0.1 # guestimate, no data

# polygenic parameters:
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.0027 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024

# years of simulations
yob_index_min = 1940
yob_index_max = 2000

yob_min = yob_index_min - 70
yob_max = yob_index_max + 10
interval = 10

# Run simulations
simulation_output = pedigree_simulations(n_simulations = 100, k=k, yob_index_min = yob_index_min, yob_index_max = yob_index_max, lambda=lambda, mean_gen_yr=mean_gen_yr, 
                                  fert_rate=fert_rate, disease_onset=disease_onset, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd,
                                  penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                                  penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                                  K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, h2_dementia=h2_dementia, life_expectancy=life_expectancy, 
                                  current_year=current_year, rg_ALSFTD=rg_ALSFTD, rg_ALSdem=rg_ALSdem, rg_FTDdem=rg_FTDdem, yob_min = yob_min, yob_max = yob_max, interval = interval, plot=plot)

results_df <- simulation_output$results_df

# Save results
write.csv(results_df, "simulations.csv", row.names = FALSE)
