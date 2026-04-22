#' ---
#' title: "simPed: a simulation scheme for pedigrees with Mendelian and complex traits"
#' output: github_document
#' author: Wouter van Rheenen
#' date: "`r format(Sys.time(), '%d %B %Y')`"
#' ---
#'

options( warn = 1 )

#' ## R code and functions.
source("src/libraries_simPed.R")
source("src/functions_simPed.R")

#' ## step-wise simulation with plots 
# Parameters for simulation
# script_path <- this.path()
# script_dir <- this.dir()
k = 3 # total of 3 generations
lambda = 2 # mean number of offspring
fert_rate = read.table("../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=TRUE)
mean_gen_yr = 30
life_expectancy = read.table("../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt", header=TRUE)
current_year = 2025
disease_onset = read.csv("../data/age_of_onset/age_of_onset.csv", header=TRUE)

DAF_common = 0.0008 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.0001 # Derived from Douglas 2024 + Gnomad 
DAF_ftd = 0.00014 # Derived from Gnomad
penetrance_ALS_common = 0.21 # read.table("../data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = 0.4 #read.table("../data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024
h2_ALS = 0.5 # additive polygenic heritability, derived from vRheenen 
K_ALS = (0.0028 * 0.9) # derived from Levison 2025, 1/375 = 0.0027. Corrected for 10% monogenic. 
penetrance_FTD_common = 0.1 # derived from Van Wijk 2024 
penetrance_FTD_patho = 0.1 # derived from Van Wijk 2024
penetrance_FTD_nonALS = 0.9 # derived from Van Wijk 2024
h2_FTD = 0.375 # additive polygenic heritability
K_FTD = (0.0013 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
penetrance_dem_common = 0.45
penetrance_dem_patho = 0.1
K_dementia = 0.418 # read.table("../data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.4 # heritabiliy is 0.5 for AD according to Gatz 2006; lowered to compensate for prevalence of e.g. vascular dementia
rg = 0.8

#' Step 1. simulate the core pedigree (with monogenic disease)
#+ step1_core_pedigree
core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=1960, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
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
color = ifelse(ped_plt2$ID %in% ped_plt1$ID, "white", "lightgray")

plot(ped_plt2, title="step 2 - simulated the in-laws", cex=1.5, carrier = carriers, fill=color)

png("step2_inlaws.png", width=800, height=400)  # Adjust size and resolution as needed
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

color = ifelse(ped_plt3$ID %in% ped_plt1$ID, "white", ifelse(ped_plt3$ID %in% ped_plt2$ID, "lightgray", "darkgray"))

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
                     penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg)

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
core_ped = add_1st_relatives(core_ped)
core_ped = add_2nd_and_3rd_relatives(core_ped)
core_ped = add_4th_to_6th_relatives(core_ped)

#core_ped_ALS = affected_pathways_ALS(core_ped)
#core_ped_FTD = affected_pathways_FTD(core_ped)
#core_ped_dem = affected_pathways_dem(core_ped)

#write.csv(core_ped_ALS, "core_ped_ALS.csv", row.names = FALSE)
#write.csv(core_ped_FTD, "core_ped_FTD.csv", row.names = FALSE)
#write.csv(core_ped_dem, "core_ped_dem.csv", row.names = FALSE)

core_ped = affected_pathways_ALS(core_ped)
core_ped = affected_pathways_FTD(core_ped)
core_ped = affected_pathways_dem(core_ped)

write.csv(core_ped, "core_ped.csv", row.names = FALSE)

#' step 6. Wrapper function
simulated_ped = sim_ped(i=100, k=2, plot=TRUE,
                        DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, disease_onset = disease_onset, mean_gen_yr=mean_gen_yr,
                        penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                        penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                        K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, 
                        h2_dementia=h2_dementia, life_expectancy=life_expectancy, current_year=current_year, rg=rg)

#' step 7. Simulate 25000 pedigrees with monogenic inheritance 

# Parameters for simulation
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table("/data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("/data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2100
disease_onset = read.csv("/data/age_of_onset/age_of_onset.csv", header = T)
DAF_common = 0.001 # Disease allele frequency, lower than reported in Van Wijk 2024
DAF_patho = 0.00015
DAF_ftd = 0.0003
penetrance_ALS_common = 0.20 # read.table("/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = 0.5 #read.table("/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = 0.0033 # read.table("/data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
penetrance_FTD_common = 0.1 # read.table("/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 
penetrance_FTD_patho = 0.1 # read.table("/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_FTD_nonALS = 0.95 # read.table("/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_FTD = 0.375 # additive polygenic heritability
K_FTD = 0.0013 # read.table("/data/lifetime_disease_risk/lifetime_FTD_risk.txt" ), header=T) # derived from Coyle-Gilchrist 2016
penetrance_dem_common = 0.5
penetrance_dem_patho = 0.2
K_dementia = 0.418 # read.table("/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.4 # heritabiliy is 0.5 for AD according to Gatz 2006; lowered to compensate for prevalence of e.g. vascular dementia
rg = 0.8
plot = FALSE # Do not plot pedigree

# Run simulations
results_df = pedigree_simulations(n_simulations = 25, k=k, lambda=lambda, mean_gen_yr=mean_gen_yr, 
                                  fert_rate=fert_rate, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd,
                                  penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                                  penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, 
                                  K_FTD=K_FTD, h2_FTD=h2_FTD, penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho=penetrance_dem_patho, K_dementia=K_dementia, h2_dementia=h2_dementia, life_expectancy=life_expectancy, 
                                  current_year=current_year, rg=rg, plot=plot)

# Save results
write.csv(results_df, "test_10000.csv", row.names = FALSE)
