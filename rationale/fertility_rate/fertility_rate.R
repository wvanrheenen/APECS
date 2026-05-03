#' ## To do's

#' make sure warnings are treated as errors
options( warn = 1 )

source("../../src/libraries_simPed.R")
source("../../src/functions_simPed.R")

## Parameters for simulation
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table("../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2025
disease_onset = read.csv(file.path(script_dir, "../../data/age_of_onset/age_of_onset.csv"), header = T)

DAF_common = 0.00075 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.00012 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0.00015 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.0027 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
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

# ## Loop for fertility rate estimates

# fert_estimates = data.frame(gen=numeric(),
#                             yob_children=numeric(),
#                             N_offspring=numeric())
# for(n in 1:1000){
#   cat("simulation:", n,"\n")
#   yob_index = sample(seq(1860 + k * mean_gen_yr,2065,1),1)
#   print(yob_index)

#   print("Initiate pedigree")
#   # 1. Initialize founder with up-to-date function arguments and allele labels
#   core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   print("Add generations")
#   # 2. Add generations with up-to-date allele labels & fertility logic
#   core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   print("Add inlaws")
#   # 3. Add "inlaw" branches (correct allele handling)
#   core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   print("Add external branches")
#   # 4. Add external branches (correct allele handling)
#   core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)   
  
#   print("Estimate fertility")
#   fert_estimates_n = data.frame(gen=numeric(),
#                                 yob_children=numeric(),
#                                 N_offspring=numeric())
#   max_gen = max(core_ped$gen)
#   if(max_gen > 0){
#     for(i in 1:max_gen){
#       women = filter(core_ped, gen == i-1 & sex == 1)
#       for(j in 1:nrow(women)){
#         II = filter(core_ped, mid == women$id[j])
#         if(nrow(II) > 0) {
#           fert_estimates_i = data.frame(
#             gen = i-1,
#             yob_children = round(mean(II$yob)),
#             N_offspring = nrow(II)
#           )
#         } else { # Later generations have a higher chance of producing no offspring, falsely increasing the mean no of offspring, and thus no yob_children. 
#           fert_estimates_i = data.frame(
#             gen = i-1,
#             yob_children = women$yob[j] + mean_gen_yr, # Expected year of birth for children; should be provided for plot
#             N_offspring = 0
#           )
#         }
#         fert_estimates_n = bind_rows(fert_estimates_n, fert_estimates_i)
#       }
#     }
#   }
# fert_estimates = bind_rows(fert_estimates, fert_estimates_n)
# }

fert_estimates <- fread("./simulated_fertility_estimate_3gen.txt", header = TRUE, sep = "\t")

fert_estimates <- fert_estimates %>%
  filter(yob_children >= 1885 & yob_children <= 2065)

fert_rate <- fert_rate %>%
  filter(year >= 1850 & year <2100)

pdf("simulated_fertility_rate.pdf", height = 3, width = 4.5)
ggplot(data = fert_estimates, aes(x = yob_children, y = N_offspring, color = "Simulated")) +
  stat_summary(fun = mean, geom = "line") + 
  geom_line(data = fert_rate, aes(x = year, y = mean_fertility, color = "Historical")) + 
  xlab("Year of offspring") +
  ylab("Mean fertility rate") +
  labs(title = "(B) Historical vs. simulated fertility rate") +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Data") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face="bold", hjust = 0.5, margin = margin(b = 5)),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )
dev.off()

# fwrite(fert_estimates, "./simulated_fertility_estimate_3gen.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
