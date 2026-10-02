#' ## To do's

#' make sure warnings are treated as errors
options( warn = 2 )

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
penetrance_ALS_common = 0.21 # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 # derived from Douglas 2024 (FUS 19%, SOD1 54%)
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.0027 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
penetrance_FTD_common = 0.1  # derived from Gao 2025
penetrance_FTD_patho = 0.1 # copied from Gao 2025
penetrance_FTD_nonALS = 0.9 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
penetrance_dem_common = 0.5 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0.1 # guestimate, no data
K_dementia = 0.418 # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024
plot = FALSE # Do not plot pedigree

## Loop for life expectancy estimates

# Initialize data frame to store life expectancy estimates
life_expectancy_estimates = data.frame(
  yob = numeric(),
  simulated_life_expectancy = numeric()
)

# # Run simulations
# for(n in 1:2500){
#   cat("simulation:", n, "\n")
#   yob_index = sample(seq(1850 + k * mean_gen_yr, 2070, 1), 1)
#   print(yob_index)
#   # Initialize and build pedigree
#   core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
#   core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)   
  
#   # Extract life expectancy data from the pedigree
#   life_expectancy_data = core_ped %>%
#     dplyr::select(yob, life_expectancy) %>%
#     group_by(yob) %>%
#     summarize(simulated_life_expectancy = mean(life_expectancy))
#   # Append to the estimates data frame
#   life_expectancy_estimates = bind_rows(life_expectancy_estimates, life_expectancy_data)
# }

life_expectancy_estimates <- read.csv("life_expectancy_comparison_data.csv")

# life_expectancy_estimates <- read.csv("life_expectancy_comparison_data.csv")
life_expectancy_estimates <- life_expectancy_estimates %>%
  filter(yob >= 1845)

# # Calculate mean simulated life expectancy for each year of birth
# mean_life_expectancy = life_expectancy_estimates %>%
#   group_by(yob) %>%
#   summarize(mean_simulated_life_expectancy = mean(simulated_life_expectancy))

# Prepare historical data
historical_data = life_expectancy %>%
  mutate(
    yob = year_of_birth,
    historical_life_expectancy = life_expectancy # Approximating life expectancy at birth, given that you reached age 15
  ) %>%
  dplyr::select(yob, historical_life_expectancy)

# Save plot to PNG
pdf("simulated_life_expectancy_comparison.pdf", height = 3, width = 3)
p <- ggplot() +
  stat_summary(data = life_expectancy_estimates, 
               aes(x = yob, y = simulated_life_expectancy, color = "Simulated"),
               fun = mean, geom = "line") +
  geom_line(data = historical_data, 
            aes(x = yob, y = historical_life_expectancy, color = "Historical")) +
  scale_x_continuous(breaks = c(1850, 1900, 1950, 2000, 2050)) +
  labs(title = "(A) Historical vs. simulated life expectancy",
       x = "Year of Birth", y = "Life Expectancy") +
  scale_color_manual(
    values = c(
      "Historical" = RColorBrewer::brewer.pal(4, "RdBu")[1],
      "Simulated" = RColorBrewer::brewer.pal(4, "RdBu")[4]
    ),
    name = "Data"
  ) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )
print(p)
dev.off()


# # # Save the data
# write.csv(life_expectancy_estimates, "life_expectancy_comparison_data.csv", row.names = FALSE)
