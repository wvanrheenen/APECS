# Load libraries and functions
source("../../src/libraries_simPed.R")
source("../../src/functions_ALS_dementia.R")

# Get the script's path
script_path <- this.path()
script_dir  <- this.dir()

#' write function to simulate pedigrees without phenotypes (to speed things up)
#' goal is to check if the observed life expectancy matches historical data that will be used as input.
#' 
#' first read in the table with fertility rate over time (for Netherlands in this instance)
# Read in data and set parameters (as in your original script)
fert_rate = read.table(file.path("../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
life_expectancy = read.table(file.path("../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)

monogenic = FALSE
k = 2
lambda = NA
DAF = 0
penetrance = 0
mean_gen_yr = 30
current_year = 2025

# Initialize data frame to store life expectancy estimates
life_expectancy_estimates = data.frame(
  yob = numeric(),
  simulated_life_expectancy = numeric()
)

# Run simulations
for(n in 1:200){
  cat("simulation:", n, "\n")
  yob_index = sample(seq(1850 + k * mean_gen_yr, 2050, 1), 1)
  print(yob_index)
  # Initialize and build pedigree
  core_ped = init_ped(monogenic=monogenic, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF=DAF, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_inlaws(core_ped, DAF=DAF, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)
  
  # Extract life expectancy data from the pedigree
  life_expectancy_data = core_ped %>%
    dplyr::select(yob, life_expectancy) %>%
    group_by(yob) %>%
    summarize(simulated_life_expectancy = mean(life_expectancy))
  # Append to the estimates data frame
  life_expectancy_estimates = bind_rows(life_expectancy_estimates, life_expectancy_data)
}

# Calculate mean simulated life expectancy for each year of birth
mean_life_expectancy = life_expectancy_estimates %>%
  group_by(yob) %>%
  summarize(mean_simulated_life_expectancy = mean(simulated_life_expectancy))

# Prepare historical data
historical_data = life_expectancy %>%
  mutate(
    yob = year_of_birth,
    historical_life_expectancy = life_expectancy # Approximating life expectancy at birth, given that you reached age 15
  ) %>%
  dplyr::select(yob, historical_life_expectancy)

# Create plot
pdf("simulated_life_expectancy_comparison.pdf", h=3, w=3)
ggplot() +
  stat_summary(data = life_expectancy_estimates, 
               aes(x = yob, y = simulated_life_expectancy, color = "Simulated"),
               fun = mean, geom = "line") +
  geom_line(data = historical_data, 
            aes(x = yob, y = historical_life_expectancy, color = "Historical")) +
  theme(legend.position = "bottom") +
  xlab("Year of Birth") +
  ylab("Life Expectancy") +
  scale_color_manual(values = wes_palette(n=2, name="Darjeeling1"),
                     name = "Data") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")
  )
dev.off()

# Save the data
write.csv(life_expectancy_estimates, "life_expectancy_comparison_data.csv", row.names = FALSE)