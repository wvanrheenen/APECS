# Load libraries and functions
source("../../src/libraries_simPed.R")
source("../../src/functions_simPed.R")

# Get the script's path
script_path <- this.path()
script_dir  <- this.dir()

# Load real-world data on age of parenthood
real_age_data = read_csv(file.path(script_dir, "../../data/age_of_parenthood/swedish_data_1891_2023.csv"))

# Parameters for simulation
fert_rate = read.table(file.path("../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
monogenic  = FALSE
k          = 2
lambda     = NA
DAF        = 0
penetrance = 0
fert_rate  = fert_rate
life_expectancy = read.table(file.path("../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025

# Initialize data frame to store simulated age at parenthood
simulated_age_data = data.frame(
  gen = numeric(),
  yob = numeric(),
  age_at_parenthood = numeric()
)

# Run simulations
for(n in 1:10000){
  cat("simulation:", n,"\n")
  yob_index = sample(seq(1850 + k * mean_gen_yr,2050,1),1)
  print(yob_index)
  print("init")
  core_ped = init_ped(monogenic=monogenic, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  print("add gen")
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF=DAF, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)
  print("add inlaws")
  core_ped = add_inlaws(core_ped, DAF=DAF, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  print("add unlinked")
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)
  
  # Estimate age at parenthood
  max_gen = max(core_ped$gen)
  if(max_gen > 0){
    for(i in 1:max_gen){
      parents = filter(core_ped, gen == i-1)
      for(j in 1:nrow(parents)){
        children = filter(core_ped, pid == parents$id[j] | mid == parents$id[j])
        if(nrow(children) > 0){
          age_at_parenthood = min(children$yob) - parents$yob[j]
          simulated_age_i = data.frame(
            gen = i-1,
            yob = parents$yob[j],
            age_at_parenthood = age_at_parenthood
          )
          simulated_age_data = bind_rows(simulated_age_data, simulated_age_i)
        }
      }
    }
  }
}

print(warnings())

# Create a filtered version of simulated_age_data
#simulated_age_data_filtered <- simulated_age_data %>%
#  filter(yob >= 1890 & yob <= 2025)

# Plot using the filtered data
pdf("simulated_age_at_parenthood.pdf", h=3, w=3)
ggplot(data=simulated_age_data, aes(x=yob, y=age_at_parenthood, color="Simulated")) +
  stat_summary(fun=mean, geom="line") + 
  geom_line(data=real_age_data, aes(x=Year, y=`Mean age at childbearing, historical`, color="Historical")) + 
  scale_x_continuous(limits = c(1890, 2025)) +
  theme(legend.position = "bottom")+
  xlab("Year of Birth") +
  ylab("mean age at parenthood") +
  scale_color_manual(values = wes_palette(n=2, name="Darjeeling1"),  
                   name = "Data") +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black")
  )
dev.off()

# Save simulated data
fwrite(simulated_age_data, "./simulated_age_at_parenthood_estimates.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")
