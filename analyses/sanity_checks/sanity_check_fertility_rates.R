#' loading libraries and functions for src folder seems a bit complicated. However, makes sure the script can be run from any location as long as the file/directory structure is maintained as in the github repository.
#' Could be improved...

# Load libraries and functions
source("../../src/libraries_simPed.R")
source("../../src/functions_ALS_dementia.R")

# Get the script's path
script_path <- this.path()
script_dir  <- this.dir()

#' write function to simulate pedigrees without phenotypes (to speed things up)
#' goal is to check if the observed pedigree size (N offspring per pair, e.g. mean fertility) match historical data that will be used as input.
#' 
#' first read in the table with fertility rate over time (for Netherlands in this instance)
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

fert_estimates = data.frame(gen=numeric(),
                            yob_children=numeric(),
                            N_offspring=numeric())
for(n in 1:1000){
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

  print("estimate fertility")
  fert_estimates_n = data.frame(gen=numeric(),
                              yob_children=numeric(),
                              N_offspring=numeric())
  max_gen = max(core_ped$gen)
  if(max_gen > 0){
    for(i in 1:max_gen){
      women = filter(core_ped, gen == i-1 & sex == 1)
      for(j in 1:nrow(women)){
        II = filter(core_ped, mid == women$id[j])
        if(nrow(II) > 0) {
          fert_estimates_i = data.frame(
            gen = i-1,
            yob_children = round(mean(II$yob)),
            N_offspring = nrow(II)
          )
        } else { # Later generations have a higher chance of producing no offspring, falsely increasing the mean no of offspring, and thus no yob_children. 
          fert_estimates_i = data.frame(
            gen = i-1,
            yob_children = women$yob[j] + mean_gen_yr, # Expected year of birth for children; should be provided for plot
            N_offspring = 0
          )
        }
        fert_estimates_n = bind_rows(fert_estimates_n, fert_estimates_i)
      }
    }
  }
fert_estimates = bind_rows(fert_estimates, fert_estimates_n)
}

print(warnings())

pdf("simulated_fertility_rate_PB.pdf", h=3, w=3)
ggplot(data=fert_estimates, aes(x=yob_children, y=N_offspring, color="Simulated")) +
  stat_summary(fun=mean, geom="line") + 
  geom_line(data=fert_rate, aes(x=year, y=mean_fertility, color="Historical")) + 
  theme(legend.position = "bottom")+
  xlab("year") +
  ylab("mean fertility rate") +
  scale_color_manual(values = wes_palette(n=2, name="Darjeeling1"),  ,
                   name = "Data") +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black")
  )
dev.off()
fwrite(fert_estimates, "./simulated_fertility_estimates_check.txt", col.names=T, row.names=F, quote=F, sep="\t")
