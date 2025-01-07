#' loading libraries and functions for src folder seems a bit complicated. However, makes sure the script can be run from any location as long as the file/directory structure is maintained as in the github repository.
#' Could be improved...

#' Get the script's path using commandArgs
script_path = normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
script_dir  = dirname(script_path)
#' load libraries and functions
source(file.path(script_dir, "/../src/libraries_simPed.R"))
source(file.path(script_dir, "/../src/functions_simPed.R"))

#' write function to simulate pedigrees without phenotypes (to speed things up)
#' goal is to check if the observed pedigree size (N offspring per pair, e.g. mean fertility) match historical data that will be used as input.
#' 
#' first read in the table with fertility rate over time (for Europe in this instance)
fert_rate = read.table(file.path(script_dir, "/../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)

monogenic  = FALSE
k          = 2
gen_yr     = 25
lambda     = NA
DAF        = 0
penetrance = 0
fert_rate  = fert_rate

fert_estimates = data.frame(gen=numeric(),
                            yob=numeric(),
                            N_offspring=numeric())
for(n in 1:1500){
  cat("simulation:", n,"\n")
  yob_index  = sample(seq(1800 + k * gen_yr,2050,1),1)
  print(yob_index)
  print("init")
  core_ped = init_ped(monogenic=monogenic, k=k, yob_index=yob_index, gen_yr=gen_yr)
  print("add gen")
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF=DAF, gen_yr=gen_yr, fert_rate=fert_rate)
  print("add inlaws")
  core_ped = add_inlaws(core_ped, DAF=DAF, gen_yr=gen_yr)
  print("add unlinked")
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF, gen_yr=gen_yr, fert_rate=fert_rate)

  print("estimate fertility")
  fert_estimates_n = data.frame(gen=numeric(),
                                yob=numeric(),
                                N_offspring=numeric())
  max_gen = max(core_ped$gen)
  if(max_gen > 0){
    for(i in 1:max_gen){
      women = filter(core_ped, gen == i-1 & sex == 1)
      for(j in 1:nrow(women)){
        II = filter(core_ped, mid == women$id[j])
        fert_estimates_i = data.frame(gen = i-1,
                                      yob = women$yob[j],
                                      N_offspring = nrow(II))
        fert_estimates_n = bind_rows(fert_estimates_n, fert_estimates_i)
      }
    }
  }
  fert_estimates = bind_rows(fert_estimates, fert_estimates_n)
}

print(warnings())

ggplot(data=fert_estimates, aes(x=yob, y=N_offspring, color="Simulated")) +
  stat_summary(fun=mean, geom="line") + 
  geom_line(data=fert_rate, aes(x=year, y=mean_fertility, color="Historical")) + 
  theme(legend.position = "bottom")+
  xlab("year") +
  ylab("mean fertility rate") +
  scale_color_manual(values = c("Simulated" = "blue", "Historical" = "red"),
                   name = "Data") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95),     # Position legend inside axes in top-right
    legend.justification = c("right", "top"), # Align legend by its top-right corner
    legend.spacing.y = unit(0.1, "lines"),
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.title = element_text(size = 10), # Adjust legend title size (optional),
    axis.text.x = element_text(color="black"),
    axis.ticks = element_line(color = "black")
  )

fwrite(fert_estimates, "./simulated_fertility_estimates.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")