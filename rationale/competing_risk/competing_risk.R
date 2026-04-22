#' ## To do's

#' make sure warnings are treated as errors
options( warn = 1 )
library(wesanderson)
library(scales)

source("../src/libraries_simPed.R")

#' function to generate realistic years between generations
generate_gen_yr = function(mean_gen_yr, min_gen_yr = (mean_gen_yr-5), max_gen_yr = (mean_gen_yr+5)) {
  gen_yr = rnorm(1, mean = mean_gen_yr, sd = (max_gen_yr - min_gen_yr) / 6)
  return(round(max(min_gen_yr, min(max_gen_yr, gen_yr))))
}

#' functions to generate life expectancies from left-skewed normal distribution
get_xi_for_mean <- function(mean_desired, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  xi <- mean_desired - omega * sqrt(2/pi) * delta
  return(xi)
}

get_omega <- function(life_expectancy, max_age = 100, k = 5.5) {
  return(k * (max_age / life_expectancy))
}

sample_life_expectancy <- function(mean_life_exp, max_age = 100, k = 5.5, alpha = -3) {
  omega <- k * (max_age / mean_life_exp)
  delta <- alpha / sqrt(1 + alpha^2)
  xi <- mean_life_exp - omega * sqrt(2/pi) * delta
  round(rsn(1, xi = xi, omega = omega, alpha = alpha))
}

#' ### Function to initiate pedigree with one founder with a mutation
init_ped = function(DAF_common, DAF_patho, DAF_ftd, k, mean_gen_yr, yob_index, life_expectancy, current_year){
  core_ped = data.frame(gen = 0, # generation
                        id = "C0", # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1, # sex, 0 = male, 1 = female
                        a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
                        a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
                        a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
                        a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
                        a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency
                        a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd))) # sample disease allele from population frequency 
  # get the birth year of the founder:
  gen_yr = generate_gen_yr(mean_gen_yr)
  core_ped$yob = yob_index - k * gen_yr
  # simulate life expectancy
  core_ped$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == core_ped$yob]
  core_ped$life_expectancy <- sample_life_expectancy(core_ped$mean_life_exp)
  # simulate age
  core_ped$age = current_year - core_ped$yob
  if(core_ped$age > core_ped$life_expectancy){
    core_ped$status = "dead"
    core_ped$age_censored = core_ped$life_expectancy
  } else {
    core_ped$status = "alive"
    core_ped$age_censored = core_ped$age
  }
  return(core_ped)
}

#' ### Function to add next generation to core pedigree
add_gen = function(df_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, fert_rate, mean_gen_yr, life_expectancy, current_year){
  g = 0
  while(g < k){
    # select individuals from youngest generation
    I1s = filter(df_ped, gen == g) # call I1 for parent 1
    # check if there are individuals in this generation, if not re-simulate generation
    if(nrow(I1s) > 0) {
      # for each individual
      for(i in 1:nrow(I1s)){
        gen_yr = generate_gen_yr(mean_gen_yr)
        # simulate partner (I2 for parent 2)
        I2 = data.frame(gen = g, 
                        id  = gsub("C", "P", I1s$id[i]),
                        pid = NA, # at this point parents are irrelevant, will be simulated in later step
                        mid = NA, # at this point parents are irrelevant, will be simulated in later step
                        sex = abs(I1s$sex[i] - 1), # opposite sex partners only...
                        a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
                        a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
                        a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
                        a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
                        a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency
                        a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency 
                        yob = I1s$yob[i])
        I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
        I2$life_expectancy <- sample_life_expectancy(I2$mean_life_exp) 
        I2$age = current_year - I2$yob
        if(I2$age > I2$life_expectancy){
          I2$status = "dead"
          I2$age_censored = I2$life_expectancy
        } else {
          I2$status = "alive"
          I2$age_censored = I2$age
        }
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rpois(1, fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
        } else { 
          n_II = rpois(1, lambda)
        }
        # make sure generation 0 gets offspring
        while(n_II == 0 & g == 0){
          if(is.na(lambda)){
            n_II = rpois(1, fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
          } else { 
            n_II = rpois(1, lambda)
          }
        }
        # create data_frame for offspring:
        IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
        colnames(IIs) = colnames(df_ped)
        j = 0
        while(j < n_II){
          j = j+1
          IIs$gen[j] = g + 1
          IIs$id[j]  = paste(I1s$id[i], j, sep="_")
          IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
          IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
          IIs$sex[j] = sample(c(0,1), 1)
          IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1)
          IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
          IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1)
          IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)          
          IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1)
          IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1)          
          IIs$yob[j] = I1s$yob[i] + gen_yr
          IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
          IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
          IIs$age[j] = current_year - IIs$yob[j]
          if(IIs$age[j] > IIs$life_expectancy[j]){
            IIs$status[j] = "dead"
            IIs$age_censored[j] = IIs$life_expectancy[j]
          } else {
            IIs$status[j] = "alive"
            IIs$age_censored[j] = IIs$age[j]
          }
        }
        if(n_II > 0){
          df_ped = bind_rows(df_ped, I2, IIs)
        } else {
          df_ped = bind_rows(df_ped, I2)
        }
      }
      g = g+1
    } else {
      g = g-1 # go back one generation and re-simulate to make sure there is any offspring, this will over-estimate pedigree size...
    }
  }
  return(df_ped)
}

#' ### Function to simulate the "inlaws"
#' These are the ancestors for those who married into this pedigree
add_inlaws = function(df_ped, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, life_expectancy, current_year){
  adj_ped = filter(df_ped, grepl("P", id))
  g = max(adj_ped$gen)
  while(! g == 0){
    IIs = filter(adj_ped, gen == g) # get individuals from offspring generation
    for(i in 1:nrow(IIs)){
      gen_yr = generate_gen_yr(mean_gen_yr) 
      # simulate parents 1 (I1)
      I1 = data.frame(gen = g-1, 
                      id  = paste0(IIs$id[i], "_m"), # "_m" suffix for mother
                      pid = NA, 
                      mid = NA, 
                      sex = 1, # mother
                      a1_common  = ifelse(IIs$a1_common[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
                      a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # non-transmitted allele sampled from population
                      a1_patho  = ifelse(IIs$a1_patho[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
                      a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # non-transmitted allele sampled from population 
                      a1_ftd  = ifelse(IIs$a1_ftd[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
                      a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # non-transmitted allele sampled from population
                      yob = IIs$yob[i] - gen_yr) # define year of birth
      I1$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I1$yob]
      I1$life_expectancy = sample_life_expectancy(I1$mean_life_exp)      
      I1$age = current_year - I1$yob
      if(I1$age > I1$life_expectancy){
        I1$status = "dead"
        I1$age_censored = I1$life_expectancy
      } else {
        I1$status = "alive"
        I1$age_censored = I1$age
      }
      # simulate father
      I2 = data.frame(gen = g-1, 
                      id = paste0(IIs$id[i], "_p"), # "_p" suffix for father
                      pid = NA, 
                      mid = NA, 
                      sex = 0,
                      a1_common = ifelse(IIs$a2_common[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
                      a2_common = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # non-transmitted allele sampled from population
                      a1_patho = ifelse(IIs$a2_patho[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
                      a2_patho = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # non-transmitted allele sampled from population                      
                      a1_ftd = ifelse(IIs$a2_ftd[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
                      a2_ftd = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # non-transmitted allele sampled from population 
                      yob = IIs$yob[i] - gen_yr) # define year of birth
      I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
      I2$life_expectancy = sample_life_expectancy(I2$mean_life_exp) 
      I2$age = current_year - I2$yob
      if(I2$age > I2$life_expectancy){
        I2$status = "dead"
        I2$age_censored = I2$life_expectancy
      } else {
        I2$status = "alive"
        I2$age_censored = I2$age
      }                
      # add parents to pedigree
      adj_ped = bind_rows(adj_ped, I1, I2)
      # update parental ID in ped
      adj_ped[adj_ped$id == IIs$id[i],"mid"] = I1$id
      adj_ped[adj_ped$id == IIs$id[i],"pid"] = I2$id
    }
    g = g - 1
  }
  df_ped = filter(df_ped, ! grepl("P", id)) %>% 
    bind_rows(., adj_ped) %>% 
    arrange(., gen)
  return(df_ped)
}

#' ### Function to simulate all external branches of the pedigree 
#' These are the branches with individuals unlinked to the core pedigree.
add_ext_branches = function(df_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, fert_rate, life_expectancy, current_year){
  g = 0
  while(g < k){
    # First, for pairs already in df_ped then select only mothers from the parental generation # sex is irrelevant in this simulation
    I1s = mutate(df_ped, id_parents = gsub("_[a-z]$", "", id)) %>% # create temporary variable ID minus _p or _m (so parents have same ID)
      group_by(id_parents) %>%
      filter(n() > 1) %>% # filter for parental couple
      ungroup() %>%
      filter(grepl("m$", id) & gen == g) # select only mothers in generation of offspring
    if(nrow(I1s > 0)){
      # for each mother
      for(i in 1:nrow(I1s)){
        gen_yr = generate_gen_yr(mean_gen_yr)
        # find father that is already in pedigree
        I2_id = gsub("_m$", "_p", I1s$id[i])
        I2 = filter(df_ped, id == I2_id)
        if(! nrow(I2 == 1)){
          print(I2)
          print(I2_id)
          stop("Non-unique IDs found!!")
        }
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rpois(1, fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"]) - 1 # minus 1 because one child has already been simulated in first round
        } else { 
          n_II = rpois(1, lambda) - 1 # minus 1 because one child has already been simulated in first round
        }
        # simulate offspring simulate to add_gen() function
        if(n_II > 0){
          IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
          colnames(IIs) = colnames(df_ped)
          j = 0
          while(j < n_II){
            j = j+1
            IIs$gen[j] = g + 1
            IIs$id[j]  = paste(I1s$id[i], j, sep="_")
            IIs$pid[j] = I2_id
            IIs$mid[j] = I1s$id[i]
            IIs$sex[j] = sample(c(0,1), 1)
            IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1) # a1 is always from mother - here, I1
            IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
            IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1) # a1 is always from mother - here, I1
            IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)
            IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1) # a1 is always from mother - here, I1
            IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1)            
            IIs$yob[j] = I1s$yob[i] + gen_yr
            IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
            IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
            IIs$age[j] = current_year - IIs$yob[j]
            if(IIs$age[j] > IIs$life_expectancy[j]){
              IIs$status[j] = "dead"
              IIs$age_censored[j] = IIs$life_expectancy[j]
            } else {
              IIs$status[j] = "alive"
              IIs$age_censored[j] = IIs$age[j]
            }
          }
          df_ped = bind_rows(df_ped, IIs)
        }
      }
    }
    
    # Second, find NF who don't have a partner yet (the offspring simulated in previous loop)
    # these have IDs ending with a digit, but do have either "p" or "m" in their IDs
    I1s = filter(df_ped, gen == g & 
                   (grepl("p", id) | grepl("m", id)) &
                   ! (grepl("p$", id) | grepl("m$", id)))
    if(nrow(I1s > 0)){
      for(i in 1:nrow(I1s)){
        gen_yr = generate_gen_yr(mean_gen_yr)
        I2 = data.frame(gen = g, 
                        id = paste0("P",I1s$id[i]),
                        pid = NA, 
                        mid = NA, 
                        sex = abs(I1s$sex[i] - 1),
                        a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)),
                        a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)),
                        a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)),
                        a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)),
                        a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)),
                        a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)),
                        yob = I1s$yob[i])
        I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
        I2$life_expectancy = sample_life_expectancy(I2$mean_life_exp)        
        I2$age = current_year - I2$yob
        if(I2$age > I2$life_expectancy){
          I2$status = "dead"
          I2$age_censored = I2$life_expectancy
        } else {
          I2$status = "alive"
          I2$age_censored = I2$age
        }                
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rpois(1, fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
        } else { 
          n_II = rpois(1, lambda)
        }
        # create data_frame for offspring:
        IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
        colnames(IIs) = colnames(df_ped)
        j = 0
        while(j < n_II){
          j = j+1
          IIs$gen[j] = g + 1
          IIs$id[j]  = paste(I1s$id[i], j, sep="_")
          IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
          IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
          IIs$sex[j] = sample(c(0,1), 1)
          IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1)
          IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
          IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1)
          IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)
          IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1)
          IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1) 
          IIs$yob[j] = I1s$yob[i] + gen_yr
          IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
          IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
          IIs$age[j] = current_year - IIs$yob[j]
          if(IIs$age[j] > IIs$life_expectancy[j]){
            IIs$status[j] = "dead"
            IIs$age_censored[j] = IIs$life_expectancy[j]
          } else {
            IIs$status[j] = "alive"
            IIs$age_censored[j] = IIs$age[j]
          }
        }
        if(n_II > 0){
          df_ped = bind_rows(df_ped, I2, IIs)
        } else {
          df_ped = bind_rows(df_ped, I2)
        }
      }
    }
    g = g+1
  }
  return(df_ped)
}

## Parameters for simulation
k = 2 # total of 4 generations
lambda = NA # mean number of offspring
fert_rate = read.table("GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2200

disease_onset <- read.csv("age_of_onset.csv", header = TRUE)

## Phenotype dependent parameters, can be set to 0
DAF_common = 0 # Disease allele frequency, lower than reported in Van Wijk 2024
DAF_patho = 0
DAF_ftd = 0

Penetrance_ALS = 0.21

## Loop for life expectancy estimates

all_dead_ages <- data.frame()  # to store all dead individuals' ages
total_individuals <- 0         # to count all individuals simulated

for(n in 1:2500){
  cat("simulation:", n, "\n")
  yob_index = sample(seq(1900 + k * mean_gen_yr, 2070, 1), 1)
  
  # Initialize and build pedigree
  core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, 
                      k=k, mean_gen_yr=mean_gen_yr, yob_index=yob_index, 
                      life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, 
                    DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, 
                    mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, 
                        DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, 
                        current_year=current_year)
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, 
                              DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, 
                              fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)
  
  # Append dead ages for this run
  dead_ages_run <- core_ped %>% 
    filter(status == "dead") %>% 
    dplyr::select(age_censored) %>% 
    rename(age = age_censored)
  all_dead_ages <- bind_rows(all_dead_ages, dead_ages_run)
  
  # Accumulate total individuals
  total_individuals <- total_individuals + nrow(core_ped)
}

# Calculate cumulative death % by age
ages_seq <- 0:max(all_dead_ages$age, na.rm=TRUE)

cum_dead_pct <- sapply(ages_seq, function(a){
  sum(all_dead_ages$age <= a) / total_individuals * 100
})

df_plot <- data.frame(
  age = ages_seq,
  cumulative_dead_percent = cum_dead_pct
)
write.csv(df_plot, "cumulative_death_percentages.csv", row.names = FALSE)
df_plot <- read.csv("cumulative_death_percentages.csv")

## Add C9 cumulative incidence
# Extract age and cumulative incidence for C9 ALS
c9_cuminc <- disease_onset %>%
  dplyr::select(age, C9_ALS) %>%
  mutate(cumulative_incidence = cumsum(C9_ALS)) %>%     # cumulative sum of incidence rates
  mutate(scaled_cuminc = cumulative_incidence / max(cumulative_incidence) * Penetrance_ALS * 100)  # scale to penetrance


# Plot
pdf("competing_risk.pdf", width=3, height=5)
ggplot() +
  geom_line(data = df_plot, aes(x = age, y = cumulative_dead_percent, color = "Mortality"), size = 1) +
  geom_line(data = c9_cuminc, aes(x = age, y = scaled_cuminc, color = "ALS Penetrance C9orf72"), size = 1) +
  scale_x_continuous(limits = c(25, 100), breaks = seq(25, 100, 25)) +
  scale_y_continuous(limits = c(0, 100), labels = percent_format(scale = 1)) +
  scale_color_manual(
    name = NULL,
    values = c("Mortality" = wes_palette("Darjeeling1")[1], "ALS Penetrance C9orf72" = wes_palette("Darjeeling1")[2])
  ) +
  guides(color = guide_legend(nrow = 2, ncol = 1, byrow = TRUE)) +  # <--- legend layout here  
  xlab("Age (years)") +
  ylab("Cumulative %") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  )
dev.off()