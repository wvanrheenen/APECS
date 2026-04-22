
#' make sure warnings are treated as errors
options( warn = 1 )
source("../src/libraries_simPed.R")
library(wesanderson)
library(scales)

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

#' ### Helper function to sample the age of disease onset
sample_onset_age <- function(freq_table, freq_col) {
  probs <- freq_table[[freq_col]]
  probs <- probs / sum(probs) # ensure sum to 1
  ages <- freq_table$age
  sample(ages, size=1, prob=probs, replace=TRUE)
}

#' ### helper functions to simulate survival after disease onset
# Sample survival time after ALS onset (in years)
sample_survival_ALS <- function() {
  # Log-normal with median ≈ 2.5, IQR ≈ 1.7–5.3
  rlnorm(1, meanlog = log(2.5), sdlog = 0.6)
}

# Sample survival time after FTD onset (in years)
sample_survival_FTD <- function() {
  # Log-normal with median ≈ 6.6, 95% CI ≈ 5.6–7.6
  rlnorm(1, meanlog = log(6.6), sdlog = 0.25)
}

#' ### Function to add genetic values for Mendelian and polygenic inheritance
add_pheno = function(df_ped, disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, 
                      penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg) {
  
  df_ped = df_ped %>%
    rowwise() %>%
    mutate(
      mendel_ALS_y1_common = ifelse(a1_common == 1, rbinom(1, 1, penetrance_ALS_common), 0),
      mendel_ALS_y2_common = ifelse(a2_common == 1, rbinom(1, 1, penetrance_ALS_common), 0),
      mendel_ALS_y1_patho  = ifelse(a1_patho == 1, rbinom(1, 1, penetrance_ALS_patho), 0),
      mendel_ALS_y2_patho  = ifelse(a2_patho == 1, rbinom(1, 1, penetrance_ALS_patho), 0), 
      
      age_mendel_ALS_Y_common = ifelse(mendel_ALS_y1_common == 1 | mendel_ALS_y2_common == 1, 
                                       sample_onset_age(disease_onset, "C9_ALS"), NA),
      age_mendel_ALS_Y_patho = ifelse(mendel_ALS_y1_patho == 1 | mendel_ALS_y2_patho == 1, 
                                      sample_onset_age(disease_onset, "C9_ALS"), NA),  # Change column if needed
      
      mendel_ALS_Y_common = ifelse(!is.na(age_mendel_ALS_Y_common) & age_mendel_ALS_Y_common <= age_censored, 1, 0),
      mendel_ALS_Y_patho  = ifelse(!is.na(age_mendel_ALS_Y_patho)  & age_mendel_ALS_Y_patho  <= age_censored, 1, 0),
      mendel_ALS_Y = ifelse(mendel_ALS_Y_common + mendel_ALS_Y_patho > 0, 1, 0)
    ) %>%
    ungroup()

  # FTD section
  df_ped = df_ped %>%
    rowwise() %>%
    mutate(
      mendel_FTD_y1_common = ifelse(a1_common == 1, rbinom(1, 1, penetrance_FTD_common), 0),
      mendel_FTD_y2_common = ifelse(a2_common == 1, rbinom(1, 1, penetrance_FTD_common), 0),
      mendel_FTD_y1_patho  = ifelse(a1_patho == 1, rbinom(1, 1, penetrance_FTD_patho), 0),
      mendel_FTD_y2_patho  = ifelse(a2_patho == 1, rbinom(1, 1, penetrance_FTD_patho), 0),
      mendel_FTD_y1_ftd    = ifelse(a1_ftd   == 1, rbinom(1, 1, penetrance_FTD_nonALS), 0),
      mendel_FTD_y2_ftd    = ifelse(a2_ftd   == 1, rbinom(1, 1, penetrance_FTD_nonALS), 0),
      age_mendel_FTD_Y_common = ifelse(mendel_FTD_y1_common == 1 | mendel_FTD_y2_common == 1, 
                                       sample_onset_age(disease_onset, "C9_FTD"), NA),
      age_mendel_FTD_Y_patho  = ifelse(mendel_FTD_y1_patho == 1 | mendel_FTD_y2_patho == 1, 
                                       sample_onset_age(disease_onset, "C9_FTD"), NA), # Change column if needed
      age_mendel_FTD_Y_ftd    = ifelse(mendel_FTD_y1_ftd == 1 | mendel_FTD_y2_ftd == 1, 
                                       sample_onset_age(disease_onset, "C9_FTD"), NA), # Or GRN/MAPT if available
      mendel_FTD_Y_common = ifelse(!is.na(age_mendel_FTD_Y_common) & age_mendel_FTD_Y_common <= age_censored, 1, 0),
      mendel_FTD_Y_patho  = ifelse(!is.na(age_mendel_FTD_Y_patho)  & age_mendel_FTD_Y_patho  <= age_censored, 1, 0),
      mendel_FTD_Y_ftd    = ifelse(!is.na(age_mendel_FTD_Y_ftd)    & age_mendel_FTD_Y_ftd    <= age_censored, 1, 0),
      mendel_FTD_Y = ifelse(mendel_FTD_Y_common + mendel_FTD_Y_patho + mendel_FTD_Y_ftd > 0, 1, 0)
    ) %>%
    ungroup()

  df_ped = df_ped %>%
    rowwise() %>%
    mutate(
      mendel_dem_y1_common = ifelse(a1_common == 1, rbinom(1, 1, penetrance_dem_common), 0),
      mendel_dem_y2_common = ifelse(a2_common == 1, rbinom(1, 1, penetrance_dem_common), 0),
      mendel_dem_y1_patho  = ifelse(a1_patho == 1, rbinom(1, 1, penetrance_dem_patho), 0),
      mendel_dem_y2_patho  = ifelse(a2_patho == 1, rbinom(1, 1, penetrance_dem_patho), 0), 
      
      age_mendel_dem_Y_common = ifelse(mendel_dem_y1_common == 1 | mendel_dem_y2_common == 1, 
                                       sample_onset_age(disease_onset, "C9_Dementia"), NA),
      age_mendel_dem_Y_patho = ifelse(mendel_dem_y1_patho == 1 | mendel_dem_y2_patho == 1, 
                                      sample_onset_age(disease_onset, "C9_Dementia"), NA),  # Change column if needed
      
      mendel_dem_Y_common = ifelse(!is.na(age_mendel_dem_Y_common) & age_mendel_dem_Y_common <= age_censored, 1, 0),
      mendel_dem_Y_patho  = ifelse(!is.na(age_mendel_dem_Y_patho)  & age_mendel_dem_Y_patho  <= age_censored, 1, 0),
      mendel_dem_Y = ifelse(mendel_dem_Y_common + mendel_dem_Y_patho > 0, 1, 0)
    ) %>%
    ungroup()


  # polygenic model:
  df_ped$G_ALS = NA
  df_ped$G_FTD = NA
  df_ped$G_dementia = NA
  
  # Sample G from multivariate normal for founders
  founders = which(is.na(df_ped$pid)) 
  # ALS and FTD, correlated
  Vg = matrix(c(h2_ALS, rg * sqrt(h2_ALS * h2_FTD), 
              rg * sqrt(h2_ALS * h2_FTD), h2_FTD), nrow=2, byrow=TRUE)
  Gs = mvrnorm(length(founders), mu = c(0, 0), Sigma = Vg)

  df_ped$G_ALS[founders] = Gs[,1]
  df_ped$G_FTD[founders] = Gs[,2]

  empirical_rg = cor(df_ped$G_ALS[founders], df_ped$G_FTD[founders])
  
  # dementia 
  df_ped$G_dementia[founders] = rnorm(length(founders), mean=0, sd=sqrt(h2_dementia))


  # loop through nonfounders and estimate G based on G of parents and h2
  for(i in 1:nrow(df_ped)){
    if(is.na(df_ped$G_ALS[i])){
      pid = df_ped$pid[i]
      mid = df_ped$mid[i]
      Gp_ALS = df_ped[which(df_ped$id == pid),]$G_ALS
      Gm_ALS = df_ped[which(df_ped$id == mid),]$G_ALS
      Gp_FTD = df_ped[which(df_ped$id == pid),]$G_FTD
      Gm_FTD = df_ped[which(df_ped$id == mid),]$G_FTD
      
      # Calculate mean and variance for offspring's G
      mean_ALS = (Gp_ALS + Gm_ALS) / 2
      var_ALS = 0.5 * h2_ALS
      mean_FTD = (Gp_FTD + Gm_FTD) / 2
      var_FTD = 0.5 * h2_FTD
      
      # Simulate correlated G for offspring
      Vg_offspring = matrix(c(var_ALS, rg * sqrt(var_ALS * var_FTD), 
                              rg * sqrt(var_ALS * var_FTD), var_FTD), nrow=2, byrow=TRUE)
      Gs_offspring = mvrnorm(1, mu = c(mean_ALS, mean_FTD), Sigma=Vg_offspring)
      
      df_ped$G_ALS[i] = Gs_offspring[1]
      df_ped$G_FTD[i] = Gs_offspring[2]
    }
    if(is.na(df_ped$G_dementia[i])){
      pid = df_ped$pid[i]
      mid = df_ped$mid[i]
      Gp_dem = df_ped[which(df_ped$id == pid),]$G_dementia
      Gm_dem = df_ped[which(df_ped$id == mid),]$G_dementia
      df_ped$G_dementia[i] = rnorm(1, mean= (Gp_dem+Gm_dem)/2, sd = sqrt(0.5*h2_dementia)) 
    }
  }
  
  # sample non-genetic value E:
  # Without genetic correlation
  df_ped$E_ALS = rnorm(nrow(df_ped), 0, sqrt(1-h2_ALS))
  df_ped$E_FTD = rnorm(nrow(df_ped), 0, sqrt(1-h2_FTD))
  df_ped$E_dementia = rnorm(nrow(df_ped), 0, sqrt(1-h2_dementia))

  # Polygenic ALS phenotype using constant (lifetime) risk
  df_ped$P_ALS = df_ped$G_ALS + df_ped$E_ALS
  df_ped$LT_ALS = -qnorm(K_ALS, 0, 1)
  df_ped$passed_ALS_LT = ifelse(df_ped$P_ALS > df_ped$LT_ALS, 1, 0)
  df_ped$age_polygenicY_ALS = mapply(
    function(passed) if (passed == 1) sample_onset_age(disease_onset, "Polygenic_ALS") else NA,
    df_ped$passed_ALS_LT
  ) 
  df_ped$polygenicY_ALS = ifelse(!is.na(df_ped$age_polygenicY_ALS) & df_ped$age_polygenicY_ALS <= df_ped$age_censored, 1, 0)


  # Polygenic FTD phenotype using constant (lifetime) risk
  df_ped$P_FTD = df_ped$G_FTD + df_ped$E_FTD
  df_ped$LT_FTD = -qnorm(K_FTD, 0, 1)
  df_ped$passed_FTD_LT = ifelse(df_ped$P_FTD > df_ped$LT_FTD, 1, 0)
  df_ped$age_polygenicY_FTD = mapply(
    function(passed) if (passed == 1) sample_onset_age(disease_onset, "C9_FTD") else NA,
    df_ped$passed_FTD_LT
  ) 
  df_ped$polygenicY_FTD = ifelse(!is.na(df_ped$age_polygenicY_FTD) & df_ped$age_polygenicY_FTD <= df_ped$age_censored, 1, 0)

  # define phenotype other dementias
  # Polygenic dementia phenotype using constant (lifetime) risk
  df_ped$P_dementia = df_ped$G_dementia + df_ped$E_dementia
  df_ped$LT_dem = -qnorm(K_dementia, 0, 1)
  df_ped$passed_dem_LT = ifelse(df_ped$P_dementia > df_ped$LT_dem, 1, 0)
  df_ped$age_polygenicY_dementia = mapply(
    function(passed) if (passed == 1) sample_onset_age(disease_onset, "Polygenic_Dementia") else NA,
    df_ped$passed_dem_LT
  ) 
  df_ped$polygenicY_dementia = ifelse(!is.na(df_ped$age_polygenicY_dementia) & df_ped$age_polygenicY_dementia <= df_ped$age_censored, 1, 0)

  # final phenotype
  df_ped$Y_ALS = ifelse(df_ped$polygenicY_ALS + df_ped$mendel_ALS_Y > 0, 1, 0)
  df_ped$Y_FTD = ifelse(df_ped$polygenicY_FTD + df_ped$mendel_FTD_Y > 0, 1, 0)
  df_ped$Y_dementia_other = ifelse(df_ped$polygenicY_dementia > 0, 1, 0)
  df_ped$Y_dementia = ifelse(df_ped$Y_dementia_other + df_ped$mendel_dem_Y + df_ped$Y_FTD > 0, 1, 0) 

  # age and year of onset if Y = 1
  # Age of onset for ALS
  df_ped$age_ALS = ifelse(df_ped$mendel_ALS_Y_common == 1, df_ped$age_mendel_ALS_Y_common,
                        ifelse(df_ped$mendel_ALS_Y_patho == 1, df_ped$age_mendel_ALS_Y_patho,
                               ifelse(df_ped$polygenicY_ALS == 1, df_ped$age_polygenicY_ALS, NA)))
  df_ped$year_onset_ALS = ifelse(df_ped$Y_ALS == 1 & !is.na(df_ped$age_ALS), df_ped$yob + df_ped$age_ALS, NA)

  # Age of onset for FTD
  df_ped$age_FTD = ifelse(df_ped$mendel_FTD_Y_common == 1, df_ped$age_mendel_FTD_Y_common,
                        ifelse(df_ped$mendel_FTD_Y_patho == 1, df_ped$age_mendel_FTD_Y_patho,
                               ifelse(df_ped$mendel_FTD_Y_ftd == 1, df_ped$age_mendel_FTD_Y_ftd,
                                      ifelse(df_ped$polygenicY_FTD == 1, df_ped$age_polygenicY_FTD, NA))))
  df_ped$year_onset_FTD = ifelse(df_ped$Y_FTD == 1 & !is.na(df_ped$age_FTD), df_ped$yob + df_ped$age_FTD, NA)

  # Age of onset for overall dementia (FTD or other)
  df_ped$age_dementia = ifelse(df_ped$Y_FTD == 1, df_ped$age_FTD,
                            ifelse(df_ped$mendel_dem_Y_common == 1, df_ped$age_mendel_dem_Y_common,
                              ifelse(df_ped$mendel_dem_Y_patho == 1, df_ped$age_mendel_dem_Y_patho,
                                ifelse(df_ped$polygenicY_dementia == 1, df_ped$age_polygenicY_dementia, NA))))
  df_ped$year_onset_dementia = ifelse(df_ped$Y_dementia == 1 & !is.na(df_ped$age_dementia), 
                                    df_ped$yob + df_ped$age_dementia, NA)

  df_ped = df_ped %>%
    rowwise() %>%
    mutate(
      post_surv_ALS = ifelse(Y_ALS == 1 & !is.na(age_ALS), sample_survival_ALS(), NA_real_),
      post_surv_FTD = ifelse(Y_FTD == 1 & !is.na(age_FTD), sample_survival_FTD(), NA_real_),
      # Compute possible censoring times
      censor_ALS = ifelse(!is.na(age_ALS) & !is.na(post_surv_ALS), age_ALS + post_surv_ALS, NA_real_),
      censor_FTD = ifelse(!is.na(age_FTD) & !is.na(post_surv_FTD), age_FTD + post_surv_FTD, NA_real_),
      # Take the minimum of all non-NA censoring times
      age_censored_updated = min(c(age_censored, censor_ALS, censor_FTD), na.rm = TRUE)
    ) %>%
    ungroup()

  df_ped$age_censored <- df_ped$age_censored_updated

  # Update phenotypes depending on the update age_censored;
  # Recalculate phenotypes after updating age_censored
  df_ped = df_ped %>%
    rowwise() %>%
    mutate(
      mendel_ALS_Y_common = ifelse(!is.na(age_mendel_ALS_Y_common) & age_mendel_ALS_Y_common <= age_censored, 1, 0),
      mendel_ALS_Y_patho  = ifelse(!is.na(age_mendel_ALS_Y_patho)  & age_mendel_ALS_Y_patho  <= age_censored, 1, 0),
      mendel_ALS_Y = ifelse(mendel_ALS_Y_common + mendel_ALS_Y_patho > 0, 1, 0),
      polygenicY_ALS = ifelse(!is.na(age_polygenicY_ALS) & age_polygenicY_ALS <= age_censored, 1, 0),
      Y_ALS = ifelse(polygenicY_ALS + mendel_ALS_Y > 0, 1, 0),

      mendel_FTD_Y_common = ifelse(!is.na(age_mendel_FTD_Y_common) & age_mendel_FTD_Y_common <= age_censored, 1, 0),
      mendel_FTD_Y_patho  = ifelse(!is.na(age_mendel_FTD_Y_patho)  & age_mendel_FTD_Y_patho  <= age_censored, 1, 0),
      mendel_FTD_Y_ftd    = ifelse(!is.na(age_mendel_FTD_Y_ftd)    & age_mendel_FTD_Y_ftd    <= age_censored, 1, 0),
      mendel_FTD_Y = ifelse(mendel_FTD_Y_common + mendel_FTD_Y_patho + mendel_FTD_Y_ftd > 0, 1, 0),
      polygenicY_FTD = ifelse(!is.na(age_polygenicY_FTD) & age_polygenicY_FTD <= age_censored, 1, 0),
      Y_FTD = ifelse(polygenicY_FTD + mendel_FTD_Y > 0, 1, 0),  
    
      mendel_dem_Y_common = ifelse(!is.na(age_mendel_dem_Y_common) & age_mendel_dem_Y_common <= age_censored, 1, 0),
      mendel_dem_Y_patho  = ifelse(!is.na(age_mendel_dem_Y_patho)  & age_mendel_dem_Y_patho  <= age_censored, 1, 0),
      mendel_dem_Y = ifelse(mendel_dem_Y_common + mendel_dem_Y_patho > 0, 1, 0),
      polygenicY_dementia = ifelse(!is.na(age_polygenicY_dementia) & age_polygenicY_dementia <= age_censored, 1, 0),
      Y_dementia_other = ifelse(polygenicY_dementia > 0, 1, 0),
      Y_dementia = ifelse(Y_dementia_other + mendel_dem_Y + Y_FTD > 0, 1, 0)
    ) %>%
    ungroup()
  
  return(df_ped)
}

## Actual analysis ## 
k = 2 # total of 4 generations
lambda = NA # mean number of offspring
fert_rate = read.table("GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2200
disease_onset = read.csv("age_of_onset.csv", header = T)
DAF_common = 1 # Disease allele frequency, lower than reported in Van Wijk 2024
DAF_patho = 0.0
DAF_ftd = 0.0
penetrance_ALS_common = 1 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = 0 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = 0.0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
penetrance_FTD_common = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 
penetrance_FTD_patho = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_FTD_nonALS = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_FTD = 0.375 # additive polygenic heritability
K_FTD = 0.0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_FTD_risk.txt" ), header=T) # derived from Coyle-Gilchrist 2016
penetrance_dem_common = 0.0
penetrance_dem_patho = 0.0
K_dementia = 0.0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.0 # heritabiliy is 0.5 for AD according to Gatz 2006; lowered to compensate for prevalence of e.g. vascular dementia
rg = 0.5


# Number of simulation runs
num_sims <- 1000

# Data frame to store age of onset from all sims
age_onset_sim <- data.frame(sim = integer(), age_ALS = numeric())

for (sim_i in 1:num_sims) {
  cat("Running simulation:", sim_i, "\n")
  
  # Initialize and simulate pedigree
  yob_index <- sample(seq(1880 + k * mean_gen_yr, 1950, 1), 1)
  core_ped <- init_ped(DAF_common, DAF_patho, DAF_ftd, k, mean_gen_yr, yob_index, life_expectancy, current_year)
  core_ped <- add_gen(core_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, fert_rate, mean_gen_yr, life_expectancy, current_year)
  core_ped <- add_inlaws(core_ped, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, life_expectancy, current_year)
  core_ped <- add_ext_branches(core_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, fert_rate, life_expectancy, current_year)
  
  # Add phenotype (ALS onset, etc.)
  core_ped <- add_pheno(core_ped, disease_onset, 
                       penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, 
                       penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD,
                       penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg)
  
  # Extract age_ALS for affected individuals (Y_ALS == 1 and non-NA age_ALS)
  affected_ages <- core_ped %>% 
    filter(Y_ALS == 1 & !is.na(age_ALS)) %>%
    dplyr::select(age_ALS)
  
  # Append to output dataframe with simulation index
  if(nrow(affected_ages) > 0) {
    age_onset_sim <- bind_rows(age_onset_sim, data.frame(sim = sim_i, age_ALS = affected_ages$age_ALS))
  }
}

# Summarize simulated age distribution
sim_age_dist <- age_onset_sim %>%
  group_by(age_ALS) %>%
  summarise(count = n()) %>%
  mutate(sim_prob = count / sum(count))

# Prepare empirical age onset distribution from input file (column "C9_ALS")
empirical_age_dist <- disease_onset %>%
  dplyr::select(age, C9_ALS) %>%
  mutate(emp_prob = C9_ALS / sum(C9_ALS)) %>%
  dplyr::select(age, emp_prob)

# Merge datasets on age for comparison
compare_df <- full_join(sim_age_dist, empirical_age_dist, by = c("age_ALS" = "age")) %>%
  mutate(sim_prob = ifelse(is.na(sim_prob), 0, sim_prob),
         emp_prob = ifelse(is.na(emp_prob), 0, emp_prob)) %>%
  rename(age = age_ALS)

# Plot comparison of empirical and simulated age of ALS onset distributions
png("simulated_monogenic_disease_onset_comparison.png", width = 8, height = 8, units = "in", res = 300)
plot_onset <- ggplot(compare_df) +
  geom_line(aes(x = age, y = emp_prob, color = "Historical"), linewidth = 1) +
  geom_line(aes(x = age, y = sim_prob, color = "Simulated"), linewidth = 1) +
  labs(title = "Age at Monogenic ALS Onset: Simulated vs Historical",
       x = "Age at Onset",
       y = "Probability",
       color = "Data") +
  scale_color_manual(values = wes_palette(n=2, name="Darjeeling1"), name = "Data") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # formats 0.04 as 4%
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
print(plot_onset)
dev.off()

write.csv(compare_df, "monogenic_onset.csv", row.names = FALSE)
compare_df <- read.csv("monogenic_onset.csv", header = T)


png("simulated_monogenic_disease_onset_1880_1950.png", height = 3, width = 3, units = "in", res = 300)

ggplot(compare_df) +
  geom_line(aes(x = age, y = emp_prob, color = "Historical"), linewidth = 1) +
  geom_line(aes(x = age, y = sim_prob, color = "Simulated"), linewidth = 1) +
  labs(x = "Age at Onset",
       y = "Probability of monogenic ALS onset",
       color = "Data",
       title = NULL) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Data") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(15, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 8),  # smaller y-axis label font size
        axis.ticks = element_line(color = "black"))

dev.off()

## Actual analysis ## 
k = 2 # total of 4 generations
lambda = NA # mean number of offspring
fert_rate = read.table("GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2200
disease_onset = read.csv("age_of_onset.csv", header = T)
DAF_common = 0.0 # Disease allele frequency, lower than reported in Van Wijk 2024
DAF_patho = 0.0
DAF_ftd = 0.0
penetrance_ALS_common = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = 0 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = 1 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
penetrance_FTD_common = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 
penetrance_FTD_patho = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_FTD_nonALS = 0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_FTD = 0.375 # additive polygenic heritability
K_FTD = 0.0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_FTD_risk.txt" ), header=T) # derived from Coyle-Gilchrist 2016
penetrance_dem_common = 0.0
penetrance_dem_patho = 0.0
K_dementia = 0.0 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.0 # heritabiliy is 0.5 for AD according to Gatz 2006; lowered to compensate for prevalence of e.g. vascular dementia
rg = 0.5


# Number of simulation runs
num_sims <- 1000

# Data frame to store age of onset from all sims
age_onset_sim <- data.frame(sim = integer(), age_ALS = numeric())

for (sim_i in 1:num_sims) {
  cat("Running simulation:", sim_i, "\n")
  
  # Initialize and simulate pedigree
  yob_index <- sample(seq(1880 + k * mean_gen_yr, 1950, 1), 1)
  core_ped <- init_ped(DAF_common, DAF_patho, DAF_ftd, k, mean_gen_yr, yob_index, life_expectancy, current_year)
  core_ped <- add_gen(core_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, fert_rate, mean_gen_yr, life_expectancy, current_year)
  core_ped <- add_inlaws(core_ped, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, life_expectancy, current_year)
  core_ped <- add_ext_branches(core_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, fert_rate, life_expectancy, current_year)
  
  # Add phenotype (ALS onset, etc.)
  core_ped <- add_pheno(core_ped, disease_onset, 
                       penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, 
                       penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD,
                       penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg)
  
  # Extract age_ALS for affected individuals (Y_ALS == 1 and non-NA age_ALS)
  affected_ages <- core_ped %>% 
    filter(Y_ALS == 1 & !is.na(age_ALS)) %>%
    dplyr::select(age_ALS)
  
  # Append to output dataframe with simulation index
  if(nrow(affected_ages) > 0) {
    age_onset_sim <- bind_rows(age_onset_sim, data.frame(sim = sim_i, age_ALS = affected_ages$age_ALS))
  }
}

# Summarize simulated age distribution
sim_age_dist <- age_onset_sim %>%
  group_by(age_ALS) %>%
  summarise(count = n()) %>%
  mutate(sim_prob = count / sum(count))

# Prepare empirical age onset distribution from input file (column "Polygenic_ALS")
empirical_age_dist <- disease_onset %>%
  dplyr::select(age, Polygenic_ALS) %>%
  mutate(emp_prob = Polygenic_ALS / sum(Polygenic_ALS)) %>%
  dplyr::select(age, emp_prob)

# Merge datasets on age for comparison
compare_df <- full_join(sim_age_dist, empirical_age_dist, by = c("age_ALS" = "age")) %>%
  mutate(sim_prob = ifelse(is.na(sim_prob), 0, sim_prob),
         emp_prob = ifelse(is.na(emp_prob), 0, emp_prob)) %>%
  rename(age = age_ALS)
  

write.csv(compare_df, "polygenic_onset.csv", row.names = FALSE)

compare_df <- read.csv("polygenic_onset.csv", header = T)

# Plot comparison of Historical and simulated age of ALS onset distributions
# Save plot as PNG

png("simulated_polygenic_disease_onset_1880_1950.png", height = 3, width = 3, units = "in", res = 300)

ggplot(compare_df) +
  geom_line(aes(x = age, y = emp_prob, color = "Historical"), linewidth = 1) +
  geom_line(aes(x = age, y = sim_prob, color = "Simulated"), linewidth = 1) +
  labs(x = "Age at Onset",
       y = "Probability of polygenic ALS onset",
       color = "Data",
       title = NULL) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Data") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(15, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 8),  # smaller y-axis label font size
        axis.ticks = element_line(color = "black"))

dev.off()