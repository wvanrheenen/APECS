#' ## To do's
#TODO: solve over-estimation of pedigree size,
#TODO: adjust lifetime ALS/FTD risk for age

#' make sure warnings are treated as errors
options( warn = 2 )

#' function to generate realistic years between generations
generate_gen_yr = function(mean_gen_yr, min_gen_yr = (mean_gen_yr-5), max_gen_yr = (mean_gen_yr+5)) {
  gen_yr = rnorm(1, mean = mean_gen_yr, sd = (max_gen_yr - min_gen_yr) / 6)
  return(round(max(min_gen_yr, min(max_gen_yr, gen_yr))))
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
  core_ped$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == (core_ped$yob)])/150))
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
        I2$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob])/150))
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
          IIs$life_expectancy[j] = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]])/150))
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
      I1$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == I1$yob])/150))
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
      I2$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob])/150))
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
            IIs$life_expectancy[j] = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]])/150))
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
        I2$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob])/150))
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
          IIs$life_expectancy[j] = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]])/150))
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

#' ### Function to add genetic values for Mendelian and polygenic inheritance
add_pheno = function(df_ped, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, 
                      penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, h2_dementia, K_dementia, rg) {
  
  # Mendelian ALS (autosomal dominant)
  df_ped$age_penetrance_ALS_common = sapply(df_ped$age_censored, function(age) {
    penetrance_ALS_common$cum_incidence_c9_ALS[which.min(abs(penetrance_ALS_common$age - age))]
  })
  df_ped$age_penetrance_ALS_patho = sapply(df_ped$age_censored, function(age) {
    penetrance_ALS_patho$cum_incidence_patho_ALS[which.min(abs(penetrance_ALS_patho$age - age))]
  })

  df_ped = df_ped %>%
  mutate(
    mendel_ALS_y1_common = ifelse(a1_common == 1, rbinom(n(), 1, age_penetrance_ALS_common), 0),
    mendel_ALS_y2_common = ifelse(a2_common == 1, rbinom(n(), 1, age_penetrance_ALS_common), 0),
    mendel_ALS_y1_patho = ifelse(a1_patho == 1, rbinom(n(), 1, age_penetrance_ALS_patho), 0),
    mendel_ALS_y2_patho = ifelse(a2_patho == 1, rbinom(n(), 1, age_penetrance_ALS_patho), 0), 
    mendel_ALS_Y_common = ifelse(mendel_ALS_y1_common + mendel_ALS_y2_common > 0, 1, 0),
    mendel_ALS_Y_patho = ifelse(mendel_ALS_y1_patho + mendel_ALS_y2_patho > 0, 1, 0),   
    mendel_ALS_Y = ifelse(mendel_ALS_y1_common + mendel_ALS_y2_common + mendel_ALS_y1_patho + mendel_ALS_y2_patho > 0, 1, 0)
  )

  # Mendelian FTD (autosomal dominant)
  df_ped$age_penetrance_FTD_common = sapply(df_ped$age_censored, function(age) {
    penetrance_FTD_common$cum_incidence_c9_FTD[which.min(abs(penetrance_FTD_common$age - age))]
  })
  df_ped$age_penetrance_FTD_patho = sapply(df_ped$age_censored, function(age) {
    penetrance_FTD_patho$cum_incidence_patho_FTD[which.min(abs(penetrance_FTD_patho$age - age))]
  })
  df_ped$age_penetrance_FTD_nonALS = sapply(df_ped$age_censored, function(age) {
    penetrance_FTD_nonALS$cum_incidence_FTD_GRNMAPT[which.min(abs(penetrance_FTD_nonALS$age - age))]
  })

  df_ped = df_ped %>%
  mutate(
    mendel_FTD_y1_common = ifelse(a1_common == 1, rbinom(n(), 1, age_penetrance_FTD_common), 0),
    mendel_FTD_y2_common = ifelse(a2_common == 1, rbinom(n(), 1, age_penetrance_FTD_common), 0),
    mendel_FTD_y1_patho = ifelse(a1_patho == 1, rbinom(n(), 1, age_penetrance_FTD_patho), 0),
    mendel_FTD_y2_patho = ifelse(a2_patho == 1, rbinom(n(), 1, age_penetrance_FTD_patho), 0),
    mendel_FTD_y1_ftd = ifelse(a1_ftd == 1, rbinom(n(), 1, age_penetrance_FTD_nonALS), 0),
    mendel_FTD_y2_ftd = ifelse(a2_ftd == 1, rbinom(n(), 1, age_penetrance_FTD_nonALS), 0),
    mendel_FTD_Y_common = ifelse(mendel_FTD_y1_common + mendel_FTD_y2_common > 0, 1, 0),  
    mendel_FTD_Y_patho = ifelse(mendel_FTD_y1_patho + mendel_FTD_y2_patho > 0, 1, 0),   
    mendel_FTD_Y_ftd = ifelse(mendel_FTD_y1_ftd + mendel_FTD_y2_ftd > 0, 1, 0),   
    mendel_FTD_Y = ifelse(mendel_FTD_y1_common + mendel_FTD_y2_common + mendel_FTD_y1_patho + mendel_FTD_y2_patho + mendel_FTD_y1_ftd + mendel_FTD_y2_ftd > 0, 1, 0)
  )

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

  # Define phenotype ALS with age-specific liability threshold
  df_ped$P_ALS = df_ped$G_ALS + df_ped$E_ALS
  df_ped$age_related_K_ALS = sapply(df_ped$age_censored, function(age) {
    K_ALS$cum_risk_adjusted[which.min(abs(K_ALS$age - age))]
  })
  df_ped$LT_ALS = -qnorm(df_ped$age_related_K_ALS, 0, 1)
  df_ped$polygenicY_ALS = ifelse(df_ped$P_ALS > df_ped$LT_ALS, 1, 0)

  # Define phenotype FTD with age-specific liability threshold
  df_ped$P_FTD = df_ped$G_FTD + df_ped$E_FTD
  df_ped$age_related_K_FTD = sapply(df_ped$age_censored, function(age) {
    K_FTD$cum_risk_adjusted[which.min(abs(K_FTD$age - age))]
  })
  df_ped$LT_FTD = -qnorm(df_ped$age_related_K_FTD, 0, 1)
  df_ped$polygenicY_FTD = ifelse(df_ped$P_FTD > df_ped$LT_FTD, 1, 0)

  # define phenotype other dementias 
  df_ped$P_dementia = df_ped$G_dementia + df_ped$E_dementia
  df_ped$age_related_K_dementia = sapply(df_ped$age_censored, function(age) {
  K_dementia$cum_incidence[which.min(abs(K_dementia$age - age))]
  })
  df_ped$LT_dem = -qnorm(df_ped$age_related_K_dementia, 0, 1) # liability threshold
  df_ped$polygenicY_dementia = ifelse(df_ped$P_dementia > df_ped$LT_dem, 1, 0)
  
  # final phenotype
  df_ped$Y_ALS = ifelse(df_ped$polygenicY_ALS + df_ped$mendel_ALS_Y > 0, 1, 0)
  df_ped$Y_FTD = ifelse(df_ped$polygenicY_FTD + df_ped$mendel_FTD_Y > 0, 1, 0)
  df_ped$Y_dementia_other = ifelse(df_ped$polygenicY_dementia > 0, 1, 0)
  df_ped$Y_dementia = ifelse(df_ped$Y_dementia_other + df_ped$Y_FTD > 0, 1, 0) 
  
  return(df_ped)
}

#' function to identify relatedness 

add_1st_relatives = function(df_ped) {
  # Initialize columns for relative counts
  df_ped$relatives_1st = 0
  df_ped$relatives_1st_als = 0
  df_ped$relatives_1st_ftd = 0
  df_ped$relatives_1st_ftd_unique = 0
  df_ped$relatives_1st_dementia = 0
  df_ped$relatives_1st_dementia_unique = 0
  
  # Loop through each individual
  for(i in 1:nrow(df_ped)) {
    id = df_ped$id[i]
    relatives_1st = 0
    relatives_1st_als = 0
    relatives_1st_ftd = 0
    relatives_1st_ftd_unique = 0
    relatives_1st_dementia = 0
    relatives_1st_dementia_unique = 0
    
    # Check for siblings
    siblings = df_ped %>% filter(pid == df_ped$pid[i] & mid == df_ped$mid[i] & id != df_ped$id[i])
    relatives_1st = relatives_1st + nrow(siblings)
    relatives_1st_als = relatives_1st_als + sum(siblings$Y_ALS == 1)
    relatives_1st_ftd = relatives_1st_ftd + sum(siblings$Y_FTD == 1)
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(siblings$Y_FTD == 1 & siblings$Y_ALS == 0)
    relatives_1st_dementia = relatives_1st_dementia + sum(siblings$Y_dementia == 1)
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(siblings$Y_dementia == 1 & siblings$Y_ALS == 0)
    
    # Check for parents
    parents = df_ped %>% filter(id %in% c(df_ped$pid[i], df_ped$mid[i]))
    relatives_1st = relatives_1st + nrow(parents)
    relatives_1st_als = relatives_1st_als + sum(parents$Y_ALS == 1)
    relatives_1st_ftd = relatives_1st_ftd + sum(parents$Y_FTD == 1)
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(parents$Y_FTD == 1 & parents$Y_ALS == 0)
    relatives_1st_dementia = relatives_1st_dementia + sum(parents$Y_dementia == 1)
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(parents$Y_dementia == 1 & parents$Y_ALS == 0)
    
    # Check for children
    children = df_ped %>% filter(pid == df_ped$id[i] | mid == df_ped$id[i])
    relatives_1st = relatives_1st + nrow(children)
    relatives_1st_als = relatives_1st_als + sum(children$Y_ALS == 1)
    relatives_1st_ftd = relatives_1st_ftd + sum(children$Y_FTD == 1)
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(children$Y_FTD == 1 & children$Y_ALS == 0)
    relatives_1st_dementia = relatives_1st_dementia + sum(children$Y_dementia == 1)
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(children$Y_dementia == 1 & children$Y_ALS == 0)

    # Update the data frame
    df_ped$relatives_1st[i] = relatives_1st
    df_ped$relatives_1st_als[i] = relatives_1st_als
    df_ped$relatives_1st_ftd[i] = relatives_1st_ftd
    df_ped$relatives_1st_ftd_unique[i] = relatives_1st_ftd_unique
    df_ped$relatives_1st_dementia[i] = relatives_1st_dementia
    df_ped$relatives_1st_dementia_unique[i] = relatives_1st_dementia_unique
  }
  return(df_ped)
}

add_2nd_and_3rd_relatives = function(df_ped) {
  # Initialize columns for relative counts
  df_ped$relatives_2nd = 0
  df_ped$relatives_2nd_als = 0
  df_ped$relatives_2nd_ftd = 0
  df_ped$relatives_2nd_ftd_unique = 0
  df_ped$relatives_2nd_dementia = 0
  df_ped$relatives_2nd_dementia_unique = 0
  df_ped$relatives_3rd = 0
  df_ped$relatives_3rd_als = 0
  df_ped$relatives_3rd_ftd = 0
  df_ped$relatives_3rd_ftd_unique = 0
  df_ped$relatives_3rd_dementia = 0
  df_ped$relatives_3rd_dementia_unique = 0
 
  # Loop through each individual
  for(i in 1:nrow(df_ped)) {
    id = df_ped$id[i]
    relatives_2nd = 0
    relatives_2nd_als = 0
    relatives_2nd_ftd = 0
    relatives_2nd_ftd_unique = 0
    relatives_2nd_dementia = 0
    relatives_2nd_dementia_unique = 0
    relatives_3rd = 0
    relatives_3rd_als = 0
    relatives_3rd_ftd = 0    
    relatives_3rd_ftd_unique = 0
    relatives_3rd_dementia = 0    
    relatives_3rd_dementia_unique = 0
    
    # Identify grandparents (second-degree relatives)
    parent_ids = c(df_ped$pid[i], df_ped$mid[i])
    parent_ids = parent_ids[!is.na(parent_ids)]
    grandparent_ids = c()
    for(parent_id in parent_ids) {
      parent_row = df_ped %>% filter(id == parent_id)
      if(nrow(parent_row) > 0) {
        grandparent_pid = parent_row$pid
        grandparent_mid = parent_row$mid
        grandparent_ids = c(grandparent_ids, grandparent_pid, grandparent_mid)
      }
    }
    grandparent_ids = grandparent_ids[!is.na(grandparent_ids)]
    grandparent_rows = df_ped %>% filter(id %in% grandparent_ids)
    if(nrow(grandparent_rows) > 0) {
      relatives_2nd = relatives_2nd + nrow(grandparent_rows)
      relatives_2nd_als = relatives_2nd_als + sum(grandparent_rows$Y_ALS == 1)
      relatives_2nd_ftd = relatives_2nd_ftd + sum(grandparent_rows$Y_FTD == 1)
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(grandparent_rows$Y_FTD == 1 & grandparent_rows$Y_ALS == 0)
      relatives_2nd_dementia = relatives_2nd_dementia + sum(grandparent_rows$Y_dementia == 1)
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(grandparent_rows$Y_dementia == 1 & grandparent_rows$Y_ALS == 0)
      
      # Identify great-grandparents (third-degree relatives)
      great_grandparent_ids = c()
      for(grandparent_id in grandparent_ids) {
        grandparent_row = df_ped %>% filter(id == grandparent_id)
        if(nrow(grandparent_row) > 0) {
          great_grandparent_pid = grandparent_row$pid
          great_grandparent_mid = grandparent_row$mid
          great_grandparent_ids = c(great_grandparent_ids, great_grandparent_pid, great_grandparent_mid)
        }
      }
      great_grandparent_ids = great_grandparent_ids[!is.na(great_grandparent_ids)]
      great_grandparent_rows = df_ped %>% filter(id %in% great_grandparent_ids)
      if(nrow(great_grandparent_rows) > 0) {
        relatives_3rd = relatives_3rd + nrow(great_grandparent_rows)
        relatives_3rd_als = relatives_3rd_als + sum(great_grandparent_rows$Y_ALS == 1)
        relatives_3rd_ftd = relatives_3rd_ftd + sum(great_grandparent_rows$Y_FTD == 1)
        relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(great_grandparent_rows$Y_FTD == 1 & great_grandparent_rows$Y_ALS == 0)
        relatives_3rd_dementia = relatives_3rd_dementia + sum(great_grandparent_rows$Y_dementia == 1)
        relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(great_grandparent_rows$Y_dementia == 1 & great_grandparent_rows$Y_ALS == 0)
      }
    }
    
    # Identify grandchildren (second-degree relatives)
    child_ids = df_ped %>% filter(pid == df_ped$id[i] | mid == df_ped$id[i]) %>% pull(id)
    grandchild_ids = c()
    for(cid in child_ids) {
      grandchild_rows = df_ped %>% filter(pid == cid | mid == cid)
      grandchild_ids = c(grandchild_ids, grandchild_rows$id)
    }
    grandchild_rows = df_ped %>% filter(id %in% grandchild_ids)
    if(nrow(grandchild_rows) > 0) {
      relatives_2nd = relatives_2nd + nrow(grandchild_rows)
      relatives_2nd_als = relatives_2nd_als + sum(grandchild_rows$Y_ALS == 1)
      relatives_2nd_ftd = relatives_2nd_ftd + sum(grandchild_rows$Y_FTD == 1)
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(grandchild_rows$Y_FTD == 1 & grandchild_rows$Y_ALS == 0)
      relatives_2nd_dementia = relatives_2nd_dementia + sum(grandchild_rows$Y_dementia == 1)
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(grandchild_rows$Y_dementia == 1 & grandchild_rows$Y_ALS == 0)
    }
    
    # Identify great-grandchildren (third-degree relatives)
    great_grandchild_ids = c()
    for(gcid in grandchild_ids) {
      great_grandchild_rows = df_ped %>% filter(pid == gcid | mid == gcid)
      great_grandchild_ids = c(great_grandchild_ids, great_grandchild_rows$id)
    }
    great_grandchild_rows = df_ped %>% filter(id %in% great_grandchild_ids)
    if(nrow(great_grandchild_rows) > 0) {
      relatives_3rd = relatives_3rd + nrow(great_grandchild_rows)
      relatives_3rd_als = relatives_3rd_als + sum(great_grandchild_rows$Y_ALS == 1)
      relatives_3rd_ftd = relatives_3rd_ftd + sum(great_grandchild_rows$Y_FTD == 1)
      relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(great_grandchild_rows$Y_FTD == 1 & great_grandchild_rows$Y_ALS == 0)
      relatives_3rd_dementia = relatives_3rd_dementia + sum(great_grandchild_rows$Y_dementia == 1)
      relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(great_grandchild_rows$Y_dementia == 1 & great_grandchild_rows$Y_ALS == 0)
    }
    
    # Identify aunts and uncles (second-degree relatives)
    aunt_uncle_ids = c()
    for(parent_id in parent_ids) {
      parent_row = df_ped %>% filter(id == parent_id)
      if(nrow(parent_row) > 0) {
        parent_pid = parent_row$pid
        parent_mid = parent_row$mid
        
        if(!is.na(parent_pid)) {
          sibling_ids_pid = df_ped %>% filter(pid == parent_pid & id != parent_row$id) %>% pull(id)
          aunt_uncle_ids = c(aunt_uncle_ids, sibling_ids_pid)
        }
        
        if(!is.na(parent_mid)) {
          sibling_ids_mid = df_ped %>% filter(mid == parent_mid & id != parent_row$id) %>% pull(id)
          aunt_uncle_ids = c(aunt_uncle_ids, sibling_ids_mid)
        }
      }
    }
    aunt_uncle_rows = df_ped %>% filter(id %in% aunt_uncle_ids)
    if(nrow(aunt_uncle_rows) > 0) {
      relatives_2nd = relatives_2nd + nrow(aunt_uncle_rows)
      relatives_2nd_als = relatives_2nd_als + sum(aunt_uncle_rows$Y_ALS == 1)
      relatives_2nd_ftd = relatives_2nd_ftd + sum(aunt_uncle_rows$Y_FTD == 1)
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(aunt_uncle_rows$Y_FTD == 1 & aunt_uncle_rows$Y_ALS == 0)
      relatives_2nd_dementia = relatives_2nd_dementia + sum(aunt_uncle_rows$Y_dementia == 1)
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(aunt_uncle_rows$Y_dementia == 1 & aunt_uncle_rows$Y_ALS == 0)
    }
    
    # Identify nephews and nieces (second-degree relatives)
    sibling_ids = c()
    if(!is.na(df_ped$pid[i])) {
      sibling_ids_pid = df_ped %>% filter(pid == df_ped$pid[i] & id != df_ped$id[i]) %>% pull(id)
      sibling_ids = c(sibling_ids, sibling_ids_pid)
    }
    if(!is.na(df_ped$mid[i])) {
      sibling_ids_mid = df_ped %>% filter(mid == df_ped$mid[i] & id != df_ped$id[i]) %>% pull(id)
      sibling_ids = c(sibling_ids, sibling_ids_mid)
    }
    
    nephew_niece_ids = c()
    for(sib in sibling_ids) {
      nephew_niece_rows = df_ped %>% filter(pid == sib | mid == sib)
      nephew_niece_ids = c(nephew_niece_ids, nephew_niece_rows$id)
    }
    nephew_niece_rows = df_ped %>% filter(id %in% nephew_niece_ids)
    if(nrow(nephew_niece_rows) > 0) {
      relatives_2nd = relatives_2nd + nrow(nephew_niece_rows)
      relatives_2nd_als = relatives_2nd_als + sum(nephew_niece_rows$Y_ALS == 1)
      relatives_2nd_ftd = relatives_2nd_ftd + sum(nephew_niece_rows$Y_FTD == 1)
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(nephew_niece_rows$Y_FTD == 1 & nephew_niece_rows$Y_ALS == 0)
      relatives_2nd_dementia = relatives_2nd_dementia + sum(nephew_niece_rows$Y_dementia == 1)
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(nephew_niece_rows$Y_dementia == 1 & nephew_niece_rows$Y_ALS == 0)

      # Identify grand-nephews and grand-nieces (third-degree relatives)
      grand_nephew_niece_ids = c()
      for(nephew_niece_id in nephew_niece_ids) {
        grand_nephew_niece_rows = df_ped %>% filter(pid == nephew_niece_id | mid == nephew_niece_id)
        grand_nephew_niece_ids = c(grand_nephew_niece_ids, grand_nephew_niece_rows$id)
      }
      grand_nephew_niece_rows = df_ped %>% filter(id %in% grand_nephew_niece_ids)
      if(nrow(grand_nephew_niece_rows) > 0) {
        relatives_3rd = relatives_3rd + nrow(grand_nephew_niece_rows)
        relatives_3rd_als = relatives_3rd_als + sum(grand_nephew_niece_rows$Y_ALS == 1)
        relatives_3rd_ftd = relatives_3rd_ftd + sum(grand_nephew_niece_rows$Y_FTD == 1)
        relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(grand_nephew_niece_rows$Y_FTD == 1 & grand_nephew_niece_rows$Y_ALS == 0)
        relatives_3rd_dementia = relatives_3rd_dementia + sum(grand_nephew_niece_rows$Y_dementia == 1)
        relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(grand_nephew_niece_rows$Y_dementia == 1 & grand_nephew_niece_rows$Y_ALS == 0)
      }
    }
    
    # Identify first cousins (third-degree relatives)
parent_ids = c(df_ped$pid[i], df_ped$mid[i])
parent_ids = parent_ids[!is.na(parent_ids)]
aunt_uncle_ids = c()
for(parent_id in parent_ids) {
  parent_row = df_ped %>% filter(id == parent_id)
  if(nrow(parent_row) > 0) {
    parent_pid = parent_row$pid
    parent_mid = parent_row$mid
    
    if(!is.na(parent_pid)) {
      sibling_ids_pid = df_ped %>% filter(pid == parent_pid & id != parent_row$id) %>% pull(id)
      aunt_uncle_ids = c(aunt_uncle_ids, sibling_ids_pid)
    }
    
    if(!is.na(parent_mid)) {
      sibling_ids_mid = df_ped %>% filter(mid == parent_mid & id != parent_row$id) %>% pull(id)
      aunt_uncle_ids = c(aunt_uncle_ids, sibling_ids_mid)
    }
  } else {
    # Handle missing parent
  }
}
first_cousin_ids = c()
for(aunt_uncle_id in aunt_uncle_ids) {
  first_cousin_rows = df_ped %>% filter(pid == aunt_uncle_id | mid == aunt_uncle_id)
  first_cousin_ids = c(first_cousin_ids, first_cousin_rows$id)
}
first_cousin_rows = df_ped %>% filter(id %in% first_cousin_ids)
if(nrow(first_cousin_rows) > 0) {
  relatives_3rd = relatives_3rd + nrow(first_cousin_rows)
  relatives_3rd_als = relatives_3rd_als + sum(first_cousin_rows$Y_ALS == 1)
  relatives_3rd_ftd = relatives_3rd_ftd + sum(first_cousin_rows$Y_FTD == 1)
  relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(first_cousin_rows$Y_FTD == 1 & first_cousin_rows$Y_ALS == 0)
  relatives_3rd_dementia = relatives_3rd_dementia + sum(first_cousin_rows$Y_dementia == 1)
  relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(first_cousin_rows$Y_dementia == 1 & first_cousin_rows$Y_ALS == 0)
}

# Update the data frame
df_ped$relatives_2nd[i] = relatives_2nd
df_ped$relatives_2nd_als[i] = relatives_2nd_als
df_ped$relatives_2nd_ftd[i] = relatives_2nd_ftd
df_ped$relatives_2nd_ftd_unique[i] = relatives_2nd_ftd_unique
df_ped$relatives_2nd_dementia[i] = relatives_2nd_dementia
df_ped$relatives_2nd_dementia_unique[i] = relatives_2nd_dementia_unique
df_ped$relatives_3rd[i] = relatives_3rd
df_ped$relatives_3rd_als[i] = relatives_3rd_als
df_ped$relatives_3rd_ftd[i] = relatives_3rd_ftd
df_ped$relatives_3rd_ftd_unique[i] = relatives_3rd_ftd_unique
df_ped$relatives_3rd_dementia[i] = relatives_3rd_dementia
df_ped$relatives_3rd_dementia_unique[i] = relatives_3rd_dementia_unique
  }
  return(df_ped)
}


#' wrapper function to simulate full pedigree
sim_ped = function(i, k, lambda=NA, yob_index=1960, mean_gen_yr, fert_rate, 
                   DAF_common, DAF_patho, DAF_ftd, penetrance_ALS_common, penetrance_ALS_patho, 
                   K_ALS, h2_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, K_FTD, h2_FTD, K_dementia, h2_dementia, rg, 
                   life_expectancy, current_year, plot=TRUE){
  core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)   
  core_ped = add_pheno(core_ped, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, h2_dementia, K_dementia, rg) 
  core_ped = add_1st_relatives(core_ped)
  core_ped = add_2nd_and_3rd_relatives(core_ped)
  if(plot){
    last_gen_affected = filter(core_ped, gen == k & (Y_ALS == 1 | Y_dementia ==1))
    if(nrow(last_gen_affected) > 0) {
    ped_plt = ped(id = core_ped$id, 
                  fid = core_ped$pid,
                  mid = core_ped$mid,
                  sex = core_ped$sex + 1,
                  isConnected = TRUE,
                  reorder = FALSE)
    carriers = filter(core_ped, a1_common + a2_common + a1_patho + a2_patho > 0)$id
    affected_ALS = filter(core_ped, Y_ALS == 1)$id
    affected_FTD = filter(core_ped, Y_FTD == 1)$id
    affected_dementia = filter(core_ped, Y_dementia_other == 1)$id
    affected_either = filter(core_ped, Y_ALS == 1 | Y_dementia == 1)$id    
    # color according to percentile from N(0,1)
    pal_ALS = colorRampPalette(c("white", "orange", "red"))(1000)
    pal_FTD = colorRampPalette(c("white", "skyblue", "blue"))(1000)
    pal_dementia = colorRampPalette(c("white", "lightgreen", "darkgreen"))(1000)    
    Pcolors_ALS = pal_ALS[ceiling(pnorm(core_ped$P_ALS)*1000)]
    Pcolors_FTD = pal_FTD[ceiling(pnorm(core_ped$P_FTD)*1000)]
    Pcolors_dementia = pal_dementia[ceiling(pnorm(core_ped$P_dementia)*1000)]    
   
    # make plot
    plot(ped_plt, title="Polygenic phenotype for ALS on liability scale", 
    carrier = carriers, fill=Pcolors_ALS, cex=0.8)
    plot(ped_plt, title="Polygenic phenotype for FTD on liability scale", 
    carrier = carriers, fill=Pcolors_FTD, cex=0.8)
    plot(ped_plt, title="Polygenic phenotype for other dementias on liability scale", 
    carrier = carriers, fill=Pcolors_dementia, cex=0.8)    
    plot(ped_plt, title="Binary ALS phenotype", 
    carrier = carriers, aff = affected_ALS, cex=0.8)
    plot(ped_plt, title="Binary FTD phenotype", 
    carrier = carriers, aff = affected_FTD, cex=0.8)
    plot(ped_plt, title="Binary dementia phenotype", 
    carrier = carriers, aff = affected_dementia, cex=0.8)
    plot(ped_plt, title="Binary phenotype for either ALS or dementia (incl FTD)", 
    carrier = carriers, aff = affected_either, cex=0.8)
  }}
  return(core_ped)
}

#' final function to run a desired number of pedigree simulations
pedigree_simulations = function(n_simulations, k=k, lambda=NA, mean_gen_yr=mean_gen_yr, fert_rate,
                                DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, penetrance_ALS_common=penetrance_ALS_common, 
                                penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, penetrance_FTD_common=penetrance_FTD_common, 
                                penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS, K_FTD=K_FTD, h2_FTD=h2_FTD, K_dementia=K_dementia, h2_dementia=h2_dementia, 
                                life_expectancy=life_expectancy, current_year=current_year, rg=rg, plot=FALSE) {
  # Initialize data frame to store results
  results_df = data.frame(
    monogenic = numeric(),
    monogenic_source = character(),
    id = character(),
    relatives_1st = numeric(),
    relatives_1st_als = numeric(),
    relatives_1st_ftd = numeric(),
    relatives_1st_ftd_unique = numeric(),
    relatives_1st_dementia = numeric(),
    relatives_1st_dementia_unique = numeric(),
    relatives_2nd = numeric(),
    relatives_2nd_als = numeric(),
    relatives_2nd_ftd = numeric(),
    relatives_2nd_ftd_unique = numeric(),
    relatives_2nd_dementia = numeric(),
    relatives_2nd_dementia_unique = numeric(),
    relatives_3rd = numeric(),
    relatives_3rd_als = numeric(),
    relatives_3rd_ftd = numeric(),
    relatives_3rd_ftd_unique = numeric(),
    relatives_3rd_dementia = numeric(),
    relatives_3rd_dementia_unique = numeric(),
    sex = numeric(),
    a1_common = numeric(),
    a2_common = numeric(),
    a1_patho = numeric(),
    a2_patho = numeric(),    
    a1_ftd = numeric(),
    a2_ftd = numeric(),
    yob = numeric(),
    age = numeric(),
    age_censored = numeric(),
    status = character(),
    mendel_ALS_y1_common = numeric(),
    mendel_ALS_y2_common = numeric(),
    mendel_FTD_y1_common = numeric(),
    mendel_FTD_y2_common = numeric(),
    mendel_ALS_y1_patho = numeric(),
    mendel_ALS_y2_patho = numeric(),
    mendel_FTD_y1_patho = numeric(),
    mendel_FTD_y2_patho = numeric(),
    mendel_FTD_y1_ftd = numeric(),
    mendel_FTD_y2_ftd = numeric(),    
    mendel_ALS_Y_common = numeric(),
    mendel_ALS_Y_patho = numeric(),
    mendel_ALS_Y = numeric(),
    mendel_FTD_Y_common = numeric(),
    mendel_FTD_Y_patho = numeric(),
    mendel_FTD_Y_ftd = numeric(),
    mendel_FTD_Y = numeric(),
    G_ALS = numeric(),
    G_FTD = numeric(),
    G_dementia = numeric(),
    E_ALS = numeric(),
    E_FTD = numeric(),
    E_dementia = numeric(),
    P_ALS = numeric(),
    P_FTD = numeric(),
    P_dementia = numeric(),
    age_related_K_dementia = numeric(),
    LT_dem = numeric(),
    polygenicY_ALS = numeric(),
    polygenicY_FTD = numeric(),
    polygenicY_dementia = numeric(),
    Y_ALS = numeric(),
    Y_FTD = numeric(),
    Y_dementia_other= numeric(),
    Y_dementia = numeric(),
    gen = numeric() 
  )
  
  print((paste("Parameters for this run of simulations are as follows:")))
  print((paste("Generations simulated:", (k+1))))
  print((paste("Fertility rate is:", max(fert_rate$mean_fertility))))
  print((paste("Max lifetime expectancy is:", max(life_expectancy$life_expectancy))))
  print((paste("Disease allele frequency common variant:", DAF_common)))
  print((paste("Cumulative ALS penetrance of common variant:", max(penetrance_ALS_common$cum_incidence_c9_ALS))))
  print((paste("Cumulative FTD penetrance of common variant:", max(penetrance_FTD_common$cum_incidence_c9_FTD))))
  print((paste("Disease allele frequency rare, pathogenic variant:", DAF_patho)))
  print((paste("Disease allele ftd specific variants:", DAF_ftd)))
  print((paste("ALS lifetime risk:", max(K_ALS$cum_risk_adjusted))))
  print((paste("FTD lifetime risk:", max(K_FTD$cum_risk_adjusted))))
  print((paste("Dementia lifetime risk:", max(K_dementia$cum_incidence, na.rm = TRUE))))
  print((paste("Heritability of ALS:", h2_ALS)))
  print((paste("Heritability of FTD:", h2_FTD)))
  print((paste("Heritability of other dementias:", h2_dementia)))
  print((paste("Genetic correlation between ALS and FTD:", rg)))


  # Create simulation metrics file 
  simulation_metrics = data.frame(
    Total_Simulations = numeric(),
    Total_Individuals = numeric(),
    Simulated_ALS_LifetimeRisk = numeric(),
    Total_ALS_Cases = numeric(),
    Total_Mendelian_ALS_Cases = numeric(),
    Total_Mendelian_ALS_Cases_common = numeric(),
    Total_Mendelian_ALS_Cases_patho = numeric(),
    Mendelian_ALS_Percentage = numeric(),
    Total_Polygenic_ALS_Cases = numeric(),
    Polygenic_ALS_Percentage = numeric(),
    Simulated_FTD_LifetimeRisk = numeric(),
    Total_FTD_Cases = numeric(),
    Total_Mendelian_FTD_Cases = numeric(),
    Total_Mendelian_FTD_Cases_common = numeric(),
    Total_Mendelian_FTD_Cases_patho = numeric(),
    Total_Mendelian_FTD_Cases_ftd = numeric(),
    Mendelian_FTD_Percentage = numeric(),
    Total_Polygenic_FTD_Cases = numeric(),
    Polygenic_FTD_Percentage = numeric(),  
    Total_Carriers = numeric(),
    Total_Penetrance_ALS = numeric(),
    Total_Penetrance_FTD = numeric(),
    Total_ALS_FTD_Cases = numeric(),
    Total_Dementia_Cases = numeric(),
    Total_other_dementia_Cases = numeric(),
    Simulated_other_dementia_LifetimeRisk = numeric(),
    LastGen_Mendelian_ALS_Cases = numeric(),
    LastGen_Mendelian_FTD_Cases = numeric(),
    LastGen_Polygenic_ALS_Cases = numeric(),
    LastGen_Polygenic_FTD_Cases = numeric(),
    LastGen_Dementia_Cases = numeric(),
    LastGen_Carriers = numeric(),
    LastGen_Sporadic = numeric(),
    LastGen_Sporadic_Mendel = numeric(),
    LastGen_Sporadic_Mendel_perc = numeric()
  )
  
  # Initialize variables to track total metrics    
  total_simulations = 0
  total_individuals = 0
  total_als_cases = 0
  total_risk_alleles_ALS = 0
  total_risk_alleles_FTD = 0
  mendelian_als_cases = 0
  mendelian_als_cases_common = 0
  mendelian_als_cases_patho = 0
  polygenic_als_cases = 0
  penetrant_risk_alleles_ALS = 0
  total_ftd_cases = 0
  mendelian_ftd_cases = 0
  mendelian_ftd_cases_common = 0
  mendelian_ftd_cases_patho = 0  
  mendelian_ftd_cases_ftd = 0  
  polygenic_ftd_cases = 0  
  penetrant_risk_alleles_FTD = 0
  total_als_ftd_cases = 0
  total_other_dementia_cases = 0
  total_dementia_cases = 0

  # Initialize last generation counters
  mendelian_ALS_individuals = 0
  polygenic_ALS_individuals = 0
  mendelian_FTD_individuals = 0
  polygenic_FTD_individuals = 0 
  dementia_individuals = 0 
  carrier_individuals = 0
  isolated_als_with_risk_alleles = 0
  sporadic_appearing_cases = 0

  # Run simulations
  while(total_simulations < n_simulations) {
    print(paste("Running simulation", total_simulations + 1, "of", n_simulations, "..."))
    yob_index = sample(1960:1975, 1)
    print(paste("Year of birth index is:", yob_index))
    
    # Run simulation with error handling
    tryCatch(
      expr = {
        ped = sim_ped(total_simulations + 1, k=k, lambda=lambda, yob_index=yob_index, mean_gen_yr=mean_gen_yr, 
                      fert_rate=fert_rate, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd,
                      penetrance_ALS_common=penetrance_ALS_common, penetrance_ALS_patho=penetrance_ALS_patho, K_ALS=K_ALS, h2_ALS=h2_ALS, 
                      penetrance_FTD_common=penetrance_FTD_common, penetrance_FTD_patho=penetrance_FTD_patho, penetrance_FTD_nonALS=penetrance_FTD_nonALS,K_FTD=K_FTD, h2_FTD=h2_FTD, 
                      K_dementia=K_dementia, h2_dementia=h2_dementia, 
                      life_expectancy=life_expectancy, current_year=current_year, rg=rg, plot=plot) 
        # Check if ped contains required columns
                    required_columns = c("id", "pid", "mid", "relatives_1st", "relatives_1st_als", "relatives_1st_ftd", 
                     "relatives_1st_dementia", "relatives_2nd", "relatives_2nd_als", "relatives_2nd_ftd", 
                     "relatives_2nd_dementia", "relatives_3rd", "relatives_3rd_als", "relatives_3rd_ftd", 
                     "relatives_3rd_dementia", "sex", "a1_common", "a2_common", "a1_patho", "a2_patho", "a1_ftd", "a2_ftd", "yob", "age", "life_expectancy", "age_censored", 
                     "status", "mendel_ALS_y1_common", "mendel_ALS_y2_common", "mendel_ALS_Y_common", "mendel_ALS_y1_patho", "mendel_ALS_y2_patho", "mendel_ALS_Y_patho",
                     "mendel_ALS_Y", "mendel_FTD_y1_common", "mendel_FTD_y2_common", "mendel_FTD_y1_patho", "mendel_FTD_y2_patho", "mendel_FTD_y1_ftd", "mendel_FTD_y2_ftd",  
                     "mendel_FTD_Y_common", "mendel_FTD_Y_patho", "mendel_FTD_Y_ftd", "mendel_FTD_Y", "G_ALS", "G_FTD", "G_dementia", 
                     "E_ALS", "E_FTD", "E_dementia", "P_ALS", "P_FTD", "P_dementia", "polygenicY_ALS", "polygenicY_FTD", "polygenicY_dementia", 
                     "Y_ALS", "Y_FTD", "Y_dementia_other", "Y_dementia", "gen")
        if(!all(required_columns %in% colnames(ped))) {
          print(paste("Simulation", total_simulations + 1, "is missing required columns. Skipping..."))
          next
        }
        
        # Update total counts
        total_individuals = total_individuals + nrow(ped)
 
        # Loop through all individuals in the pedigree
        for(i in 1:nrow(ped)) {
          row = ped[i,]
          
          # Check for ALS cases
          if(row$Y_ALS == 1) {
            total_als_cases = total_als_cases + 1
            
            # Check for Mendelian ALS inheritance
            if(row$mendel_ALS_Y == 1) {
              mendelian_als_cases = mendelian_als_cases + 1
              if(row$mendel_ALS_Y_common == 1) {
                mendelian_als_cases_common = mendelian_als_cases_common + 1
              } else if(row$mendel_ALS_Y_patho == 1) {
                mendelian_als_cases_patho = mendelian_als_cases_patho + 1
              }

            } else if(row$polygenicY_ALS == 1) {
              polygenic_als_cases = polygenic_als_cases + 1
            }
          }
          # Check for FTD cases
          if(row$Y_FTD == 1) {
            total_ftd_cases = total_ftd_cases + 1
            
            # Check for Mendelian ALS inheritance
            if(row$mendel_FTD_Y == 1) {
              mendelian_ftd_cases = mendelian_ftd_cases + 1
              if(row$mendel_FTD_Y_common == 1) {
                mendelian_ftd_cases_common = mendelian_ftd_cases_common + 1
              } else if(row$mendel_FTD_Y_patho == 1) {
                mendelian_ftd_cases_patho = mendelian_ftd_cases_patho + 1
              } else if(row$mendel_FTD_Y_ftd == 1) {
                mendelian_ftd_cases_ftd = mendelian_ftd_cases_ftd + 1
              }
            } else if(row$polygenicY_FTD == 1) {
              polygenic_ftd_cases = polygenic_ftd_cases + 1
            }
          }

          # Check for ALS/FTD co-occurrence
          if(row$Y_ALS == 1 & row$Y_FTD == 1) {
            total_als_ftd_cases = total_als_ftd_cases + 1
          }
          
          # Check for other dementia cases 
          if(row$Y_dementia == 1) {
            total_dementia_cases = total_dementia_cases+1
            if(row$Y_dementia_other == 1) { 
              total_other_dementia_cases = total_other_dementia_cases+1
            }
          }

          # Check for risk alleles
          if(row$a1_common == 1 | row$a2_common == 1 | row$a1_patho == 1 | row$a2_patho == 1) {
            total_risk_alleles_ALS = total_risk_alleles_ALS + 1
            if(row$mendel_ALS_Y == 1) {
              penetrant_risk_alleles_ALS = penetrant_risk_alleles_ALS + 1
            }
          }
          if(row$a1_common == 1 | row$a2_common == 1 | row$a1_patho == 1 | row$a2_patho == 1) {
            total_risk_alleles_FTD = total_risk_alleles_FTD + 1
            if(row$mendel_FTD_Y == 1) {
              penetrant_risk_alleles_FTD = penetrant_risk_alleles_FTD + 1
            }
          }
        }
             print("Finished loop through individuals")
 
        # Increment total simulations counter
        total_simulations = total_simulations + 1
        
        # Filter for affected phenotype in the last generation
        last_gen_cases = ped %>% filter(gen == k & (Y_ALS))
        
        if(nrow(last_gen_cases) == 0) {
          print(paste("Simulation", total_simulations, "has no ALS cases in the last generation. Skipping..."))
          next
        }
        print(str(last_gen_cases))
        # Check for monogenic inheritance in any form
        monogenic_inheritance = ifelse(any(ped$mendel_ALS_Y == 1 | ped$mendel_FTD_Y == 1), TRUE, FALSE)
        
        # Check for affected individuals in the last generation
        mendelian_ALS_last_gen = ped %>% filter(gen == k & mendel_ALS_Y == 1)
        mendelian_ALS_individuals = mendelian_ALS_individuals + nrow(mendelian_ALS_last_gen)

        mendelian_FTD_last_gen = ped %>% filter(gen == k & mendel_FTD_Y == 1)
        mendelian_FTD_individuals = mendelian_FTD_individuals + nrow(mendelian_FTD_last_gen)

        # Separate counts for polygenic ALS and FTD
        polygenic_ALS_last_gen = ped %>% filter(gen == k & polygenicY_ALS == 1)
        polygenic_ALS_individuals = polygenic_ALS_individuals + nrow(polygenic_ALS_last_gen)

        polygenic_FTD_last_gen = ped %>% filter(gen == k & polygenicY_FTD == 1)
        polygenic_FTD_individuals = polygenic_FTD_individuals + nrow(polygenic_FTD_last_gen)    

        # Counts for other dementias and any dementia (=FTD+others) in last gen
        dementia_last_gen = ped %>% filter(gen == k & Y_dementia == 1)
        dementia_individuals = dementia_individuals + nrow(dementia_last_gen)

        # Check for carriers in the last generation who are not affected
        carriers_last_gen = ped %>% filter(gen == k & (a1_common == 1 | a2_common == 1 | a1_patho == 1 | a2_patho == 1) & (Y_ALS == 0 & Y_FTD == 0))
        carrier_individuals = carrier_individuals + nrow(carriers_last_gen)
        
        # Select one affected individual per simulation
        selected_case = last_gen_cases[startsWith(last_gen_cases$id, "C"), ][1,]
        if(is.na(selected_case$id) || nrow(selected_case) == 0) {
          print(paste("Simulation", total_simulations, "has no ALS cases in the founder pedigree. Skipping..."))
          next
        }

        # Extract relevant information
        id = paste0("sim", total_simulations, "_", selected_case$id)
        
        relatives_1st = selected_case$relatives_1st
        relatives_1st_als = selected_case$relatives_1st_als
        relatives_1st_ftd = selected_case$relatives_1st_ftd
        relatives_1st_ftd_unique = selected_case$relatives_1st_ftd_unique
        relatives_1st_dementia = selected_case$relatives_1st_dementia
        relatives_1st_dementia_unique = selected_case$relatives_1st_dementia_unique
        relatives_2nd = selected_case$relatives_2nd
        relatives_2nd_als = selected_case$relatives_2nd_als
        relatives_2nd_ftd = selected_case$relatives_2nd_ftd
        relatives_2nd_ftd_unique = selected_case$relatives_2nd_ftd_unique
        relatives_2nd_dementia = selected_case$relatives_2nd_dementia
        relatives_2nd_dementia_unique = selected_case$relatives_2nd_dementia_unique
        relatives_3rd = selected_case$relatives_3rd
        relatives_3rd_als = selected_case$relatives_3rd_als
        relatives_3rd_ftd = selected_case$relatives_3rd_ftd
        relatives_3rd_ftd_unique = selected_case$relatives_3rd_ftd_unique
        relatives_3rd_dementia = selected_case$relatives_3rd_dementia
        relatives_3rd_dementia_unique = selected_case$relatives_3rd_dementia_unique
        sex = selected_case$sex
        a1_common = selected_case$a1_common
        a2_common = selected_case$a2_common
        a1_patho = selected_case$a1_patho
        a2_patho = selected_case$a2_patho
        a1_ftd = selected_case$a1_ftd
        a2_ftd = selected_case$a2_ftd        
        yob = selected_case$yob
        age = selected_case$age
        age_censored = selected_case$age_censored
        status = selected_case$status
        mendel_ALS_y1_common = selected_case$mendel_ALS_y1_common
        mendel_ALS_y2_common = selected_case$mendel_ALS_y2_common
        mendel_ALS_y1_patho = selected_case$mendel_ALS_y1_patho
        mendel_ALS_y2_patho = selected_case$mendel_ALS_y2_patho        
        mendel_ALS_Y_common = selected_case$mendel_ALS_Y_common
        mendel_ALS_Y_patho = selected_case$mendel_ALS_Y_patho
        mendel_ALS_Y = selected_case$mendel_ALS_Y
        mendel_FTD_y1_common = selected_case$mendel_FTD_y1_common
        mendel_FTD_y2_common = selected_case$mendel_FTD_y2_common
        mendel_FTD_Y_common = selected_case$mendel_FTD_Y_common
        mendel_FTD_y1_patho = selected_case$mendel_FTD_y1_patho
        mendel_FTD_y2_patho = selected_case$mendel_FTD_y2_patho
        mendel_FTD_y1_ftd = selected_case$mendel_FTD_y1_ftd
        mendel_FTD_y2_ftd = selected_case$mendel_FTD_y2_ftd
        mendel_FTD_Y_patho = selected_case$mendel_FTD_Y_patho 
        mendel_FTD_Y_ftd = selected_case$mendel_FTD_Y_ftd 
        mendel_FTD_Y = selected_case$mendel_FTD_Y
        G_ALS = selected_case$G_ALS
        G_FTD = selected_case$G_FTD
        G_dementia = selected_case$G_dementia
        E_ALS = selected_case$E_ALS
        E_FTD = selected_case$E_FTD
        E_dementia = selected_case$E_dementia
        P_ALS = selected_case$P_ALS
        P_FTD = selected_case$P_FTD
        P_dementia = selected_case$P_dementia
        age_related_K_dementia = selected_case$age_related_K_dementia
        LT_dem = selected_case$LT_dem
        polygenicY_ALS = selected_case$polygenicY_ALS
        polygenicY_FTD = selected_case$polygenicY_FTD
        polygenicY_dementia = selected_case$polygenicY_dementia
        Y_ALS = selected_case$Y_ALS
        Y_FTD = selected_case$Y_FTD
        Y_dementia_other = selected_case$Y_dementia_other
        Y_dementia = selected_case$Y_dementia
        gen = selected_case$gen

        # Check for monogenic inheritance source
        if(monogenic_inheritance) {
          if(any(ped$id == "C0" & (ped$a1_common == 1 | ped$a2_common == 1 | ped$a1_patho == 1 | ped$a2_patho == 1))) {
            monogenic_source = "Disease allele came from founder (C0)"
          } else {
            monogenic_source = "Disease allele came from extended branches/inlaws"
          }
        } else {
          monogenic_source = "No monogenic affected invididuals in pedigree"
        }
        
        # Check for seemingly sporadic pedigrees with carriers in relatives
        if (!is.na(selected_case$id) && nrow(selected_case) > 0) {
        # Check if the selected case is the only one with ALS, FTD, or dementia
        is_isolated <- relatives_1st_als == 0 && relatives_2nd_als == 0 && relatives_3rd_als == 0
        if (is_isolated) {
            sporadic_appearing_cases <- sporadic_appearing_cases + 1
        }
        # Check if there are relatives with risk alleles
        relatives_with_risk_alleles <- any((ped$a1_common == 1 | ped$a2_common == 1 | 
                                    ped$a1_patho == 1 | ped$a2_patho == 1))
        
        # Increment counter if conditions are met 
        if (mendel_ALS_Y == 1 && is_isolated && relatives_with_risk_alleles) {
          isolated_als_with_risk_alleles <- isolated_als_with_risk_alleles + 1
        } 
        } else {
          print("No valid selected_case found, skipping sporadic pedigree check.")
        }

        # Check for missing values
        if(is.na(id) | is.na(sex) | is.na(a1_common) | is.na(a2_common) | is.na(a1_patho) | is.na(a2_patho) | is.na(yob) | 
           is.na(mendel_ALS_y1_common) | is.na(mendel_ALS_y2_common) | is.na(mendel_ALS_Y_common) | 
           is.na(mendel_ALS_y1_patho) | is.na(mendel_ALS_y2_patho) | is.na(mendel_ALS_Y_patho) | 
           is.na(mendel_FTD_y1_patho) | is.na(mendel_FTD_y2_patho) | is.na(mendel_FTD_Y_patho) | 
           is.na(mendel_FTD_y1_common) | is.na(mendel_FTD_y2_common) | is.na(mendel_FTD_Y_common) | is.na(mendel_FTD_Y_ftd) | 
           is.na(G_ALS) | is.na(G_FTD) | is.na(G_dementia) | 
           is.na(E_ALS) | is.na(E_FTD) | is.na(E_dementia) | 
           is.na(P_ALS) | is.na(P_FTD) | is.na(P_dementia) | 
           is.na(polygenicY_ALS) | is.na(polygenicY_FTD) | is.na(polygenicY_dementia) | 
           is.na(Y_ALS) | is.na(Y_FTD) | is.na(Y_dementia_other) | is.na(Y_dementia) | 
           is.na(gen) | 
           is.na(relatives_1st) | is.na(relatives_1st_als) | is.na(relatives_1st_ftd) | is.na(relatives_1st_dementia) | 
           is.na(relatives_2nd) | is.na(relatives_2nd_als) | is.na(relatives_2nd_ftd) | is.na(relatives_2nd_dementia) | 
           is.na(relatives_3rd) | is.na(relatives_3rd_als) | is.na(relatives_3rd_ftd) | is.na(relatives_3rd_dementia)) {
          print(paste("Simulation", total_simulations, "has missing values for individual", id, ". Skipping..."))
          next
        }
      print("Adding to results data frame...")
        
        # Add to results data frame
        results_df = rbind(results_df, data.frame(
          monogenic = ifelse(monogenic_inheritance, 1, 0),
          monogenic_source = monogenic_source,
          id = id,
          relatives_1st = relatives_1st,
          relatives_1st_als = relatives_1st_als,
          relatives_1st_ftd = relatives_1st_ftd,
          relatives_1st_ftd_unique = relatives_1st_ftd_unique,
          relatives_1st_dementia = relatives_1st_dementia,
          relatives_1st_dementia_unique = relatives_1st_dementia_unique,
          relatives_2nd = relatives_2nd,
          relatives_2nd_als = relatives_2nd_als,
          relatives_2nd_ftd = relatives_2nd_ftd,
          relatives_2nd_ftd_unique = relatives_2nd_ftd_unique,
          relatives_2nd_dementia = relatives_2nd_dementia,
          relatives_2nd_dementia_unique = relatives_2nd_dementia_unique,
          relatives_3rd = relatives_3rd,
          relatives_3rd_als = relatives_3rd_als,
          relatives_3rd_ftd = relatives_3rd_ftd,
          relatives_3rd_ftd_unique = relatives_3rd_ftd_unique,
          relatives_3rd_dementia = relatives_3rd_dementia,
          relatives_3rd_dementia_unique = relatives_3rd_dementia_unique,
          sex = sex,
          a1_common = a1_common,
          a2_common = a2_common,
          a1_patho = a1_patho,
          a2_patho = a2_patho,          
          a1_ftd = a1_ftd,
          a2_ftd = a2_ftd,
          yob = yob,
          age = age,
          age_censored = age_censored,
          status = status,
          mendel_ALS_y1_common = mendel_ALS_y1_common,
          mendel_ALS_y2_common = mendel_ALS_y2_common,
          mendel_ALS_y1_patho = mendel_ALS_y1_patho,
          mendel_ALS_y2_patho = mendel_ALS_y2_patho,      
          mendel_ALS_Y_common = mendel_ALS_Y_common,
          mendel_ALS_Y_patho = mendel_ALS_Y_patho,              
          mendel_ALS_Y = mendel_ALS_Y,
          mendel_FTD_y1_common = mendel_FTD_y1_common,
          mendel_FTD_y2_common = mendel_FTD_y2_common,
          mendel_FTD_y1_patho = mendel_FTD_y1_patho,
          mendel_FTD_y2_patho = mendel_FTD_y2_patho,      
          mendel_FTD_y1_ftd = mendel_FTD_y1_ftd,
          mendel_FTD_y2_ftd = mendel_FTD_y2_ftd,    
          mendel_FTD_Y_common = mendel_FTD_Y_common,
          mendel_FTD_Y_patho = mendel_FTD_Y_patho,   
          mendel_FTD_Y_ftd = mendel_FTD_Y_ftd,
          mendel_FTD_Y = mendel_FTD_Y,
          G_ALS = G_ALS,
          G_FTD = G_FTD,
          G_dementia = G_dementia,
          E_ALS = E_ALS,
          E_FTD = E_FTD,
          E_dementia = E_dementia,
          P_ALS = P_ALS,
          P_FTD = P_FTD,
          P_dementia = P_dementia,
          age_related_K_dementia = age_related_K_dementia,
          LT_dem = LT_dem,
          polygenicY_ALS = polygenicY_ALS,
          polygenicY_FTD = polygenicY_FTD,
          polygenicY_dementia = polygenicY_dementia,
          Y_ALS = Y_ALS,
          Y_FTD = Y_FTD,
          Y_dementia_other = Y_dementia_other,
          Y_dementia = Y_dementia,
          gen = gen
        ))
      },
      error = function(e) {
        print(paste("Simulation", total_simulations + 1, "failed with error:", e$message))
      }
    )
  }

  # Calculate metrics
  simulated_lifetime_als_risk = total_als_cases / total_individuals
  simulated_lifetime_ftd_risk = total_ftd_cases / total_individuals
  simulated_lifetime_otherdementia_risk = total_other_dementia_cases / total_individuals
  mendelian_als_percentage = mendelian_als_cases / total_als_cases
  polygenic_als_percentage = polygenic_als_cases / total_als_cases
  mendelian_ftd_percentage = mendelian_ftd_cases / total_ftd_cases
  polygenic_ftd_percentage = polygenic_ftd_cases / total_ftd_cases
  penetrance_percentage_ALS = penetrant_risk_alleles_ALS / total_risk_alleles_ALS
  penetrance_percentage_FTD = penetrant_risk_alleles_FTD / total_risk_alleles_FTD
  perc_seemingly_sporadic <- if(sporadic_appearing_cases > 0) {
    isolated_als_with_risk_alleles / sporadic_appearing_cases
  } else {
    0
  }
  
  # Save simulation metrics to a CSV file
  simulation_metrics = data.frame(
  Total_Simulations = total_simulations,
  Total_Individuals = total_individuals,
  Simulated_ALS_LifetimeRisk = simulated_lifetime_als_risk,
  Total_ALS_Cases = total_als_cases,
  Total_Mendelian_ALS_Cases = mendelian_als_cases,
  Total_Mendelian_ALS_Cases_common = mendelian_als_cases_common,
  Total_Mendelian_ALS_Cases_patho = mendelian_als_cases_patho,
  Mendelian_ALS_Percentage = mendelian_als_percentage,
  Total_Polygenic_ALS_Cases = polygenic_als_cases,
  Polygenic_ALS_Percentage = polygenic_als_percentage,
  Simulated_FTD_LifetimeRisk = simulated_lifetime_ftd_risk,
  Total_FTD_Cases = total_ftd_cases,
  Total_Mendelian_FTD_Cases = mendelian_ftd_cases,
  Total_Mendelian_FTD_Cases_common = mendelian_ftd_cases_common,
  Total_Mendelian_FTD_Cases_patho = mendelian_ftd_cases_patho,
  Total_Mendelian_FTD_Cases_ftd = mendelian_ftd_cases_ftd,
  Mendelian_FTD_Percentage = mendelian_ftd_percentage,
  Total_Polygenic_FTD_Cases = polygenic_ftd_cases,
  Polygenic_FTD_Percentage = polygenic_ftd_percentage,
  Total_Carriers = total_risk_alleles_ALS,
  Total_Penetrance_ALS = penetrance_percentage_ALS,
  Total_Penetrance_FTD = penetrance_percentage_FTD,
  Total_ALS_FTD_Cases = total_als_ftd_cases, 
  Total_Dementia_Cases = total_dementia_cases,
  Simulated_other_dementia_LifetimeRisk = simulated_lifetime_otherdementia_risk,
  Total_other_dementia_Cases = total_other_dementia_cases,
  LastGen_Mendelian_ALS_Cases = mendelian_ALS_individuals,
  LastGen_Polygenic_ALS_Cases = polygenic_ALS_individuals,
  LastGen_Mendelian_FTD_Cases = mendelian_FTD_individuals,
  LastGen_Polygenic_FTD_Cases = polygenic_FTD_individuals,
  LastGen_Dementia_Cases = dementia_individuals,
  LastGen_Carriers = carrier_individuals,
  LastGen_Sporadic = sporadic_appearing_cases,
  LastGen_Sporadic_Mendel = isolated_als_with_risk_alleles,
  LastGen_Sporadic_Mendel_perc = perc_seemingly_sporadic
  )

  # Print metrics
  print(paste("Statistics of all simulations:"))
  print(paste("Total simulations completed:", total_simulations))
  print(paste("Total number of individuals simulated:", total_individuals))
  print(paste("Total number of ALS cases:", total_als_cases))
  print(paste("Total number of Mendelian ALS individuals:", mendelian_als_cases))
  print(paste("Total number of common Mendelian ALS individuals:", mendelian_als_cases_common))
  print(paste("Total number of patho Mendelian ALS individuals:", mendelian_als_cases_patho))
  print(paste("Total number of polygenic ALS individuals:", polygenic_als_cases))
  print(paste("Simulated lifetime risk of ALS:", round(simulated_lifetime_als_risk * 100, 2), "%"))
  print(paste("Percentage of ALS cases due to Mendelian inheritance:", round(mendelian_als_percentage * 100, 2), "%"))
  print(paste("Percentage of ALS cases due to polygenic inheritance:", round(polygenic_als_percentage * 100, 2), "%"))
  print(paste("Disease penetrance of ALS among individuals with a risk allele:", round(penetrance_percentage_ALS * 100, 2), "%"))

  print(paste("Total number of FTD cases:", total_ftd_cases))
  print(paste("Total number of Mendelian FTD individuals:", mendelian_ftd_cases))
  print(paste("Total number of common Mendelian FTD individuals:", mendelian_ftd_cases_common))
  print(paste("Total number of patho Mendelian FTD individuals:", mendelian_ftd_cases_patho))
  print(paste("Total number of GRN/MAPT Mendelian FTD individuals:", mendelian_ftd_cases_ftd))
  print(paste("Total number of polygenic FTD individuals:", polygenic_ftd_cases))
  print(paste("Simulated lifetime risk of FTD:", round(simulated_lifetime_ftd_risk * 100, 2), "%"))
  print(paste("Percentage of FTD cases due to Mendelian inheritance:", round(mendelian_ftd_percentage * 100, 2), "%"))
  print(paste("Percentage of FTD cases due to polygenic inheritance:", round(polygenic_ftd_percentage * 100, 2), "%"))
  print(paste("Disease penetrance of FTD among individuals with a risk allele:", round(penetrance_percentage_FTD * 100, 2), "%"))
  print(paste("Total number of cases with both ALS and FTD:", total_als_ftd_cases))

  print(paste("Total number of dementia (FTD+others) cases:", total_dementia_cases))
  print(paste("Total number of other dementia cases:", total_other_dementia_cases))
  print(paste("Simulated lifetime risk of other dementia:", round(simulated_lifetime_otherdementia_risk * 100, 2), "%"))

  print(paste("Statistics of last generation/index patients:"))
  print(paste("Index patients with Mendelian ALS:", mendelian_ALS_individuals))
  print(paste("Index patients with polygenic ALS:", polygenic_ALS_individuals))
  print(paste("Index patients with Mendelian FTD:", mendelian_FTD_individuals))
  print(paste("Index patients with polygenic FTD:", polygenic_FTD_individuals))
  print(paste("Index patients with any dementia (FTD+others):", dementia_individuals))
  print(paste("Carrier individuals in the last generation:", carrier_individuals))

  print(paste("We calculate our models with one individual per pedigree, which involves:"))
  print(paste("Mendelian ALS:", sum(results_df$mendel_ALS_Y)))
  print(paste("Polygenic ALS:", sum(results_df$polygenicY_ALS)))
  print(paste("Mendelian FTD:", sum(results_df$mendel_FTD_Y)))
  print(paste("Polygenic FTD:", sum(results_df$polygenicY_FTD)))
  print(paste("People with both ALS and FTD:", sum(results_df$Y_ALS == 1 & results_df$Y_FTD == 1)))
  print(paste("People with any (FTD+others) dementia:", sum(results_df$Y_dementia == 1)))
  print(paste("ALS index patients seemingly sporadic:", sporadic_appearing_cases))
  print(paste("Isolated Mendelian cases, seemingly sporadic, with carriers in family:", isolated_als_with_risk_alleles))
  print(paste("Percentage isolated Mendelian cases, among seemingly sporadic cases:", perc_seemingly_sporadic))
  return(list(results_df = results_df, metrics = simulation_metrics))
}