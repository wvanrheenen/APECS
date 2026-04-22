#' ## To do's
#TODO: solve over-estimation of pedigree size,

#' make sure warnings are treated as errors
options( warn = 2 )

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
  core_offspring_created <- FALSE   
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
        # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
        } else { 
          n_II = rnbinom(1, size = 5, mu = lambda)
        }
        # Force offspring for at least one C*_1 individual; this will over-estimate pedigree size...
        while(n_II == 0 & (g < (k-1) | (grepl("^C", I1s$id[i]) & grepl("_1$", I1s$id[i])))) { # (g < (k-1); because the re-simulating one generation back was only a problem in for the gens that were not (k-1))
          n_II = 1 # minimal number of offspring needed to continue core pedigree; resampling has risk of generating >1 offspring when original there was none 
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
        # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"]) - 1 # minus 1 because one child has already been simulated in first round
        } else { 
          n_II = rnbinom(1, size = 5, mu = lambda) - 1 # minus 1 because one child has already been simulated in first round
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
        # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        if(is.na(lambda)){
          n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
        } else { 
          n_II = rnbinom(1, size = 5, mu = lambda)
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
  rlnorm(1, meanlog = log(100), sdlog = 0.6)
}

# Sample survival time after FTD onset (in years)
sample_survival_FTD <- function() {
  # Log-normal with median ≈ 6.6, 95% CI ≈ 5.6–7.6
  rlnorm(1, meanlog = log(100), sdlog = 0.25)
}

#' ### Function to add genetic values for Mendelian and polygenic inheritance
add_pheno = function(df_ped, disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, 
                      penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, 
                      rg_ALSFTD, rg_ALSdem, rg_FTDdem) {
  
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
  Vg = matrix(c(
    h2_ALS,
    rg_ALSFTD * sqrt(h2_ALS * h2_FTD),
    rg_ALSdem * sqrt(h2_ALS * h2_dementia),

    rg_ALSFTD * sqrt(h2_ALS * h2_FTD),
    h2_FTD,
    rg_FTDdem * sqrt(h2_FTD * h2_dementia),

    rg_ALSdem * sqrt(h2_ALS * h2_dementia),
    rg_FTDdem * sqrt(h2_FTD * h2_dementia),
    h2_dementia
  ), nrow = 3, byrow = TRUE)

  Gs = mvrnorm(length(founders), mu = c(0, 0, 0), Sigma = Vg)

  df_ped$G_ALS[founders]       = Gs[, 1]
  df_ped$G_FTD[founders]       = Gs[, 2]
  df_ped$G_dementia[founders]  = Gs[, 3]

  empirical_rg_ALSFTD = cor(df_ped$G_ALS[founders], df_ped$G_FTD[founders])
  empirical_rg_ALSdem = cor(df_ped$G_ALS[founders], df_ped$G_dementia[founders])
  empirical_rg_FTDdem = cor(df_ped$G_FTD[founders], df_ped$G_dementia[founders])
  

  # loop through nonfounders and estimate G based on G of parents and h2
  for(i in 1:nrow(df_ped)){
    if(is.na(df_ped$G_ALS[i])){
      pid = df_ped$pid[i]
      mid = df_ped$mid[i]
      Gp_ALS = df_ped[which(df_ped$id == pid),]$G_ALS
      Gm_ALS = df_ped[which(df_ped$id == mid),]$G_ALS
      Gp_FTD = df_ped[which(df_ped$id == pid),]$G_FTD
      Gm_FTD = df_ped[which(df_ped$id == mid),]$G_FTD
      Gp_dem = df_ped[which(df_ped$id == pid),]$G_dementia
      Gm_dem = df_ped[which(df_ped$id == mid),]$G_dementia      
      
      # Calculate mean and variance for offspring's G
      mean_ALS = (Gp_ALS + Gm_ALS) / 2
      var_ALS = 0.5 * h2_ALS
      mean_FTD = (Gp_FTD + Gm_FTD) / 2
      var_FTD = 0.5 * h2_FTD
      mean_dem = (Gp_dem + Gm_dem) / 2
      var_dem = 0.5 * h2_dementia

      # Simulate correlated G for offspring
      Vg_offspring = matrix(c(
        var_ALS,
        rg_ALSFTD * sqrt(var_ALS * var_FTD),
        rg_ALSdem  * sqrt(var_ALS * var_dem),

        rg_ALSFTD * sqrt(var_ALS * var_FTD),
        var_FTD,
        rg_FTDdem * sqrt(var_FTD * var_dem),

        rg_ALSdem  * sqrt(var_ALS * var_dem),
        rg_FTDdem * sqrt(var_FTD * var_dem),
        var_dem
      ), nrow = 3, byrow = TRUE)

      Gs_offspring = mvrnorm(1,
        mu    = c(mean_ALS, mean_FTD, mean_dem),
        Sigma = Vg_offspring)
      
      df_ped$G_ALS[i] = Gs_offspring[1]
      df_ped$G_FTD[i] = Gs_offspring[2]
      df_ped$G_dementia[i] = Gs_offspring[3]
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

#' function to identify relatedness of monogenic 

add_1st_relatives = function(df_ped, k) {
  # Initialize columns for relative counts
  df_ped$relatives_1st = 0
  df_ped$relatives_1st_als = 0
  df_ped$relatives_1st_als_monogenic = 0
  df_ped$relatives_1st_als_polygenic = 0
  df_ped$relatives_1st_ftd = 0
  df_ped$relatives_1st_ftd_unique = 0
  df_ped$relatives_1st_dementia = 0
  df_ped$relatives_1st_dementia_unique = 0
  
  # Filter to get ONLY index patients: last generation, ALS affected, founder pedigree
  index_patients <- df_ped %>%
    filter(gen == k, Y_ALS == 1, str_starts(id, "C"))
  
  if (nrow(index_patients) == 0) {
    message("No ALS index patients in founder pedigree (gen == k, Y_ALS == 1, id starts with 'C'). Skipping...")
    return(df_ped)
  }
  index_rows <- match(index_patients$id, df_ped$id)

  # Loop through each individual
  for(i in index_rows) {
    id = df_ped$id[i]
    relatives_1st = 0
    relatives_1st_als = 0
    relatives_1st_als_monogenic = 0
    relatives_1st_als_polygenic = 0
    relatives_1st_ftd = 0
    relatives_1st_ftd_unique = 0
    relatives_1st_dementia = 0
    relatives_1st_dementia_unique = 0
    
    proband_onset = df_ped$year_onset_ALS[i]

    # Check for siblings
    siblings = df_ped %>% filter(pid == df_ped$pid[i] & mid == df_ped$mid[i] & id != df_ped$id[i])
    relatives_1st = relatives_1st + nrow(siblings)
    relatives_1st_als = relatives_1st_als + sum(siblings$Y_ALS == 1 & siblings$year_onset_ALS <= proband_onset & !is.na(siblings$year_onset_ALS))
    relatives_1st_als_monogenic = relatives_1st_als_monogenic + sum(siblings$mendel_ALS_Y == 1 & siblings$year_onset_ALS <= proband_onset & !is.na(siblings$year_onset_ALS))
    relatives_1st_als_polygenic = relatives_1st_als_polygenic + sum(siblings$polygenicY_ALS == 1 & siblings$year_onset_ALS <= proband_onset & !is.na(siblings$year_onset_ALS))
    relatives_1st_ftd = relatives_1st_ftd + sum(siblings$Y_FTD == 1 & siblings$year_onset_FTD <= proband_onset & !is.na(siblings$year_onset_FTD))
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(siblings$Y_FTD == 1 & siblings$Y_ALS == 0 & siblings$year_onset_FTD <= proband_onset & !is.na(siblings$year_onset_FTD))
    relatives_1st_dementia = relatives_1st_dementia + sum(siblings$Y_dementia == 1 & siblings$year_onset_dementia <= proband_onset & !is.na(siblings$year_onset_dementia))
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(siblings$Y_dementia == 1 & siblings$Y_ALS == 0 & siblings$year_onset_dementia <= proband_onset & !is.na(siblings$year_onset_dementia))

    # Check for parents
    parents = df_ped %>% filter(id %in% c(df_ped$pid[i], df_ped$mid[i]))
    relatives_1st = relatives_1st + nrow(parents)
    relatives_1st_als = relatives_1st_als + sum(parents$Y_ALS == 1 & parents$year_onset_ALS <= proband_onset & !is.na(parents$year_onset_ALS))
    relatives_1st_als_monogenic = relatives_1st_als_monogenic + sum(parents$mendel_ALS_Y == 1 & parents$year_onset_ALS <= proband_onset & !is.na(parents$year_onset_ALS))
    relatives_1st_als_polygenic = relatives_1st_als_polygenic + sum(parents$polygenicY_ALS == 1 & parents$year_onset_ALS <= proband_onset & !is.na(parents$year_onset_ALS))
    relatives_1st_ftd = relatives_1st_ftd + sum(parents$Y_FTD == 1 & parents$year_onset_FTD <= proband_onset & !is.na(parents$year_onset_FTD))
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(parents$Y_FTD == 1 & parents$Y_ALS == 0 & parents$year_onset_FTD <= proband_onset & !is.na(parents$year_onset_FTD))
    relatives_1st_dementia = relatives_1st_dementia + sum(parents$Y_dementia == 1 & parents$year_onset_dementia <= proband_onset & !is.na(parents$year_onset_dementia))
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(parents$Y_dementia == 1 & parents$Y_ALS == 0 & parents$year_onset_dementia <= proband_onset & !is.na(parents$year_onset_dementia))
    
    # Check for children
    children = df_ped %>% filter(pid == df_ped$id[i] | mid == df_ped$id[i])
    relatives_1st = relatives_1st + nrow(children)
    relatives_1st_als = relatives_1st_als + sum(children$Y_ALS == 1 & children$year_onset_ALS <= proband_onset & !is.na(children$year_onset_ALS))
    relatives_1st_als_monogenic = relatives_1st_als_monogenic + sum(children$mendel_ALS_Y == 1 & children$year_onset_ALS <= proband_onset & !is.na(children$year_onset_ALS))
    relatives_1st_als_polygenic = relatives_1st_als_polygenic + sum(children$polygenicY_ALS == 1 & children$year_onset_ALS <= proband_onset & !is.na(children$year_onset_ALS))
    relatives_1st_ftd = relatives_1st_ftd + sum(children$Y_FTD == 1 & children$year_onset_FTD <= proband_onset & !is.na(children$year_onset_FTD))
    relatives_1st_ftd_unique = relatives_1st_ftd_unique + sum(children$Y_FTD == 1 & children$Y_ALS == 0 & children$year_onset_FTD <= proband_onset & !is.na(children$year_onset_FTD))
    relatives_1st_dementia = relatives_1st_dementia + sum(children$Y_dementia == 1 & children$year_onset_dementia <= proband_onset & !is.na(children$year_onset_dementia))
    relatives_1st_dementia_unique = relatives_1st_dementia_unique + sum(children$Y_dementia == 1 & children$Y_ALS == 0 & children$year_onset_dementia <= proband_onset & !is.na(children$year_onset_dementia))

    # Update the data frame
    df_ped$relatives_1st[i] = relatives_1st
    df_ped$relatives_1st_als[i] = relatives_1st_als
    df_ped$relatives_1st_als_monogenic[i] = relatives_1st_als_monogenic
    df_ped$relatives_1st_als_polygenic[i] = relatives_1st_als_polygenic  
    df_ped$relatives_1st_ftd[i] = relatives_1st_ftd
    df_ped$relatives_1st_ftd_unique[i] = relatives_1st_ftd_unique
    df_ped$relatives_1st_dementia[i] = relatives_1st_dementia
    df_ped$relatives_1st_dementia_unique[i] = relatives_1st_dementia_unique
  }
  return(df_ped)
}

add_2nd_and_3rd_relatives = function(df_ped, k) {
  # Initialize columns for relative counts
  df_ped$relatives_2nd = 0
  df_ped$relatives_2nd_als = 0
  df_ped$relatives_2nd_als_monogenic = 0
  df_ped$relatives_2nd_als_polygenic = 0
  df_ped$relatives_2nd_ftd = 0
  df_ped$relatives_2nd_ftd_unique = 0
  df_ped$relatives_2nd_dementia = 0
  df_ped$relatives_2nd_dementia_unique = 0
  df_ped$relatives_3rd = 0
  df_ped$relatives_3rd_als = 0
  df_ped$relatives_3rd_als_monogenic = 0
  df_ped$relatives_3rd_als_polygenic = 0
  df_ped$relatives_3rd_ftd = 0
  df_ped$relatives_3rd_ftd_unique = 0
  df_ped$relatives_3rd_dementia = 0
  df_ped$relatives_3rd_dementia_unique = 0
 
  # Filter to get ONLY index patients: last generation, ALS affected, founder pedigree
  index_patients <- df_ped %>%
    filter(gen == k, Y_ALS == 1, str_starts(id, "C"))
  
  if (nrow(index_patients) == 0) {
    message("No ALS index patients in founder pedigree (gen == k, Y_ALS == 1, id starts with 'C'). Skipping...")
    return(df_ped)
  }
  index_rows <- match(index_patients$id, df_ped$id)

  # Loop through each individual
  for(i in index_rows) {
    id = df_ped$id[i]
    relatives_2nd = 0
    relatives_2nd_als = 0
    relatives_2nd_als_monogenic = 0
    relatives_2nd_als_polygenic = 0
    relatives_2nd_ftd = 0
    relatives_2nd_ftd_unique = 0
    relatives_2nd_dementia = 0
    relatives_2nd_dementia_unique = 0
    relatives_3rd = 0
    relatives_3rd_als = 0
    relatives_3rd_als_monogenic = 0
    relatives_3rd_als_polygenic = 0
    relatives_3rd_ftd = 0    
    relatives_3rd_ftd_unique = 0
    relatives_3rd_dementia = 0    
    relatives_3rd_dementia_unique = 0
    
    proband_onset = df_ped$year_onset_ALS[i] 
    
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
      relatives_2nd_als = relatives_2nd_als + sum(grandparent_rows$Y_ALS == 1 & grandparent_rows$year_onset_ALS <= proband_onset & !is.na(grandparent_rows$year_onset_ALS))
      relatives_2nd_als_monogenic = relatives_2nd_als_monogenic + sum(grandparent_rows$mendel_ALS_Y == 1 & grandparent_rows$year_onset_ALS <= proband_onset & !is.na(grandparent_rows$year_onset_ALS))
      relatives_2nd_als_polygenic = relatives_2nd_als_polygenic + sum(grandparent_rows$polygenicY_ALS == 1 & grandparent_rows$year_onset_ALS <= proband_onset & !is.na(grandparent_rows$year_onset_ALS))
      relatives_2nd_ftd = relatives_2nd_ftd + sum(grandparent_rows$Y_FTD == 1 & grandparent_rows$year_onset_FTD <= proband_onset & !is.na(grandparent_rows$year_onset_FTD))
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(grandparent_rows$Y_FTD == 1 & grandparent_rows$Y_ALS == 0 & grandparent_rows$year_onset_FTD <= proband_onset & !is.na(grandparent_rows$year_onset_FTD))
      relatives_2nd_dementia = relatives_2nd_dementia + sum(grandparent_rows$Y_dementia == 1 & grandparent_rows$year_onset_dementia <= proband_onset & !is.na(grandparent_rows$year_onset_dementia))
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(grandparent_rows$Y_dementia == 1 & grandparent_rows$Y_ALS == 0 & grandparent_rows$year_onset_dementia <= proband_onset & !is.na(grandparent_rows$year_onset_dementia))
    
      # After grandparent identification, before great-grandparents
      grand_uncle_aunt_ids = c()

      for (grandparent_id in grandparent_ids) {
        # Find grandparents' row
        grandparent_row = df_ped %>%
          filter(id == grandparent_id)
        
        if (nrow(grandparent_row) > 0) {
          grand_uncle_aunt_pid = grandparent_row$pid
          grand_uncle_aunt_mid = grandparent_row$mid
          
          sibling_rows = df_ped %>%
            filter(pid == grand_uncle_aunt_pid & mid == grand_uncle_aunt_mid & id != grandparent_id)
          
          # Collect their IDs
          grand_uncle_aunt_ids = c(grand_uncle_aunt_ids, sibling_rows$id)
        }
      }

      grand_uncle_aunt_rows = df_ped %>% filter(id %in% grand_uncle_aunt_ids)

      if(nrow(grand_uncle_aunt_rows) > 0) {
        relatives_3rd = relatives_3rd + nrow(grand_uncle_aunt_rows)
        relatives_3rd_als = relatives_3rd_als + sum(grand_uncle_aunt_rows$Y_ALS == 1 & grand_uncle_aunt_rows$year_onset_ALS <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_ALS))
        relatives_3rd_als_monogenic = relatives_3rd_als_monogenic + sum(grand_uncle_aunt_rows$mendel_ALS_Y == 1 & grand_uncle_aunt_rows$year_onset_ALS <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_ALS))
        relatives_3rd_als_polygenic = relatives_3rd_als_polygenic + sum(grand_uncle_aunt_rows$polygenicY_ALS == 1 & grand_uncle_aunt_rows$year_onset_ALS <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_ALS))
        relatives_3rd_ftd = relatives_3rd_ftd + sum(grand_uncle_aunt_rows$Y_FTD == 1 & grand_uncle_aunt_rows$year_onset_FTD <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_FTD))
        relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(grand_uncle_aunt_rows$Y_FTD == 1 & grand_uncle_aunt_rows$Y_ALS == 0 & grand_uncle_aunt_rows$year_onset_FTD <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_FTD))
        relatives_3rd_dementia = relatives_3rd_dementia + sum(grand_uncle_aunt_rows$Y_dementia == 1 & grand_uncle_aunt_rows$year_onset_dementia <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_dementia))
        relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(grand_uncle_aunt_rows$Y_dementia == 1 & grand_uncle_aunt_rows$Y_ALS == 0 & grand_uncle_aunt_rows$year_onset_dementia <= proband_onset & !is.na(grand_uncle_aunt_rows$year_onset_dementia))
      }

 
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
        relatives_3rd_als = relatives_3rd_als + sum(great_grandparent_rows$Y_ALS == 1 & great_grandparent_rows$year_onset_ALS <= proband_onset & !is.na(great_grandparent_rows$year_onset_ALS))
        relatives_3rd_als_monogenic = relatives_3rd_als_monogenic + sum(great_grandparent_rows$mendel_ALS_Y == 1 & great_grandparent_rows$year_onset_ALS <= proband_onset & !is.na(great_grandparent_rows$year_onset_ALS))
        relatives_3rd_als_polygenic = relatives_3rd_als_polygenic + sum(great_grandparent_rows$polygenicY_ALS == 1 & great_grandparent_rows$year_onset_ALS <= proband_onset & !is.na(great_grandparent_rows$year_onset_ALS))
        relatives_3rd_ftd = relatives_3rd_ftd + sum(great_grandparent_rows$Y_FTD == 1 & great_grandparent_rows$year_onset_FTD <= proband_onset & !is.na(great_grandparent_rows$year_onset_FTD))
        relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(great_grandparent_rows$Y_FTD == 1 & great_grandparent_rows$Y_ALS == 0 & great_grandparent_rows$year_onset_FTD <= proband_onset & !is.na(great_grandparent_rows$year_onset_FTD))
        relatives_3rd_dementia = relatives_3rd_dementia + sum(great_grandparent_rows$Y_dementia == 1 & great_grandparent_rows$year_onset_dementia <= proband_onset & !is.na(great_grandparent_rows$year_onset_dementia))
        relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(great_grandparent_rows$Y_dementia == 1 & great_grandparent_rows$Y_ALS == 0 & great_grandparent_rows$year_onset_dementia <= proband_onset & !is.na(great_grandparent_rows$year_onset_dementia))
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
      relatives_2nd_als = relatives_2nd_als + sum(grandchild_rows$Y_ALS == 1 & grandchild_rows$year_onset_ALS <= proband_onset & !is.na(grandchild_rows$year_onset_ALS))
      relatives_2nd_als_monogenic = relatives_2nd_als_monogenic + sum(grandchild_rows$mendel_ALS_Y == 1 & grandchild_rows$year_onset_ALS <= proband_onset & !is.na(grandchild_rows$year_onset_ALS))
      relatives_2nd_als_polygenic = relatives_2nd_als_polygenic + sum(grandchild_rows$polygenicY_ALS == 1 & grandchild_rows$year_onset_ALS <= proband_onset & !is.na(grandchild_rows$year_onset_ALS))
      relatives_2nd_ftd = relatives_2nd_ftd + sum(grandchild_rows$Y_FTD == 1 & grandchild_rows$year_onset_FTD <= proband_onset & !is.na(grandchild_rows$year_onset_FTD))
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(grandchild_rows$Y_FTD == 1 & grandchild_rows$Y_ALS == 0 & grandchild_rows$year_onset_FTD <= proband_onset & !is.na(grandchild_rows$year_onset_FTD))
      relatives_2nd_dementia = relatives_2nd_dementia + sum(grandchild_rows$Y_dementia == 1 & grandchild_rows$year_onset_dementia <= proband_onset & !is.na(grandchild_rows$year_onset_dementia))
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(grandchild_rows$Y_dementia == 1 & grandchild_rows$Y_ALS == 0 & grandchild_rows$year_onset_dementia <= proband_onset & !is.na(grandchild_rows$year_onset_dementia))
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
      relatives_3rd_als = relatives_3rd_als + sum(great_grandchild_rows$Y_ALS == 1 & great_grandchild_rows$year_onset_ALS <= proband_onset & !is.na(great_grandchild_rows$year_onset_ALS))
      relatives_3rd_als_monogenic = relatives_3rd_als_monogenic + sum(great_grandchild_rows$mendel_ALS_Y == 1 & great_grandchild_rows$year_onset_ALS <= proband_onset & !is.na(great_grandchild_rows$year_onset_ALS))
      relatives_3rd_als_polygenic = relatives_3rd_als_polygenic + sum(great_grandchild_rows$polygenicY_ALS == 1 & great_grandchild_rows$year_onset_ALS <= proband_onset & !is.na(great_grandchild_rows$year_onset_ALS))
      relatives_3rd_ftd = relatives_3rd_ftd + sum(great_grandchild_rows$Y_FTD == 1 & great_grandchild_rows$year_onset_FTD <= proband_onset & !is.na(great_grandchild_rows$year_onset_FTD))
      relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(great_grandchild_rows$Y_FTD == 1 & great_grandchild_rows$Y_ALS == 0 & great_grandchild_rows$year_onset_FTD <= proband_onset & !is.na(great_grandchild_rows$year_onset_FTD))
      relatives_3rd_dementia = relatives_3rd_dementia + sum(great_grandchild_rows$Y_dementia == 1 & great_grandchild_rows$year_onset_dementia <= proband_onset & !is.na(great_grandchild_rows$year_onset_dementia))
      relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(great_grandchild_rows$Y_dementia == 1 & great_grandchild_rows$Y_ALS == 0 & great_grandchild_rows$year_onset_dementia <= proband_onset & !is.na(great_grandchild_rows$year_onset_dementia))
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
      relatives_2nd_als = relatives_2nd_als + sum(aunt_uncle_rows$Y_ALS == 1 & aunt_uncle_rows$year_onset_ALS <= proband_onset & !is.na(aunt_uncle_rows$year_onset_ALS))
      relatives_2nd_als_monogenic = relatives_2nd_als_monogenic + sum(aunt_uncle_rows$mendel_ALS_Y == 1 & aunt_uncle_rows$year_onset_ALS <= proband_onset & !is.na(aunt_uncle_rows$year_onset_ALS))
      relatives_2nd_als_polygenic = relatives_2nd_als_polygenic + sum(aunt_uncle_rows$polygenicY_ALS == 1 & aunt_uncle_rows$year_onset_ALS <= proband_onset & !is.na(aunt_uncle_rows$year_onset_ALS))
      relatives_2nd_ftd = relatives_2nd_ftd + sum(aunt_uncle_rows$Y_FTD == 1 & aunt_uncle_rows$year_onset_FTD <= proband_onset & !is.na(aunt_uncle_rows$year_onset_FTD))
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(aunt_uncle_rows$Y_FTD == 1 & aunt_uncle_rows$Y_ALS == 0 & aunt_uncle_rows$year_onset_FTD <= proband_onset & !is.na(aunt_uncle_rows$year_onset_FTD))
      relatives_2nd_dementia = relatives_2nd_dementia + sum(aunt_uncle_rows$Y_dementia == 1 & aunt_uncle_rows$year_onset_dementia <= proband_onset & !is.na(aunt_uncle_rows$year_onset_dementia))
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(aunt_uncle_rows$Y_dementia == 1 & aunt_uncle_rows$Y_ALS == 0 & aunt_uncle_rows$year_onset_dementia <= proband_onset & !is.na(aunt_uncle_rows$year_onset_dementia))
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
      relatives_2nd_als = relatives_2nd_als + sum(nephew_niece_rows$Y_ALS == 1 & nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(nephew_niece_rows$year_onset_ALS))
      relatives_2nd_als_monogenic = relatives_2nd_als_monogenic + sum(nephew_niece_rows$mendel_ALS_Y == 1 & nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(nephew_niece_rows$year_onset_ALS))
      relatives_2nd_als_polygenic = relatives_2nd_als_polygenic + sum(nephew_niece_rows$polygenicY_ALS == 1 & nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(nephew_niece_rows$year_onset_ALS))
      relatives_2nd_ftd = relatives_2nd_ftd + sum(nephew_niece_rows$Y_FTD == 1 & nephew_niece_rows$year_onset_FTD <= proband_onset & !is.na(nephew_niece_rows$year_onset_FTD))
      relatives_2nd_ftd_unique = relatives_2nd_ftd_unique + sum(nephew_niece_rows$Y_FTD == 1 & nephew_niece_rows$Y_ALS == 0 & nephew_niece_rows$year_onset_FTD <= proband_onset & !is.na(nephew_niece_rows$year_onset_FTD))
      relatives_2nd_dementia = relatives_2nd_dementia + sum(nephew_niece_rows$Y_dementia == 1 & nephew_niece_rows$year_onset_dementia <= proband_onset & !is.na(nephew_niece_rows$year_onset_dementia))
      relatives_2nd_dementia_unique = relatives_2nd_dementia_unique + sum(nephew_niece_rows$Y_dementia == 1 & nephew_niece_rows$Y_ALS == 0 & nephew_niece_rows$year_onset_dementia <= proband_onset & !is.na(nephew_niece_rows$year_onset_dementia))
    
      # Identify grand-nephews and grand-nieces (third-degree relatives)
      grand_nephew_niece_ids = c()
      for(nephew_niece_id in nephew_niece_ids) {
        grand_nephew_niece_rows = df_ped %>% filter(pid == nephew_niece_id | mid == nephew_niece_id)
        grand_nephew_niece_ids = c(grand_nephew_niece_ids, grand_nephew_niece_rows$id)
      }
      grand_nephew_niece_rows = df_ped %>% filter(id %in% grand_nephew_niece_ids)
      if(nrow(grand_nephew_niece_rows) > 0) {
        relatives_3rd = relatives_3rd + nrow(grand_nephew_niece_rows)
        relatives_3rd_als = relatives_3rd_als + sum(grand_nephew_niece_rows$Y_ALS == 1 & grand_nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_ALS))
        relatives_3rd_als_monogenic = relatives_3rd_als_monogenic + sum(grand_nephew_niece_rows$mendel_ALS_Y == 1 & grand_nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_ALS))
        relatives_3rd_als_polygenic = relatives_3rd_als_polygenic + sum(grand_nephew_niece_rows$polygenicY_ALS == 1 & grand_nephew_niece_rows$year_onset_ALS <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_ALS))
        relatives_3rd_ftd = relatives_3rd_ftd + sum(grand_nephew_niece_rows$Y_FTD == 1 & grand_nephew_niece_rows$year_onset_FTD <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_FTD))
        relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(grand_nephew_niece_rows$Y_FTD == 1 & grand_nephew_niece_rows$Y_ALS == 0 & grand_nephew_niece_rows$year_onset_FTD <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_FTD))
        relatives_3rd_dementia = relatives_3rd_dementia + sum(grand_nephew_niece_rows$Y_dementia == 1 & grand_nephew_niece_rows$year_onset_dementia <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_dementia))
        relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(grand_nephew_niece_rows$Y_dementia == 1 & grand_nephew_niece_rows$Y_ALS == 0 & grand_nephew_niece_rows$year_onset_dementia <= proband_onset & !is.na(grand_nephew_niece_rows$year_onset_dementia))
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
      relatives_3rd_als = relatives_3rd_als + sum(first_cousin_rows$Y_ALS == 1 & first_cousin_rows$year_onset_ALS <= proband_onset & !is.na(first_cousin_rows$year_onset_ALS))
      relatives_3rd_als_monogenic = relatives_3rd_als_monogenic + sum(first_cousin_rows$mendel_ALS_Y == 1 & first_cousin_rows$year_onset_ALS <= proband_onset & !is.na(first_cousin_rows$year_onset_ALS))
      relatives_3rd_als_polygenic = relatives_3rd_als_polygenic + sum(first_cousin_rows$polygenicY_ALS == 1 & first_cousin_rows$year_onset_ALS <= proband_onset & !is.na(first_cousin_rows$year_onset_ALS))
      relatives_3rd_ftd = relatives_3rd_ftd + sum(first_cousin_rows$Y_FTD == 1 & first_cousin_rows$year_onset_FTD <= proband_onset & !is.na(first_cousin_rows$year_onset_FTD))
      relatives_3rd_ftd_unique = relatives_3rd_ftd_unique + sum(first_cousin_rows$Y_FTD == 1 & first_cousin_rows$Y_ALS == 0 & first_cousin_rows$year_onset_FTD <= proband_onset & !is.na(first_cousin_rows$year_onset_FTD))
      relatives_3rd_dementia = relatives_3rd_dementia + sum(first_cousin_rows$Y_dementia == 1 & first_cousin_rows$year_onset_dementia <= proband_onset & !is.na(first_cousin_rows$year_onset_dementia))
      relatives_3rd_dementia_unique = relatives_3rd_dementia_unique + sum(first_cousin_rows$Y_dementia == 1 & first_cousin_rows$Y_ALS == 0 & first_cousin_rows$year_onset_dementia <= proband_onset & !is.na(first_cousin_rows$year_onset_dementia))
    }


    # Update the data frame
    df_ped$relatives_2nd[i] = relatives_2nd
    df_ped$relatives_2nd_als[i] = relatives_2nd_als
    df_ped$relatives_2nd_als_monogenic[i] = relatives_2nd_als_monogenic
    df_ped$relatives_2nd_als_polygenic[i] = relatives_2nd_als_polygenic
    df_ped$relatives_2nd_ftd[i] = relatives_2nd_ftd
    df_ped$relatives_2nd_ftd_unique[i] = relatives_2nd_ftd_unique
    df_ped$relatives_2nd_dementia[i] = relatives_2nd_dementia
    df_ped$relatives_2nd_dementia_unique[i] = relatives_2nd_dementia_unique
    df_ped$relatives_3rd[i] = relatives_3rd
    df_ped$relatives_3rd_als[i] = relatives_3rd_als
    df_ped$relatives_3rd_als_monogenic[i] = relatives_3rd_als_monogenic
    df_ped$relatives_3rd_als_polygenic[i] = relatives_3rd_als_polygenic
    df_ped$relatives_3rd_ftd[i] = relatives_3rd_ftd
    df_ped$relatives_3rd_ftd_unique[i] = relatives_3rd_ftd_unique
    df_ped$relatives_3rd_dementia[i] = relatives_3rd_dementia
    df_ped$relatives_3rd_dementia_unique[i] = relatives_3rd_dementia_unique
    }
    return(df_ped)
  }


#' wrapper function to simulate full pedigree
sim_ped = function(i, k, lambda=NA, yob_index=1960, mean_gen_yr, fert_rate, disease_onset,
                   DAF_common, DAF_patho, DAF_ftd, penetrance_ALS_common, penetrance_ALS_patho, 
                   K_ALS, h2_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, K_FTD, h2_FTD, 
                   penetrance_dem_common, penetrance_dem_patho, K_dementia, h2_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem, 
                   life_expectancy, current_year, plot=FALSE){
  core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)   
  core_ped = add_pheno(core_ped, disease_onset=disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem) 
  core_ped = add_1st_relatives(core_ped, k)
  core_ped = add_2nd_and_3rd_relatives(core_ped, k)

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
pedigree_simulations = function(
  n_simulations, k, yob_index_min, yob_index_max, lambda = NA, mean_gen_yr, fert_rate, disease_onset,
  DAF_common, DAF_patho, DAF_ftd, penetrance_ALS_common, penetrance_ALS_patho, 
  K_ALS, h2_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, 
  K_FTD, h2_FTD, penetrance_dem_common, penetrance_dem_patho, K_dementia, h2_dementia, 
  life_expectancy, current_year, rg_ALSFTD, rg_ALSdem, rg_FTDdem, plot = FALSE, yob_min = 1975, yob_max = 2025, interval = 10
) {
  library(dplyr)
  
  # Define yob intervals
  breaks = seq(yob_min, yob_max, by = interval)
  labels = paste(head(breaks, -1), tail(breaks, -1), sep = "-")
  
  # Initialize counters
  interval_metrics = data.frame(
    yob_interval = labels,
    n_individuals = 0,
    n_ALS = 0,
    n_mendelian_ALS = 0,
    n_mendelian_c9_ALS = 0,
    n_mendelian_SOD1FUS_ALS = 0,
    n_polygenic_ALS = 0,
    n_FTD = 0,
    n_mendelian_FTD = 0,
    n_c9_FTD = 0,
    n_SOD1FUS_FTD = 0,
    n_ftd_FTD = 0,
    n_polygenic_FTD = 0,
    n_dementia = 0,
    n_c9_dementia = 0,
    n_other_dementia = 0,
    n_ALS_FTD_cooccurrence = 0,
    n_c9_alleles = 0,
    n_SOD1FUS_alleles = 0,
    stringsAsFactors = FALSE
  )
  results_df = data.frame() # For index patients

  for (sim in seq_len(n_simulations)) {
    yob_index = sample(yob_index_min:yob_index_max, 1)
    ped = tryCatch(
      sim_ped(sim, k = k, lambda = lambda, yob_index = yob_index, mean_gen_yr = mean_gen_yr,
              fert_rate = fert_rate, disease_onset = disease_onset,
              DAF_common = DAF_common, DAF_patho = DAF_patho, DAF_ftd = DAF_ftd,
              penetrance_ALS_common = penetrance_ALS_common, penetrance_ALS_patho = penetrance_ALS_patho,
              K_ALS = K_ALS, h2_ALS = h2_ALS,
              penetrance_FTD_common = penetrance_FTD_common, penetrance_FTD_patho = penetrance_FTD_patho,
              penetrance_FTD_nonALS = penetrance_FTD_nonALS, K_FTD = K_FTD, h2_FTD = h2_FTD,
              penetrance_dem_common = penetrance_dem_common, penetrance_dem_patho = penetrance_dem_patho,
              K_dementia = K_dementia, h2_dementia = h2_dementia,
              life_expectancy = life_expectancy, current_year = current_year, rg_ALSFTD = rg_ALSFTD, rg_ALSdem = rg_ALSdem, rg_FTDdem = rg_FTDdem, plot = plot),
      error = function(e) { message(sprintf("Simulation %d failed: %s", sim, e$message)); NULL }
    )
    if (is.null(ped)) next
    
    # Assign yob_interval for all individuals
    ped$yob_interval = cut(ped$yob, breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
    
    # 1. All individuals: count by interval
    ped_all = ped[!is.na(ped$yob_interval), ]
    if (nrow(ped_all) > 0) {
      interval_counts = ped_all %>%
        group_by(yob_interval) %>%
        summarize(
          n_individuals = n(),
          n_ALS = sum(Y_ALS == 1, na.rm = TRUE),
          n_mendelian_ALS = sum(mendel_ALS_Y == 1, na.rm = TRUE),
          n_mendelian_c9_ALS = sum(mendel_ALS_Y_common == 1, na.rm = TRUE),
          n_mendelian_SOD1FUS_ALS = sum(mendel_ALS_Y_patho == 1, na.rm = TRUE),
          n_polygenic_ALS = sum(polygenicY_ALS == 1, na.rm = TRUE),
          n_FTD = sum(Y_FTD == 1, na.rm = TRUE),
          n_mendelian_FTD = sum(mendel_FTD_Y == 1, na.rm = TRUE),
          n_c9_FTD = sum(mendel_FTD_Y_common == 1, na.rm = TRUE),
          n_SOD1FUS_FTD = sum(mendel_FTD_Y_patho == 1, na.rm = TRUE),
          n_ftd_FTD = sum(mendel_FTD_Y_ftd == 1, na.rm = TRUE),
          n_polygenic_FTD = sum(polygenicY_FTD == 1, na.rm = TRUE),
          n_dementia = sum(Y_dementia == 1, na.rm = TRUE),
          n_c9_dementia = sum(mendel_dem_Y_common == 1, na.rm = TRUE),
          n_other_dementia = sum(Y_dementia_other == 1, na.rm = TRUE),
          n_ALS_FTD_cooccurrence = sum(Y_ALS == 1 & Y_FTD == 1, na.rm = TRUE),
          n_c9_alleles = sum(a1_common == 1 | a2_common == 1, na.rm = TRUE),
          n_SOD1FUS_alleles = sum(a1_patho == 1 | a2_patho == 1, na.rm = TRUE)
        ) %>% ungroup()
      for (i in seq_len(nrow(interval_counts))) {
        idx = which(interval_metrics$yob_interval == interval_counts$yob_interval[i])
        for (col in names(interval_counts)[-1]) {
          interval_metrics[idx, col] = interval_metrics[idx, col] + interval_counts[[col]][i]
        }
      }
    }
    
    # 2. Index patient row for results_df + phenocopy tracking  
    case_rows = ped %>% 
      filter(gen == k, Y_ALS == 1, startsWith(id, "C"), year_onset_ALS > (current_year - 15))

    if (nrow(case_rows) == 0) {
      message(sprintf("Simulation %d has no recent-onset ALS index patients in the founder pedigree. Skipping...", sim))
      next
    }

    # Process ALL matching index patients (no head(1))
    for (j in 1:nrow(case_rows)) {
      case_row = case_rows[j, ]
      
      # Phenocopy = monogenic INDEX patient + any RELATED polygenic ALS in pedigree
      has_polygenic_ALS_1st = case_row$relatives_1st_als_polygenic > 0
      has_polygenic_ALS_2nd = case_row$relatives_2nd_als_polygenic > 0
      has_polygenic_ALS_3rd = case_row$relatives_3rd_als_polygenic > 0

      has_monogenic_index = case_row$mendel_ALS_Y == 1
      has_phenocopy_1st = has_monogenic_index & has_polygenic_ALS_1st
      has_phenocopy_2nd = has_monogenic_index & has_polygenic_ALS_2nd
      has_phenocopy_3rd = has_monogenic_index & has_polygenic_ALS_3rd

      message(sprintf("Simulation %d has %d recent-onset ALS index patient(s) in founder pedigree. Adding to dataframe...", sim, nrow(case_rows)))

      # Add simulation identifier and phenocopy flag to index patient row
      case_row$sim_id = sim
      case_row$phenocopy = ifelse(
        has_phenocopy_1st | has_phenocopy_2nd | has_phenocopy_3rd, 1, 0
      )
      case_row$phenocopy_1st = ifelse(has_phenocopy_1st, 1, 0)
      case_row$phenocopy_2nd = ifelse(has_phenocopy_2nd, 1, 0)
      case_row$phenocopy_3rd = ifelse(has_phenocopy_3rd, 1, 0)
      case_row$id = paste0("sim", sim, "_", case_row$id)
      results_df = bind_rows(results_df, case_row)
    }
  }
  
  return(list(
    results_df = results_df,
    interval_metrics_combined = interval_metrics
  ))
}

