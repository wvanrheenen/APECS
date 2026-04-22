
#' make sure warnings are treated as errors
# options( warn = 2 )

source("../../src/libraries_simPed.R")

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


#' Sample offspring conditional on parent's sibship size (rho correlation)
sample_offspring_correlated <- function(parent_sibsize, rho = 0.15, fert_mu, theta = 5) {
  cond_mu <- fert_mu + rho * (parent_sibsize - fert_mu)
  cond_mu <- max(0.1, cond_mu)  # Floor to avoid negatives/zeros prematurely
  
  rnbinom(1, size = theta, mu = cond_mu)
}

#' ### Function to initiate pedigree with one founder with a mutation
init_ped = function(k, mean_gen_yr, yob_index, life_expectancy, current_year){
  core_ped = data.frame(gen = 0, # generation
                        id = "C0", # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1) # sex, 0 = male, 1 = female
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
add_gen = function(df_ped, lambda, k, fert_rate, mean_gen_yr, life_expectancy, current_year, rho = 0.15){
  g = 0
  if(!"fert_rate" %in% names(df_ped)) df_ped$fert_rate <- NA_real_
  if(!"n_offspring" %in% names(df_ped)) df_ped$n_offspring <- NA_real_
  if(!"parent_sibsize" %in% names(df_ped)) df_ped$parent_sibsize <- NA_real_
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
        # CORRELATED OFFSPRING SAMPLING
        birth_yr <- I1s$yob[i] + gen_yr
        fert_mu <- if(is.na(lambda)){
          fert_rate[fert_rate$year == birth_yr, "mean_fertility"]
        } else { lambda }
        
        # Independent for founders (gen 0), correlated for later generations
        if(g == 0){
          # Founders: pure NB
          n_II <- rnbinom(1, size = 5, mu = fert_mu)
          while(n_II == 0 & g == 0){
            n_II <- rnbinom(1, size = 5, mu = fert_mu)
          }
        } else {
          # Later gens: condition on THEIR parent's sibship size (i.e. grandmother's n_offspring)
          parent_sibsize <- I1s$parent_sibsize[i]
          if(is.na(parent_sibsize)) parent_sibsize <- fert_mu  # Fallback
          n_II <- sample_offspring_correlated(parent_sibsize, rho, fert_mu)
        }
        
        print(paste("Gen", g, "id", I1s$id[i], ": fert_mu =", round(fert_mu,2), 
                    "parent_sibsize =", I1s$parent_sibsize[i], "n_II =", n_II))
        
        # TRACK (your fixed version)
        i1_idx <- which(df_ped$id == I1s$id[i]) 
        df_ped$fert_rate[i1_idx] <- fert_mu
        df_ped$n_offspring[i1_idx] <- n_II
        I2$fert_rate <- fert_mu
        I2$n_offspring <- n_II

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
      # Move the parent_sibsize assignment INSIDE the if(n_II > 0) block
      if(n_II > 0){
        IIs$parent_sibsize <- n_II  # Only when IIs has rows!
        df_ped <- bind_rows(df_ped, I2, IIs)
      } else {
        df_ped <- bind_rows(df_ped, I2)  # No IIs, no problem
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
add_inlaws = function(df_ped, mean_gen_yr, life_expectancy, current_year){
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
add_ext_branches = function(df_ped, lambda, k, mean_gen_yr, fert_rate, life_expectancy, current_year, rho = 0.15){
  # Initialize tracking columns
  if(!"fert_rate" %in% names(df_ped)) df_ped$fert_rate <- NA_real_
  if(!"n_offspring" %in% names(df_ped)) df_ped$n_offspring <- NA_real_
  if(!"parent_sibsize" %in% names(df_ped)) df_ped$parent_sibsize <- NA_real_

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
        # CORRELATED SAMPLING (minus 1 for existing child)
        birth_yr <- I1s$yob[i] + gen_yr
        fert_mu <- if(is.na(lambda)){
          fert_rate[fert_rate$year == birth_yr, "mean_fertility"]
        } else { lambda }
        
        if(is.na(I1s$parent_sibsize[i]) || g == 0){
          # Independent for founders
          n_II <- rnbinom(1, size = 5, mu = fert_mu) - 1 # minus 1 since already a kid
        } else {
          # Correlated
          n_II <- sample_offspring_correlated(I1s$parent_sibsize[i], rho, fert_mu) - 1 # minus 1 since already a kid
        }
        n_II <- max(0, n_II)  # Ensure non-negative

        # Track for parents
        i1_idx <- which(df_ped$id == I1s$id[i]) 
        if(is.na(df_ped$n_offspring[i1_idx])){
          total_offspring <- 1 + n_II  # Start counting from existing child
        } else {
          total_offspring <- df_ped$n_offspring[i1_idx] + n_II
        }
        df_ped$n_offspring[i1_idx] <- total_offspring
        df_ped$fert_rate[i1_idx] <- fert_mu
        i2_idx <- which(df_ped$id == I2_id)         
        df_ped$n_offspring[i2_idx] <- total_offspring
        df_ped$fert_rate[i2_idx] <- fert_mu

        existing_child_id <- gsub("_[mp]$", "", I1s$id[i])  # P0_1 from P0_1_m
        existing_child_idx <- which(df_ped$id == existing_child_id)
        if(length(existing_child_idx) > 0){
          df_ped$parent_sibsize[existing_child_idx] <- total_offspring
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
          IIs$parent_sibsize <- total_offspring  # Full sibship!
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

        # CORRELATED SAMPLING
        birth_yr <- I1s$yob[i] + gen_yr
        fert_mu <- if(is.na(lambda)){
          fert_rate[fert_rate$year == birth_yr, "mean_fertility"]
        } else { lambda }
        
        if(is.na(I1s$parent_sibsize[i]) || g == 0){
          n_II <- rnbinom(1, size = 5, mu = fert_mu)
        } else {
          n_II <- sample_offspring_correlated(I1s$parent_sibsize[i], rho, fert_mu)
        }

        I2 = data.frame(gen = g, 
                        id = paste0("P",I1s$id[i]),
                        pid = NA, 
                        mid = NA, 
                        sex = abs(I1s$sex[i] - 1),
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

        # Track for both parents
        i1_idx <- which(df_ped$id == I1s$id[i])      
        df_ped$n_offspring[i1_idx] <- n_II
        df_ped$fert_rate[i1_idx] <- fert_mu
        I2$fert_rate <- fert_mu
        I2$n_offspring <- n_II

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
          IIs$parent_sibsize <- n_II
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

# Parameters for simulation
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table("GM_fertility_rate_Netherlands_1800_2100.txt", header=T)
mean_gen_yr = 30
life_expectancy = read.table("OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)
current_year = 2200
plot = FALSE # Do not plot pedigree
rho = 0.15

fert_estimates = data.frame(gen=numeric(),
                            yob_children=numeric(),
                            N_offspring=numeric())
# Initialize collectors OUTSIDE the n-loop (already have fert_estimates)
cor_gen_all <- data.frame(  # NEW: across all sims
  sim = integer(),
  parent_gen = integer(),
  parent_n_offspring = numeric(),
  child_gen = integer(),
  child_n_offspring = numeric()
)

for(n in 1:10){
  cat("simulation:", n,"\n")
  yob_index = sample(seq(1850 + k * mean_gen_yr,2070,1),1)
  print(paste("Selected YOB:", yob_index, ""))
  print(paste("no of gens:", k, "")
  )
  print("Initiate pedigree")
  # 1. Initialize founder with up-to-date function arguments and allele labels
  core_ped <- init_ped(
    k = k,
    mean_gen_yr = mean_gen_yr,
    yob_index = yob_index,
    life_expectancy = life_expectancy,
    current_year = current_year
  )
  print("Add generations")
  # 2. Add generations with up-to-date allele labels & fertility logic
  core_ped <- add_gen(
    df_ped = core_ped,
    lambda = lambda,
    k = k,
    fert_rate = fert_rate,
    mean_gen_yr = mean_gen_yr,
    life_expectancy = life_expectancy,
    current_year = current_year,
    rho = rho
  )

  print("Add inlaws")
  # 3. Add "inlaw" branches (correct allele handling)
  core_ped <- add_inlaws(
    df_ped = core_ped,
    mean_gen_yr = mean_gen_yr,
    life_expectancy = life_expectancy,
    current_year = current_year
  )

  print("Add external branches")
  # 4. Add external branches (correct allele handling)
  core_ped <- add_ext_branches(
    df_ped = core_ped,
    lambda = lambda,
    k = k,
    mean_gen_yr = mean_gen_yr,
    fert_rate = fert_rate,
    life_expectancy = life_expectancy,
    current_year = current_year,
    rho = rho
  )
# write.table(core_ped, "core_ped.tsv", sep="\t", row.names=FALSE, quote=FALSE)

  # INSIDE the n-loop, after add_ext_branches():
  print("Extract correlations")
  cor_gen_n <- data.frame(  # Per simulation
    parent_gen = integer(),
    parent_n_offspring = numeric(),
    child_gen = integer(),
    child_n_offspring = numeric()
  )
  max_gen = max(core_ped$gen)
  if(max_gen > 0){
    for(g in 0:(max_gen-1)) {
      parents = filter(core_ped, gen == g & !is.na(n_offspring) & n_offspring > 0)
      for(i in 1:nrow(parents)) {
        parent_id <- parents$id[i]
        parent_n_off <- parents$n_offspring[i]
        
        # Find kids
        kids <- filter(core_ped, pid == parent_id | mid == parent_id)
        
        if(nrow(kids) > 0) {
          # Kids exist - use their n_offspring (may be NA=0)
          for(n in 1:nrow(kids)) {
            child_n_off <- ifelse(is.na(kids$n_offspring[n]), 0, kids$n_offspring[k])
            cor_gen_n <- bind_rows(cor_gen_n, data.frame(
              parent_gen = g,
              parent_n_offspring = parent_n_off,
              child_gen = g+1,
              child_n_offspring = child_n_off
            ))
          }
        } else {
          # NO kids → record 0 offspring
          cor_gen_n <- bind_rows(cor_gen_n, data.frame(
            parent_gen = g,
            parent_n_offspring = parent_n_off,
            child_gen = g+1,
            child_n_offspring = 0
          ))
        }
      }
    }
  }


  # Add simulation ID and combine
  cor_gen_n$sim <- n
  cor_gen_all <- bind_rows(cor_gen_all, cor_gen_n)

  print("Estimate fertility")
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

# fert_estimates <- fread("./simulated_fertility_estimates_check.txt", header = TRUE, sep = "\t")
fert_rate_plot <- fert_rate %>% filter(year >= 1850 & year <= 2100)

pdf("simulated_fertility_rate_smooth.pdf", height = 3, width = 3)  # Wider for CI
fert_estimates$group <- "Simulated"
ggplot() +
  # Simulated data: smooth trend + 95% CI ribbon
  geom_smooth(data = fert_estimates, 
              aes(x = yob_children, y = N_offspring, color = "Simulated"), 
              method = "loess", span = 0.1, se = TRUE, linewidth = 1) +
  
  # Historical data: exact line
  geom_line(data = fert_rate_plot, 
            aes(x = year, y = mean_fertility, color = "Historical"), 
            linewidth = 1) +
  
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Data") +
  xlab("Year") + ylab("Mean fertility rate") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()

pdf("simulated_fertility_rate_mean.pdf", height = 3, width = 3)  # Wider for CI
ggplot(data = fert_estimates, aes(x = yob_children, y = N_offspring, color = "Simulated")) +
  stat_summary(fun = mean, geom = "line") + 
  # Historical data: exact line
  geom_line(data = fert_rate_plot, 
            aes(x = year, y = mean_fertility, color = "Historical"), 
            linewidth = 1) +
  
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Data") +
  xlab("Year") + ylab("Mean fertility rate") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

dev.off()


# Histogram of simulated years of birth
pdf("simulated_yob_histogram.pdf", height = 3, width = 6)
ggplot(fert_estimates, aes(x = yob_children)) +
  geom_histogram(bins = 30, fill = wes_palette(n = 1, name = "Darjeeling1")[1], 
                 color = "black", alpha = 0.8, linewidth = 0.3) +
  geom_density(color = wes_palette(n = 1, name = "Darjeeling1")[2], 
               linewidth = 1, alpha = 0.7) +
  labs(title = "Distribution of Simulated Years of Birth",
       subtitle = paste("N =", nrow(fert_estimates), "individuals"),
       x = "Year of Birth", y = "Count") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

# Optional: Print summary stats
cat("YOB Summary:\n")
print(summary(fert_estimates$yob_children))
cat("N individuals:", nrow(fert_estimates), "\n")

# Final correlation analysis
cor_summary <- cor_gen_all %>%
  group_by(parent_gen, child_gen) %>%
  summarise(
    n_pairs = n(),
    cor = round(cor(parent_n_offspring, child_n_offspring), 3),
    .groups = "drop"
  )

print(cor_summary)
fwrite(cor_gen_all, "cor_gen_all_sims.txt", sep="\t")
fwrite(cor_summary, "cor_gen_summary.txt", sep="\t")

# Plot
pdf("intergen_correlations.pdf", height=4, width=8)
ggplot(cor_gen_all, aes(x=parent_n_offspring, y=child_n_offspring)) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", se=TRUE) +
  facet_wrap(~paste("Gen", parent_gen, "→", child_gen)) +
  ggtitle("Intergenerational Fertility Correlation (All Sims)")
dev.off()

fwrite(fert_estimates, "./simulated_fertility_estimates_check.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
