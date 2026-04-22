#' loading libraries and functions for src folder seems a bit complicated. However, makes sure the script can be run from any location as long as the file/directory structure is maintained as in the github repository.
#' Could be improved...


### Helper functions ### 

#' make sure warnings are treated as errors
options( warn = 2 )

source("../../src/libraries_simPed.R")

#' ### Function to initiate pedigree with one founder with a mutation
init_ped = function(DAF_common){
  ped = data.frame(gen = 0, # generation
                        id = "C0", # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1, # sex, 0 = male, 1 = female
                        a1_common  = 1, # sample disease allele from population frequency
                        a2_common  = 0) # sample disease allele from population frequency
  return(ped)
}

#' ### Function to add next generation to core pedigree
add_gen = function(df_ped, lambda, k, DAF_common, fert_rate){
  g = 0
  while(g < k){
    # select individuals from youngest generation
    I1s = filter(df_ped, gen == g) # call I1 for parent 1
    # check if there are individuals in this generation, if not re-simulate generation
    if(nrow(I1s) > 0) {
      # for each individual
      for(i in 1:nrow(I1s)){
        # simulate partner (I2 for parent 2)
        I2 = data.frame(gen = g, 
                        id  = gsub("C", "P", I1s$id[i]),
                        pid = NA, # at this point parents are irrelevant, will be simulated in later step
                        mid = NA, # at this point parents are irrelevant, will be simulated in later step
                        sex = abs(I1s$sex[i] - 1), # opposite sex partners only...
                        a1_common  = 0, # sample disease allele from population frequency
                        a2_common  = 0 # sample disease allele from population frequency
                        )
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        n_II = rpois(1, lambda)        # make sure generation 0 gets offspring
        while(n_II == 0 & g == 0){
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
add_inlaws = function(df_ped, DAF_common){
  adj_ped = filter(df_ped, grepl("P", id))
  g = max(adj_ped$gen)
  while(! g == 0){
    IIs = filter(adj_ped, gen == g) # get individuals from offspring generation
    for(i in 1:nrow(IIs)){
      # simulate parents 1 (I1)
      I1 = data.frame(gen = g-1, 
                      id  = paste0(IIs$id[i], "_m"), # "_m" suffix for mother
                      pid = NA, 
                      mid = NA, 
                      sex = 1, # mother
                      a1_common  = ifelse(IIs$a1_common[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
                      a2_common  = 0) # non-transmitted allele sampled from population
      # simulate father
      I2 = data.frame(gen = g-1, 
                      id = paste0(IIs$id[i], "_p"), # "_p" suffix for father
                      pid = NA, 
                      mid = NA, 
                      sex = 0,
                      a1_common = ifelse(IIs$a2_common[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
                      a2_common = 0) # non-transmitted allele sampled from population
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


add_ext_branches = function(df_ped, lambda, k, DAF_common, fert_rate){
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
        # find father that is already in pedigree
        I2_id = gsub("_m$", "_p", I1s$id[i])
        I2 = filter(df_ped, id == I2_id)
        if(! nrow(I2 == 1)){
          print(I2)
          print(I2_id)
          stop("Non-unique IDs found!!")
        }
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        n_II = rpois(1, lambda) - 1 # minus 1 because one child has already been simulated in first round
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
            df_ped = bind_rows(df_ped, IIs)
          }  
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
        I2 = data.frame(gen = g, 
                        id = paste0("P",I1s$id[i]),
                        pid = NA, 
                        mid = NA, 
                        sex = abs(I1s$sex[i] - 1),
                        a1_common  = 0,
                        a2_common  = 0)         
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
        n_II = rpois(1, lambda)
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
add_pheno = function(df_ped, h2_ALS, K_ALS, h2_FTD, K_FTD, rg) {
  
  # polygenic model:
  df_ped$G_ALS = NA
  df_ped$G_FTD = NA
  
  # Sample G from multivariate normal for founders
  founders = which(is.na(df_ped$pid)) 
  # ALS and FTD, correlated
  Vg = matrix(c(h2_ALS, rg * sqrt(h2_ALS * h2_FTD), 
              rg * sqrt(h2_ALS * h2_FTD), h2_FTD), nrow=2, byrow=TRUE)
  Gs = mvrnorm(length(founders), mu = c(0, 0), Sigma = Vg)

  df_ped$G_ALS[founders] = Gs[,1]
  df_ped$G_FTD[founders] = Gs[,2]

  empirical_rg = cor(df_ped$G_ALS[founders], df_ped$G_FTD[founders])
  
  # loop through nonfounders and estimate G based on G of parents and h2
  for(i in 1:nrow(df_ped)){
    if(is.na(df_ped$G_ALS[i])){
      pid = df_ped$pid[i]
      mid = df_ped$mid[i]
      Gp_ALS = df_ped[which(df_ped$id == pid),]$G_ALS[1]
      Gm_ALS = df_ped[which(df_ped$id == mid),]$G_ALS[1]
      Gp_FTD = df_ped[which(df_ped$id == pid),]$G_FTD[1]
      Gm_FTD = df_ped[which(df_ped$id == mid),]$G_FTD[1]
      
      # Calculate mean and variance for offspring's G
      mean_ALS = (Gp_ALS + Gm_ALS) / 2
      var_ALS = 0.5 * h2_ALS
      mean_FTD = (Gp_FTD + Gm_FTD) / 2
      var_FTD = 0.5 * h2_FTD
      
      # Simulate correlated G for offspring
      Vg_offspring = matrix(c(var_ALS, rg * sqrt(var_ALS * var_FTD), 
                              rg * sqrt(var_ALS * var_FTD), var_FTD), nrow=2, byrow=TRUE)
      # cat(sprintf("  mean_ALS = %0.2f\n", mean_ALS))
      # cat(sprintf("  mean_FTD= %0.2f\n", mean_FTD))                 
      # cat(sprintf("  var_ALS = %0.2f\n", var_ALS))
      # cat(sprintf("  var_FTD= %0.2f\n", var_FTD))    
      # cat(sprintf("  rg = %0.2f\n", rg))    
      # cat("Vg_offspring:")
      # print(Vg_offspring)           
      Gs_offspring = mvrnorm(1, mu = c(mean_ALS, mean_FTD), Sigma=Vg_offspring)
      
      df_ped$G_ALS[i] = Gs_offspring[1]
      df_ped$G_FTD[i] = Gs_offspring[2]
    }
  }
  
  # sample non-genetic value E:
  # Without genetic correlation
  df_ped$E_ALS = rnorm(nrow(df_ped), 0, sqrt(1-h2_ALS))
  df_ped$E_FTD = rnorm(nrow(df_ped), 0, sqrt(1-h2_FTD))

  # Polygenic ALS phenotype using constant (lifetime) risk
  df_ped$P_ALS = df_ped$G_ALS + df_ped$E_ALS
  df_ped$LT_ALS = -qnorm(K_ALS, 0, 1)
  df_ped$polygenicY_ALS = ifelse(df_ped$P_ALS > df_ped$LT_ALS, 1, 0)

  # Polygenic FTD phenotype using constant (lifetime) risk
  df_ped$P_FTD = df_ped$G_FTD + df_ped$E_FTD
  df_ped$LT_FTD = -qnorm(K_FTD, 0, 1)
  df_ped$polygenicY_FTD = ifelse(df_ped$P_FTD > df_ped$LT_FTD, 1, 0)
 
  # final phenotype
  df_ped$Y_ALS = ifelse(df_ped$polygenicY_ALS > 0, 1, 0)
  df_ped$Y_FTD = ifelse(df_ped$polygenicY_FTD > 0, 1, 0)
  
  return(df_ped)
}

# Calculate the theoretical first-degree recurrence risk (lambda1)
lambda1_theory <- function(K, h2) {
  T  <- -qnorm(K)
  z  <- dnorm(T)
  i  <- z / K
  # Heritability must be between 0 and 1
  # Equation for first-degree relative (n = 1)
  numer  <- T - h2 * i / 2
  denom  <- sqrt(1 - h2^2 / 4)
  K1     <- 1 - pnorm(numer / denom)
  lambda1 <- K1 / K
  return(lambda1)
}

lambda2_theory <- function(K, h2) {
  T  <- -qnorm(K)
  z  <- dnorm(T)
  i  <- z / K
  
  # Equation for second-degree relative (n=2)
    # numerator inside pnorm argument
  numer  <- T - h2 * i / 4
    # denominator inside pnorm argument
  denom  <- sqrt(1 - h2^2 / 16)
  K2     <- 1 - pnorm(numer / denom)
  lambda2 <- K2 / K
  return(lambda2)
}

#' argument parser
parser <- ArgumentParser()
parser$add_argument("-h2_als", "--heritability_als", type="character", required=TRUE,
                    help="Comma-separated list of ALS heritability values, e.g. 0.2,0.4,0.6,0.8")
parser$add_argument("-K_als", "--lifetimerisk_als", type="character", required=TRUE,
                    help="Comma-separated list of ALS lifetime risk values, e.g. 0.2,0.4,0.6,0.8")
parser$add_argument("-h2_ftd", "--heritability_ftd", type="character", required=TRUE,
                    help="Comma-separated list of FTD heritability values, e.g. 0.2,0.4,0.6,0.8")
parser$add_argument("-K_ftd", "--lifetimerisk_ftd", type="character", required=TRUE,
                    help="Comma-separated list of FTD lifetime risk values, e.g. 0.2,0.4,0.6,0.8")
parser$add_argument("-rg", "--genetic_correlation", type="character", required=TRUE,
                    help="Comma-separated list of ALS~FTD genetic correlation values, e.g. 0.2,0.4,0.6,0.8,0,10")            
parser$add_argument("-l", "--lambda", type="character", required=TRUE,
                    help="Comma-separated list of lambda values, e.g. 1,2,3,6")
parser$add_argument("-k", "--generations", type="character", required=TRUE,
                    help="Comma-separated list of generations after founder, e.g. 2,3 [default %(default)s]")
parser$add_argument("-P", "--peds", type="numeric", default=500,
                    help="Number of pedigrees per simulation [default %(default)s]")
parser$add_argument("-o", "--out", type="character", default="polygenic_inheritance.out",
                    help="Output filename [default %(default)s]")
parser$add_argument("-t", "--times", type="numeric", default=1,
                    help="Number of times to repeat the simulation per parameter set [default %(default)s]")                    
args <- parser$parse_args()

# Parse penetrance and lambda vectors
h2_als_vals <- as.numeric(strsplit(args$heritability_als, ",")[[1]])
K_als_vals  <- as.numeric(strsplit(args$lifetimerisk_als, ",")[[1]])
h2_ftd_vals <- as.numeric(strsplit(args$heritability_ftd, ",")[[1]])
K_ftd_vals  <- as.numeric(strsplit(args$lifetimerisk_ftd, ",")[[1]])
rg_vals     <- as.numeric(strsplit(args$genetic_correlation, ",")[[1]])
lambda_vals <- as.numeric(strsplit(args$lambda, ",")[[1]])
k_vals      <- as.numeric(strsplit(args$generations, ",")[[1]])
n_peds      <- args$peds
n_times     <- args$times
DAF_common  <- 0

#' Compare theory and simulations for polygenic trait parameters
sim_polygenic_all <- function(N, peds, k, lambda, K_als, h2_als, K_ftd, h2_ftd, rg) {
  obs_K  <- rep(NA, N)
  lambda1_emp <- rep(NA, N)
  lambda1_theor <- rep(NA, N)
  lambda2_emp <- rep(NA, N)
  lambda2_theor <- rep(NA, N)
  
  for(x in 1:N){
    polygenic_pedigrees <- vector("list", peds)

    # Simulate pedigrees once
    for(h in 1:peds){
      polygenic_pedigrees[[h]] <- init_ped(DAF_common = 0) %>%
        add_gen(lambda = lambda, k = k, DAF_common = 0, fert_rate = NULL) %>%
        add_inlaws(DAF_common = 0) %>%
        add_ext_branches(DAF_common = 0, lambda = lambda, k = k, fert_rate=NULL) %>%
        add_pheno(h2_ALS = h2_als, K_ALS = K_als,
                  h2_FTD = h2_ftd, K_FTD = K_ftd, rg = rg)
    }

    # Initialize counters for lambda1 & lambda2 metrics
    population_N   <- 0
    affected_N     <- 0
    first_degree_N <- 0
    first_degree_Y1 <- 0
    second_degree_N <- 0
    second_degree_Y1 <- 0

    for(h in seq_along(polygenic_pedigrees)){
      pedigree       <- polygenic_pedigrees[[h]]
      population_N   <- population_N + nrow(pedigree)
      affected       <- dplyr::filter(pedigree, Y_ALS == 1)
      affected_N     <- affected_N + nrow(affected)

      for(affected_id in affected$id){
        # Offspring of affected
        offspring    <- dplyr::filter(pedigree, mid == affected_id | pid == affected_id)
        
        # Parents of affected
        individual   <- dplyr::filter(pedigree, id == affected_id)
        pid_val      <- individual$pid[1]  # extract scalar
        mid_val      <- individual$mid[1]  # extract scalar

        parent_ids   <- c(pid_val, mid_val)
        parent_ids   <- parent_ids[!is.na(parent_ids)]
        parents      <- dplyr::filter(pedigree, id %in% parent_ids)

        # Siblings: individuals sharing both parents but excluding self
        sibs <- dplyr::filter(pedigree, 
                              pid == pid_val & 
                              mid == mid_val & 
                              id != affected_id)
        
        # Update counts: sum of parents + offspring + siblings
        first_degree_N  <- first_degree_N + nrow(offspring) + nrow(parents) + nrow(sibs)
        first_degree_Y1 <- first_degree_Y1 + sum(offspring$Y_ALS) + sum(parents$Y_ALS) + sum(sibs$Y_ALS)

        ## Second degree risk
        individual <- dplyr::filter(pedigree, id == affected_id)
        # Using your add_2nd_and_3rd_relatives_als() logic, implement similar here:
        
        # GRANDPARENTS
        parent_ids <- c(individual$mid[1], individual$pid[1])
        grandparent_ids <- c()
        for(pid in parent_ids){
          if(!is.na(pid)){
            p_row <- dplyr::filter(pedigree, id == pid)
            grandparent_ids <- c(grandparent_ids, p_row$mid, p_row$pid)
          }
        }
        grandparent_ids <- grandparent_ids[!is.na(grandparent_ids)]
        
        # AUNTS AND UNCLES
        aunt_uncle_ids <- c()
        for(pid in parent_ids){
          if(!is.na(pid)){
            p_row <- dplyr::filter(pedigree, id == pid)
            if(nrow(p_row) > 0){
              grand_pid <- p_row$mid[1]
              grand_mid <- p_row$pid[1]
              
              siblings1 <- dplyr::filter(pedigree, mid == grand_mid & pid == grand_pid & id != pid)
              aunt_uncle_ids <- c(aunt_uncle_ids, siblings1$id)
            }
          }
        }
        
        # GRANDCHILDREN
        children_ids <- c(dplyr::filter(pedigree, mid == affected_id | pid == affected_id)$id)
        grandchildren_ids <- c()
        for(child_id in children_ids){
          grandchildren_ids <- c(grandchildren_ids, dplyr::filter(pedigree, mid == child_id | pid == child_id)$id)
        }
        
        # NAMES of 2nd degree relatives (combine all)
        second_degree_ids <- unique(c(grandparent_ids, aunt_uncle_ids, grandchildren_ids))
        second_degree_rows <- dplyr::filter(pedigree, id %in% second_degree_ids)
        
        if(nrow(second_degree_rows) > 0){
          second_degree_N <- second_degree_N + nrow(second_degree_rows)
          second_degree_Y1 <- second_degree_Y1 + sum(second_degree_rows$Y_ALS)
        }
      }
    }
    
    # Calculate prevalence and empirical lambdas
    prevalence_pop <- affected_N / population_N
    
    lambda1_emp[x] <- first_degree_Y1 / max(1, first_degree_N) / prevalence_pop
    lambda2_emp[x] <- second_degree_Y1 / max(1, second_degree_N) / prevalence_pop
    obs_K[x] <- prevalence_pop
    
    lambda1_theor[x] <- lambda1_theory(K_als, h2_als)
    lambda2_theor[x] <- lambda2_theory(K_als, h2_als)

    cat("SAMPLE", x,
        ": lambda1_emp =", round(lambda1_emp[x], 3),
        ", lambda1_theor =", round(lambda1_theor[x], 3),
        ", lambda2_emp =", round(lambda2_emp[x], 3),
        ", lambda2_theor =", round(lambda2_theor[x], 3), "\n")
  }

  return(list(
    K = obs_K,
    lambda1_empirical = lambda1_emp,
    lambda1_theory = lambda1_theor,
    lambda2_empirical = lambda2_emp,
    lambda2_theory = lambda2_theor
  ))
}

# Grid wrapper calling sim_polygenic for each param combo
results <- bind_rows(
  lapply(h2_als_vals, function(h2a) {
    cat(sprintf("Running simulations for h2_als = %0.2f\n", h2a))
    bind_rows(lapply(K_als_vals, function(Ka) {
      cat(sprintf("  K_ALS = %0.3f\n", Ka))
      bind_rows(lapply(h2_ftd_vals, function(h2f) {
        bind_rows(lapply(K_ftd_vals, function(Kf) {
          bind_rows(lapply(rg_vals, function(rg) {
            bind_rows(lapply(lambda_vals, function(lam) {
              bind_rows(lapply(k_vals, function(gen) {
                cat(sprintf("    Lambda = %0.1f, Gen = %d, rg = %0.2f\n", lam, gen, rg))
                
                # Call the combined simulation function once
                sim_results <- sim_polygenic_all(
                  N = n_times,
                  peds = n_peds,
                  k = gen,
                  lambda = lam,
                  K_als = Ka,
                  h2_als = h2a,
                  K_ftd = Kf,
                  h2_ftd = h2f,
                  rg = rg
                )
                
                # Organize replicate-wise results as a tibble
                tibble(
                  h2_als = h2a,
                  K_als = Ka,
                  h2_ftd = h2f,
                  K_ftd = Kf,
                  rg = rg,
                  lambda = lam,
                  gen = gen,
                  replicate = seq_len(n_times),
                  K_empirical = sim_results$K,
                  lambda1_empirical = sim_results$lambda1_empirical,
                  lambda1_theory = sim_results$lambda1_theory,
                  lambda2_empirical = sim_results$lambda2_empirical,
                  lambda2_theory = sim_results$lambda2_theory
                )
              }))
            }))
          }))
        }))
      }))
    }))
  })
)
cat("All simulations complete.\n")

# Save output
write.table(results, file = args$out, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Results written to", args$out, "\n")