## This is a script to provide a rationale for the way the simulations work under a polygenic disease assumption. 
## Use this is as a rationale and evidence that inheritance patterns follow actual polygenic inheritance 


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
                      a1_common  = 0, # assume mother always transmits a1 for coding convenience
                      a2_common  = 0) # non-transmitted allele sampled from population
      # simulate father
      I2 = data.frame(gen = g-1, 
                      id = paste0(IIs$id[i], "_p"), # "_p" suffix for father
                      pid = NA, 
                      mid = NA, 
                      sex = 0,
                      a1_common = 0, # assume father always transmits a2 for coding convenience
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

#' compare theory and simulations for polygenic trait parameters
#' Compare theory and simulations for polygenic trait parameters
sim_polygenic = function(N, peds, k, lambda, K_als, h2_als, K_ftd, h2_ftd, rg) {
  obs_h2 <- rep(NA, N)
  obs_K  <- rep(NA, N)
  
  for(x in 1:N){
    polygenic_pedigrees <- vector("list", peds)
    
    # --- simulate multiple pedigrees ---
    for(h in 1:peds){
      polygenic_pedigrees[[h]] <- init_ped(DAF_common = 0) %>%
        add_gen(lambda = lambda, k = k, DAF_common = 0, fert_rate = NULL) %>%
        add_inlaws(DAF_common = 0) %>%
        add_pheno(h2_ALS = h2_als, K_ALS = K_als,
                  h2_FTD = h2_ftd, K_FTD = K_ftd, rg = rg)
    }
    
    # --- summary counters ---
    population_N   <- 0
    affected_N     <- 0
    offspring_N    <- 0
    offspring_Y1   <- 0
    
    for(h in seq_along(polygenic_pedigrees)){
      pedigree       <- polygenic_pedigrees[[h]]
      population_N   <- population_N + nrow(pedigree)
      affected       <- dplyr::filter(pedigree, Y_ALS == 1)
      affected_N     <- affected_N + nrow(affected)
      for(affected_id in affected$id){
        offspring    <- dplyr::filter(pedigree, mid == affected_id | pid == affected_id)
        offspring_N  <- offspring_N + nrow(offspring)
        offspring_Y1 <- offspring_Y1 + sum(offspring$Y_ALS)
      }
    }
    
    # --- empirical prevalence across all pedigrees ---
    obs_K[x] <- affected_N / population_N
    
    # --- Falconer-style recurrence heritability ---
    LT  <- -qnorm(K_als, 0, 1)
    Z   <- dnorm(LT)
    i   <- Z / K_als
    
    Kr  <- offspring_Y1 / max(1, offspring_N)  # guard div/0
    if(Kr > 0 && Kr < 1) {
      LTr <- -qnorm(Kr, 0, 1)
      Zr  <- dnorm(LTr)
      aR  <- 0.5
      inside <- 1 - (1 - LT/i) * (LT^2 - LTr^2)
      if(inside < 0) inside <- 0   # clamp numerical artifacts
      obs_h2[x] <- (LT - LTr * sqrt(inside)) / (aR * (i + (i - LT) * LTr^2))
    } else {
      obs_h2[x] <- NA
    }
    
    cat("SAMPLE", x, ": h2 =", round(obs_h2[x], 3),
        ", obs_K =", round(obs_K[x], 3),
        "(input K =", K_als, ")\n")
  }
  
  return(list(h2 = obs_h2, K = obs_K))
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
                # Run n_times replicates per parameter set
                sim_results <- sim_polygenic(
                  N = n_times,
                  peds = n_peds,
                  k = gen,
                  lambda = lam,
                  K_als = Ka,
                  h2_als = h2a,
                  K_ftd = Kf,       # scalar lifetime risk for FTD (single value)
                  h2_ftd = h2f,     # scalar heritability for FTD (single value)
                  rg = rg           # scalar genetic correlation
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
                  h2_empirical = sim_results$h2,
                  K_empirical  = sim_results$K
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