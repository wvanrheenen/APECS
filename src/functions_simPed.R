#' ## To do's
#TODO: recode lambda to vector with length of k to vary over generations
#TODO: calculate penetrance per pedigree and check if this works
#TODO: Pleiotropy and genetic correlation: extend to two correlated traits

#' ### Function to initiate pedigree with one founder with a mutation
init_ped = function(monogenic=TRUE){
  
  core_ped = data.frame(gen = 0, # generation
                        id = "C0", # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1, # sex, 0 = male, 1 = female
                        a1 = 0) # first allele at disease locus
  if(monogenic){
    core_ped$a2 = 1 # disease allele at second allele at disease locus
  }else{
    core_ped$a2 = 0 # second allele at disease locus   
  }
  return(core_ped)
}

#' ### Function to simulate a next generation
add_gen = function(df_ped, lambda, k, DAF){
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
                        id = gsub("C", "P", I1s$id[i]),
                        pid = NA, # at this point parents are irrelevant, will be simulated in later step
                        mid = NA, # at this point parents are irrelevant, will be simulated in later step
                        sex = abs(I1s$sex[i] - 1), # opposite sex partners only...
                        a1 = sample(c(0,1), 1, prob=c(1-DAF, DAF)), # sample disease allele from population frequency
                        a2 = sample(c(0,1), 1, prob=c(1-DAF, DAF))) # sample disease allele from population frequency
        # simulate number of offspring from Poisson distribution with mean lambda (as defined by general pedigree parameters)
        n_II = rpois(1,lambda)
        # make sure generation 0 gets offspring
        while(n_II == 0 & g == 0){
          n_II = rpois(1,lambda)
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
          IIs$mid[j] = ifelse(I1s$sex[i] == 0, I2$id, I1s$id[i])
          IIs$sex[j] = sample(c(0,1), 1)
          IIs$a1[j]  = sample(c(I1s$a1[i], I1s$a2[i]), 1)
          IIs$a2[j]  = sample(c(I2$a1, I2$a2), 1)
        }
        if(n_II > 0){
          df_ped = bind_rows(df_ped, I2, IIs)
        }
      }
      g = g+1
    } else {
      g = g-1
    }
  }
  return(df_ped)
}

#' ### Function to simulate the "inlaws"
#' These are the ancestors for those who married into this pedigree
add_inlaws = function(df_ped, DAF){
  adj_ped = filter(df_ped, grepl("P", id))
  g = max(adj_ped$gen)
  while(! g == 0){
    IIs = filter(adj_ped, gen == g) # get individuals from offspring generation
    for(i in 1:nrow(IIs)){
      # simulate parents 1 (I1)
      I1 = data.frame(gen = g-1, 
                      id = paste0(IIs$id[i], "_m"), # "_m" suffix for mother
                      pid = NA, 
                      mid = NA, 
                      sex = 1, # mother
                      a1 = ifelse(IIs$a1[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
                      a2 = sample(c(0,1), 1, prob=c(1-DAF, DAF))) # non-transmitted allele sampled from population
      # simulate father
      I2 = data.frame(gen = g-1, 
                      id = paste0(IIs$id[i], "_p"), # "_p" suffix for father
                      pid = NA, 
                      mid = NA, 
                      sex = 0,
                      a1 = ifelse(IIs$a2[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
                      a2 = sample(c(0,1), 1, prob=c(1-DAF, DAF))) # non-transmitted allele sampled from population
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
add_ext_branches = function(df_ped, lambda, k, DAF){
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
        # simulate number of offspring from Poisson distribution with mean lambda
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
            IIs$pid[j] = I2$id
            IIs$mid[j] = I1s$id[i]
            IIs$sex[j] = sample(c(0,1), 1)
            IIs$a1[j]  = sample(c(I1s$a1[i], I1s$a2[i]), 1) # a1 is always from mother - here, I1
            IIs$a2[j]  = sample(c(I2$a1, I2$a2), 1)
          }
          df_ped = bind_rows(df_ped, IIs)
        }
      }
    }
    
    # Second, find NF who don't have a partner yet (offspring simulated in previous loop)
    # these have IDs ending with a digit, but do have either "p" or "m" in their IDs
    I1s = filter(df_ped, gen == g & 
                   (grepl("p", id) | grepl("m", id)) &
                   ! (grepl("p$", id) | grepl("m$", id)))
    if(nrow(I1s > 0)){
      for(i in 1:nrow(I1s)){
        # simulate partner
        # cat("simulate partner...")
        I2 = data.frame(gen = g, 
                        id = paste0("P",I1s$id[i]),
                        pid = NA, 
                        mid = NA, 
                        sex = abs(I1s$sex[i] - 1),
                        a1 = sample(c(0,1), 1, prob=c(1-DAF, DAF)),
                        a2 = sample(c(0,1), 1, prob=c(1-DAF, DAF)))
        # simulate number of offspring from Poisson distribution with mean 2
        n_II = rpois(1,lambda)
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
          IIs$a1[j]  = sample(c(I1s$a1[i], I1s$a2[i]), 1)
          IIs$a2[j]  = sample(c(I2$a1, I2$a2), 1)
        }
        if(n_II > 0){
          df_ped = bind_rows(df_ped, I2, IIs)
        }
      }
    }
    g = g+1
  }
  return(df_ped)
}

#' ### Function to add genetic values for Mendelian and polygenic inheritance
add_pheno = function(df_ped, penetrance, h2, K){
  
  # Mendelian (autosomal dominant)
  df_ped = mutate(df_ped, mendel_y1 = ifelse(a1 == 1, sample(c(0,1), 1, prob=c(1-penetrance, penetrance)), 0))
  df_ped = mutate(df_ped, mendel_y2 = ifelse(a2 == 1, sample(c(0,1), 1, prob=c(1-penetrance, penetrance)), 0))
  df_ped = mutate(df_ped, mendelY = ifelse(mendel_y1 + mendel_y2 > 0, 1, 0))
  
  # polygenic model:
  df_ped$G = NA
  # sample G from N(0,h2) for founders
  founders = which(is.na(df_ped$pid)) # select founders in pedigree for which no parental information is available.
  df_ped$G[founders] = rnorm(length(founders), mean=0, sd=sqrt(h2))
  # loop through nonfounders and estimate G based on G of parents and h2
  # for offspring, E(G) = (Gm+Gp)/2
  # the sampling variance is 0.5 * h2 (see Lynch and Walsh)
  for(i in 1:nrow(df_ped)){
    if(is.na(df_ped$G[i])){
      pid = df_ped$pid[i]
      mid = df_ped$mid[i]
      Gp = df_ped[which(df_ped$id == pid),]$G
      Gm = df_ped[which(df_ped$id == mid),]$G
      df_ped$G[i] = rnorm(1, mean= (Gp+Gm)/2, sd = sqrt(0.5*h2)) 
    }
  }
  # sample non-genetic value E:
  df_ped$E = rnorm(nrow(df_ped), 0, sqrt(1-h2))
  # define phenotype
  df_ped$P = df_ped$G + df_ped$E
  LT = -qnorm(K, 0, 1) # liability threshold
  df_ped$polygenicY = ifelse(df_ped$P > LT, 1, 0)
  
  # final phenotype
  df_ped$Y = ifelse(df_ped$polygenicY + df_ped$mendelY > 0, 1, 0)
  
  return(df_ped)
}

#' ### Alternative function to add genetic values for Mendelian and polygenic inheritance
#' This function samples G from multivariate normal distribution.
#' This adds genetic values for polygenic inheritance FAST for smaller pedigrees
add_pheno_small = function(df_ped, penetrance, h2, K){
  
  # Mendelian (autosomal dominant)
  df_ped = mutate(df_ped, mendel_y1 = ifelse(a1 == 1, sample(c(0,1), 1, prob=c(1-penetrance, penetrance)), 0))
  df_ped = mutate(df_ped, mendel_y2 = ifelse(a2 == 1, sample(c(0,1), 1, prob=c(1-penetrance, penetrance)), 0))
  df_ped = mutate(df_ped, mendelY = ifelse(mendel_y1 + mendel_y2 > 0, 1, 0))
  
  # polygenic model:
  # get kinship coefficients from pedigree
  km = kinship(ped(id = df_ped$id, 
                   fid = df_ped$pid,
                   mid = df_ped$mid,
                   sex = df_ped$sex + 1,
                   isConnected = TRUE))
  # define variance covariance matrix of G in pedigree (GRM)
  Vg = km*2*h2
  # sample G from multivariate normal distribution
  N_ind = nrow(df_ped)
  # df_ped$G = mvrnorm(1, mu = rep(0, N_ind), Sigma = Vg) # SLOWER - TRADITIONAL - ALTERNATIVE
  df_ped$G = rmvn(1, mu = rep(0, N_ind), sigma = Vg)[1,]
  # sample non-genetic value E (no shared environment):
  df_ped$E = rnorm(nrow(df_ped), 0, sqrt(1-h2))
  # define phenotype using liability threshold
  df_ped$P = df_ped$G + df_ped$E
  LT = -qnorm(K, 0, 1) # liability threshold
  df_ped$polygenicY = ifelse(df_ped$P > LT, 1, 0)
  
  # final phenotype
  df_ped$Y = ifelse(df_ped$polygenicY + df_ped$mendelY > 0, 1, 0)
  
  return(df_ped)
}

#' wrapper function to simulate full pedigree
sim_ped = function(i=numeric, 
                   k=numeric(), lambda=numeric(), 
                   monogenic=TRUE, DAF=numeric(), penetrance=numeric(), 
                   K=numeric(), h2=numeric(),
                   small=TRUE,
                   plot=TRUE){
  core_ped = init_ped(monogenic=monogenic)
  core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF=DAF)
  core_ped = add_inlaws(core_ped, DAF=DAF)
  core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF)
  if(small) {
    core_ped = add_pheno_small(core_ped, penetrance, h2, K) # traditional ALTERNATIVE
  } else {
    core_ped = add_pheno(core_ped, penetrance, h2, K) 
  }
  if(plot){
    ped_plt = ped(id = core_ped$id, 
                  fid = core_ped$pid,
                  mid = core_ped$mid,
                  sex = core_ped$sex + 1,
                  isConnected = TRUE,
                  reorder = FALSE)
    carriers = filter(core_ped, a1 + a2 > 0)$id
    affected = filter(core_ped, Y == 1)$id
    # color according to percentile from N(0,1)
    pal = colorRampPalette(c("white", "orange", "red"))(1000) 
    Pcolors = pal[ceiling(pnorm(core_ped$P)*1000)]
    # make plot
    pdf(paste0("PED",i,"_by_polygenic_P.pdf"), w=30, h=12)
    plot(ped_plt, title="step 4 - simulated polygenic phenotype on liability scale", carrier = carriers, fill=Pcolors)
    dev.off()
    pdf(paste0("PED",i,"_by_phenotype.pdf"), w=30, h=12)
    plot(ped_plt, title="step 4 - colored by binary phenotype", carrier = carriers, aff = affected)
    dev.off()
  }
  return(core_ped)
}
