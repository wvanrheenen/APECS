# SimPlex

A framework to simulate ALS and ALS-associated disease under a monogenic/Mendelian and polygenic/complex disease model.

================

Wouter van Rheenen, Paul Beele
23 april 2026

## General outline of simulation scheme

### Simulating family members:

The pedigree consists of core pedigree (C) with one single founder this
is initiated by function `init_ped()`. The birthyear of each founder is selected
by taking the birthyear of the final `k`th generation and sampling `k` times
from a normal distibution with mean = 30, truncated at 25 to 35 years back. 

Then offspring generation will be simulated using the `add_gen()`
function, IDs start with “C” and the number reflects order of offspring
(C0_1 for oldest, C0_2 for second child). In the third generation offspring
of the oldest is denoted as C0_1\_\[0-9\] and second C0_2\_\[0-9\] etc.
etc... This way, all individuals in this lineage IDs can be traced to
founder. 

The `add_gen()` function starts with simulating the spouse. The
number of offspring is sampled from a negative binomial distribution with mean
number of offspring $\lambda =$ `lambda` or, if available, mean fertility provided
by demographic data, based on the birthyear of the offspring. To guarantee 
offspring of the "core pedigree" (i.e. individuals in the lineage whose 
ID starts with 'C'), we force at least one of offspring for the eldest child
in each generation of the core pedigree member. 

The spouses married into this pedigree are denoted with “P\*” and have same code as partner 
(P0_1 for partner of first child in generation 1). The number of generations
to be added is defined by `k`, so for `k=2` a three-generation pedigree is
simulated.

When all generations of core pedigree are simulated, we simulate the
in-laws e.g. the ancestors of the married in spouses (P\*) using the
`add_inlaws()` function. To complete the pedigree, we simulate offspring 
of the inlaws using `add_ext_branches()`. These individuals are unlinked
to the core pedigree, e.g. sibs, cousins, etc. of the spouses married 
into this pedigree.

### Simulating demographic data:
For each individual, we sample a life expectancy from a left-skewed
normal distribution, with the mean life expectancy based on demographic
life expectancies for an individuals year of birth. 

At this point in the simulation, we first assess if the if the individual
is either alive or dead based on whether the `current year` of the simulation
lies before or after the year where an individual is simulated to decease,
based on their life expectancy. 

### Genetics and phenotypes

#### Mendelian disease alleles

We simulated carriership of pathogenic mutations based no disease allele
frequencies. We simulated an ALS-FTD common (C9orf72-like) and rare but more
pathogenic (FUS/SOD1-like) disease allele, and an FTD-specific (GRN/MAPT-like)
disease allele. In the `init_ped` function, the probability of carrying 
the disease allele is defined by the disease-allele frequency (`DAF`) 
and determined via Bernoulli sampling. 

For spouses in the `add_gen` function, the probability of carrying the 
disease allele is also determined via Bernoulli sampling. When generating offspring
in the `add_gen` function, this disease allele is autosomal, so probability 
of transmission is 0.5, sampled from either parent's disease alleles. 
In the `add_inlaws` function, when sampling parents, one of their alleles 
is based on the already simulated offsprings disease alleles. 

If an individual carries the disease allele, we used the cumulative, lifetime 
penetrance of a disease to sample whether the individual developed disease. The
age at which an individual developed disease and handling of competing mortality 
is handled in the following paragraphs. 

#### Polygenic model

The general polygenic model is defined by $P = G + E$ where the
phenotypic value ($P$) is the sum of a genetic value ($G$) and
non-genetic value ($E$). The total phenotypic variance is standardized
to 1, with heritability ($h^2$) defining the proportion of variance explained
by genetic effects. The residual variance is defined as $1 - h^2$.

For founders, genetic values for correlated traits (ALS, FTD, and dementia)
are sampled jointly from a multivariate normal distribution: 
$\mathbf{G}_{founder} \sim N(\mathbf{0}, \mathbf{V}_g)$

where $\mathbf{V}_g$ is the covariance matrix constructed from trait 
heritabilities and their pairwise genetic correlations ($rg$). 

For non-founders, genetic values are simulated conditional on parental values 
($G_p, G_m$). The offspring mean is the mid-parent value, and the variance is 
$0.5 \times h^2$ per trait. Offspring values are drawn from a multivariate 
normal distribution to maintain specified trait correlations:
$\mathbf{G}_{offspring} \sim N\left(\frac{\mathbf{G}_p + \mathbf{G}_m}{2}, \mathbf{V}_{g, offspring}\right)$
where $\mathbf{V}_{g, offspring}$ uses the same correlation structure but 
scaled to reflect the reduction in variance due to segregation ($0.5 \times h^2$). 
Environmental components are subsequently assigned as $E \sim N(0, 1 - h^2)$.

Disease status is determined by the liability threshold model. For each trait,
the individual liability $P$ is compared to a lifetime-risk threshold ($t$), defined as:
$$
t = -\text{qnorm}(K)
$$
where $K$ is the population lifetime risk. If $P > t$, the individual is assigned 
an age at onset—if this occurs before their age at censoring, they are classified as affected. 


Limitations to polygenic model:
- Shared environment between relatives is not modeled
- There is no assortative mating
- Heritability is additive, there is no epistasis/dominance 

\## R code and functions.

``` r
source("src/libraries_simPed.R")
source("src/functions_simPed.R")
```

### Function to initiate pedigree with one founder with a mutation

``` r
print(init_ped)
```

init_ped = function(DAF_common, DAF_patho, DAF_ftd, k, mean_gen_yr, yob_index, life_expectancy, current_year){
    #   core_ped = data.frame(gen = 0, # generation
    #                         id = "C0", # individual id - founder
    #                         pid = as.character(NA), # parental id
    #                         mid = as.character(NA), # maternal id
    #                         sex = 1, # sex, 0 = male, 1 = female
    #                         a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
    #                         a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
    #                         a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
    #                         a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
    #                         a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency
    #                         a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd))) # sample disease allele from population frequency 
    #   # get the birth year of the founder:
    #   gen_yr = generate_gen_yr(mean_gen_yr)
    #   core_ped$yob = yob_index - k * gen_yr
    #   # simulate life expectancy
    #   core_ped$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == core_ped$yob]
    #   core_ped$life_expectancy <- sample_life_expectancy(core_ped$mean_life_exp)
    #   # simulate age
    #   core_ped$age = current_year - core_ped$yob
    #   if(core_ped$age > core_ped$life_expectancy){
    #     core_ped$status = "dead"
    #     core_ped$age_censored = core_ped$life_expectancy
    #   } else {
    #     core_ped$status = "alive"
    #     core_ped$age_censored = core_ped$age
    #   }
    #   return(core_ped)
    # }

### Function to simulate a next generation

``` r
print(add_gen)
```

add_gen = function(df_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, fert_rate, mean_gen_yr, life_expectancy, current_year){
    #   g = 0
    #   core_lineage <- "C0"
    #   core_offspring_created <- FALSE   
    #   while(g < k){
    #     # select individuals from youngest generation
    #     I1s = filter(df_ped, gen == g) # call I1 for parent 1
    #     # check if there are individuals in this generation, if not re-simulate generation
    #     if(nrow(I1s) > 0) {
    #       # for each individual
    #       for(i in 1:nrow(I1s)){
    #         gen_yr = generate_gen_yr(mean_gen_yr)
    #         # simulate partner (I2 for parent 2)
    #         I2 = data.frame(gen = g, 
    #                         id  = gsub("C", "P", I1s$id[i]),
    #                         pid = NA, # at this point parents are irrelevant, will be simulated in later step
    #                         mid = NA, # at this point parents are irrelevant, will be simulated in later step
    #                         sex = abs(I1s$sex[i] - 1), # opposite sex partners only...
    #                         a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
    #                         a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # sample disease allele from population frequency
    #                         a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
    #                         a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # sample disease allele from population frequency
    #                         a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency
    #                         a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # sample disease allele from population frequency 
    #                         yob = I1s$yob[i])
    #         I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
    #         I2$life_expectancy <- sample_life_expectancy(I2$mean_life_exp) 
    #         I2$age = current_year - I2$yob
    #         if(I2$age > I2$life_expectancy){
    #           I2$status = "dead"
    #           I2$age_censored = I2$life_expectancy
    #         } else {
    #           I2$status = "alive"
    #           I2$age_censored = I2$age
    #         }
    #         # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
    #         if(is.na(lambda)){
    #           n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
    #         } else { 
    #           n_II = rnbinom(1, size = 5, mu = lambda)
    #         }
    #         n_II = min(n_II, 15) # Cap at 15
    #         # Force offspring for at least one C*_1 individual; this will over-estimate pedigree size...
    #         while(n_II == 0 & (g == 0 | (grepl("^C0(_1)+$", I1s$id[i])))) { # (g < (k-1); because the re-simulating one generation back was only a problem in for the gens that were not (k-1))
    #             message(sprintf("FORCE: %s (gen=%d) had n_II=%d, forcing 1 offspring", 
    #                   I1s$id[i], g, n_II))      
    #             n_II = 1 # minimal number of offspring needed to continue core pedigree; resampling has risk of generating >1 offspring when original there was none 
    #         }
    #         # create data_frame for offspring:
    #         IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
    #         colnames(IIs) = colnames(df_ped)
    #         j = 0
    #         while(j < n_II){
    #           j = j+1
    #           IIs$gen[j] = g + 1
    #           IIs$id[j]  = paste(I1s$id[i], j, sep="_")
    #           IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
    #           IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
    #           IIs$sex[j] = sample(c(0,1), 1)
    #           IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1)
    #           IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
    #           IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1)
    #           IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)          
    #           IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1)
    #           IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1)          
    #           IIs$yob[j] = I1s$yob[i] + gen_yr
    #           IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
    #           IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
    #           IIs$age[j] = current_year - IIs$yob[j]
    #           if(IIs$age[j] > IIs$life_expectancy[j]){
    #             IIs$status[j] = "dead"
    #             IIs$age_censored[j] = IIs$life_expectancy[j]
    #           } else {
    #             IIs$status[j] = "alive"
    #             IIs$age_censored[j] = IIs$age[j]
    #           }
    #         }
    #         if(n_II > 0){
    #           df_ped = bind_rows(df_ped, I2, IIs)
    #         } else {
    #           df_ped = bind_rows(df_ped, I2)
    #         }
    #       }
    #       g = g+1
    #     } else {
    #       g = g-1 # go back one generation and re-simulate to make sure there is any offspring, this will over-estimate pedigree size...
    #     }
    #   }
    #   return(df_ped)
    # }

### Function to simulate the “inlaws”

These are the ancestors for those who married into this pedigree

``` r
print(add_inlaws)
```

    # add_inlaws = function(df_ped, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, life_expectancy, current_year){
    #   adj_ped = filter(df_ped, grepl("P", id))
    #   g = max(adj_ped$gen)
    #   while(! g == 0){
    #     IIs = filter(adj_ped, gen == g) # get individuals from offspring generation
    #     for(i in 1:nrow(IIs)){
    #       gen_yr = generate_gen_yr(mean_gen_yr) 
    #       # simulate parents 1 (I1)
    #       I1 = data.frame(gen = g-1, 
    #                       id  = paste0(IIs$id[i], "_m"), # "_m" suffix for mother
    #                       pid = NA, 
    #                       mid = NA, 
    #                       sex = 1, # mother
    #                       a1_common  = ifelse(IIs$a1_common[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
    #                       a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # non-transmitted allele sampled from population
    #                       a1_patho  = ifelse(IIs$a1_patho[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
    #                       a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # non-transmitted allele sampled from population 
    #                       a1_ftd  = ifelse(IIs$a1_ftd[i] == 1, 1, 0), # assume mother always transmits a1 for coding convenience
    #                       a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # non-transmitted allele sampled from population
    #                       yob = IIs$yob[i] - gen_yr) # define year of birth
    #       I1$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I1$yob]
    #       I1$life_expectancy = sample_life_expectancy(I1$mean_life_exp)      
    #       I1$age = current_year - I1$yob
    #       if(I1$age > I1$life_expectancy){
    #         I1$status = "dead"
    #         I1$age_censored = I1$life_expectancy
    #       } else {
    #         I1$status = "alive"
    #         I1$age_censored = I1$age
    #       }
    #       # simulate father
    #       I2 = data.frame(gen = g-1, 
    #                       id = paste0(IIs$id[i], "_p"), # "_p" suffix for father
    #                       pid = NA, 
    #                       mid = NA, 
    #                       sex = 0,
    #                       a1_common = ifelse(IIs$a2_common[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
    #                       a2_common = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)), # non-transmitted allele sampled from population
    #                       a1_patho = ifelse(IIs$a2_patho[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
    #                       a2_patho = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)), # non-transmitted allele sampled from population                      
    #                       a1_ftd = ifelse(IIs$a2_ftd[i] == 1, 1, 0), # assume father always transmits a2 for coding convenience
    #                       a2_ftd = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)), # non-transmitted allele sampled from population 
    #                       yob = IIs$yob[i] - gen_yr) # define year of birth
    #       I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
    #       I2$life_expectancy = sample_life_expectancy(I2$mean_life_exp) 
    #       I2$age = current_year - I2$yob
    #       if(I2$age > I2$life_expectancy){
    #         I2$status = "dead"
    #         I2$age_censored = I2$life_expectancy
    #       } else {
    #         I2$status = "alive"
    #         I2$age_censored = I2$age
    #       }                
    #       # add parents to pedigree
    #       adj_ped = bind_rows(adj_ped, I1, I2)
    #       # update parental ID in ped
    #       adj_ped[adj_ped$id == IIs$id[i],"mid"] = I1$id
    #       adj_ped[adj_ped$id == IIs$id[i],"pid"] = I2$id
    #     }
    #     g = g - 1
    #   }
    #   df_ped = filter(df_ped, ! grepl("P", id)) %>% 
    #     bind_rows(., adj_ped) %>% 
    #     arrange(., gen)
    #   return(df_ped)
    # }

### Function to add genetic values for Mendelian and polygenic inheritance

``` r
print(add_pheno)
```

add_ext_branches = function(df_ped, lambda, k, DAF_common, DAF_patho, DAF_ftd, mean_gen_yr, fert_rate, life_expectancy, current_year){
    #   g = 0
    #   while(g < k){
    #     # First, for pairs already in df_ped then select only mothers from the parental generation # sex is irrelevant in this simulation
    #     I1s = mutate(df_ped, id_parents = gsub("_[a-z]$", "", id)) %>% # create temporary variable ID minus _p or _m (so parents have same ID)
    #       group_by(id_parents) %>%
    #       filter(n() > 1) %>% # filter for parental couple
    #       ungroup() %>%
    #       filter(grepl("m$", id) & gen == g) # select only mothers in generation of offspring
    #     if(nrow(I1s > 0)){
    #       # for each mother
    #       for(i in 1:nrow(I1s)){
    #         gen_yr = generate_gen_yr(mean_gen_yr)
    #         # find father that is already in pedigree
    #         I2_id = gsub("_m$", "_p", I1s$id[i])
    #         I2 = filter(df_ped, id == I2_id)
    #         if(! nrow(I2 == 1)){
    #           print(I2)
    #           print(I2_id)
    #           stop("Non-unique IDs found!!")
    #         }
    #         # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
    #         if(is.na(lambda)){
    #           n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"]) - 1 # minus 1 because one child has already been simulated in first round
    #         } else { 
    #           n_II = rnbinom(1, size = 5, mu = lambda) - 1 # minus 1 because one child has already been simulated in first round
    #         }
    #         n_II = min(n_II, 15) # Cap at 15
    #         # simulate offspring simulate to add_gen() function
    #         if(n_II > 0){
    #           IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
    #           colnames(IIs) = colnames(df_ped)
    #           j = 0
    #           while(j < n_II){
    #             j = j+1
    #             IIs$gen[j] = g + 1
    #             IIs$id[j]  = paste(I1s$id[i], j, sep="_")
    #             IIs$pid[j] = I2_id
    #             IIs$mid[j] = I1s$id[i]
    #             IIs$sex[j] = sample(c(0,1), 1)
    #             IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1) # a1 is always from mother - here, I1
    #             IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
    #             IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1) # a1 is always from mother - here, I1
    #             IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)
    #             IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1) # a1 is always from mother - here, I1
    #             IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1)            
    #             IIs$yob[j] = I1s$yob[i] + gen_yr
    #             IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
    #             IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
    #             IIs$age[j] = current_year - IIs$yob[j]
    #             if(IIs$age[j] > IIs$life_expectancy[j]){
    #               IIs$status[j] = "dead"
    #               IIs$age_censored[j] = IIs$life_expectancy[j]
    #             } else {
    #               IIs$status[j] = "alive"
    #               IIs$age_censored[j] = IIs$age[j]
    #             }
    #           }
    #           df_ped = bind_rows(df_ped, IIs)
    #         }
    #       }
    #     }
        
    #     # Second, find NF who don't have a partner yet (the offspring simulated in previous loop)
    #     # these have IDs ending with a digit, but do have either "p" or "m" in their IDs
    #     I1s = filter(df_ped, gen == g & 
    #                    (grepl("p", id) | grepl("m", id)) &
    #                    ! (grepl("p$", id) | grepl("m$", id)))
    #     if(nrow(I1s > 0)){
    #       for(i in 1:nrow(I1s)){
    #         gen_yr = generate_gen_yr(mean_gen_yr)
    #         I2 = data.frame(gen = g, 
    #                         id = paste0("P",I1s$id[i]),
    #                         pid = NA, 
    #                         mid = NA, 
    #                         sex = abs(I1s$sex[i] - 1),
    #                         a1_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)),
    #                         a2_common  = sample(c(0,1), 1, prob=c(1-DAF_common, DAF_common)),
    #                         a1_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)),
    #                         a2_patho  = sample(c(0,1), 1, prob=c(1-DAF_patho, DAF_patho)),
    #                         a1_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)),
    #                         a2_ftd  = sample(c(0,1), 1, prob=c(1-DAF_ftd, DAF_ftd)),
    #                         yob = I1s$yob[i])
    #         I2$mean_life_exp = life_expectancy$life_expectancy[life_expectancy$year_of_birth == I2$yob]
    #         I2$life_expectancy = sample_life_expectancy(I2$mean_life_exp)        
    #         I2$age = current_year - I2$yob
    #         if(I2$age > I2$life_expectancy){
    #           I2$status = "dead"
    #           I2$age_censored = I2$life_expectancy
    #         } else {
    #           I2$status = "alive"
    #           I2$age_censored = I2$age
    #         }                
    #         # simulate number of offspring from negative binomial distribution with mean lambda (as defined by general pedigree parameters, or obtained from fertility rate and birthyear)
    #         if(is.na(lambda)){
    #           n_II = rnbinom(1, size = 5, mu = fert_rate[fert_rate$year == (I1s$yob[i] + gen_yr), "mean_fertility"])
    #         } else { 
    #           n_II = rnbinom(1, size = 5, mu = lambda)
    #         }
    #         n_II = min(n_II, 15) # Cap at 15
    #         # create data_frame for offspring:
    #         IIs = as.data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
    #         colnames(IIs) = colnames(df_ped)
    #         j = 0
    #         while(j < n_II){
    #           j = j+1
    #           IIs$gen[j] = g + 1
    #           IIs$id[j]  = paste(I1s$id[i], j, sep="_")
    #           IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
    #           IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
    #           IIs$sex[j] = sample(c(0,1), 1)
    #           IIs$a1_common[j]  = sample(c(I1s$a1_common[i], I1s$a2_common[i]), 1)
    #           IIs$a2_common[j]  = sample(c(I2$a1_common, I2$a2_common), 1)
    #           IIs$a1_patho[j]  = sample(c(I1s$a1_patho[i], I1s$a2_patho[i]), 1)
    #           IIs$a2_patho[j]  = sample(c(I2$a1_patho, I2$a2_patho), 1)
    #           IIs$a1_ftd[j]  = sample(c(I1s$a1_ftd[i], I1s$a2_ftd[i]), 1)
    #           IIs$a2_ftd[j]  = sample(c(I2$a1_ftd, I2$a2_ftd), 1) 
    #           IIs$yob[j] = I1s$yob[i] + gen_yr
    #           IIs$mean_life_exp[j] = life_expectancy$life_expectancy[life_expectancy$year_of_birth == IIs$yob[j]]
    #           IIs$life_expectancy[j] = sample_life_expectancy(IIs$mean_life_exp[j])
    #           IIs$age[j] = current_year - IIs$yob[j]
    #           if(IIs$age[j] > IIs$life_expectancy[j]){
    #             IIs$status[j] = "dead"
    #             IIs$age_censored[j] = IIs$life_expectancy[j]
    #           } else {
    #             IIs$status[j] = "alive"
    #             IIs$age_censored[j] = IIs$age[j]
    #           }
    #         }
    #         if(n_II > 0){
    #           df_ped = bind_rows(df_ped, I2, IIs)
    #         } else {
    #           df_ped = bind_rows(df_ped, I2)
    #         }
    #       }
    #     }
    #     g = g+1
    #   }
    #   return(df_ped)
    # }

wrapper function to simulate full pedigree

``` r
print(sim_ped)
```

sim_ped = function(i, k, lambda=NA, yob_index=1960, mean_gen_yr, fert_rate, disease_onset,
    #                    DAF_common, DAF_patho, DAF_ftd, penetrance_ALS_common, penetrance_ALS_patho, 
    #                    K_ALS, h2_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, K_FTD, h2_FTD, 
    #                    penetrance_dem_common, penetrance_dem_patho, K_dementia, h2_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem, 
    #                    life_expectancy, current_year, plot=FALSE){
    #   core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=yob_index, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
    #   core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
    #   core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
    #   core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, fert_rate=fert_rate, life_expectancy=life_expectancy, current_year=current_year)   
    #   core_ped = add_pheno(core_ped, disease_onset=disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem) 
    #   core_ped = add_1st_relatives(core_ped, k)
    #   core_ped = add_2nd_and_3rd_relatives(core_ped, k)
    #   core_ped = affected_pathways_ALS(core_ped)
    #   core_ped = affected_pathways_FTD(core_ped)
    #   core_ped = affected_pathways_dem(core_ped)

    #   if(plot){
    #     last_gen_affected = filter(core_ped, gen == k & (Y_ALS == 1 | Y_dementia ==1))
    #     if(nrow(last_gen_affected) > 0) {
    #     ped_plt = ped(id = core_ped$id, 
    #                   fid = core_ped$pid,
    #                   mid = core_ped$mid,
    #                   sex = core_ped$sex + 1,
    #                   isConnected = TRUE,
    #                   reorder = FALSE)
    #     carriers = filter(core_ped, a1_common + a2_common + a1_patho + a2_patho > 0)$id
    #     affected_ALS = filter(core_ped, Y_ALS == 1)$id
    #     affected_FTD = filter(core_ped, Y_FTD == 1)$id
    #     affected_dementia = filter(core_ped, Y_dementia_other == 1)$id
    #     affected_either = filter(core_ped, Y_ALS == 1 | Y_dementia == 1)$id    
    #     # color according to percentile from N(0,1)
    #     pal_ALS = colorRampPalette(c("white", "orange", "red"))(1000)
    #     pal_FTD = colorRampPalette(c("white", "skyblue", "blue"))(1000)
    #     pal_dementia = colorRampPalette(c("white", "lightgreen", "darkgreen"))(1000)    
    #     Pcolors_ALS = pal_ALS[ceiling(pnorm(core_ped$P_ALS)*1000)]
    #     Pcolors_FTD = pal_FTD[ceiling(pnorm(core_ped$P_FTD)*1000)]
    #     Pcolors_dementia = pal_dementia[ceiling(pnorm(core_ped$P_dementia)*1000)]    
      
    #     # make plot
    #     plot(ped_plt, title="Polygenic phenotype for ALS on liability scale", 
    #     carrier = carriers, fill=Pcolors_ALS, cex=0.8)
    #     plot(ped_plt, title="Polygenic phenotype for FTD on liability scale", 
    #     carrier = carriers, fill=Pcolors_FTD, cex=0.8)
    #     plot(ped_plt, title="Polygenic phenotype for other dementias on liability scale", 
    #     carrier = carriers, fill=Pcolors_dementia, cex=0.8)    
    #     plot(ped_plt, title="Binary ALS phenotype", 
    #     carrier = carriers, aff = affected_ALS, cex=0.8)
    #     plot(ped_plt, title="Binary FTD phenotype", 
    #     carrier = carriers, aff = affected_FTD, cex=0.8)
    #     plot(ped_plt, title="Binary dementia phenotype", 
    #     carrier = carriers, aff = affected_dementia, cex=0.8)
    #     plot(ped_plt, title="Binary phenotype for either ALS or dementia (incl FTD)", 
    #     carrier = carriers, aff = affected_either, cex=0.8)
    #   }}
    #   return(core_ped)
    # }

## step-wise simulation with plots

step 0. Define parameters:

``` r
# Demographic parameters
script_path <- this.path()
script_dir <- this.dir()
k = 2 # total of 3 generations
lambda = NA # mean number of offspring
fert_rate = read.table(file.path(script_dir, "../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt"), header=T)
mean_gen_yr = 30
life_expectancy = read.table(file.path(script_dir, "../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025
disease_onset = read.csv(file.path(script_dir, "../data/age_of_onset/age_of_onset.csv"), header = T)

# Mendelian parameters
DAF_common = 0.00075 # Disease allele frequency, closer to Douglas 2024 than Van Wijk 2024
DAF_patho = 0.00012 # Derived from Douglas 2024 + Gnomad (SOD1+FUS)
DAF_ftd = 0.00015 # Derived from Gnomad (GRN+MAPT)
penetrance_ALS_common = 0.21 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024 (24%), Gao 2025
penetrance_ALS_patho = 0.50 #read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Douglas 2024 (FUS 19%, SOD1 54%)
penetrance_FTD_common = 0.1  # derived from Gao 2025
penetrance_FTD_patho = 0.1 # copied from Gao 2025
penetrance_FTD_nonALS = 0.9 
penetrance_dem_common = 0.5 # derived from Gao 2025, bit lower
penetrance_dem_patho = 0.1 # guestimate, no data

# polygenic parameters:
h2_ALS = 0.45 # additive polygenic heritability, derived from vRheenen; Russel 2019; lowered for polygenic proportion
K_ALS = (0.0027 * 0.85) # derived from Ryan 2019, 1/375 = 0.0027. Corrected for 15% monogenic. 
h2_FTD = 0.45 # additive polygenic heritability, Dijkstra 2025
K_FTD = (0.00134 * 0.75) # derived from Coyle-Gilchrist 2016, 1 in 750 = 0.0013, corrected for 25% monogenic. 
K_dementia = 0.418 # read.table(file.path(script_dir, "/data/lifetime_disease_risk/lifetime_dementia_risk.txt" ), header=T) # derived from Fang 2025 
h2_dementia = 0.55 # Dijkstra 2025; lowered to compensate for prevalence of e.g. vascular dementia
rg_ALSFTD = 0.6 # vRheenen 2021
rg_ALSdem = 0.25 # vRheenen 2021, Wainberg 2023, Chen 2024
rg_FTDdem = 0.35 # vRheenen 2021, Chen 2024
```

Step 1. simulate the core pedigree (with monogenic disease)

``` r
core_ped = init_ped(DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, k=k, yob_index=2010, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
# run the loop to add generations of offspring
core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
ped_plt1 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
color = "white"

png("step1_core_pedigree.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt1,
  title="Step 1 - Core Pedigree",
  cex=1.1,
  carrier = carriers,
  fill=color
)
dev.off()

```

![](intro_simPed_files/figure-gfm/step1_core_pedigree-1.png)<!-- -->

step 2. Simulate the inlaws, ancestors of the spouses married into this
pedigree

``` r
core_ped = add_inlaws(core_ped, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, mean_gen_yr=mean_gen_yr, life_expectancy=life_expectancy, current_year=current_year)
ped_plt2 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
color = ifelse(ped_plt2$ID %in% ped_plt1$ID, "white", "lightgray")

plot(ped_plt2, title="Step 2 - simulated the in-laws", cex=1.5, carrier = carriers, fill=color)

png("step2_inlaws.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt2,
  title="Step 2 - In-laws",
  cex=1.1,
  carrier = carriers,
  fill=color
)
dev.off()
```

![](intro_simPed_files/figure-gfm/step2_inlaws-1.png)<!-- -->

``` r
# step 3. Simulate external branches, unlinked to founder
```

``` r
core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF_common=DAF_common, DAF_patho=DAF_patho, DAF_ftd=DAF_ftd, fert_rate=fert_rate, mean_gen_yr=mean_gen_yr, 
                            life_expectancy=life_expectancy, current_year=current_year)
ped_plt3 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id

color = ifelse(ped_plt3$ID %in% ped_plt1$ID, "white", ifelse(ped_plt3$ID %in% ped_plt2$ID, "lightgray", "darkgray"))

png("step3_external_branches.png", width=400, height=400)  # Adjust size and resolution as needed
plot(
  ped_plt3,
  title="Step 3 - External branches",
  cex=1.0,
  carrier = carriers,
  fill=color
)
dev.off()
```

![](intro_simPed_files/figure-gfm/step3_unlinked-1.png)<!-- -->

``` r
# step 4. Simulate phenotypes
```

``` r
core_ped = add_pheno(core_ped, disease_onset, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, penetrance_FTD_common, 
                     penetrance_FTD_patho, penetrance_FTD_nonALS, h2_FTD, K_FTD, penetrance_dem_common, penetrance_dem_patho, h2_dementia, K_dementia, rg_ALSFTD, rg_ALSdem, rg_FTDdem)

carriers_ALS_FTD = filter(core_ped, (a1_common + a2_common + a1_patho + a2_patho) > 0)$id
affected_ALS = filter(core_ped, Y_ALS == 1)$id
affected_FTD = filter(core_ped, Y_FTD == 1)$id
affected_ALS_FTD = filter(core_ped, Y_ALS == 1 | Y_FTD == 1)$id
affected_dementia = filter(core_ped, Y_dementia_other == 1)$id
affected_either = filter(core_ped, Y_ALS == 1 | Y_dementia == 1)$id
status = filter(core_ped, status == "dead")$id
# make plots
png("plot_ALS.png", width=800, height=600)  
plot(ped_plt3, title="Binary ALS phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_ALS, cex=1.0)
dev.off()
png("plot_FTD.png", width=800, height=600)  
plot(ped_plt3, title="Binary ALS+FTD phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_ALS_FTD, cex=1.0)
dev.off()
png("plot_dem.png", width=800, height=600) 
plot(ped_plt3, title="Binary ALS+dementia phenotype", 
     carrier = carriers_ALS_FTD, aff = affected_either, cex=1.0)
dev.off()
```

![](intro_simPed_files/figure-gfm/step4_pheno-2.png)<!-- -->
