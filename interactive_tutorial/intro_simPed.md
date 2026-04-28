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
$$ \mathbf{G}_{founder} \sim N(\mathbf{0}, \mathbf{V}_g) $$

where $\mathbf{V}_g$ is the covariance matrix constructed from trait 
heritabilities and their pairwise genetic correlations ($rg$). 

For non-founders, genetic values are simulated conditional on parental values 
($G_p, G_m$). The offspring mean is the mid-parent value, and the variance is 
$0.5 \times h^2$ per trait. Offspring values are drawn from a multivariate 
normal distribution to maintain specified trait correlations:
$$
\mathbf{G}_{offspring} \sim N\left(\frac{\mathbf{G}_p + \mathbf{G}_m}{2}, \mathbf{V}_{g, offspring}\right)
$$
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

    ## function (monogenic = TRUE) 
    ## {
    ##     core_ped = data.frame(gen = 0, id = "C0", pid = as.character(NA), 
    ##         mid = as.character(NA), sex = 1, a1 = 0)
    ##     if (monogenic) {
    ##         core_ped$a2 = 1
    ##     }
    ##     else {
    ##         core_ped$a2 = 0
    ##     }
    ##     return(core_ped)
    ## }

### Function to simulate a next generation

``` r
print(add_gen)
```

    ## function (df_ped, lambda, k, DAF) 
    ## {
    ##     g = 0
    ##     while (g < k) {
    ##         I1s = filter(df_ped, gen == g)
    ##         if (nrow(I1s) > 0) {
    ##             for (i in 1:nrow(I1s)) {
    ##                 I2 = data.frame(gen = g, id = gsub("C", "P", 
    ##                   I1s$id[i]), pid = NA, mid = NA, sex = abs(I1s$sex[i] - 
    ##                   1), a1 = sample(c(0, 1), 1, prob = c(1 - DAF, 
    ##                   DAF)), a2 = sample(c(0, 1), 1, prob = c(1 - 
    ##                   DAF, DAF)))
    ##                 n_II = rpois(1, lambda)
    ##                 while (n_II == 0 & g == 0) {
    ##                   n_II = rpois(1, lambda)
    ##                 }
    ##                 IIs = as.data.frame(matrix(NA, nrow = n_II, ncol = ncol(df_ped)))
    ##                 colnames(IIs) = colnames(df_ped)
    ##                 j = 0
    ##                 while (j < n_II) {
    ##                   j = j + 1
    ##                   IIs$gen[j] = g + 1
    ##                   IIs$id[j] = paste(I1s$id[i], j, sep = "_")
    ##                   IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], 
    ##                     I2$id)
    ##                   IIs$mid[j] = ifelse(I1s$sex[i] == 0, I2$id, 
    ##                     I1s$id[i])
    ##                   IIs$sex[j] = sample(c(0, 1), 1)
    ##                   IIs$a1[j] = sample(c(I1s$a1[i], I1s$a2[i]), 
    ##                     1)
    ##                   IIs$a2[j] = sample(c(I2$a1, I2$a2), 1)
    ##                 }
    ##                 if (n_II > 0) {
    ##                   df_ped = bind_rows(df_ped, I2, IIs)
    ##                 }
    ##             }
    ##             g = g + 1
    ##         }
    ##         else {
    ##             g = g - 1
    ##         }
    ##     }
    ##     return(df_ped)
    ## }

### Function to simulate the “inlaws”

These are the ancestors for those who married into this pedigree

``` r
print(add_inlaws)
```

    ## function (df_ped, DAF) 
    ## {
    ##     adj_ped = filter(df_ped, grepl("P", id))
    ##     g = max(adj_ped$gen)
    ##     while (!g == 0) {
    ##         IIs = filter(adj_ped, gen == g)
    ##         for (i in 1:nrow(IIs)) {
    ##             I1 = data.frame(gen = g - 1, id = paste0(IIs$id[i], 
    ##                 "_m"), pid = NA, mid = NA, sex = 1, a1 = ifelse(IIs$a1[i] == 
    ##                 1, 1, 0), a2 = sample(c(0, 1), 1, prob = c(1 - 
    ##                 DAF, DAF)))
    ##             I2 = data.frame(gen = g - 1, id = paste0(IIs$id[i], 
    ##                 "_p"), pid = NA, mid = NA, sex = 0, a1 = ifelse(IIs$a2[i] == 
    ##                 1, 1, 0), a2 = sample(c(0, 1), 1, prob = c(1 - 
    ##                 DAF, DAF)))
    ##             adj_ped = bind_rows(adj_ped, I1, I2)
    ##             adj_ped[adj_ped$id == IIs$id[i], "mid"] = I1$id
    ##             adj_ped[adj_ped$id == IIs$id[i], "pid"] = I2$id
    ##         }
    ##         g = g - 1
    ##     }
    ##     df_ped = filter(df_ped, !grepl("P", id)) %>% bind_rows(., 
    ##         adj_ped) %>% arrange(., gen)
    ##     return(df_ped)
    ## }

### Function to simulate all external branches of the pedigree

These are the branches with individuals unlinked to the core pedigree.

``` r
print(add_ext_branches)
```

    ## function (df_ped, lambda, k, DAF) 
    ## {
    ##     g = 0
    ##     while (g < k) {
    ##         I1s = mutate(df_ped, id_parents = gsub("_[a-z]$", "", 
    ##             id)) %>% group_by(id_parents) %>% filter(n() > 1) %>% 
    ##             ungroup() %>% filter(grepl("m$", id) & gen == g)
    ##         if (nrow(I1s > 0)) {
    ##             for (i in 1:nrow(I1s)) {
    ##                 I2_id = gsub("_m$", "_p", I1s$id[i])
    ##                 I2 = filter(df_ped, id == I2_id)
    ##                 n_II = rpois(1, lambda) - 1
    ##                 if (n_II > 0) {
    ##                   IIs = as.data.frame(matrix(NA, nrow = n_II, 
    ##                     ncol = ncol(df_ped)))
    ##                   colnames(IIs) = colnames(df_ped)
    ##                   j = 0
    ##                   while (j < n_II) {
    ##                     j = j + 1
    ##                     IIs$gen[j] = g + 1
    ##                     IIs$id[j] = paste(I1s$id[i], j, sep = "_")
    ##                     IIs$pid[j] = I2$id
    ##                     IIs$mid[j] = I1s$id[i]
    ##                     IIs$sex[j] = sample(c(0, 1), 1)
    ##                     IIs$a1[j] = sample(c(I1s$a1[i], I1s$a2[i]), 
    ##                       1)
    ##                     IIs$a2[j] = sample(c(I2$a1, I2$a2), 1)
    ##                   }
    ##                   df_ped = bind_rows(df_ped, IIs)
    ##                 }
    ##             }
    ##         }
    ##         I1s = filter(df_ped, gen == g & (grepl("p", id) | grepl("m", 
    ##             id)) & !(grepl("p$", id) | grepl("m$", id)))
    ##         if (nrow(I1s > 0)) {
    ##             for (i in 1:nrow(I1s)) {
    ##                 I2 = data.frame(gen = g, id = paste0("P", I1s$id[i]), 
    ##                   pid = NA, mid = NA, sex = abs(I1s$sex[i] - 
    ##                     1), a1 = sample(c(0, 1), 1, prob = c(1 - 
    ##                     DAF, DAF)), a2 = sample(c(0, 1), 1, prob = c(1 - 
    ##                     DAF, DAF)))
    ##                 n_II = rpois(1, lambda)
    ##                 IIs = as.data.frame(matrix(NA, nrow = n_II, ncol = ncol(df_ped)))
    ##                 colnames(IIs) = colnames(df_ped)
    ##                 j = 0
    ##                 while (j < n_II) {
    ##                   j = j + 1
    ##                   IIs$gen[j] = g + 1
    ##                   IIs$id[j] = paste(I1s$id[i], j, sep = "_")
    ##                   IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], 
    ##                     I2$id)
    ##                   IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], 
    ##                     I2$id)
    ##                   IIs$sex[j] = sample(c(0, 1), 1)
    ##                   IIs$a1[j] = sample(c(I1s$a1[i], I1s$a2[i]), 
    ##                     1)
    ##                   IIs$a2[j] = sample(c(I2$a1, I2$a2), 1)
    ##                 }
    ##                 if (n_II > 0) {
    ##                   df_ped = bind_rows(df_ped, I2, IIs)
    ##                 }
    ##             }
    ##         }
    ##         g = g + 1
    ##     }
    ##     return(df_ped)
    ## }

### Function to add genetic values for Mendelian and polygenic inheritance

``` r
print(add_pheno)
```

    ## function (df_ped, penetrance, h2, K) 
    ## {
    ##     df_ped = mutate(df_ped, mendel_y1 = ifelse(a1 == 1, sample(c(0, 
    ##         1), 1, prob = c(1 - penetrance, penetrance)), 0))
    ##     df_ped = mutate(df_ped, mendel_y2 = ifelse(a2 == 1, sample(c(0, 
    ##         1), 1, prob = c(1 - penetrance, penetrance)), 0))
    ##     df_ped = mutate(df_ped, mendelY = ifelse(mendel_y1 + mendel_y2 > 
    ##         0, 1, 0))
    ##     df_ped$G = NA
    ##     founders = which(is.na(df_ped$pid))
    ##     df_ped$G[founders] = rnorm(length(founders), mean = 0, sd = sqrt(h2))
    ##     for (i in 1:nrow(df_ped)) {
    ##         if (is.na(df_ped$G[i])) {
    ##             pid = df_ped$pid[i]
    ##             mid = df_ped$mid[i]
    ##             Gp = df_ped[which(df_ped$id == pid), ]$G
    ##             Gm = df_ped[which(df_ped$id == mid), ]$G
    ##             df_ped$G[i] = rnorm(1, mean = (Gp + Gm)/2, sd = sqrt(0.5 * 
    ##                 h2))
    ##         }
    ##     }
    ##     df_ped$E = rnorm(nrow(df_ped), 0, sqrt(1 - h2))
    ##     df_ped$P = df_ped$G + df_ped$E
    ##     LT = -qnorm(K, 0, 1)
    ##     df_ped$polygenicY = ifelse(df_ped$P > LT, 1, 0)
    ##     df_ped$Y = ifelse(df_ped$polygenicY + df_ped$mendelY > 0, 
    ##         1, 0)
    ##     return(df_ped)
    ## }

### Alternative function to add genetic values for Mendelian and polygenic inheritance

This function samples G from multivariate normal distribution. This adds
genetic values for polygenic inheritance FAST for smaller pedigrees

``` r
print(add_pheno_small)
```

    ## function (df_ped, penetrance, h2, K) 
    ## {
    ##     df_ped = mutate(df_ped, mendel_y1 = ifelse(a1 == 1, sample(c(0, 
    ##         1), 1, prob = c(1 - penetrance, penetrance)), 0))
    ##     df_ped = mutate(df_ped, mendel_y2 = ifelse(a2 == 1, sample(c(0, 
    ##         1), 1, prob = c(1 - penetrance, penetrance)), 0))
    ##     df_ped = mutate(df_ped, mendelY = ifelse(mendel_y1 + mendel_y2 > 
    ##         0, 1, 0))
    ##     km = kinship(ped(id = df_ped$id, fid = df_ped$pid, mid = df_ped$mid, 
    ##         sex = df_ped$sex + 1, isConnected = TRUE))
    ##     Vg = km * 2 * h2
    ##     N_ind = nrow(df_ped)
    ##     df_ped$G = rmvn(1, mu = rep(0, N_ind), sigma = Vg)[1, ]
    ##     df_ped$E = rnorm(nrow(df_ped), 0, sqrt(1 - h2))
    ##     df_ped$P = df_ped$G + df_ped$E
    ##     LT = -qnorm(K, 0, 1)
    ##     df_ped$polygenicY = ifelse(df_ped$P > LT, 1, 0)
    ##     df_ped$Y = ifelse(df_ped$polygenicY + df_ped$mendelY > 0, 
    ##         1, 0)
    ##     return(df_ped)
    ## }

wrapper function to simulate full pedigree

``` r
print(sim_ped)
```

    ## function (i = numeric, k = numeric(), lambda = numeric(), monogenic = TRUE, 
    ##     DAF = numeric(), penetrance = numeric(), K = numeric(), h2 = numeric(), 
    ##     small = TRUE, plot = TRUE) 
    ## {
    ##     core_ped = init_ped(monogenic = monogenic)
    ##     core_ped = add_gen(core_ped, lambda = lambda, k = k, DAF = DAF)
    ##     core_ped = add_inlaws(core_ped, DAF = DAF)
    ##     core_ped = add_ext_branches(core_ped, lambda = lambda, k = k, 
    ##         DAF = DAF)
    ##     if (small) {
    ##         core_ped = add_pheno_small(core_ped, penetrance, h2, 
    ##             K)
    ##     }
    ##     else {
    ##         core_ped = add_pheno(core_ped, penetrance, h2, K)
    ##     }
    ##     if (plot) {
    ##         ped_plt = ped(id = core_ped$id, fid = core_ped$pid, mid = core_ped$mid, 
    ##             sex = core_ped$sex + 1, isConnected = TRUE, reorder = FALSE)
    ##         carriers = filter(core_ped, a1 + a2 > 0)$id
    ##         affected = filter(core_ped, Y == 1)$id
    ##         pal = colorRampPalette(c("white", "orange", "red"))(1000)
    ##         Pcolors = pal[ceiling(pnorm(core_ped$P) * 1000)]
    ##         pdf(paste0("PED", i, "_by_polygenic_P.pdf"), w = 30, 
    ##             h = 12)
    ##         plot(ped_plt, title = "step 4 - simulated polygenic phenotype on liability scale", 
    ##             carrier = carriers, fill = Pcolors)
    ##         dev.off()
    ##         pdf(paste0("PED", i, "_by_phenotype.pdf"), w = 30, h = 12)
    ##         plot(ped_plt, title = "step 4 - colored by binary phenotype", 
    ##             carrier = carriers, aff = affected)
    ##         dev.off()
    ##     }
    ##     return(core_ped)
    ## }

## step-wise simulation with plots

step 0. Define parameters:

``` r
# pedigree parameters
k = 2 # total of 4 generations
lambda = 2 # mean number of offspring
# Mendelian parameters
DAF = 0.002 # disease allele frequency
penetrance = 0.9 # penetrance
# polygenic parameters:
h2 = 0.8 # additive polygenic heritability
K = 0.05 # life-time risk
```

Step 1. simulate the core pedigree (with monogenic disease)

``` r
core_ped = init_ped(monogenic=TRUE)
# run the loop to add generations of offspring
core_ped = add_gen(core_ped, lambda=lambda, k=k, DAF=DAF)
ped_plt1 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, a1 + a2 > 0)$id
color = "red"
plot(ped_plt1, title="step 1 - simulated core pedigree", cex=0.8, carrier = carriers, fill=color)
```

![](intro_simPed_files/figure-gfm/step1_core_pedigree-1.png)<!-- -->

step 2. Simulate the inlaws, ancestors of the spouses married into this
pedigree

``` r
core_ped = add_inlaws(core_ped, DAF=DAF)
ped_plt2 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, a1 + a2 > 0)$id
color = ifelse(ped_plt2$ID %in% ped_plt1$ID, "red", "orange")
plot(ped_plt2, title="step 2 - simulated the in-laws", cex=0.8, carrier = carriers, fill=color)
```

![](intro_simPed_files/figure-gfm/step2_inlaws-1.png)<!-- -->

``` r
# step 3. Simulate external branches, unlinked to founder
```

``` r
core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF)
ped_plt3 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, a1 + a2 > 0)$id
color = ifelse(ped_plt3$ID %in% ped_plt1$ID, "red", ifelse(ped_plt3$ID %in% ped_plt2$ID, "orange", "purple"))
plot(ped_plt3, title="step 3 - simulated external branches to pedigree", carrier = carriers, fill=color, cex=0.8)
```

![](intro_simPed_files/figure-gfm/step3_unlinked-1.png)<!-- -->

``` r
# step 4. Simulate phenotypes
```

``` r
core_ped = add_pheno_small(core_ped, penetrance, h2, K)
carriers = filter(core_ped, a1 + a2 > 0)$id
affected = filter(core_ped, Y == 1)$id
# define pallette with 100 colors from ramp
pal = colorRampPalette(c("white", "orange", "red"))(1000)
# color according to percentile from N(0,1)
Pcolors = pal[ceiling(pnorm(core_ped$P)*1000)]
# make plot
plot(ped_plt3, title="step 4 - simulated polygenic phenotype on liability scale", carrier = carriers, fill=Pcolors, cex=0.8)
```

![](intro_simPed_files/figure-gfm/step4_pheno-1.png)<!-- -->

``` r
plot(ped_plt3, title="step 4 - colored by binary phenotype", carrier = carriers, aff = affected, cex=0.8)
```

![](intro_simPed_files/figure-gfm/step4_pheno-2.png)<!-- -->
