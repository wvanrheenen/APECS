#' ---
#' title: "simPed: a simulation scheme for pedigrees with Mendelian and complex traits"
#' output: github_document
#' author: Wouter van Rheenen
#' date: "`r format(Sys.time(), '%d %B %Y')`"
#' ---
#'

#' ## General outline of simulation scheme
#' 
#' ### Simulating family members:
#' The pedigree consists of core pedigree (C) with one single founder this is 
#' initiated by function `init_ped()`. 
#' 
#' Then offspring generation will be simulated using the `add_gen()`
#' function, IDs start with "C" and number reflects order of offspring (C0_0 for oldest, C0_1 for second child)
#' in third generation offspring of oldest is denoted as C0_0_[0-9] and second 
#' C0_1_[0-9] etc. etc.. Therefore, all individuals in this lineage IDs can be 
#' traced to founder. The `add_gen()` function starts with simulating the spouse.
#' The number of offspring is samles from a Poisson distribution with mean
#' number of offspring $\lambda =$ `lambda`. The spouses married into this pedigree are denoted with 
#' "P*" and have same code as partner (P_0_0 for partner of firs child in 
#' generation 2). The number of generations to be added is defined by `k`, so for
#' `k=2` a three-generation pedigree is simulated.
#' 
#' When all generations of core pedigree are simulated, we simulate the in-laws 
#' e.g. the ancestors of the married in spouses (P*) using the `add_inlaws()`
#' function. 
#' 
#' To complete the pedigree, we simulate offspring of the inlaws using `add_ext_branches()`. 
#' These individuals are unlinked to the core pedigree, e.g. sibs, cousins, etc. 
#' of the spouses married in to this pedigree.
#' 
#' ### Genetics and phenotypes
#' #### Mendelian disease alleles
#' The founder of the core pedigree can be defined to carry a pathogenic mutation,
#' when `init_ped(monogenic=TRUE)` is used. For now, this disease allele is autosomal,
#' so probability of transmission is 0.5. For spouses married into this pedigree,
#' the probability of carrying the disease allele is defined by the disease-allele
#' frequency (`DAF`). Note, disease alleles introduced by spouses can be transmitted
#' in the core pedigree as well. This can be prevented by setting `DAF=0`.
#' 
#' The probability of developing the disease is defined by the disease-allele
#' penetrance `penetrance`. Therefore, the disease is modeled to be autosomal dominant with or without reduced penetrance.
#' 
#' #### Polygenic model
#' The general polygenic model is defined by $P = G + E$ where the phenotypic value ($P$)
#' is the sum of a genetic value ($G$) and non-genetic value ($E$). The phenotypic
#' value $P$ is standardized to zero mean and unit variance and thus follows $N(0,1)$.
#' The heritabiliy ($h^2$), defined as the proportion of phenotypic variance ($V_p$) explained by
#' additive genetic variance. The non-genetic component has variance $1 - h^2$ by definition.
#' 
#' For a founder $G$ can be drawn from $N(0,h^2)$. For offspring of parents $p$ and $m$
#' with genetic values of $G_p$ and $G_m$ respectively $G_{offspring} ~ N(\frac{G_p + G_m}{2},\frac{1}{2}h^2)$ 
#' of small pedigree G can be simulated using a multivariate normal distribution.
#' This requires looping through all offspring which can be slow. Alternatively,
#' for simulating a large number of small pedigrees up to four generations, G can
#' be simulating from a multivariate normal distribution (`MASS::mvrnorm()` or 
#' `mvnfast::rmvn()`) where the variance-covariance matrix is defined as
#' $2\textbf{K}h^2$ where $\textbf{K}$ is the kinship matrix of the pedigree `ped`
#' obtained through `ribd::kinship(ped)`. Once $G$ is simulated, $E$ is assigned
#' from $N(0,1-h^2)$. Once $P$ is simulated, disease status is defined by the
#' liability threshold model, where the threhold $t$ is defined such that $\Phi_p$, 
#' the area under the tail of the standard normal distribution from $t$ is 
#' the population life-time risk $K$. Note, $K$ is defined by the user. Subsequently,
#' when for one individual $P > t$ disease status is defined as affected.
#' 
#' Equation for heritability sanity check:
#' 
#' Limitations to polygenic model:
#' 
#' * Shared environment between relatives is not modeled
#' 
#' * There is no assortative mating
#' 
#' * Heritability is additive, there is no epistasis/dominance


#' ## R code and functions.
source("src/libraries_simPed.R")
source("src/functions_simPed.R")

#' ### Function to initiate pedigree with one founder with a mutation
print(init_ped)

#' ### Function to simulate a next generation
print(add_gen)

#' ### Function to simulate the "inlaws"
#' These are the ancestors for those who married into this pedigree
print(add_inlaws)

#' ### Function to simulate all external branches of the pedigree 
#' These are the branches with individuals unlinked to the core pedigree.
print(add_ext_branches)

#' ### Function to add genetic values for Mendelian and polygenic inheritance
print(add_pheno)

#' ### Alternative function to add genetic values for Mendelian and polygenic inheritance
#' This function samples G from multivariate normal distribution.
#' This adds genetic values for polygenic inheritance FAST for smaller pedigrees
print(add_pheno_small)

#' wrapper function to simulate full pedigree
print(sim_ped)

#' ## step-wise simulation with plots 
#' step 0. Define parameters:
# pedigree parameters
k = 2 # total of 4 generations
lambda = 2 # mean number of offspring
# Mendelian parameters
DAF = 0.002 # disease allele frequency
penetrance = 0.9 # penetrance
# polygenic parameters:
h2 = 0.8 # additive polygenic heritability
K = 0.05 # life-time risk

#' Step 1. simulate the core pedigree (with monogenic disease)
#+ step1_core_pedigree
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

#' step 2. Simulate the inlaws, ancestors of the spouses married into this pedigree
#+ step2_inlaws
core_ped = add_inlaws(core_ped, DAF=DAF)
ped_plt2 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, a1 + a2 > 0)$id
color = ifelse(ped_plt2$ID %in% ped_plt1$ID, "red", "orange")
plot(ped_plt2, title="step 2 - simulated the in-laws", cex=0.8, carrier = carriers, fill=color)

#' step 3. Simulate external branches, unlinked to founder
#+ step3_unlinked
core_ped = add_ext_branches(core_ped, lambda=lambda, k=k, DAF=DAF)
ped_plt3 = ped(id = core_ped$id,
               fid = core_ped$pid,
               mid = core_ped$mid,
               sex = core_ped$sex + 1,
               isConnected = TRUE)
carriers = filter(core_ped, a1 + a2 > 0)$id
color = ifelse(ped_plt3$ID %in% ped_plt1$ID, "red", ifelse(ped_plt3$ID %in% ped_plt2$ID, "orange", "purple"))
plot(ped_plt3, title="step 3 - simulated external branches to pedigree", carrier = carriers, fill=color, cex=0.8)

#' step 4. Simulate phenotypes
#+ step4_pheno
core_ped = add_pheno_small(core_ped, penetrance, h2, K)
carriers = filter(core_ped, a1 + a2 > 0)$id
affected = filter(core_ped, Y == 1)$id
# define pallette with 100 colors from ramp
pal = colorRampPalette(c("white", "orange", "red"))(1000)
# color according to percentile from N(0,1)
Pcolors = pal[ceiling(pnorm(core_ped$P)*1000)]
# make plot
plot(ped_plt3, title="step 4 - simulated polygenic phenotype on liability scale", carrier = carriers, fill=Pcolors, cex=0.8)
plot(ped_plt3, title="step 4 - colored by binary phenotype", carrier = carriers, aff = affected, cex=0.8)











