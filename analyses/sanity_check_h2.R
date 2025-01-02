source("src/libraries_simPed.R")
source("src/functions_simPed.R")

#' argument parser
parser = ArgumentParser()

parser$add_argument("-h2", "--heritability", type="numeric", default=0.8,
                    dest="h2", help="true heritability [default %(default)s]")
parser$add_argument("-K", "--prevalence", type="numeric", default=0.1,
                    dest="K", help="Life-time risk [default %(default)s]")
parser$add_argument("-l", "--lambda", type="numeric", default=2,
                    dest="lambda", help="average number of offspring [default %(default)s]")
parser$add_argument("-k", "--generations", type="numeric", default=2,
                    dest="k", help="offspring generations, total number of generations = k+1 [default %(default)s]")
parser$add_argument("-N", "--n-sim", type="numeric", default=25,
                    dest="n_sim", help="number of simulations [default %(default)s]")
parser$add_argument("-p", "--n-ped", type="numeric", default=1000,
                    dest="n_ped", help="number of pedigrees per simulation [default %(default)s]")
parser$add_argument("-o", "--out", type="character", required=TRUE,
                    dest="out", help="name of the output file [REQUIRED]")
args = parser$parse_args()


#' compare theory and simulations for polygenic trait parameters
sim_h2 = function(N=as.numeric(), peds=as.numeric(), k=as.numeric(), lambda=as.numeric(), K=as.numeric(), h2=as.numeric(), small=F){
  obs_h2 = rep(NA, N)
  for(x in 1:N){
    polygenic_pedigrees = list()
    for(h in 1:peds){
      # cat("simulate pedigree", h, "\n")
      polygenic_pedigrees[[h]] = sim_ped(h, k=k, lambda=lambda, monogenic=FALSE, DAF=0, penetrance=0, K=K, h2=h2, small=small, plot=F)
    }
    population_N = 0
    affected_N = 0
    offspring_N = 0
    offspring_Y1 = 0
    for(h in 1:length(polygenic_pedigrees)){
      pedigree = polygenic_pedigrees[[h]]
      population_N = population_N + nrow(pedigree)
      affected = filter(pedigree, Y == 1)
      affected_N = affected_N + nrow(affected)
      for(affected_id in affected$id){
        offspring = filter(pedigree, mid == affected_id | pid == affected_id)
        offspring_N = offspring_N + nrow(offspring)
        offspring_Y1 = offspring_Y1 + sum(offspring$Y)
      }
    }

    LT  = -qnorm(K, 0, 1)
    Z  = dnorm(LT)
    i  = Z / K
    
    Kr = offspring_Y1 / offspring_N
    LTr = -qnorm(Kr, 0, 1)
    Zr = dnorm(LTr)
    ir = Zr / Kr

    aR = 0.5
    
    obs_h2[x] = (LT - LTr * sqrt(1-(1-LT/i)*(LT^2 - LTr^2))) / (aR*(i + (i-LT) * LTr^2))
    cat("SAMPLE", x, ":", obs_h2[x], "\n")
  }
  return(obs_h2)
}

results = sim_h2(N=args$n_sim, peds=args$n_ped, k=args$k, lambda=args$l, K=args$K, h2=args$h2)

write.table(results, args$out, col.names=F, row.names=F, quote=F, sep="\t")

# # heritabilities:
# h2s = c(0.2, 0.4, 0.6, 0.8)
# Ks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
# lambdas = c(1, 1.5, 2, 3)
# ks = c(1, 2, 3, 4)
