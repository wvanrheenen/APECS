source("src/libraries_simPed.R")
source("src/functions_simPed.R")

#' compare theory and simulations for polygenic trait parameters
sim_h2 = function(N=100, peds=as.numeric(), k=as.numeric(), lambda=as.numeric(), K=as.numeric(), h2=as.numeric(), small=T){
  obs_h2 = rep(NA, N)
  for(x in 1:N){
    cat("POPULATION", x, "\n")
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
  }
  return(obs_h2)
}

# heritabilities:
h2s = c(0.2, 0.4, 0.6, 0.8)

# by K
Ks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
sim_by_K = list(expand.grid(h2s, Ks))
n_sim = 25
n_param = nrow(sim_by_K[[1]])
names(sim_by_K[[1]]) <- c("h2", "K")
sim_by_K[[2]] = list()
for(i in 1:nrow(sim_by_K[[1]])){
  cat("simulation", i, "/", nrow(sim_by_K[[1]]), "\n")
  sim_by_K[[2]][[i]] = sim_h2(N=n_sim, peds=1000, k=2, lambda=2, K=sim_by_K[[1]]$K[i], h2=sim_by_K[[1]]$h2[i], small=T)
}
df_plt1 = as.data.frame(matrix(NA, ncol=3, nrow=n_sim*n_param))
colnames(df_plt1) = c("h2", "k", "sim_h2")
i = 0
for(j in 1:n_sim){
  for(k in 1:n_param){
    i = i+1 
    df_plt1[i,] = c(sim_by_K[[1]]$h2[[k]], sim_by_K[[1]]$K[[k]], sim_by_K[[2]][[k]][j])
  }
}
write.table(df_plt1, file="analyses/sim_by_K.txt", col.names=T, row.names=F, quote=F, sep="\t")
  
# by lambda
lambdas = c(1, 1.3, 2, 3)
sim_by_L = list(expand.grid(h2s, lambdas))
n_sim = 25
n_param = nrow(sim_by_L[[1]])
names(sim_by_L[[1]]) <- c("h2", "L")
sim_by_L[[2]] = list()
for(i in 1:nrow(sim_by_L[[1]])){
  cat("simulation", i, "/", nrow(sim_by_L[[1]]), "\n")
  sim_by_L[[2]][[i]] = sim_h2(N=n_sim, peds=1000, k=2, lambda=sim_by_L[[1]]$L[i], K=0.1, h2=sim_by_L[[1]]$h2[i], small=T)
}
df_plt2 = as.data.frame(matrix(NA, ncol=3, nrow=n_sim*n_param))
colnames(df_plt2) = c("h2", "k", "sim_h2")
i = 0
for(j in 1:n_sim){
  for(k in 1:n_param){
    i = i+1 
    df_plt2[i,] = c(sim_by_L[[1]]$h2[[k]], sim_by_L[[1]]$L[[k]], sim_by_L[[2]][[k]][j])
  }
}
write.table(df_plt2, file="analyses/sim_by_L.txt", col.names=T, row.names=F, quote=F, sep="\t")

# by k
ks = c(1, 2, 3, 4)
sim_by_k = list(expand.grid(h2s, ks))
n_sim = 25
n_param = nrow(sim_by_k[[1]])
names(sim_by_k[[1]]) <- c("h2", "k")
sim_by_k[[2]] = list()
for(i in 1:nrow(sim_by_k[[1]])){
  cat("simulation", i, "/", nrow(sim_by_k[[1]]), "\n")
  sim_by_k[[2]][[i]] = sim_h2(N=n_sim, peds=1000, k=sim_by_k[[1]]$k[i], lambda=2, K=0.1, h2=sim_by_k[[1]]$h2[i], small=F)
}
df_plt3 = as.data.frame(matrix(NA, ncol=3, nrow=n_sim*n_param))
colnames(df_plt3) = c("h2", "k", "sim_h2")
i = 0
for(j in 1:n_sim){
  for(k in 1:n_param){
    i = i+1 
    df_plt3[i,] = c(sim_by_k[[1]]$h2[[k]], sim_by_k[[1]]$k[[k]], sim_by_k[[2]][[k]][j])
  }
}
write.table(df_plt3, file="analyses/sim_by_k.txt", col.names=T, row.names=F, quote=F, sep="\t")

