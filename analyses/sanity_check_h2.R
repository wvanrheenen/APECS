## compare theory and simulations for polygenic trait parameters
# sim_h2 = function(N=100, peds=as.numeric(), k=as.numeric(), lambda=as.numeric(), K=as.numeric(), h2=as.numeric(), small=T){
#   obs_h2 = rep(NA, N)
#   for(i in 1:N){
#     cat("POPULATION", i, "\n")
#     polygenic_pedigrees = list()
#     for(h in 1:peds){
#       # cat("simulate pedigree", h, "\n")
#       polygenic_pedigrees[[h]] = sim_ped(h, k=k, lambda=lambda, monogenic=FALSE, DAF=0, penetrance=0, K=K, h2=h2, small=small, plot=F)
#     }
#     population_N = 0
#     affected_N = 0
#     offspring_N = 0
#     offspring_Y1 = 0
#     for(h in 1:length(polygenic_pedigrees)){
#       pedigree = polygenic_pedigrees[[h]]
#       population_N = population_N + nrow(pedigree)
#       affected = filter(pedigree, Y == 1)
#       affected_N = affected_N + nrow(affected)
#       for(affected_id in affected$id){
#         offspring = filter(pedigree, mid == affected_id | pid == affected_id)
#         offspring_N = offspring_N + nrow(offspring)
#         offspring_Y1 = offspring_Y1 + sum(offspring$Y)
#       }
#     }
#     # obs_K  = affected_N / population_N
#     obs_K  = K
#     obs_Kr = offspring_Y1 / offspring_N
#     obs_T  = -qnorm(obs_K, 0, 1)
#     obs_Tr = -qnorm(obs_Kr, 0, 1)
#     obs_Z  = dnorm(obs_T)
#     obs_Zr = dnorm(obs_Tr)
#     obs_i  = obs_Z / obs_K
#     obs_ir = obs_Zr / obs_Kr
#     aR = 0.5
#     obs_h2[i] = (obs_T - obs_Tr * sqrt(1-(1-obs_T/obs_i)*(obs_T^2 - obs_Tr^2))) / (aR*(obs_i + (obs_i-obs_T) * obs_Tr^2))
#   }
#   obs_h2_hat = mean(obs_h2)
#   obs_h2_se = sd(obs_h2) / sqrt(N)
#   return(list(obs_h2, obs_K, obs_Kr, obs_h2_hat, obs_h2_se))
# }
# t1 = Sys.time()
# s1 = sim_h2(N=100, peds=100, k=2, lambda=2, K=0.15, h2=0.7, small=F)
# t2 = Sys.time()
# s2 = sim_h2(N=100, peds=100, k=2, lambda=2, K=0.15, h2=0.7, small=T)
# t3 = Sys.time()
# print(t2 - t1)
# print(t3 - t2)
