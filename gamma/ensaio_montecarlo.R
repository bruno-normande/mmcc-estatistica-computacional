## Essaio Monte Carlo
#   Distribuição Gamma
#   author: Bruno Normande

library(lattice)
library(boot)
source("gamma_estimators.R")

R = 100
N = c(100, 1000, 10000, 100000)
K = c(1, 2, 3, 5, 9)
BOOT_S = 200

nrows = R*length(K)*length(N)
each_estim = R*length(K)*length(N)

estimatives = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))
estimatives2 = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))
estimatives3 = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))


i = 1
count_eq = 0
count_dif = 0
# para cada shape k
for(k in 1:length(K)){
  # para cada tamanho da amostragem
  for(n in 1:length(N)){
    print(sprintf("Shape: %d N: %d", K[k], N[n]))
    for(r in 1:R){
      dist = rgamma(N[n], K[k]) # rate = scale = 1 por default      
      k_hat = max_ver(dist, 1:N[n])
      b = boot(dist, max_ver, BOOT_S)
      k_til = 2*k_hat - mean(b$t[,1])
      
      k_hat2 = first_moment(dist, 1:N[n])
      b = boot(dist, first_moment, BOOT_S)
      k_til2 = 2*k_hat - mean(b$t[,1])
      
      k_hat3 = second_cent_moment(dist, 1:N[n])
      b = boot(dist, second_cent_moment, BOOT_S)
      k_til3 = 2*k_hat - mean(b$t[,1])
      
      estimatives[i,] = c(N[n],K[k], k_hat, k_til)
      estimatives2[i,] = c(N[n],K[k], k_hat2, k_til2)
      estimatives3[i,] = c(N[n],K[k], k_hat3, k_til3)
      
      i = i + 1
    }
  }
}