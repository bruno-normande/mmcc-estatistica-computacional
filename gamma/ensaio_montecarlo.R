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
nrows_bias_eqm = length(N)*length(K)

estimatives = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))
bias_eqm_1 = data.frame(N = rep(0, nrows_bias_eqm),
                        K = rep(0, nrows_bias_eqm),
                        Bias = rep(0, nrows_bias_eqm),
                        EQM = rep(0, nrows_bias_eqm),
                        Bias_til = rep(0, nrows_bias_eqm),
                        EQM_til = rep(0, nrows_bias_eqm))

estimatives2 = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))
bias_eqm_2 = data.frame(N = rep(0, nrows_bias_eqm),
                        K = rep(0, nrows_bias_eqm),
                        Bias = rep(0, nrows_bias_eqm),
                        EQM = rep(0, nrows_bias_eqm),
                        Bias_til = rep(0, nrows_bias_eqm),
                        EQM_til = rep(0, nrows_bias_eqm))

estimatives3 = data.frame(N = rep(0,nrows),
                         K = rep(0,nrows),
                         k_hat = rep(0, nrows),
                         k_til = rep(0, nrows))
bias_eqm_3 = data.frame(N = rep(0, nrows_bias_eqm),
                        K = rep(0, nrows_bias_eqm),
                        Bias = rep(0, nrows_bias_eqm),
                        EQM = rep(0, nrows_bias_eqm),
                        Bias_til = rep(0, nrows_bias_eqm),
                        EQM_til = rep(0, nrows_bias_eqm))

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

get_bias_eqm = function(est, k)
{
  bias = mean(est) - k
  EQM = bias^2 + var(est)
  return(c(bias,EQM))
}

counter = 1
# Retirando Viés e EQM
# para cada shape k
for(k in 1:length(K))
{
  # para cada tamanho da amostragem
  the_k = K[k]
  for(n in 1:length(N))
  {
    the_n = N[n]
    print(sprintf("%d -- %d", the_k, the_n))
    # para k_hat
    est = estimatives[estimatives$N == the_n & estimatives$K == the_k,]
    bias_eqm = get_bias_eqm(est$k_hat, the_k)
    bias_eqm_til = get_bias_eqm(est$k_til, the_k)
    bias_eqm_1[counter,] = c(the_n, the_k, bias_eqm, bias_eqm_til)
    
    # para k_hat2
    est = estimatives2[estimatives2$N == the_n & estimatives2$K == the_k,]
    bias_eqm = get_bias_eqm(est$k_hat, the_k)
    bias_eqm_til = get_bias_eqm(est$k_til, the_k)
    bias_eqm_2[counter,] = c(the_n, the_k, bias_eqm, bias_eqm_til)
    
    # para k_hat3
    est = estimatives3[estimatives3$N == the_n & estimatives3$K == the_k,]
    bias_eqm = get_bias_eqm(est$k_hat, the_k)
    bias_eqm_til = get_bias_eqm(est$k_til, the_k)
    bias_eqm_3[counter,] = c(the_n, the_k, bias_eqm, bias_eqm_til)
    
    counter = counter + 1
  }
}

get_row = function(x)
{
  return(sprintf("%d & %d & %s & %s\\", x[1], x[2], 
                 (abs(x[3]) > abs(x[5])), 
                  (x[4] > x[6])))
}

print("Maxima Verossimilhança")
table_rows = apply(bias_eqm_1,1,get_row)
for(r in table_rows){print(r)}

print("Primeiro Momento")
table_rows = apply(bias_eqm_2,1,get_row)
for(r in table_rows){print(r)}

print("Segundo Momento Central")
table_rows = apply(bias_eqm_3,1,get_row)
for(r in table_rows){print(r)}
