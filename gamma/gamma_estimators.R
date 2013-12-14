## Esse arquivo contém os estimadores para gamma

# Estimador pelo primeiro momento
# Primeiro momento = E[X] = mean(x)
# E[X] = shape/rate, rate = 1
# shape_hat = mean(x)
first_moment = function(x, d)
{
  return(mean(x[d]))
}

## Estimador pela segundo momento central
# var[X] = shape/rate^2
# rate = 1
# shape_hat = var(x)
seg_cent_moment = function(x, d)
{
  return(var(x[d]))
}

# temos que maximizar:
# log(p(D|a,ˆb)) = n*(a−1)*mean(log(x)) − n*log(Γ(a)) − n*a*log(mean(x))+ n*a*log(a)−n*a
# log(p(D|a,ˆb)) >= n*(a−1)*mean(log(x)) − n*log(Γ(a)) - n*a*log(mean(x))
#                               + n*(1 + log(a_0))*(a - a_0) + n*a_0*log(a_0) - n*a
# Usando "generalized Newton" fazer a seguinte aproximação
# log(p(D|a,ˆb)) ≈ c0 + c1*a + c2*log(a)
# 1/a = 1/a_0 + (mean(log(x)) - log(mean(x)) + log(a_0) - digamma(a_0)) / 
#                                                         a_0^2*(1/a_0 - trigamma(a_0))
max_ver = function(x, d){
  EPSILON = 0.00001
  dist_sample = x[d]
  # inicializando a
  log_meanx = log(mean(dist_sample))
  mean_logx = mean(log(dist_sample))
  a_hat = 0.5/(log_meanx - mean_logx)
  
  while(digamma(a_hat) != mean_logx - log_meanx + log(a_hat)){
    aux = mean_logx - log_meanx + log(a_hat) - digamma(a_hat)
    aux = aux/(a_hat^2 * (1/a_hat - trigamma(a_hat)))
    aux = aux + 1/a_hat
    aux = 1/aux
    if(abs(a_hat - aux) < EPSILON){
      break
    }
    a_hat = aux
  }
  
  return(a_hat)
}