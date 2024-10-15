library(RTMB)
library(LaMa)


# Simulating data ---------------------------------------------------------

set.seed(123)
n = 5200
z = rnorm(n) # white noise covariate

N = 2
omega = matrix(c(0,1,1,0),2,2)
mu = c(1, 5)
sigma = c(1, 1)

beta = matrix(c(log(5), log(8), 0.2, -0.2), nrow = 2)
Lambda = exp(cbind(1,z) %*% t(beta))

n_tilde = ceiling(n / (1 + mean(apply(Lambda, 1, mean)))*1.5) # crude estimate of distinct # of stays
s = rep(NA, n_tilde)
C = x = rep(NA, 3*n)
s[1] = 1
for(t in 2:n_tilde){
  s[t] = sample(1:N, size=1, prob = omega[s[t-1],])
}
times = rpois(1, lambda = Lambda[1, s[1]]) + 1
C[1:times] = rep(s[1], times)
x[1:times] = rnorm(times, mu[s[1]], sigma[s[1]])
currentInd = times
t = 2
while(currentInd <= n){
  times = rpois(1, lambda = Lambda[currentInd + 1, s[t]]) + 1
  C[currentInd + 1:times] = rep(s[t], times)
  x[currentInd + 1:times] = rnorm(times, mu[s[t]], sigma[s[t]])
  
  currentInd = currentInd + times
  t = t+1
}
x = x[1:n]
z = z[1:n]



# Functions ---------------------------------------------------------------

mllk = function(par){
  getAll(par, dat)
  
  ## parameter transformation
  sigma = exp(logsigma); REPORT(sigma); REPORT(mu)
  if(N > 2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(logitomega)),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  
  ## computing all dwell-time distributions (pmf and cdf)
  Lambda = exp(cbind(1,Z) %*% t(beta))
  dm = lapply(1:N, function(j) sapply(1:agsizes[j]-1, dpois, lambda = Lambda[,j]))
  
  # Gamma = tpm_ihsmm(omega, dm)
  # delta = stationary_sparse(Gamma[[1]])
  
  ## calculating all state-dependent densities
  allprobs = matrix(NA, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind,j] = dnorm(x[ind], mu[j], sigma[j])
  # allprobsL = t(apply(allprobs, 1, rep, times = agsizes))
  
  ## forward algo
  # - forward_sparse(delta, Gamma, allprobs[-(1:(max(agsizes)-1)),], agsizes = agsizes)
  - forward_ihsmm(dm, omega, allprobs)
}


# Fitting model -----------------------------------------------------------

par = list(mu = mu, logsigma = log(sigma), beta = beta)
dat = list(x = x, N = 2, Z = z, agsizes = c(25,30))

obj = MakeADFun(mllk, par)
opt = nlminb(obj$par, obj$fn, obj$gr)

opt$par

