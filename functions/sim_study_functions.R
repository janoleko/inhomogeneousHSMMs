

# Function for simulating data --------------------------------------------

sim_phsmm = function(n, beta, omega, stateparams, L=24, K=1, seed = 123){
  # set.seed(seed)
  tod = rep(1:L, n)
  
  mu = stateparams$mu
  sigma = stateparams$sigma
  kappa = stateparams$kappa
  
  N = dim(omega)[1]
  Z = LaMa::trigBasisExp(1:L, L, degree = K)
  Lambda = exp(cbind(1, Z)%*%t(beta))
  n_tilde = ceiling(n / (1+mean(apply(Lambda, 1, mean)))*1.5) # crude estimate of distinct # of stays
  
  s = rep(NA, n_tilde)
  C = step = angle = rep(NA, 3*n)
  
  s[1] = sample(1:N, size=1)
  for(t in 2:n_tilde){
    s[t] = sample(1:N, size=1, prob = omega[s[t-1],])
  }
  
  times = rpois(1, lambda = Lambda[tod[1],s[1]]) + 1
  C[1:times] = rep(s[1], times)
  step[1:times] = rgamma(times, shape = mu[s[1]]^2/sigma[s[1]]^2, scale = sigma[s[1]]^2/mu[s[1]])
  angle[1:times] = CircStats::rvm(times, 0, kappa[s[1]])
  
  currentInd = times
  
  for(t in 2:n_tilde){
    # times = rpois(1, lambda = Lambda[tod[currentInd],s[t]]) + 1
    times = rpois(1, lambda = Lambda[tod[currentInd + 1], s[t]]) + 1
    C[currentInd + 1:times] = rep(s[t], times)
    step[currentInd + 1:times] = rgamma(times, shape = mu[s[t]]^2/sigma[s[t]]^2, scale = sigma[s[t]]^2/mu[s[t]])
    angle[currentInd + 1:times] = CircStats::rvm(times, 0, kappa[s[t]])
    
    currentInd = currentInd + times
  }
  
  X = data.frame(step = step, angle = angle, C = C, tod = tod[1:length(step)])
  return(X[1:n,])
}


# beta = matrix(c(log(c(5,8,7)), -0.2, 0.2, -0.5, 0.3, -0.2, -0.2), nrow = 3)
# omega = matrix(c(0, 0.7, 0.3, 
#                  0.2, 0, 0.8,
#                  0.5, 0.5, 0), nrow = 3, byrow = TRUE)
# 
# color = c("orange", "deepskyblue", "seagreen2")
# Ztod = cbind(1, LaMa::trigBasisExp(seq(1,24,length=200), 24))
# dM = exp(Ztod%*%t(beta))
# plot(dM[,1], type = "l", ylim = c(0,15), col = color[1])
# lines(dM[,2], type = "l", col = color[2])
# lines(dM[,3], type = "l", col = color[3])
# 
# lambda = colMeans(dM)
# 
# stateparams = list(
#   mu = c(10, 100, 700),
#   sigma = c(10, 80, 500),
#   kappa = c(0.2, 1, 2.5)
# )
# 
# data = sim_phsmm(1000, beta, omega, stateparams)


# Functions for fitting models --------------------------------------------

# HMM
mllkHMM = function(theta.star, X, Z, ind, N=3, L=24, K=1){
  n = nrow(X)
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  kappa = exp(theta.star[2*N + 1:N])
  beta = matrix(theta.star[3*N + 1:(N*(N-1)*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
  Gamma = LaMa::tpm_p(tod=1:L, L=L, beta=beta, degree=K, Z=Z)
  delta = LaMa::stationary_p(Gamma, t = X$tod[1])
  allprobs = matrix(NA, nrow = n, ncol = N)
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  -LaMa::forward_p(delta, Gamma, allprobs, X$tod)
}

fitHMM = function(data, stateparams, N=3, L=24, K=1){
  Z = LaMa::trigBasisExp(1:L, L=L, degree=K)
  ind = which(!is.na(data$step) & !is.na(data$angle))
  
  thetainit = c(log(c(stateparams$mu,
                      stateparams$sigma,
                      stateparams$kappa)), rep(-2,N*(N-1)), rep(0, 2*N*(N-1)*K))
  
  t1 = Sys.time()
  mod = tryCatch(
    nlm(mllkHMM, thetainit, X = data, Z = Z, ind = ind, N=N, L=L, K=1,
        iterlim = 500, print.level = 0),
    error = function(e) e
  )
  est_time = Sys.time()
  
  if(!methods::is(mod, "error")){
    theta.star = mod$estimate
    mu = exp(theta.star[1:N])
    sigma = exp(theta.star[N+1:N])
    kappa = exp(theta.star[2*N + 1:N])
    beta = matrix(theta.star[3*N + 1:(N*(N-1)*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
    Gamma = LaMa::tpm_p(tod=1:L, L=L, beta=beta, degree=K, Z=Z)
    Delta = LaMa::stationary_p(Gamma)
    
    return(
      list(llk = -mod$minimum, time = est_time, beta=beta, Gamma=Gamma, Delta=Delta, mu=mu, sigma=sigma, kappa=kappa)
    )
  } else{
    return("No convergence")
  }
}


# HSMMs
## homogeneous HSMM

mllkHSMM = function(theta.star, X, ind, N=3, agsizes){
  n = nrow(X)
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  kappa = exp(theta.star[2*N + 1:N])
  lambda = exp(theta.star[3*N + 1:N])
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  dm = list()
  for(j in 1:N){ dm[[j]] = dpois(1:agsizes[j]-1, lambda[j]) }
  Gamma = LaMa::tpm_hsmm(omega, dm)
  delta = LaMa::stationary(Gamma)
  
  allprobs = matrix(NA, nrow = n, ncol = N)
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  -LaMa::forward_s(delta, Gamma, allprobs, agsizes)
}

fitHSMM = function(data, lambda, stateparams, agsizes=rep(20,N), N=3, L=24, K=1){
  ind = which(!is.na(data$step) & !is.na(data$angle))
  
  thetainit = c(log(c(stateparams$mu,
                      stateparams$sigma,
                      stateparams$kappa,
                      lambda)), 
                rep(0,N*(N-2))
  )
  
  t1 = Sys.time()
  mod = tryCatch(
    nlm(mllkHSMM, thetainit, X=data, ind=ind, N=N, agsizes=agsizes,
        iterlim = 500, print.level = 0),
    error = function(e) e
  )
  est_time = Sys.time()-t1
  
  if(!methods::is(mod, "error")){
    theta.star = mod$estimate
    mu = exp(theta.star[1:N])
    sigma = exp(theta.star[N+1:N])
    kappa = exp(theta.star[2*N + 1:N])
    lambda = exp(theta.star[3*N + 1:N])
    if(N>2){ # only needed if N>2
      omega = matrix(0,N,N)
      omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[4*N+1:(N*(N-2))])),N,N-1)))
      omega = t(omega)/apply(omega,2,sum)
    }else{ omega = matrix(c(0,1,1,0),2,2) }
    for(j in 1:N){ dm[[j]] = dpois(1:agsizes[j]-1, lambda[j]) }
    Gamma = LaMa::tpm_hsmm(omega, dm)
    delta = LaMa::stationary(Gamma)
    
    return(
      list(llk = -mod$minimum, time = est_time, lambda=lambda, omega=omega, mu=mu, sigma=sigma, kappa=kappa, Gamma=Gamma, delta=delta)
    )
  } else{
    return("No convergence")
  }
}


## inhomogeneous HSMM
mllkpHSMM = function(theta.star, X, Z, ind, N=3, L=24, K=1, agsizes){
  n = nrow(X)
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  kappa = exp(theta.star[2*N + 1:N])
  beta = matrix(theta.star[3*N + 1:(N*(1+2*K))], nrow = N, ncol = 1+2*K)
  Lambda = exp(cbind(1,Z)%*%t(beta))
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[3*N+N*(1+2*K)+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  dm = list()
  for(j in 1:N){ 
    dm[[j]] = sapply(1:agsizes[j]-1, dpois, lambda = Lambda[,j])
  }
  Gamma = LaMa::tpm_phsmm(omega, dm)
  delta = LaMa::stationary_p(Gamma, t = X$tod[1])
  
  allprobs = matrix(NA, nrow = n, ncol = N)
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta, Gamma, allprobs_large, X$tod)
}

fitpHSMM = function(data, beta, stateparams, agsizes=rep(30,N), N=3, L=24, K=1, stepmax = 20){
  Z = LaMa::trigBasisExp(1:L-1, L=L, degree=K)
  # indexshift is because in LaMa gamma^(t) = Pr(S_t | S_t-1) while in this paper gamma^(t) = Pr(S_t+1 | S_t)
  ind = which(!is.na(data$step) & !is.na(data$angle))
  
  thetainit = c(log(c(stateparams$mu,
                      stateparams$sigma,
                      stateparams$kappa)),
                as.numeric(beta), 
                rep(0,N*(N-2))
  )
  
  t1 = Sys.time()
  mod = tryCatch(
    nlm(mllkpHSMM, thetainit, X=data, Z=Z, ind=ind, N=N, L=L, K=K, agsizes=agsizes,
        iterlim = 300, print.level = 2, stepmax=stepmax),
    error = function(e) e
  )
  est_time = Sys.time()-t1
  
  if(!methods::is(mod, "error")){
    theta.star = mod$estimate
    mu = exp(theta.star[1:N])
    sigma = exp(theta.star[N+1:N])
    kappa = exp(theta.star[2*N + 1:N])
    beta = matrix(theta.star[3*N + 1:(N*(1+2*K))], nrow = N, ncol = 1+2*K)
    Lambda = exp(cbind(1,Z)%*%t(beta))
    if(N>2){ # only needed if N>2
      omega = matrix(0,N,N)
      omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[3*N+N*(1+2*K)+1:(N*(N-2))])),N,N-1)))
      omega = t(omega)/apply(omega,2,sum)
    }else{ omega = matrix(c(0,1,1,0),2,2) }
    dm = list()
    for(j in 1:N){ 
      dm[[j]] = sapply(1:agsizes[j]-1, dpois, lambda = Lambda[,j])
    }
    Gamma = LaMa::tpm_phsmm(omega, dm)
    delta = LaMa::stationary_p(Gamma, t = data$tod[1])
    
    return(
      list(llk = -mod$minimum, time = est_time, beta=beta, Lambda=Lambda, 
           omega=omega, mu=mu, sigma=sigma, kappa=kappa, 
           Gamma=Gamma, delta=delta, theta.star = theta.star)
    )
  } else{
    return("No convergence")
  }
}

