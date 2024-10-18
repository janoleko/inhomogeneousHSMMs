## Likelihood functions for simulation study
# This file contains three different likelihood functions for the simulation study.
# The main difference to the likelihoods for the applications are that Poisson dwell times are used
# and the incorporation of RTMB::DataEval(), such that the observations can be changed without rebuilding the tape.

## periodically inhomogeneous HMM
nllpHMM = function(par){
  getAll(par, .GlobalEnv$dat)
  
  # data eval allows for changing the data without rebuilding the tape
  getStep = function(x) obs$step
  getAngle = function(x) obs$angle
  step = DataEval(getStep, rep(advector(1), 0))
  angle = DataEval(getAngle, rep(advector(1), 0))
  ###############################################
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  Gamma = tpm_g(Z, beta, ad = TRUE)
  delta = stationary_p(Gamma, t = tod[1], ad = TRUE)
  allprobs = matrix(NA, nrow = length(step), ncol = N)
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], 0, kappa = kappa[j])
  }
  -forward_g(delta, Gamma[,,tod], allprobs, tod, ad = TRUE)
}

## homogeneous HSMM
nllHSMM = function(par){
  getAll(par, .GlobalEnv$dat)
  
  # data eval allows for changing the data without rebuilding the tape
  getStep = function(x) obs$step
  getAngle = function(x) obs$angle
  step = DataEval(getStep, rep(advector(1), 0))
  angle = DataEval(getAngle, rep(advector(1), 0))
  ###############################################
  
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  lambda = exp(loglambda); REPORT(lambda) # dwell time mean
  
  # dwell-time pmfs and embedded tpm
  dm = lapply(1:N, function(i) dpois(1:agsizes[i]-1, lambda[i])); REPORT(dm)
  omega = if (N == 2) tpm_emb() else tpm_emb(logitomega) # case distinction: 2 states -> no parameters to vary
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], 0, kappa[j])
  }
  
  # forward algorithm for HSMMs
  -forward_hsmm(dm, omega, allprobs) # stationary start by default
}

## periodically inhomogeneous HSMM
nllpHSMM = function(par){
  getAll(par, .GlobalEnv$dat)
  
  # data eval allows for changing the data without rebuilding the tape
  getStep = function(x) obs$step
  getAngle = function(x) obs$angle
  step = DataEval(getStep, rep(advector(1), 0))
  angle = DataEval(getAngle, rep(advector(1), 0))
  ###############################################
  
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  Lambda = exp(Z %*% t(beta)); REPORT(Lambda) # dwell time means
  REPORT(beta)

  # dwell-time pmfs and embedded tpm
  dm = lapply(1:N, function(i) sapply(1:agsizes[i]-1, dpois, lambda = Lambda[,i])); REPORT(dm)
  omega = if (N == 2) tpm_emb() else tpm_emb(logitomega) # case distinction: 2 states -> no parameters to vary
  REPORT(omega)
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], 0, kappa[j])
  }
  
  # forward algorithm for periodically inhomogeneous HSMMs
  -forward_phsmm(dm, omega, allprobs, tod) # periodically stationary start by default
}
