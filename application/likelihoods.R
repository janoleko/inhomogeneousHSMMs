## this file contains the four likelihood functions used in the muskox analysis.
# These are:
# - nllHMM: likelihood for homogeneous HMM
# - nllpHMM: likelihood for periodically inhomogeneous HMM
# - nllHSMM: likelihood for homogeneous HSMM
# - nllpHSMM: likelihood for periodically inhomogeneous HSMM

nllHMM = function(par){
  getAll(par, dat) # making parameters and data accessible without $
  
  # parameter transformations for unconstrained optimisation
  mu = exp(logmu); REPORT(mu) # reporting with RTMB for easy access later
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa); REPORT(mu.turn)
  
  # transition probability matrix and initial distribution
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], mu.turn[j], kappa[j])
  }
  
  # forward algorithm
  -forward(delta, Gamma, allprobs, ad = TRUE)
}

## likelihood for periodic HMM
nllpHMM = function(par){
  getAll(par, dat)
  
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa); REPORT(mu.turn)
  
  # transition probability matrix (array) and initial distribution
  Gamma_p = tpm_g(Z, beta, ad = TRUE); REPORT(Gamma_p) # 24 slices
  delta = stationary_p(Gamma_p, ad = TRUE, t = tod[1]) # periodically stationary distribution at t1
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], mu.turn[j], kappa[j])
  }
  
  # forward algorithm
  -forward_g(delta, Gamma_p[,,tod], allprobs, ad = TRUE)
}

## likelihood for homogeneous HSMM
nllHSMM = function(par){
  getAll(par, dat)
  
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa); REPORT(mu.turn)
  # negative binomial distribution parametrised in terms of mean mu and phi = 1/size
  mu_dwell = exp(logmu_dwell); REPORT(mu_dwell) # dwell time mean
  phi_dwell = exp(logphi_dwell); REPORT(phi_dwell) # dwell time dispersion
  size = 1 / phi_dwell
  
  # dwell-time pmfs and embedded tpm
  dm = lapply(1:N, function(i) dnbinom(1:agsizes[i]-1, size[i], size[i]/(size[i] + mu_dwell[i])))
  REPORT(dm)
  omega = if (N == 2) tpm_emb() else tpm_emb(logitomega) # case distinction: 2 states -> no parameters to vary
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], mu.turn[j], kappa[j])
  }
  
  # forward algorithm for HSMMs
  -forward_hsmm(dm, omega, allprobs) # stationary start by default
}

## likelihood for periodically inhomogeneous HSMM
nllpHSMM = function(par){
  getAll(par, dat)
  
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa); REPORT(mu.turn)
  
  # dwell-time pmfs and embedded tpm
  # negative binomial distribution parametrised in terms of mean mu and phi = 1/size
  Mu_dwell = exp(Z_Mu %*% t(betaMu)); REPORT(Mu_dwell) # dwell time means
  Phi_dwell = exp(Z_Phi %*% t(betaPhi)); REPORT(Phi_dwell) # dwell time dispersions
  Size = 1 / Phi_dwell
  
  # dwell-time pmfs and embedded tpm
  dm = lapply(1:N, function(i) sapply(1:agsizes[i]-1, dnbinom, Size[,i], Size[,i]/(Size[,i] + Mu_dwell[,i])))
  REPORT(dm)
  omega = if (N == 2) tpm_emb() else tpm_emb_g(Z_omega, beta_omega) # case distinction: 2 states -> no parameters to vary
  REPORT(omega)

  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])*
      dvm(angle[ind], mu.turn[j], kappa[j])
  }
  
  # forward algorithm for periodically inhomogeneous HSMMs
  -forward_phsmm(dm, omega, allprobs, tod) # periodically stationary start by default
}