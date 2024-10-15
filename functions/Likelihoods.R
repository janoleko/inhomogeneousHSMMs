
# this file contains the negative log-likelihood functions for the HMMs and HSMMs
# as well as the auxiliary functions for parameter transformation


# HMM likelihoods ---------------------------------------------------------

### homogeneous
partransform = function(theta.star, N=3){
  p = list()
  p$Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
  p$delta = LaMa::stationary(p$Gamma)
  p$mu = exp(theta.star[N*(N-1)+1:N])
  p$sigma = exp(theta.star[N*(N-1)+N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(N-1)+2*N+1:N])
  return(p)
}
mllk = function(theta.star, X, N=3, ptransform){
  p = ptransform(theta.star, N)
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  -LaMa::forward(p$delta, p$Gamma, allprobs)
}

### inhomogeneous
partransform_p = function(theta.star, N=3, L=24, K, t1){
  p = list()
  p$beta = matrix(theta.star[1:(N*(N-1)*(1+2*K))], nrow=N*(N-1), ncol=1+2*K)
  p$Gamma = LaMa::tpm_p(tod=1:L, L=L, beta=p$beta, degree=K)
  p$delta = LaMa::stationary_p(p$Gamma, t1)
  p$mu = exp(theta.star[(1+2*K)*N*(N-1)+1:N])
  p$sigma = exp(theta.star[(1+2*K)*N*(N-1)+N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[(1+2*K)*N*(N-1)+2*N+1:N])
  return(p)
}
mllk_p = function(theta.star, X, N=3, L=24, K=1, ptransform_p){
  p = ptransform_p(theta.star, N, L, K, X$tod[1])
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  -LaMa::forward_p(p$delta, p$Gamma, allprobs, X$tod)
}


# HSMM likelihoods --------------------------------------------------------

### homogeneous
partransform_s = function(theta.star, N=3, agsizes){
  p = list()
  # parameter transform: dwell-time distributions
  p$mu_dwell = exp(theta.star[1:N])
  p$phi_dwell = exp(theta.star[N+1:N])
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[2*N+1:N])
  p$sigma = exp(theta.star[3*N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[4*N+1:N])
  # parametr transform: omega
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[5*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  p$omega = omega
  # dwell-time distributions and approximating tpm
  dm = list()
  for(j in 1:N){
    dm[[j]] = dnbinom(1:agsizes[j]-1, mu = p$mu_dwell[j], size = 1/p$phi_dwell[j])
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_hsmm(omega, dm)
  p$delta = LaMa::stationary(p$Gamma)
  return(p)
}
mllk_s = function(theta.star, X, N=3, agsizes){
  p = partransform_s(theta.star, N, agsizes)
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  -LaMa::forward_s(p$delta, p$Gamma, allprobs, agsizes)
}

### inhomogeneous
partransform_sp = function(theta.star, N=3, L=24, K, agsizes, t1){
  p = list()
  # parameter transform: dwell-time distributions
  p$beta_mu = matrix(theta.star[1:(N*(1+2*K))], nrow = N)
  p$phi_dwell = exp(theta.star[N*(1+2*K)+1:N])
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[N*(1+2*K)+N+1:N])
  p$sigma = exp(theta.star[N*(1+2*K)+2*N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(1+2*K)+3*N+1:N])
  # parametr transform: omega
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
                        exp(theta.star[N*(1+2*K)+4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  p$omega = omega
  # dwell-time distributions and approximating tpm
  Z = LaMa::trigBasisExp(tod=1:L, L=L, degree = K)
  p$Mu_dwell = exp(cbind(1,Z)%*%t(p$beta_mu))
  dm = list()
  for(j in 1:N){
    dm[[j]] = matrix(nrow=L, ncol = agsizes[j])
    for(t in 1:L){
      dm[[j]][t,] = dnbinom(1:agsizes[j]-1, mu = p$Mu_dwell[t,j], size = 1/p$phi_dwell[j])
    }
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_phsmm(omega, dm)
  p$delta = LaMa::stationary_p(p$Gamma, t = t1)
  return(p)
}
mllk_sp = function(theta.star, X, N=3, L=24, K, agsizes){
  p = partransform_sp(theta.star=theta.star, N=N, L=L, K=K, agsizes=agsizes, t1 = X$tod[1])
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta=p$delta, Gamma=p$Gamma, allprobs=allprobs_large, tod=X$tod)
}


## inhomogeneous also in dispersion parameter of dwell-time distributions
partransform_sp2 = function(theta.star, N=3, L=24, K=c(1,1), agsizes, t1){
  p = list()
  # parameter transform: dwell-time distributions
  p$beta_Mu = matrix(theta.star[1:(N*(1+2*K[1]))], nrow = N)
  p$beta_Phi = matrix(theta.star[N*(1+2*K[1])+1:(N*(1+2*K[2]))], nrow = N)
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[N*(2+2*(K[1]+K[2]))+1:N])
  p$sigma = exp(theta.star[N*(2+2*(K[1]+K[2]))+N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(2+2*(K[1]+K[2]))+2*N+1:N])
  # parametr transform: omega
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
              exp(theta.star[N*(2+2*(K[1]+K[2]))+3*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  p$omega = omega
  # dwell-time distributions and approximating tpm
  
  Z1 = LaMa::trigBasisExp(tod=1:L, L=L, degree=K[1])
  Z2 = LaMa::trigBasisExp(tod=1:L, L=L, degree=K[2])
  p$Mu_dwell = exp(cbind(1,Z1)%*%t(p$beta_Mu))
  p$Phi_dwell = exp(cbind(1,Z2)%*%t(p$beta_Phi))
  dm = list()
  for(j in 1:N){
    dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = p$Mu_dwell[,j], 
                     size = 1/p$Phi_dwell[,j])
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_phsmm(omega, dm)
  p$delta = LaMa::stationary_p(p$Gamma, t = t1)
  return(p)
}
mllk_sp2 = function(theta.star, X, N=3, L=24, K=c(1,1), agsizes){
  p = partransform_sp2(theta.star, N=N, L=L, K=K, agsizes=agsizes, t1=X$tod[1])
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(p$delta, p$Gamma, allprobs_large, X$tod)
}

### inhomogeneous in dwell time mean and omega
partransform_sp3 = function(theta.star, N=3, L=24, K, agsizes, t1){
  p = list()
  # parameter transform: dwell-time distributions
  p$beta_mu = matrix(theta.star[1:(N*(1+2*K))], nrow = N)
  p$phi_dwell = exp(theta.star[N*(1+2*K)+1:N])
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[N*(1+2*K)+N+1:N])
  p$sigma = exp(theta.star[N*(1+2*K)+2*N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(1+2*K)+3*N+1:N])
  p$beta_omega = matrix(theta.star[N*(1+2*K)+4*N+1:(N*(N-2)*(1+2*K))], nrow = N)
  Z = LaMa::trigBasisExp(tod=1:L, L=L, degree = K)
  p$expEta_omega = exp(cbind(1, Z)%*%t(p$beta_omega))
  # parameter transform: inhomogeneous omega
  if(N>2){ # only needed if N>2
    p$omega = array(dim = c(N,N,L))
    for(t in 1:L){
      omega = matrix(0,N,N)
      omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),p$expEta_omega[t,]),N,N-1)))
      omega = t(omega)/apply(omega,2,sum)
      p$omega[,,t] = omega
    }
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  # dwell-time distributions and approximating tpm
  p$Mu_dwell = exp(cbind(1,Z)%*%t(p$beta_mu))
  dm = list()
  for(j in 1:N){
    dm[[j]] = matrix(nrow=L, ncol = agsizes[j])
    for(t in 1:L){
      dm[[j]][t,] = dnbinom(1:agsizes[j]-1, mu = p$Mu_dwell[t,j], size = 1/p$phi_dwell[j])
    }
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_phsmm(p$omega, dm)
  p$delta = LaMa::stationary_p(p$Gamma, t = t1)
  return(p)
}
mllk_sp3 = function(theta.star, X, N=3, L=24, K, agsizes){
  p = partransform_sp3(theta.star=theta.star, N=N, L=L, K=K, agsizes=agsizes, t1 = X$tod[1])
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta=p$delta, Gamma=p$Gamma, allprobs=allprobs_large, tod=X$tod)
}
### inhomogeneous in dwell time mean and omega
partransform_sp3 = function(theta.star, N=3, L=24, K, agsizes, t1){
  p = list()
  # parameter transform: dwell-time distributions
  p$beta_mu = matrix(theta.star[1:(N*(1+2*K))], nrow = N)
  p$phi_dwell = exp(theta.star[N*(1+2*K)+1:N])
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[N*(1+2*K)+N+1:N])
  p$sigma = exp(theta.star[N*(1+2*K)+2*N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(1+2*K)+3*N+1:N])
  p$beta_omega = matrix(theta.star[N*(1+2*K)+4*N+1:(N*(N-2)*(1+2*K))], nrow = N)
  Z = LaMa::trigBasisExp(tod=1:L, L=L, degree = K)
  p$expEta_omega = exp(cbind(1, Z)%*%t(p$beta_omega))
  # parameter transform: inhomogeneous omega
  if(N>2){ # only needed if N>2
    p$omega = array(dim = c(N,N,L))
    for(t in 1:L){
      omega = matrix(0,N,N)
      omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),p$expEta_omega[t,]),N,N-1)))
      omega = t(omega)/apply(omega,2,sum)
      p$omega[,,t] = omega
    }
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  # dwell-time distributions and approximating tpm
  p$Mu_dwell = exp(cbind(1,Z)%*%t(p$beta_mu))
  dm = list()
  for(j in 1:N){
    dm[[j]] = matrix(nrow=L, ncol = agsizes[j])
    for(t in 1:L){
      dm[[j]][t,] = dnbinom(1:agsizes[j]-1, mu = p$Mu_dwell[t,j], size = 1/p$phi_dwell[j])
    }
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_phsmm(p$omega, dm)
  p$delta = LaMa::stationary_p(p$Gamma, t = t1)
  return(p)
}
mllk_sp3 = function(theta.star, X, N=3, L=24, K, agsizes){
  p = partransform_sp3(theta.star=theta.star, N=N, L=L, K=K, agsizes=agsizes, t1 = X$tod[1])
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta=p$delta, Gamma=p$Gamma, allprobs=allprobs_large, tod=X$tod)
}


### inhomogeneous in mean, dispersion and omega
partransform_sp4 = function(theta.star, N=3, L=24, K=c(1,1,1), agsizes, t1){
  p = list()
  # parameter transform: dwell-time distributions
  p$beta_Mu = matrix(theta.star[1:(N*(1+2*K[1]))], nrow = N)
  p$beta_Phi = matrix(theta.star[N*(1+2*K[1])+1:(N*(1+2*K[2]))], nrow = N)
  # parameter transform: state-dep. distributions
  p$mu = exp(theta.star[N*(2+2*(K[1]+K[2]))+1:N])
  p$sigma = exp(theta.star[N*(2+2*(K[1]+K[2]))+N+1:N])
  p$mu.turn = c(pi,0,0)
  p$kappa = exp(theta.star[N*(2+2*(K[1]+K[2]))+2*N+1:N])
  
  p$beta_omega = matrix(theta.star[N*(2+2*(K[1]+K[2]))+3*N + 1:(N*(N-2)*(1+2*K[3]))], nrow = N)
  Z1 = LaMa::trigBasisExp(tod=1:L, L=L, degree = K[1])
  Z2 = LaMa::trigBasisExp(tod=1:L, L=L, degree = K[2])
  Z3 = LaMa::trigBasisExp(tod=1:L, L=L, degree = K[3])
  p$expEta_omega = exp(cbind(1, Z3)%*%t(p$beta_omega))
  # parameter transform: inhomogeneous omega
  if(N>2){ # only needed if N>2
    p$omega = array(dim = c(N,N,L))
    for(t in 1:L){
      omega = matrix(0,N,N)
      omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),p$expEta_omega[t,]),N,N-1)))
      omega = t(omega)/apply(omega,2,sum)
      p$omega[,,t] = omega
    }
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  # dwell-time distributions and approximating tpm
  p$Mu_dwell = exp(cbind(1,Z1)%*%t(p$beta_Mu))
  p$Phi_dwell = exp(cbind(1,Z2)%*%t(p$beta_Phi))
  dm = list()
  for(j in 1:N){
    dm[[j]] = matrix(nrow=L, ncol = agsizes[j])
    for(t in 1:L){
      dm[[j]][t,] = dnbinom(1:agsizes[j]-1, mu = p$Mu_dwell[t,j], size = 1/p$Phi_dwell[t,j])
    }
  }
  p$dm = dm
  p$Gamma = LaMa::tpm_phsmm(p$omega, dm)
  p$delta = LaMa::stationary_p(p$Gamma, t = t1)
  return(p)
}
mllk_sp4 = function(theta.star, X, N=3, L=24, K=c(1,1,1), agsizes){
  p = partransform_sp4(theta.star=theta.star, N=N, L=L, K=K, agsizes=agsizes, t1 = X$tod[1])
  # allprobs matrix
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=p$mu[j]^2/p$sigma[j]^2,scale=p$sigma[j]^2/p$mu[j])*
      CircStats::dvm(X$angle[ind], p$mu.turn[j], p$kappa[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
  
  -LaMa::forward_p(delta=p$delta, Gamma=p$Gamma, allprobs=allprobs_large, tod=X$tod)
}

