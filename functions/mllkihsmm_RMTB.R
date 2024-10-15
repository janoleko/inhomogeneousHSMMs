## inhomogeneous HSMM

# dat = list(step, angle, Z, N, agsizes)
mllkiHSMM = function(par){
  getAll(par, dat)
  n = length(x)
  ## parameter transformation
  REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)
  if(N > 2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(logitomega)),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  ## computing all dwell-time distributions (pmf and cdf)
  Lambda = exp(cbind(1,Z) %*% t(beta))
  dm = lapply(1:N, function(j) sapply(1:agsizes[j]-1, dpois, lambda = Lambda[,j]))
  Fm = lapply(1:N, function(j) cbind(0, t(apply(dm[[j]][,-agsizes[j]], 1, cumsum))))
  ## calculating all state-dependent densities
  allprobs = matrix(NA, nrow = n, ncol = N)
  ind = which(!is.na(x))
  for(j in 1:N){
    allprobs[ind,j] = dnorm(x[ind], mu[j], sigma[j])
  }
  ## forward algorithm
  # inflating initial distribution
  largedelta = rep(0, sum(agsizes))
  largedelta[c(1, cumsum(agsizes)[-N])] = delta
  # only start after largest agsize of #observations
  foo = largedelta * rep(allprobs[max(agsizes),], times = agsizes)
  sumfoo = sum(foo)
  phi = foo / sumfoo
  l = log(sumfoo)
  for(t in (max(agsizes)+1):n){
    lNA = is.na(l)
    
    # extracting the dwell-time probabilities for this t.p.m. s
    dmt = lapply(1:N, function(j) diag(dm[[j]][t:(t-agsizes[j]+1),]))
    Fmt = lapply(1:N, function(j) diag(Fm[[j]][t:(t-agsizes[j]+1),]))
    Gamma = tpm_hsmm2(omega, dmt, Fmt)
    foo = (phi %*% Gamma) * rep(allprobs[t,], times = agsizes)
    sumfoo = sum(foo)
    phi = foo / sumfoo
    l = l + log(sumfoo)
    
    if(!lNA & is.na(l)) print(t)
    
  }
  -l
}

max2 = function(x,y){
  (x + y + abs(x - y)) / 2
}

## faster approximating-tpm constructer
tpm_hsmm2 <- function(omega, dm, Fm = NULL, eps = 1e-20) {
  Nv = sapply(dm, length)
  N = length(Nv)
  
  # Pre-allocate G matrix
  total_cols = sum(Nv)
  G = matrix(0, total_cols, total_cols)
  
  # compute cdf if not provided
  if(is.null(Fm)) Fm = lapply(1:N, function(i) cumsum(c(0, dm[[i]][-Nv[i]])))
  
  row_start = 1  # Track row start index for G
  for (i in 1:N) {
    Ni = Nv[i]
    # ci = ifelse(abs(1 - Fm[[i]]) > eps, dm[[i]] / (1 - Fm[[i]]), 1)
    # cim = pmax(1 - ci, 0)
    ci = dm[[i]] / (1 - Fm[[i]] + eps)
    cim = max2(1 - ci, 0)
    
    Gi = matrix(0, Ni, total_cols)  # Pre-allocate the block for Gi
    
    col_start = 1  # Track column start index for Gi
    for (j in 1:N) {
      Nj = Nv[j]
      
      if (i == j) {
        if (Ni == 1) {
          Gi[1, (col_start + Nj - 1)] = cim
        } else {
          diag_block = diag(cim[-Ni], Ni - 1, Ni - 1)
          Gi[1:(Ni - 1), (col_start + 1):(col_start + Nj - 1)] = diag_block
          Gi[Ni, (col_start + Nj - 1)] = cim[Ni]
        }
      } else {
        if (Ni == 1) {
          Gi[1, col_start:(col_start + Nj - 1)] = c(omega[i, j] * ci, rep(0, Nj - 1))
        } else {
          Gi[, col_start] = omega[i, j] * ci
        }
      }
      
      col_start = col_start + Nj
    }
    
    G[row_start:(row_start + Ni - 1), ] = Gi
    row_start = row_start + Ni
  }
  G
}

tpm_phsmm2 = function(omega, dm, eps = 1e-10){
  L = nrow(dm[[1]]) # period length
  N = length(dm) # number of states
  Ni = sapply(dm, ncol) # state aggregate sizes
  bigN = sum(Ni) # overall number of approximating states
  
  ## computing all cdfs
  Fm = lapply(1:N, function(i) cbind(0, t(apply(dm[[i]][,-Ni[i]], 1, cumsum))))
  
  ## if omega is just matrix, turn into array
  if(is.matrix(omega)){
    omega = array(omega, dim = c(N, N, L))
  }
  
  # for each timepoint, use tpm_hsmm applied to shifted pmfs and cdfs
  Gamma = array(NaN, dim = c(bigN, bigN, L))
  for(t in 1:L){
    dmt = Fmt = vector("list", N)
    
    for(i in 1:N){
      ind = (t:(t - Ni[i] + 1)) %% L
      ind[ind == 0] = L
      
      dmt[[i]] = diag(dm[[i]][ind,])
      Fmt[[i]] = diag(Fm[[i]][ind,])
    }
    # apply tpm_hsmm() to shifted pmfs and pdfs
    Gamma[,,t] = tpm_hsmm2(omega[,,t], dmt, Fm = Fmt, eps = eps)
  }
  Gamma
}

tpm_phsmm3 <- function(omega, dm, eps = 1e-10) {
  L = nrow(dm[[1]])  # period length
  N = length(dm)     # number of states
  Ni = sapply(dm, ncol)  # state aggregate sizes
  bigN = sum(Ni)     # total number of approximating states
  
  # Precompute all CDFs (vectorized cumulative sum)
  Fm = lapply(1:N, function(i) {
    temp = dm[[i]][, -Ni[i], drop = FALSE]
    cbind(0, t(apply(temp, 1, cumsum)))
  })
  
  # Convert omega to an array if it's a matrix
  if (is.matrix(omega)) {
    omega = array(omega, dim = c(N, N, L))
  }
  
  # Preallocate the Gamma array
  Gamma = array(NaN, dim = c(bigN, bigN, L))
  
  # Precompute the index shifts for each state
  state_indices = lapply(Ni, function(ni) {
    sapply(1:ni, function(k) (1:L - k + L) %% L + 1)
  })
  
  # Loop over each time point
  for (t in 1:L) {
    dmt = Fmt = vector("list", N)  # Preallocate lists for shifted dm and Fm
    
    for (i in 1:N) {
      # Use precomputed indices for cyclic shifts
      ind = state_indices[[i]][, t]
      
      # Extract shifted rows from dm and Fm using diagonal matrices
      dmt[[i]] = diag(dm[[i]][ind, , drop = FALSE])
      Fmt[[i]] = diag(Fm[[i]][ind, , drop = FALSE])
    }
    
    # Apply tpm_hsmm2 to the shifted pmfs and cdfs
    Gamma[,,t] = tpm_hsmm2(omega[,,t], dmt, Fm = Fmt, eps = eps)
  }
  
  return(Gamma)
}


## simulate data
set.seed(123)
n = 1100
z = rnorm(n) # white noise covariate

N = 2
omega = matrix(c(0,1,1,0),2,2)
mu = c(2, 5)
sigma = c(1, 1)

beta = matrix(c(log(5), log(8), 0.3, -0.3), nrow = 2)
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
  x[currentInd + 1:times] = rnorm(times, mu[s[1]], sigma[s[1]])

  currentInd = currentInd + times
  t = t+1
}

x = x[1:n]
z = z[1:n]

par = list(mu = mu, logsigma = log(sigma), logitdelta = 0,
           beta = beta)

theta.star = unlist(par)

dat = list(x = x, N = 2, Z = z, agsizes = c(30,30))

library(RTMB)

obj = MakeADFun(mllkiHSMM, par)
opt = nlminb(obj$par, obj$fn, obj$gr)


