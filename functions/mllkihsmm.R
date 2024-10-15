## inhomogeneous HSMM
mllkiHSMM = function(theta.star, X, Z, N=2, agsizes){
  n = nrow(X); maxagsize = max(agsizes)
  ## parameter assignment
  mu = theta.star[1:N]
  sigma = exp(theta.star[N+1:N])
  p = ncol(Z)
  delta = c(1, exp(theta.star[2*N + 1:(N-1)]))
  delta = delta / sum(delta)
  beta = matrix(theta.star[2*N + N-1 + 1:(N*p)], nrow = N, ncol = p)
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[2*N+N*p+N-1+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  ## computing all dwell-time distributions (pmf and cdf)
  Lambda = exp(Z %*% t(beta))
  dm = lapply(1:N, function(j) sapply(1:agsizes[j]-1, dpois, lambda = Lambda[,j]))
  Fm = lapply(1:N, function(j) cbind(0, t(apply(dm[[j]][,-agsizes[j]], 1, cumsum))))
  ## calculating all state-dependent densities
  allprobs = matrix(NA, nrow = n, ncol = N)
  ind = which(!is.na(x))
  for(j in 1:N){
    allprobs[ind,j] = dnorm(X$x[ind], mu[j], sigma[j])
  }
  allprobs_large = t(apply(allprobs, 1, rep, times = agsizes))
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
    # extracting the dwell-time probabilities for this t.p.m. s
    dmt = lapply(1:N, function(j) diag(dm[[j]][t:(t-agsizes[j]+1),]))
    Fmt = lapply(1:N, function(j) diag(Fm[[j]][t:(t-agsizes[j]+1),]))
    Gamma = tpm_hsmm2(omega, dmt, Fmt)
    foo = (phi %*% Gamma) * rep(allprobs[t,], times = agsizes)
    sumfoo = sum(foo)
    phi = foo / sumfoo
    l = l + log(sumfoo)
  }
  -l
}

## faster approximating-tpm constructer
tpm_hsmm2 <- function(omega, dm, Fm = NULL, eps = 1e-10) {
  Nv = sapply(dm, length)
  N = length(Nv)
  
  # Pre-allocate G matrix
  total_cols = sum(Nv)
  G = matrix(0, total_cols, total_cols)
  
  if(is.null(Fm)){
    Fm = lapply(1:N, function(x) cumsum(c(0, dm[[i]][-Nv[i]])))
  }
  
  row_start = 1  # Track row start index for G
  for (i in 1:N) {
    Ni = Nv[i]
    Fi = Fm[[i]]
    ci = ifelse(abs(1 - Fi) > eps, dm[[i]] / (1 - Fi), 1)
    cim = pmax(1 - ci, 0)
    
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


mod = nlm(mllkiHSMM, theta.star, X = X, Z = Z, agsizes = agsizes,
          iterlim = 1000, print.level = 2)
