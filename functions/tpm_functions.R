# differentiable max function
max2 = function(x,y){
  (x + y + abs(x - y)) / 2
}

# function for building the HSMM-approximating tpm
tpm_hsmm <- function(omega, dm, 
                     Fm = NULL, sparse = TRUE, eps = 1e-20) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
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
  if(sparse){
    G = methods::as(G, "sparseMatrix")
  }
  G
}

# function for building a list of sparse approximating tpms for inhomogeneous HSMMs
tpm_ihsmm = function(omega, dm, 
                     eps = 1e-20){
  n = nrow(dm[[1]]) # length of timeseries
  N = length(dm) # number of states
  Ni = sapply(dm, ncol) # state aggregate sizes
  bigN = sum(Ni) # overall number of approximating states
  maxNi = max(Ni) # maximum state aggregate size
  
  ## computing all cdfs
  Fm = lapply(1:N, function(i) cbind(0, t(apply(dm[[i]][,-Ni[i]], 1, cumsum))))
  
  ## if omega is just matrix, turn into array
  if(is.matrix(omega)){
    omega = array(omega, dim = c(N, N, n))
  }
  
  # for each timepoint, use tpm_hsmm applied to shifted pmfs and cdfs
  Gamma = vector("list", n - maxNi + 1)
  for(k in 1:(n - maxNi + 1)){
    t = k + maxNi - 1
    dmt = Fmt = vector("list", N)
    
    for(i in 1:N){
      ind = t:(t - Ni[i] + 1)
      
      dmt[[i]] = diag(dm[[i]][ind,])
      Fmt[[i]] = diag(Fm[[i]][ind,])
    }
    # apply tpm_hsmm() to shifted pmfs and pdfs
    Gamma[[k]] = tpm_hsmm(omega[,,t], dmt, Fm = Fmt, sparse = TRUE, eps = eps)
  }
  Gamma
}

# function for building a list of sparse approximating tpms for periodically inhomogeneous HSMMs
tpm_phsmm = function(omega, dm, 
                     eps = 1e-10){
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
  Gamma = vector("list", L)
  
  for(t in 1:L){
    dmt = Fmt = vector("list", N)
    
    for(i in 1:N){
      ind = (t:(t - Ni[i] + 1)) %% L
      ind[ind == 0] = L
      
      dmt[[i]] = diag(dm[[i]][ind,])
      Fmt[[i]] = diag(Fm[[i]][ind,])
    }
    # apply tpm_hsmm() to shifted pmfs and pdfs
    Gamma[[t]] = tpm_hsmm(omega[,,t], dmt, Fm = Fmt, sparse = TRUE, eps = eps)
  }
  Gamma
}

stationary_sparse = function(Gamma) 
{
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  N = dim(Gamma)[1]
  
  if(methods::is(Gamma, "sparseMatrix")) {
    Gamma = as.matrix(Gamma)
  }
  delta = Matrix::solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  delta
}

stationary_p_sparse = function(Gamma, t = NULL){
  # overloading assignment operators, currently necessary
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  L = length(Gamma) # cycle length
  N = dim(Gamma[[1]])[1] # tpm dimension
  
  if(is.null(t)) {
    Delta = matrix(NaN, nrow = L, ncol = N)
    
    GammaT = Gamma[[1]]
    for(k in 2:L) GammaT = GammaT %*% Gamma[[k]]
    
    Delta[1,] = stationary_sparse(GammaT)
    
    for(t in 2:L){
      Delta[t,] = as.numeric(Delta[t-1,, drop = FALSE] %*% Gamma[[t-1]])
    }
    return(Delta)
  } else{
    GammaT = Gamma[[t]]
    if(t < L){
      for(k in (t+1):L) GammaT = GammaT %*% Gamma[[k]]
      if(t > 1){
        for(k in 1:(t-1)) GammaT = GammaT %*% Gamma[[k]]
      }
    } else if(t == L){
      for(k in 1:(L-1)) GammaT = GammaT %*% Gamma[[k]]
    }
    delta = stationary_sparse(GammaT)
    return(delta)
  }
}


