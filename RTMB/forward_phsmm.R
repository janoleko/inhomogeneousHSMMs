forward_phsmm <- function(dm, omega, allprobs, agsizes, tod,
                          trackID = NULL, delta = NULL, eps = 1e-20, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  L = length(unique(tod)) # cycle length -> number of unique tpms
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_p_sparse(Gamma_sparse, t = tod[1])
        delta_sparse = methods::as(t(delta), "sparseMatrix")
      } else{ # if delta is provided, "stuff out" with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                            x = delta, dims = c(1, M))
      }
      
      ## report quantities for state decoding
      if(report) {
        RTMB::REPORT(delta_sparse)
        RTMB::REPORT(Gamma_sparse)
        RTMB::REPORT(allprobs)
      }
      
      # forward algorithm
      # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
      foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(is.vector(dm[[1]])){
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
          delta_sparse = Matrix::Matrix(delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = delta, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          Gamma_array = simplify2array(lapply(Gamma_sparse, as.matrix))
          delta = LaMa::stationary_p(Gamma_array, t = tod[1])
          delta_sparse = methods::as(t(delta), "sparseMatrix")
        } else{
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                              x = delta, dims = c(1, M))
        }
      }
      
      ## forward algorithm
      foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){ # dm varies
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
      }
      
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps) # still only L unique tpms
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = LaMa::calc_trackInd(trackID)
        Delta = stationary_p_sparse(Gamma_sparse) # calculate all periodically stationaries -> not very expensive because of recursion
        Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
        Delta_sparse = methods::as(Delta, "sparseMatrix") 
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
      
      ## second case: dm does not vary, only omega varies -> easier
    } else if(is.vector(dm[[1]])){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = LaMa::calc_trackInd(trackID)
          Gamma_array = simplify2array(lapply(Gamma_sparse, as.matrix))
          Delta = LaMa::stationary_p(Gamma_array) # calculate all periodically stationaries -> not very expensive because of recursion
          Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
          Delta_sparse = methods::as(Delta, "sparseMatrix") 
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
    }
  }
  
  l
}



# Testing -----------------------------------------------------------------

N = 3
n = 100
k = 10
L = 24
tod = rep(1:L, ceiling(n*(k+1) / L))[1:(n*k)]
trackID = rep(1:k, each = n)
allprobs = matrix(runif(n*k*N), nrow = n*k, ncol = N)

agsizes = c(10, 20, 25)
lambda = c(5, 7, 8)
dm = lapply(1:3, function(i) sapply(1:agsizes[i]-1, dpois, lambda = rep(lambda[i], L)))
# dm = lapply(1:3, function(i) dpois(1:agsizes[i]-1, lambda[i]))


omega = matrix(c(0, 0.8, 0.2,
                 0.1, 0, 0.9,
                 0.4, 0.6, 0), nrow = 3, ncol = 3, byrow = TRUE)

omega = array(omega, dim = c(N,N,L))
delta = c(0.3, 0.2, 0.5)

# no tracks
## stationary
forward_phsmm(dm, omega, allprobs, agsizes, tod)

## not stationary
forward_phsmm(dm, omega, allprobs, agsizes, tod, delta = delta)
## fixed dm
dm = lapply(1:3, function(i) dpois(1:agsizes[i]-1, lambda[i]))
forward_phsmm(dm, omega, allprobs, agsizes, tod)
forward_phsmm(dm, omega, allprobs, agsizes, tod, delta = delta)
# fixed omega
omega = omega[,,1]
forward_phsmm(dm, omega, allprobs, agsizes, tod)
forward_phsmm(dm, omega, allprobs, agsizes, tod, delta = delta)

# tracks, stationary
## dm varying, omega constant
forward_phsmm(dm, omega, allprobs, agsizes, tod, trackID = trackID)

## dm constant, omega varying
dm = lapply(1:3, function(i) dpois(1:agsizes[i]-1, lambda[i]))
omega = array(omega, dim = c(N,N,L))
forward_phsmm(dm, omega, allprobs, agsizes, tod, trackID = trackID)

## both varying
dm = lapply(1:3, function(i) sapply(1:agsizes[i]-1, dpois, lambda = rep(lambda[i], L)))
forward_phsmm(dm, omega, allprobs, agsizes, tod, trackID = trackID)

# tracks, not stationary, one initial distribution
## dm and omega constant
forward_hsmm(dm, omega, allprobs, agsizes, trackID = trackID, delta = delta)

## dm varying, omega constant
forward_hsmm(dml, omega, allprobs, agsizes, trackID = trackID, delta = delta)

## dm constant, omega varying
forward_hsmm(dm, omegaA, allprobs, agsizes, trackID = trackID, delta = delta)

## both varying
forward_hsmm(dml, omegaA, allprobs, agsizes, trackID = trackID, delta = delta)

# tracks, not stationary, different initial distributions
Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)

## dm and omega constant
forward_hsmm(dm, omega, allprobs, agsizes, trackID = trackID, delta = Delta)

## dm varying, omega constant
forward_hsmm(dml, omega, allprobs, agsizes, trackID = trackID, delta = Delta)

## dm constant, omega varying
forward_hsmm(dm, omegaA, allprobs, agsizes, trackID = trackID, delta = Delta)

## both varying
forward_hsmm(dml, omegaA, allprobs, agsizes, trackID = trackID, delta = Delta)
