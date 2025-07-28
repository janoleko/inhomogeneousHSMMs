forward_ihsmm <- function(dm, omega, allprobs, agsizes,
                          trackID = NULL, delta = NULL, eps = 1e-20, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  maxag = max(agsizes) # maximum number of states in the approximating HMM
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_ihsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_sparse(Gamma_sparse[[1]])
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
      foo = delta_sparse %*% diag(rep(allprobs[maxag,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in (maxag + 1):nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        foo = phi %*% Gamma_sparse[[t-maxag]] %*% diag(rep(allprobs[t,], times = agsizes))
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
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          delta = stationary_sparse(Gamma_sparse[[1]]) # build k deltas based on Gammas
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
        foo = phi %*% Gamma_sparse[[t-1]] %*% diag(rep(allprobs[t,], times = agsizes))
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
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(is.matrix(omega)){ # fixed omega
        omega = array(omega, dim = c(dim(omega), n))
      } else if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
      }
      
      # nested list of sparse approximating tpms
      Gamma_sparse = lapply(1:k, function(i){
        ind = which(trackID == uID[i])
        this_dm = lapply(1:N, function(j) dm[[j]][ind,])
        tpm_ihsmm(omega[,,ind], this_dm, eps = eps)
      })
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = LaMa::calc_trackInd(trackID)
        Delta = t(sapply(1:k, function(i) stationary_sparse(Gamma_sparse[[i]][[1]]))) # build k deltas based on Gammas
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
        Gamma_i = Gamma_sparse[[i]]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[maxag],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in (maxag + 1):length(ind)) {
          foo = phi %*% Gamma_i[[t-maxag]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
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
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = LaMa::calc_trackInd(trackID)
          Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[trackInd[k]]]))) # build k deltas based on Gammas
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
        Gamma_i = Gamma_sparse[ind]
        
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_i[[t-1]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
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
trackID = rep(1:k, each = n)
allprobs = matrix(runif(n*k*N), nrow = n*k, ncol = N)

agsizes = c(10, 20, 25)
lambda = c(5, 7, 8)
dm = lapply(1:3, function(i) sapply(1:agsizes[i]-1, dpois, lambda = rep(lambda[i], n*k)))
# dm = lapply(1:3, function(i) dpois(1:agsizes[i]-1, lambda[i]))


omega = matrix(c(0, 0.8, 0.2,
                 0.1, 0, 0.9,
                 0.4, 0.6, 0), nrow = 3, ncol = 3, byrow = TRUE)

omega = array(omega, dim = c(N,N,n*k))
delta = c(0.3, 0.2, 0.5)

# no tracks
## stationary
system.time(
  forward_ihsmm(dm, omega, allprobs, agsizes)
)
## not stationary
forward_ihsmm(dm, omega, allprobs, agsizes, delta = delta)
## fixed dm
dm = lapply(1:3, function(i) dpois(1:agsizes[i]-1, lambda[i]))
forward_ihsmm(dm, omega, allprobs, agsizes)
forward_ihsmm(dm, omega, allprobs, agsizes, delta = delta)
# fixed omega
omega = omega[,,1]
forward_ihsmm(dm, omega, allprobs, agsizes)
forward_ihsmm(dm, omega, allprobs, agsizes, delta = delta)

# tracks, stationary
## dm varying, omega constant
system.time(
  forward_ihsmm(dm, omega, allprobs, agsizes, trackID = trackID)
)


## dm constant, omega varying
forward_ihsmm(dm, omega, allprobs, agsizes, trackID = trackID)

## both varying
forward_hsmm(dml, omegaA, allprobs, agsizes, trackID = trackID)


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
