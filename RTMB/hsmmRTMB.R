mllk = function(par){
  getAll(par, dat)
  
  ## parameter transformation
  sigma = exp(logsigma); REPORT(sigma); REPORT(mu)
  lambda = exp(loglambda); REPORT(lambda)
  if(N > 2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(logitomega)),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  
  ## computing all dwell-time distributions (pmf and cdf)
  dm = lapply(1:N, function(j) dpois(1:agsizes[j]-1, lambda = lambda[j]))
  # Gamma = tpm_hsmm2(omega, dm, sparse = sparse)
  # delta = stationary_sparse(Gamma)
  
  ## calculating all state-dependent densities
  allprobs = matrix(NA, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind,j] = dnorm(x[ind], mu[j], sigma[j])
  # allprobs = t(apply(allprobs, 1, rep, times = agsizes))
  
  ## forward algo
  # -forward_sparse(delta, Gamma, allprobs)
  - forward_hsmm(dm, omega, allprobs, agsizes)
}


# Simulate from basis HSMM ------------------------------------------------

n = 1100
N = 2
omega = matrix(c(0,1,1,0),2,2)
mu = c(1, 5)
sigma = c(1, 1)
lambda = c(5, 8)
n_tilde = ceiling(n / (1 + mean(lambda))*1.5)

set.seed(123)
s = rep(NA, n_tilde)
C = x = rep(NA, 2*n)
s[1] = 1
for(t in 2:n_tilde){
  s[t] = sample(1:N, size=1, prob = omega[s[t-1],])
}
times = rpois(1, lambda = lambda[s[1]]) + 1
C[1:times] = rep(s[1], times)
x[1:times] = rnorm(times, mu[s[1]], sigma[s[1]])
currentInd = times
t = 2
while(currentInd <= n){
  times = rpois(1, lambda = lambda[s[t]]) + 1
  C[currentInd + 1:times] = rep(s[t], times)
  x[currentInd + 1:times] = rnorm(times, mu[s[t]], sigma[s[t]])
  currentInd = currentInd + times
  t = t+1
}
x = x[1:n]


par = list(mu = mu, logsigma = log(sigma), loglambda = log(lambda))
dat = list(x = x, N = 2, Z = z, agsizes = c(40,40), sparse = TRUE)

obj = MakeADFun(mllk, par)

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr) 
)

mod = obj$report()
mod$mu
mod$sigma
mod$lambda



# Sparse implementation ---------------------------------------------------

library(Matrix)

tpm_hsmm2 <- function(omega, dm, Fm = NULL, sparse = FALSE, eps = 1e-20) {
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
  if(sparse){
    G = methods::as(G, "sparseMatrix")
  }
  G
}

## faster approximating-tpm constructer with sparse matrix
# tpm_hsmm_sparse <- function(omega, dm, Fm = NULL, eps = 1e-20) {
#   Nv = sapply(dm, length)
#   N = length(Nv)
#   
#   # compute cdf if not provided
#   if (is.null(Fm)) Fm = lapply(1:N, function(i) cumsum(c(0, dm[[i]][-Nv[i]])))
#   
#   row_idx = c()  # To store row indices of non-zero elements
#   col_idx = c()  # To store column indices of non-zero elements
#   values = c()   # To store the corresponding values of non-zero elements
#   
#   row_start = 1  # Track row start index for G
#   for (i in 1:N) {
#     Ni = Nv[i]
#     ci = dm[[i]] / (1 - Fm[[i]] + eps)
#     cim = max2(1 - ci, 0)
#     
#     col_start = 1  # Track column start index for Gi
#     for (j in 1:N) {
#       Nj = Nv[j]
#       
#       if (i == j) {
#         if (Ni == 1) {
#           # Only one element to update for diagonal block
#           row_idx = c(row_idx, row_start)
#           col_idx = c(col_idx, col_start + Nj - 1)
#           values = c(values, cim)
#         } else {
#           # Diagonal block for larger sub-matrix
#           diag_block = diag(cim[-Ni], Ni - 1, Ni - 1)
#           for (k in 1:(Ni - 1)) {
#             for (l in 1:(Ni - 1)) {
#               if (diag_block[k, l] != 0) {
#                 row_idx = c(row_idx, row_start + k - 1)
#                 col_idx = c(col_idx, col_start + l)
#                 values = c(values, diag_block[k, l])
#               }
#             }
#           }
#           # Last element in the diagonal block
#           row_idx = c(row_idx, row_start + Ni - 1)
#           col_idx = c(col_idx, col_start + Nj - 1)
#           values = c(values, cim[Ni])
#         }
#       } else {
#         if (Ni == 1) {
#           row_idx = c(row_idx, row_start)
#           col_idx = c(col_idx, col_start:(col_start + Nj - 1))
#           values = c(values, c(omega[i, j] * ci, rep(0, Nj - 1)))
#         } else {
#           row_idx = c(row_idx, row_start:(row_start + Ni - 1))
#           col_idx = c(col_idx, rep(col_start, Ni))
#           values = c(values, omega[i, j] * ci)
#         }
#       }
#       
#       col_start = col_start + Nj
#     }
#     
#     row_start = row_start + Ni
#   }
#   
#   # Construct sparse matrix using the row indices, column indices, and values
#   G_sparse = Matrix::sparseMatrix(i = row_idx, j = col_idx, x = values, dims = c(sum(Nv), sum(Nv)))
#   
#   return(G_sparse)
# }

Gamma1 = tpm_hsmm2(omega, dm)
Gamma2 = tpm_hsmm_sparse(omega, dm)

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
  
  names(delta) = paste("state", 1:N)
  delta
}

forward_sparse = function(delta, Gamma, allprobs){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
  
  if(is(Gamma, "sparseMatrix")){
    if(!is(allprobs, "sparseMatrix")){
      allprobs = as(allprobs, "sparseMatrix")
    }
    
    delta = methods::as(delta, "sparseMatrix")
  }
    
  # forward algorithm
  foo = delta %*% diag(allprobs[1,])
  sumfoo = sum(foo)
  phi = foo / sumfoo
  l = log(sumfoo)
  
  for(t in 2:nrow(allprobs)) {
    foo = phi %*% Gamma %*% diag(allprobs[t,])
    sumfoo = sum(foo)
    phi = foo / sumfoo
    l = l + log(sumfoo)
  }
  
  l
}


