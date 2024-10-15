ddwell_t = function(
    r, # vector of dwell times to compute
    t, # time point at which the state is entered
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3] # cycle length
  ind = (t+(1:max(r))-1)%%L
  ind[which(ind==0)] = L
  gamma_ii = Gamma[state,state,ind]
  pmf = c(1, cumprod(gamma_ii)[-length(ind)])*(1-gamma_ii)
  return(pmf[r])
}

ddwell = function(
    r, # vector of dwell times to compute
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3]
  N = dim(Gamma)[1]
  delta = matrix(NA, L, N) # calculate p-stationary
  GammaT=Gamma[,,1]
  for (k in 2:L){ GammaT=GammaT%*%Gamma[,,k] }
  delta[1,] = solve(t(diag(N)-GammaT+1), rep(1,N))
  for(k in 2:L){ delta[k,] = delta[k-1,]%*%Gamma[,,k-1] }
  weights=numeric(L) # calculate weights
  weights[1] = sum(delta[L,-state] * Gamma[-state,state, L])
  for (k in 2:L){ weights[k] = sum(delta[k-1,-state] * Gamma[-state,state, k-1]) }
  weights = weights / sum(weights)
  pmfs_weighted = matrix(NA, L, max(r))# calculate all weighted d_i^t's
  for(k in 1:L){ pmfs_weighted[k,] = weights[k] * ddwell_t(1:max(r), k, state, Gamma) }
  pmf = apply(pmfs_weighted, 2, sum)
  return(pmf[r])
}

ddwell_hsmm = function(r, # scalar or vector of dwell times for which to compute probs
                       state, # which state to compute
                       dm, # list of L pre-calculated dwell-time distributions
                       Gamma, # list of L pre-calculated periodic state-aggregate Gamma matrices
                       aggr_sizes # vector of state-aggregate sizes
){
  L = dim(Gamma)[3]
  largeN = sum(aggr_sizes)
  N = length(aggr_sizes)
  I_minus = c(0,cumsum(aggr_sizes)[-N])[state]+1 # index of state to transition to
  # index set for sum
  if(state > 1){aggr_ind = cumsum(aggr_sizes)[state-1]+1:aggr_sizes[state]
  } else{aggr_ind = 1:aggr_sizes[1]}
  I = setdiff(1:sum(aggr_sizes), aggr_ind)
  delta = LaMa::stationary_p(Gamma)
  weights = numeric(L)
  weights[1] = sum(delta[L,I] * Gamma[I,I_minus,L])
  for (t in 2:L){ 
    weights[t] = sum(delta[t-1, I]*Gamma[I,I_minus,t-1]) 
  }
  weights = weights/sum(weights)
  pmfs = matrix(NA, L, length(r))
  for(t in 1:L){
    # pmfs[t,] = dm[[t]][[state]][r]*weights[t] 
    pmfs[t,] = dm[[state]][t,r]*weights[t] 
  }
  pmf = apply(pmfs, 2, sum)
  return(pmf)
}

empirical_pmf = function(states, agsizes){
  N = length(unique(states))
  d = rle(states)
  dtimes = data.frame(values=d$values, lengths=d$lengths)
  emppmf = list()
  for(j in 1:N){
    emppmf[[j]] = prop.table(table(factor(dtimes[which(dtimes$values == j),2], levels = 1:agsizes[j])))
  }
  return(emppmf)
}
