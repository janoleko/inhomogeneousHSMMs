##### Part 2: Consistency #####


# Functions and libraries -------------------------------------------------

source("./functions/sim_study_functions.R") # simulating
source("./simulation_study/likelihoods.R") # likelihoods

library(doParallel)
library(foreach)
library(LaMa)
library(RTMB)


# Parameters and hyperparamters -------------------------------------------

## Parameters
beta = matrix(c(log(c(8,7,6)), -0.2, 0.2, -0.6, 0.3, -0.2, -0.4), nrow = 3)
omega = matrix(c(0, 0.7, 0.3, 
                 0.2, 0, 0.8,
                 0.5, 0.5, 0), nrow = 3, byrow = TRUE)
obsparams = list(
  mu = c(20, 200, 800),
  sigma = c(20, 150, 500),
  kappa = c(0.2, 1, 2.5)
)

## Hyperparameter: aggregate sizes
todseq = seq(0, 24, length=200)
Z = cbind(1, trigBasisExp(todseq, 24))
dM = exp(Z %*% t(beta))

# finding maximum mean dwell time for each state
maxlambda = apply(dM, 2, max)
# calculating aggreagte sizes based on 99.5 percent quantile
(agsizes = ceiling((qpois(0.995, maxlambda)+1) * 1.1))


# Simulation --------------------------------------------------------------

# vector of different time series lengths T
nobs = c(1000, 2000, 5000, 10000)
# number of runs for each length
nruns = 500

Data = list()
# simulating 500 data sets of maximum length
set.seed(123)
for(i in 1:nruns){
  cat("\n", i)
  Data[[i]] = sim_phsmm(max(nobs), beta, omega, obsparams)
}

# initial parameter
par = list(
  beta = beta, 
  logitomega = rep(0,3),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# Prepare array to store results
Betas = array(dim = c(3, 3, nruns, length(nobs)))

# Fitting HSMMs to data sets of increasing lengths
for(k in 1:length(nobs)){
  cat("Number of observations:", nobs[k], "\n")
  
  # Define data
  dat = list(
    tod = Data[[1]]$tod[1:nobs[k]],
    N = 3,
    agsizes = agsizes,
    Z = cbind(1, trigBasisExp(1:24))
  )
  
  # Initial obs object
  obs = list(step = Data[[1]]$step[1:nobs[k]], angle = Data[[1]]$angle[1:nobs[k]])
  
  # Create objective function once for each nObs
  obj = MakeADFun(nllpHSMM, par, silent = TRUE)
  
  # Create a cluster of cores
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # Parallel execution using foreach
  result = foreach(i = 1:nruns,
                   .packages = c("LaMa", "RTMB"), # export packages to workers
                   .errorhandling = "pass") %dopar% {
                     
                     # Update the obs and dat environment for each iteration
                     assign("obs", list(step = Data[[i]]$step[1:nobs[k]], angle = Data[[i]]$angle[1:nobs[k]]), envir = .GlobalEnv)
                     assign("dat", list(tod = Data[[i]]$tod[1:nobs[k]], N = 3, agsizes = agsizes, Z = cbind(1, trigBasisExp(1:24))), envir = .GlobalEnv)
                     
                     # Fit the model, this will automatically use the new obs because DataEval() in nllpHSMM
                     opt = nlminb(obj$par, obj$fn, obj$gr)
                     
                     # Get the model report
                     mod = obj$report()
                     
                     # Return the beta values
                     mod$beta
                   }
  
  # Stop the cluster after execution
  stopCluster(cl)
  
  # Store results in Betas
  for(i in 1:nruns) Betas[,,i,k] = result[[i]]
}

saveRDS(Betas, file = "./simulation_study/Results/Part2_consistency/Betas.rds")




Betas = array(dim = c(3, 3, nruns, length(nobs)))

# fitting HSMMs to data sets of increasing lengths
for(k in 1:length(nobs)){
  cat("Scenario", k, "\n")
  cat("Number of observations:", nobs[k], "\n")
  
  # defining observations
  obs = list(step = Data[[1]]$step[1:nobs[1]], angle = Data[[1]]$angle[1:nobs[1]])
  
  # defining other stuff
  dat = list(tod = Data[[1]]$tod[1:nobs[[1]]],
             N = 3, agsizes = agsizes, Z = cbind(1, trigBasisExp(1:24)))
  
  cat("Creating objective function...\n")
  obj = MakeADFun(nllpHSMM, par, silent = TRUE)
  
  # for(i in 1:nruns){
  #   cat(i, "\n")
  #   obs = list(step = Data[[i]]$step[1:nobs[k]], angle = Data[[i]]$angle[1:nobs[k]])
  #   
  #   
  #   opt = nlminb(obj$par, obj$fn, obj$gr)
  #   
  #   mod = obj$report()
  #   Betas[,,i,k] = mod$beta
  # }

  # Create a cluster of cores
  cl = makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  result = foreach(i = 1:nruns,
                   .packages = c("LaMa", "RTMB"),
                   .export = c("obj", "dat"),
                   .errorhandling = "pass") %dopar% {

    obs = list(step = Data[[i]]$step[1:nobs[k]], angle = Data[[i]]$angle[1:nobs[k]])

    opt = nlminb(obj$par, obj$fn, obj$gr)

    mod = obj$report()

    mod$beta
  }
  stopCluster(cl)
  
  saveRDS(result,
          file = paste0("./simulation_study/simulation_results/consistency/results_", nobs[k], ".rds"))
}








# function for applying
one_rep = function(data, beta, omega, stateparams, n, agsizes, stepmax){
  mod = fitpHSMM(data=data[1:n, ], beta=beta, stateparams=stateparams, 
                 agsizes=agsizes, stepmax=stepmax)
  return(mod)
}


# defining colors
color = c("orange", "deepskyblue", "seagreen2")

# looking at dwell time mean as a function of the time of day
todseq = seq(0,24,length=200)
Ztod = cbind(1, LaMa::trigBasisExp(todseq, 24))
dM = exp(Ztod%*%t(beta))
plot(todseq, dM[,1], type = "l", ylim = c(0,15), col = color[1], lwd=2, bty="n", 
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
lines(todseq, dM[,2], type = "l", col = color[2], lwd=2)
lines(todseq, dM[,3], type = "l", col = color[3], lwd=2)
axis(1, at = seq(0,24,by=4), labels=seq(0,24,by=4))


# finding maximum mean dwell time for each state
maxlambda = apply(dM, 2, max)
# calculating aggreagte sizes based on 99.5 percent quantile
agsizes = ceiling((qpois(0.995, maxlambda)+1)*1.1) 

# checking state dependent distributions
par(mfrow = c(1,1))
delta = c(sum(data$C==1), sum(data$C==2), sum(data$C==3))/nrow(data)
hist(data$step, prob = T, breaks = 50)
curve(delta[1]*dgamma(x, shape=mu[1]^2/sigma[1]^2, scale=sigma[1]^2/mu[1]),add=T, lwd=2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape=mu[2]^2/sigma[2]^2, scale=sigma[2]^2/mu[2]),add=T, lwd=2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape=mu[3]^2/sigma[3]^2, scale=sigma[3]^2/mu[3]),add=T, lwd=2, col = color[3], n = 500)
curve(delta[1]*dgamma(x, shape=mu[1]^2/sigma[1]^2, scale=sigma[1]^2/mu[1])+
        delta[2]*dgamma(x, shape=mu[2]^2/sigma[2]^2, scale=sigma[2]^2/mu[2])+
        delta[3]*dgamma(x, shape=mu[3]^2/sigma[3]^2, scale=sigma[3]^2/mu[3]), add=T, lwd=2, lty=2, n=500)
hist(data$angle, prob = T, breaks = 50)
curve(delta[1]*CircStats::dvm(x, 0, kappa[1]), add=T, lwd=2, col = color[1], n = 500)
curve(delta[2]*CircStats::dvm(x, 0, kappa[2]), add=T, lwd=2, col = color[2], n = 500)
curve(delta[3]*CircStats::dvm(x, 0, kappa[3]), add=T, lwd=2, col = color[3], n = 500)
curve(delta[1]*CircStats::dvm(x, 0, kappa[1])+
        delta[2]*CircStats::dvm(x, 0, kappa[2])+
        delta[3]*CircStats::dvm(x, 0, kappa[3]), add=T, lwd=2, lty=2, n=500)



# Simulation --------------------------------------------------------------

# vector of differnt time series lengths T
nobs = c(1000, 2000, 5000, 10000)
# number of runs for each length
nruns = 500

Data = list()
# simulating 500 data sets of maximum length
set.seed(123)
for(i in 1:nruns){
  cat("\n", i)
  Data[[i]] = sim_phsmm(max(nobs), beta, omega, stateparams)
}

# fitting HSMMs to data sets of increasing lengths
for(k in 1:length(nobs)){
  cat("\nScenario", k)
  cat("\nNumber of observations:", nobs[k])

  results = parallel::mclapply(Data, FUN = one_rep,
                    beta=beta, omega=omega, stateparams=stateparams,
                    n=nobs[k], agsizes=agsizes, stepmax = 10,
                    mc.cores = parallel::detectCores()-2)
  saveRDS(results,
          file = paste0("./simulation_study/simulation_results/consistency/results_", nobs[k], ".rds"))
}



# Visualizing results -----------------------------------------------------

## Coefficients of trigonometric predictor
N=3; K=1
Betas = array(dim = c(N, 1+2*K, 500, 4))

for(k in 1:length(nobs)){
  res = readRDS(paste0("./simulation_study/simulation_results/consistency/results_", nobs[k], ".rds"))
  for(i in 1:nruns){
    if(!is.character(res[[i]])){
      Betas[,,i,k] = res[[i]]$beta
    }
  }
}

# for(state in 1:3){
#   #pdf(file = paste0("./figures/simulation/consistency_state", state, ".pdf"), width=8, height=5)
#   par(mfrow = c(3,4), mar = c(1.3,4.3,2.5,0.5)+0.1)
#   ylims = apply(Betas[state,,,1], 1, quantile, probs = c(0.0025, 0.9975), na.rm = TRUE)
#   for(p in 1:3){
#     for(k in 1:4){
#       if(p==1){
#         main = paste0("T = ", nobs[k])
#       } else{ main = "" }
#       if(k==1){
#         ylab = bquote({beta^(.(state))} [.(p-1)])
#       } else{ ylab = ""}
#       boxplot(Betas[state, p,,k], ylim = ylims[,p], ylab = ylab, 
#               col = "gray95", main=main)
#       abline(h = beta[state,p], col = color[state], lwd = 2)
#       # abline(h = mean(Betas[state, p,,k], na.rm = TRUE), col = "deepskyblue")
#     }
#   }
#   #dev.off()
# }

# plotting distribution of dwell-time parameters for different T's

## boxplots
for(state in 1:3){
  #pdf(file = paste0("./figures/simulation/consistency_state", state, ".pdf"), width=7.5, height=2.5)
  par(mfrow = c(1,3), mar = c(4,4.5,1,2)+0.1)
  # ylims = apply(Betas[state,,,1], 1, quantile, probs = c(0.001, 0.999), na.rm = TRUE)
  for(p in 1:3){
    ylims = apply(Betas[state,p,,], 2, range, na.rm = TRUE)
    ylim = c(min(ylims[1,]), max(ylims[2,]))
    b = Betas[state, p,,]
    B = data.frame(beta = as.numeric(b), nobs = rep(nobs, each=500))
    B$nobs = as.factor(B$nobs)
    boxplot(B$beta~B$nobs, ylim = ylim, 
            ylab = bquote({beta^(.(state))} [.(p-1)]), xlab = "number of observations",
            col = "gray95", main=main, frame = FALSE)
    abline(h = beta[state,p], col = scales::alpha(color[state],0.8), lwd = 1.5)
    # abline(h = mean(Betas[state, p,,k], na.rm = TRUE), col = "deepskyblue")
  }
  #dev.off()
}

## histograms
for(state in 1:3){
  # pdf(file = paste0("./figures/simulation/consistency_state_hist", state, ".pdf"), width=10, height=7)
  par(mfrow = c(3,4), mar = c(4.7,4,4,0.2)+0.1)
  xlims = apply(Betas[state,,,1], 1, quantile, probs = c(0.0025, 0.9975), na.rm = TRUE)
  for(p in 1:3){
    for(k in 1:4){
      if(p==1){
        main = paste0("T = ", nobs[k])
      } else{ 
        main = "" 
        }
      if(k==1){
        xlab = bquote({beta^(.(state))} [.(p-1)])
        ylab = "density"
      } else{ ylab = ""}
      hist(Betas[state, p,, k], xlab = xlab, bor = "white", main = main, prob = TRUE, ylab = ylab,
           xlim = xlims[,p], breaks = seq(xlims[1,p]-0.2, xlims[2,p]+0.2, length=30))
      lines(density(Betas[state, p,, k]), col = color[state], lwd = 1, lty = 2)
      abline(v = beta[state,p], col = color[state], lwd = 2)
    }
  }
  # dev.off()
}
