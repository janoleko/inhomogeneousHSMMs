##### Part 3: Model misspecification #####


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

# get maximum mean for each state and calculate aggreagte sizes
maxlambda = apply(dM, 2, max)
agsizes = ceiling(qpois(0.995, maxlambda) + 2)

# color
color = c("orange", "deepskyblue", "seagreen2")


### Simulation 1: one really long data set ----------------------------------

nObs = 1e6
set.seed(123)
data = sim_phsmm(nObs, beta, omega, obsparams)

## plotting state-dependent distributions
mu = obsparams$mu
sigma = obsparams$sigma
kappa = obsparams$kappa
delta = c(sum(data$C==1), sum(data$C==2), sum(data$C==3))/nrow(data)

# pdf("./figures/simulation/state_dependent.pdf", width=8, height = 4)
par(mfrow = c(1,2), xpd = F)

hist(data$step, prob = T, breaks = 100, bor = "white", main = "", 
     ylim = c(0, 0.0025), xlim = c(0,2500), ylab = "density", xlab = "turning angle")
for(j in 1:3) curve(delta[1]*dgamma2(x, mu[j], sigma[j]), add=T, lwd=2, col = color[j], n = 500)
curve(delta[1]*dgamma2(x, mu[1], sigma[1])+
        delta[2]*dgamma2(x, mu[2], sigma[2])+
        delta[3]*dgamma2(x, mu[3], sigma[3]), add=T, lwd=2, lty=2, n=500)

hist(data$angle, prob = T, breaks = 50, bor = "white", main = "", ylab = "density", xlab = "turning angle")
for(j in 1:3) curve(delta[j]*dvm(x, 0, kappa[j]), add = T, lwd = 2, col = color[j], n = 500)
curve(delta[1]*dvm(x, 0, kappa[1])+delta[2]*dvm(x, 0, kappa[2])+delta[3]*dvm(x, 0, kappa[3]),
      add=T, lwd=2, lty=2, n=500)
# dev.off()


# Model fitting -----------------------------------------------------------

### Model 1: true periodic HSMM ###

# initial parameter
par1 = list(
  beta = beta, 
  logitomega = rep(0,3),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# hyperparmeters and data
dat = list(
  tod = data$tod,
  N = 3,
  agsizes = agsizes,
  Z = cbind(1, trigBasisExp(1:24))
)

# steps and turns
obs = list(step = data$step, angle = data$angle)

# Create objective function
### CAUTION: Memory intensive ###
obj1 = MakeADFun(nllpHSMM, par1)

# Fit the model
opt1 = nlminb(obj$par, obj$fn, obj$gr)

# reporting
mod1 = obj1$report()
mod1$par = opt1$par
mod1$Gamma = simplify2array(lapply(mod1$Gamma_sparse, as.matrix))

# saving fitted model object
saveRDS(mod1, file = "./simulation_study/Results/Part3_misspecification/pHSMM.rds")


### Model 2: periodic HMM ###

# initial parameter
par2 = list(
  beta = matrix(c(rep(-2,6), rep(0,12)), nrow = 6),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# hyperparmeters and data
dat = list(
  tod = data$tod,
  N = 3,
  Z = cbind(1, trigBasisExp(1:24))
)

# steps and turns
obs = list(step = data$step, angle = data$angle)

# Create objective function
### CAUTION: Memory intensive ###
obj2 = MakeADFun(nllpHMM, par2)

# Fit the model
opt2 = nlminb(obj2$par, obj2$fn, obj2$gr)

# reporting
mod2 = obj2$report()
mod2$par = opt2$par

# saving fitted model object
saveRDS(mod2, file = "./simulation_study/Results/Part3_misspecification/pHMM.rds")


### Model 3: homogeneous HSMM ###

# initial parameter
par3 = list(
  loglambda = log(c(8,7,6)),
  logitomega = rep(0,3),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# hyperparmeters and data
dat = list(
  N = 3,
  agsizes = agsizes
)

# steps and turns
obs = list(step = data$step, angle = data$angle)

# Create objective function (this takes some time here)
obj3 = MakeADFun(nllHSMM, par3)

# Fit the model
opt3 = nlminb(obj3$par, obj3$fn, obj3$gr)

# reporting
mod3 = obj3$report()
mod3$par = opt3$par
mod3$Gamma = as.matrix(mod3$Gamma_sparse)

# saving fitted model object
saveRDS(mod3, file = "./simulation_study/Results/Part3_misspecification/HSMM.rds")


# Comparing dwell-time distributions --------------------------------------

# functions needed for overall dwell-time distribution
source("./functions/dwell_functions.R")

# load fitted models
mod1 = readRDS("./simulation_study/Results/Part3_misspecification/pHSMM.rds")
mod2 = readRDS("./simulation_study/Results/Part3_misspecification/pHMM.rds")
mod3 = readRDS("./simulation_study/Results/Part3_misspecification/HSMM.rds")

### overall dwell-time distributions

# pHSMM
mod1$pmf = lapply(1:3, function(j) ddwell_hsmm(1:agsizes[j], j, mod1$dm, mod1$Gamma, agsizes))

# pHMM
mod2$pmf = lapply(1:3, function(j) ddwell(1:agsizes[j], j, mod2$Gamma_p))

# HSMM
mod3$pmf = mod3$dm


### decoding states

# pHSMM
bigStates1 = viterbi_p(as.matrix(mod1$delta), mod1$Gamma, mod1$allprobs[,rep(1:3, times = agsizes)], data$tod[1:nrow(mod1$allprobs)])
mod1$states = rep(NA, nObs)
for(j in 1:3){
  mod1$states[which(bigStates1 %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
}

# pHMM
mod2$states = viterbi_p(mod2$delta, mod2$Gamma, mod2$allprobs, data$tod)

# HSMM
bigStates3 = viterbi(as.vector(mod3$delta), mod3$Gamma, mod3$allprobs[,rep(1:3, times = agsizes)])
mod3$states = rep(NA, nObs)
for(j in 1:3){
  mod3$states[which(bigStates3 %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
}

### "empirical" pmfs of decoded states

# pHSMM
mod1$emppmf = empirical_pmf(mod1$states, agsizes)

# pHMM
mod2$emppmf = empirical_pmf(mod2$states, agsizes)

# HSMM
mod3$emppmf = empirical_pmf(mod3$states, agsizes)



### plotting model-implied overall dwell-time distributions against viterbi decoded

modnames = c("periodic HSMM", "periodic HMM", "homogeneous HSMM")
mods = list(mod1, mod2, mod3)

# pdf("./figures/simulation/dwell-time_distribution.pdf", width = 7, height = 5)

N=3

m = matrix(c(1,1,1,2:10), nrow = 4, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.3, 1, 1, 1))
par(mar = c(0.5,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,10), ylim = c(1,10))

legend(x = 1.7, y = 9.5, inset = c(0.3,0), pch = c(NA, NA, 19), lwd = c(2, 1, NA), 
       lty = c(1, 1, NA), col = c("gray", "#00000030", "#00000030"),
       legend = c("model-implied", "true empirical", "Viterbi-decoded"), bty = "n", horiz = TRUE,
       text.width = c(2, 1.8, 0))


# color = c("orange", "deepskyblue", "seagreen2")
statenames = paste("state", 1:3)

emppmf_true = empirical_pmf(data$C, agsizes)

par(mar = c(5.5,4,2,1))

for(mod in 1:3){
  if(mod == 3){
    par(mar = c(5,4,3,0.5)+0.1)
  }
  for(j in 1:N){
    if(mod==1){
      main = statenames[j]
    } else{main = ""}
    if(mod == 3){
      xlab = "dwell time (hours)"
    } else{ xlab = ""}
    if(j == 1){
      ylab = "probabilities"
      par(mar = c(4,4,2,0.5)+0.1)
    } else{ylab = ""
    mar = c(2.5,3,2,0.5)+0.1
    }
    
    plot(1:agsizes[j], mods[[mod]]$pmf[[j]], type = "h", lwd = 2, col = color[j],
         ylim = c(0, 0.2), xlim = c(0,20),
         xlab = xlab, ylab = ylab, bty = "n", main = main)
    points(1:agsizes[j], mods[[mod]]$emppmf[[j]], pch = 19, col = "#00000030")
    lines(1:agsizes[j], emppmf_true[[j]], lwd = 1, lty = 1, col = "#00000030")
    
    if(j == 3){
      legend("topright", legend = modnames[mod], bty = "n")
    }
  }
}

# dev.off()



# Simulation 2: Decoding accuracies with confidence intervals -------------

## fit each model 500 time for accuracy CI
library(foreach)
library(doParallel)

nruns = 500

# simulate 100 data sets
set.seed(123)
Data = vector("list", nruns)
for(i in 1:nruns){
  print(i)
  Data[[i]] = sim_phsmm(1e4, beta, omega, obsparams)
} 


### Model 1: true periodic HSMM ###

# initial parameter
par1 = list(
  beta = beta, 
  logitomega = rep(0,3),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# Define data and hyperparams
dat = list(
  tod = Data[[1]]$tod,
  N = 3,
  agsizes = agsizes,
  Z = cbind(1, trigBasisExp(1:24))
)

# Initial obs object
obs = list(step = Data[[1]]$step, angle = Data[[1]]$angle)

# Create objective function once
obj = MakeADFun(nllpHSMM, par1, silent = TRUE)

# # Create a cluster of cores
# cl = makeCluster(detectCores())
# registerDoParallel(cl)
# 
# # Parallel execution using foreach
# acc1 = foreach(i = 1:nruns,
#                .packages = c("LaMa", "RTMB"), # export packages to workers
#                .errorhandling = "pass") %dopar% {
#                    
#                # Update the obs and dat environment for each iteration
#                assign("obs", list(step = Data[[i]]$step, angle = Data[[i]]$angle), envir = .GlobalEnv)
#                assign("dat", list(tod = Data[[i]]$tod, N = 3, agsizes = agsizes, Z = cbind(1, trigBasisExp(1:24))), envir = .GlobalEnv)
#                
#                # Fit the model, this will automatically use the new obs because DataEval() in nllpHSMM
#                opt = nlminb(obj$par, obj$fn, obj$gr)
#                
#                # Get the model report
#                mod = obj$report()
#                mod$Gamma = simplify2array(lapply(mod$Gamma_sparse, as.matrix))
#                
#                # decode states
#                rawstates_vit = viterbi_p(as.vector(mod$delta), mod$Gamma, mod$allprobs[,rep(1:3, times = agsizes)], Data[[i]]$tod)
#                probs = stateprobs_p(as.vector(mod$delta), mod$Gamma, mod$allprobs[,rep(1:3, times = agsizes)], Data[[i]]$tod)
#                rawstates_loc = apply(probs, 1, which.max)
#                if(is.list(rawstates_loc)){
#                  rawstates_loc[which(sapply(rawstates_loc, length) == 0)] = NA
#                  rawstates_loc = unlist(rawstates_loc)
#                }
#                states_vit = states_loc = rep(NA, length(rawstates_vit))
#                for(j in 1:3){
#                  states_vit[which(rawstates_vit %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
#                  states_loc[which(rawstates_loc %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
#                }
#                acc_vit = sum(states_vit == Data[[i]]$C) / length(states_vit)
#                acc_loc = sum(states_loc == Data[[i]]$C, na.rm = TRUE) / length(states_loc)
#                
#                list(acc_vit = acc_vit, acc_loc = acc_loc)
#                }
# stopCluster(cl)
# 
# saveRDS(acc1, "./simulation_study/Results/Part3_misspecification/acc_phsmm500.rds")


### Model 2: periodic HMM ###

# initial parameter
par2 = list(
  beta = matrix(c(rep(-2, 6), rep(0, 12)), nrow = 6), 
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# Define data and hyperparams
dat = list(
  tod = Data[[1]]$tod,
  N = 3,
  Z = cbind(1, trigBasisExp(1:24))
)

# Initial obs object
obs = list(step = Data[[1]]$step, angle = Data[[1]]$angle)

# Create objective function once
obj = MakeADFun(nllpHMM, par2, silent = FALSE)

# # Create a cluster of cores
# cl = makeCluster(detectCores())
# registerDoParallel(cl)
# 
# # Parallel execution using foreach
# acc2 = foreach(i = 1:nruns,
#                .packages = c("LaMa", "RTMB"), # export packages to workers
#                .errorhandling = "pass") %dopar% {
#                  
#                  # Update the obs and dat environment for each iteration
#                  assign("obs", list(step = Data[[i]]$step, angle = Data[[i]]$angle), envir = .GlobalEnv)
#                  assign("dat", list(tod = Data[[i]]$tod, N = 3, Z = cbind(1, trigBasisExp(1:24))), envir = .GlobalEnv)
#                  
#                  # Fit the model, this will automatically use the new obs because DataEval() in nllpHSMM
#                  opt = nlminb(obj$par, obj$fn, obj$gr)
#                  
#                  # Get the model report
#                  mod = obj$report()
#                  
#                  # decode states
#                  states_vit = viterbi_p(as.vector(mod$delta), mod$Gamma_p, mod$allprobs, Data[[i]]$tod)
#                  probs = stateprobs_p(as.vector(mod$delta), mod$Gamma_p, mod$allprobs, Data[[i]]$tod)
#                  states_loc = apply(probs, 1, which.max)
#                  if(is.list(states_loc)){
#                    states_loc[which(sapply(states_loc, length) == 0)] = NA
#                    states_loc = unlist(states_loc)
#                  }
#                  acc_vit = sum(states_vit == Data[[i]]$C) / length(states_vit)
#                  acc_loc = sum(states_loc == Data[[i]]$C, na.rm = TRUE) / length(states_loc)
#                  
#                  list(acc_vit = acc_vit, acc_loc = acc_loc)
#                }
# stopCluster(cl)
# 
# saveRDS(acc2, "./simulation_study/Results/Part3_misspecification/acc_phmm500.rds")


### Model 3: homogeneous HSMM ###

# initial parameter
par3 = list(
  loglambda = log(c(8,7,6)), 
  logitomega = rep(0,3),
  logmu = log(obsparams$mu), 
  logsigma = log(obsparams$sigma), 
  logkappa = log(obsparams$kappa)
)

# Define data and hyperparams
dat = list(
  N = 3,
  agsizes = agsizes
)

# Initial obs object
obs = list(step = Data[[1]]$step, angle = Data[[1]]$angle)

# Create objective function once
obj = MakeADFun(nllHSMM, par3, silent = TRUE)

# # Create a cluster of cores
# cl = makeCluster(detectCores())
# registerDoParallel(cl)
# 
# # Parallel execution using foreach
# acc3 = foreach(i = 1:nruns,
#                .packages = c("LaMa", "RTMB"), # export packages to workers
#                .errorhandling = "pass") %dopar% {
#                  
#                  # Update the obs and dat environment for each iteration
#                  assign("obs", list(step = Data[[i]]$step, angle = Data[[i]]$angle), envir = .GlobalEnv)
#                  assign("dat", list(N = 3, agsizes = agsizes), envir = .GlobalEnv)
#                  
#                  # Fit the model, this will automatically use the new obs because DataEval() in nllpHSMM
#                  opt = nlminb(obj$par, obj$fn, obj$gr)
#                  
#                  # Get the model report
#                  mod = obj$report()
#                  mod$Gamma = as.matrix(mod$Gamma_sparse)
# 
#                  # decode states
#                  rawstates_vit = viterbi(as.vector(mod$delta), mod$Gamma, mod$allprobs[,rep(1:3, times = agsizes)])
#                  probs = stateprobs(as.vector(mod$delta), mod$Gamma, mod$allprobs[,rep(1:3, times = agsizes)])
#                  rawstates_loc = apply(probs, 1, which.max)
#                  if(is.list(rawstates_loc)){
#                    rawstates_loc[which(sapply(rawstates_loc, length) == 0)] = NA
#                    rawstates_loc = unlist(rawstates_loc)
#                  }
#                  states_vit = states_loc = rep(NA, length(rawstates_vit))
#                  for(j in 1:3){
#                    states_vit[which(rawstates_vit %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
#                    states_loc[which(rawstates_loc %in% c(c(0,cumsum(agsizes)[-3])[j]+1:agsizes[j]))] = j
#                  }
#                  acc_vit = sum(states_vit == Data[[i]]$C) / length(states_vit)
#                  acc_loc = sum(states_loc == Data[[i]]$C, na.rm = TRUE) / length(states_loc)
#                  
#                  list(acc_vit = acc_vit, acc_loc = acc_loc)
#                }
# stopCluster(cl)
# 
# saveRDS(acc3, "./simulation_study/Results/Part3_misspecification/acc_hsmm500.rds")



# Results -----------------------------------------------------------------

acc_phsmm = readRDS("./simulation_study/Results/Part3_misspecification/acc_phsmm500.rds")
acc_phmm = readRDS("./simulation_study/Results/Part3_misspecification/acc_phmm500.rds")
acc_hsmm = readRDS("./simulation_study/Results/Part3_misspecification/acc_hsmm500.rds")


phsmm_acc = phmm_acc = hsmm_acc = data.frame(vit = numeric(nruns), loc = numeric(nruns))
for(i in 1:nruns){
  phsmm_acc$vit[i] = acc_phsmm[[i]]$acc_vit
  phsmm_acc$loc[i] = acc_phsmm[[i]]$acc_loc
  
  phmm_acc$vit[i] = acc_phmm[[i]]$acc_vit
  phmm_acc$loc[i] = acc_phmm[[i]]$acc_loc
  
  hsmm_acc$vit[i] = acc_hsmm[[i]]$acc_vit
  hsmm_acc$loc[i] = acc_hsmm[[i]]$acc_loc
}

phsmmCI_vit = quantile(phsmm_acc$vit, c(0.025, 0.975))
phmmCI_vit = quantile(phmm_acc$vit, c(0.025, 0.975))
hsmmCI_vit = quantile(hsmm_acc$vit, c(0.025, 0.975))


adjust = 1.5 # KDE bandwith adjustment

# pdf("./figures/simulation/accuracy.pdf", width = 7, height = 3.6)

m = matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.2, 1))
par(mar = c(0.5,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,10), ylim = c(0,1))

legend(x = 0.6, y = 1, lwd = 2, text.width = c(2, 1.8, 0),
       col = color,
       legend = c("periodic HSMM", "periodic HMM", "homogeneous HSMM"), bty = "n", horiz = TRUE)

par(mar = c(5,4,4,2) + 0.1)
plot(NA, xlim = c(0.94, 0.98), bty = "n", ylim = c(0, 150), 
     main = "Viterbi decoding", xlab = "accuracy", ylab = "density")
lines(density(phsmm_acc$vit, adjust = adjust), col = color[1], lwd = 2)
lines(density(phmm_acc$vit, adjust = adjust), col = color[2], lwd = 2)
lines(density(hsmm_acc$vit, adjust = adjust), col = color[3], lwd = 2)

plot(NA, xlim = c(0.94, 0.98), bty = "n", ylim = c(0, 150), 
     main = "Local decoding", xlab = "accuracy", ylab = "density")

lines(density(phsmm_acc$loc, adjust = adjust), col = color[1], lwd = 2)
lines(density(phmm_acc$loc, adjust = adjust), col = color[2], lwd = 2)
lines(density(hsmm_acc$loc, adjust = adjust), col = color[3], lwd = 2)

# dev.off()

## averages
# viterbi
round(mean(phsmm_acc$vit)*100, 2)
round(mean(phmm_acc$vit)*100, 2)
round(mean(hsmm_acc$vit)*100, 2)
# local
round(mean(phsmm_acc$loc)*100, 2)
round(mean(phmm_acc$loc)*100, 2)
round(mean(hsmm_acc$loc)*100, 2)

## confidence intervals
# viterbi
round(quantile(phsmm_acc$vit, c(0.025, 0.975))*100, 2)
round(quantile(phmm_acc$vit, c(0.025, 0.975))*100, 2)
round(quantile(hsmm_acc$vit, c(0.025, 0.975))*100, 2)
# local
round(quantile(phsmm_acc$loc, c(0.025, 0.975))*100, 2)
round(quantile(phmm_acc$loc, c(0.025, 0.975))*100, 2)
round(quantile(hsmm_acc$loc, c(0.025, 0.975))*100, 2)