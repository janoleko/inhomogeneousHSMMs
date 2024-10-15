
# functions
source("./functions/sim_study_functions.R")

# libraries
library(parallel)
library(LaMa)


# Part 2: aggregate sizes -------------------------------------------------

# function for applying
one_rep = function(data, beta, omega, stateparams, agsizes, stepmax){
  if(is.data.frame(data)){
    mod = fitpHSMM(data=data, beta=beta, stateparams=stateparams, 
                   agsizes=agsizes, stepmax=stepmax)
    
    return(mod)
  } else{return("no data")}
}

# parameter for sim model
beta = matrix(c(log(c(8,7,6)), -0.2, 0.2, -0.6, 0.3, -0.2, -0.4), nrow = 3)
omega = matrix(c(0, 0.7, 0.3, 
                 0.2, 0, 0.8,
                 0.5, 0.5, 0), nrow = 3, byrow = TRUE)
stateparams = list(
  mu = c(20, 200, 800),
  sigma = c(20, 150, 500),
  kappa = c(0.2, 1, 2.5)
)

# trigonometric basis expansion
Z = cbind(1, LaMa::trigBasisExp(1:24, 24))
# calculating all dwell time means
dM = exp(Z%*%t(beta))

# get maximum mean for each state
maxlambda = apply(dM, 2, max)

# compute aggregate sizes relativ to maximum quantile for each state
factors = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)
Agsizes = matrix(NA, nrow = length(factors), ncol = 3)
for(k in 1:length(factors)){
  Agsizes[k,] = ceiling((qpois(0.99, maxlambda)+1) * factors[k]) 
}


# Simulation --------------------------------------------------------------

# number of data sets for each scenario
nruns = 500

Data = list()
set.seed(123)
# simulating 500 data sets of length T = 5000
# for(i in 1:nruns){
#   cat("\n", i)
#   Data[[i]] = sim_phsmm(5000, beta, omega, stateparams)
# }

# fitting HSMMs with increasing aggregate sizes (as defined by factor)
# for(k in 1:length(factors)){
#   cat("\nScenario", k)
#   cat("\nAggregate factor:", factors[k])
#   
#   results = parallel::mclapply(Data, FUN = one_rep,
#                                beta=beta, omega=omega, stateparams=stateparams, 
#                                agsizes = Agsizes[k,], stepmax = 5,
#                                mc.cores = parallel::detectCores()-1)
#                                #mc.cores = 1)
#   saveRDS(results, 
#           file = paste0("./simulation_study/simulation_results/aggregate_size/results_", factors[k], ".rds"))
# }


# Analyzing results -------------------------------------------------------

# defining colors
color = c("orange", "deepskyblue", "seagreen2")


## Dwell time mean coefficients

Betas = array(dim = c(3, 3, 500, 9))

for(k in 1:length(factors)){
  res = readRDS(paste0("./simulation_study/simulation_results/aggregate_size/results_", factors[k], ".rds"))
  for(i in 1:nruns){
    if(!is.character(res[[i]])){
      Betas[,,i,k] = res[[i]]$beta
    }
  }
}


# for(state in 1:3){
#   pdf(file = paste0("./figures/simulation/aggregate_size", state, ".pdf"), width=8, height=5)
#   par(mfrow = c(3,5), mar = c(1.3,4.3,2.5,0.5)+0.1)
#   ylims = apply(Betas[state,,,9], 1, quantile, probs = c(0.0005, 0.9995), na.rm = TRUE)
#   for(p in 1:3){
#     for(k_bar in 1:5){
#       k = k_bar*2-1
#       if(p==1){
#         main = paste0("factor = ", factors[k])
#       } else{ main = "" }
#       if(k==1){
#         ylab = bquote({beta^(.(state))} [.(p-1)])
#       } else{ ylab = ""}
#       boxplot(Betas[state, p,,k], ylim = ylims[,p], ylab = ylab, 
#               col = "gray95", main=main, frame = FALSE)
#       abline(h = beta[state,p], col = color[state], lwd = 2)
#       # abline(h = mean(Betas[state, p,,k], na.rm = TRUE), col = "deepskyblue")
#     }
#   }
#   dev.off()
# }


# plotting distribution of dwell-time parameters for different aggregate sizes
for(state in 1:3){
  #pdf(file = paste0("./figures/simulation/aggregate_size_state", state, ".pdf"), width=7.5, height=2.5)
  par(mfrow = c(1,3), mar = c(4,4.5,1,2)+0.1)
  for(p in 1:3){
      ylims = apply(Betas[state,p,,], 2, range, na.rm = TRUE)
      ylim = c(min(ylims[1,]), max(ylims[2,]))
      # ylim = quantile(Betas[state,p,,9], probs = c(0.0001, 0.9999), na.rm = TRUE)
      b = Betas[state, p,,]
      B = data.frame(beta = as.numeric(b), factor = rep(factors, each=500))
      B$factor = as.factor(B$factor)
      B = subset(B, factor %in% c(0.5, 0.6, 0.7, 0.9, 1.3))
      B$factor = factor(B$factor, levels = c(0.5, 0.6, 0.7, 0.9, 1.3))
      boxplot(B$beta~B$factor, ylim = ylim, 
              ylab = bquote({beta^(.(state))} [.(p-1)]), xlab = "factor",
              col = "gray95", main="", frame = FALSE)
      abline(h = beta[state,p], col = scales::alpha(color[state],0.8), lwd = 1.5)
      # abline(h = mean(Betas[state, p,,k], na.rm = TRUE), col = "deepskyblue")
  }
  # dev.off()
}

# severe bias for 0.5
# 0.7 already okay

## histograms
for(state in 1:3){
  pdf(file = paste0("./figures/simulation/aggregate_size_state", state, "_hist.pdf"), width=10, height=7)
  par(mfrow = c(3,5), mar = c(4.5,4,3.5,1.1)+0.1)
  xlims = apply(Betas[state,,,1], 1, quantile, probs = c(0.0025, 0.9975), na.rm = TRUE)
  for(p in 1:3){
    for(k in c(1:3, 5, 9)){
      if(p==1){
        main = paste0("factor = ", factors[k])
      } else{ 
        main = "" 
      }
      if(k==1){
        xlab = bquote({beta^(.(state))} [.(p-1)])
        ylab = "density"
      } else{ ylab = ""}
      hist(Betas[state, p,, k], xlab = xlab, bor = "white", main = main, prob = TRUE, ylab = ylab,
           xlim = xlims[,p], breaks = seq(xlims[1,p]-0.5, xlims[2,p]+0.5, length=70))
      lines(density(Betas[state, p,, k]), col = color[state], lwd = 1, lty = 2)
      abline(v = beta[state,p], col = color[state], lwd = 2)
    }
  }
  dev.off()
}




## Likelihoods
LLks = Time = matrix(NA, 500, 9)

for(k in 1:length(factors)){
  res = readRDS(paste0("./simulation_study/simulation_results/aggregate_size/results_", factors[k], ".rds"))
  for(i in 1:nruns){
    if(!is.character(res[[i]])){
      LLks[i,k] = res[[i]]$llk
      Time[i,k] = res[[i]]$time
    }
  }
}

meanLLks = apply(LLks, 2, mean)

pdf("./figures/simulation/llks.pdf", width=6, height=4)
par(mfrow = c(1,1), mar = c(5,4,4,2)+0.1)
plot(factors, meanLLks, ylim = c(-39260, -39248), xaxt = "n",
     xlab = "factor", ylab = "log-likelihood", bty = "n", type = "b", lwd = 2)
axis(1, at = seq(0.5, 1.3, by = 0.2), labels = seq(0.5, 1.3, by = 0.2))
dev.off()

## same but with boxplots
# llks = data.frame(llks = as.numeric(LLks), factors = rep(factors, each = 500))
# llks$factors = factor(llks$factors, levels = factors)
# boxplot(llks$llks~llks$factors)


Time[Time<10] = Time[Time<10]*60
meanTime = apply(Time, 2, mean)

plot(factors, meanTime, xaxt = "n",
     xlab = "factor", ylab = "estimation time (seconds)", bty = "n", type = "b", lwd = 2)
axis(1, at = seq(0.5, 1.3, by = 0.2), labels = seq(0.5, 1.3, by = 0.2))


