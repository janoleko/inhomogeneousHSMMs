##### Part 4: Practical roadmap #####
## tradeoff between approximation and computation time


# Functions and libraries -------------------------------------------------

library(LaMa)

source("./functions/dwell_functions.R")


# Simulate data -----------------------------------------------------------

set.seed(1234)
n = 11000

# simulating AR1 covariate
phi = 0.9
zmean = 0
sigma = 0.1
z = rep(NA, n)
# stationary distribution of AR(1)
z[1] = rnorm(1, zmean, sqrt(sigma^2 / (1 - phi^2)))
for(t in 2:n) z[t] = zmean + phi * (z[t-1] - zmean) + rnorm(1, sd = 0.1)
Z = cbind(1, z) # design matrix

N = 2
omega = matrix(c(0,1,1,0), 2, 2)
mu = c(1, 5)
sigma = c(1, 2)

color = c("orange", "deepskyblue")
curve(dnorm(x, mu[1], sigma[1]), xlim = c(-3, 12), lwd = 2, col = color[1], bty = "n", ylab = "density", n = 500)
curve(dnorm(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = color[2], n = 500)

beta = matrix(c(log(2), 0,
                log(5), 2), nrow = 2, byrow = TRUE)
# relationship in state 2 -> maximum dwell time for x = 0
curve(exp(beta[2,2] + beta[2,2] * x), xlim = c(-1, 1), lwd = 2, col = color[2], n = 500, bty = "n")

Lambda = exp(Z %*% t(beta))
summary(Lambda[,2])
# most of the time rather short state durations, but sometimes bursts close to 20

n_tilde = ceiling(n / (1 + mean(apply(Lambda, 1, mean))) * 5) # crude estimate of distinct # of stays
s = rep(NA, n_tilde)
C = x = rep(NA, 3*n)
s[1] = 1
for(t in 2:n_tilde) s[t] = sample(1:N, size=1, prob = omega[s[t-1],])
times = rpois(1, lambda = Lambda[1, s[1]]) + 1
C[1:times] = rep(s[1], times)
x[1:times] = rnorm(times, mu[s[1]], sigma[s[1]])
currentInd = times
t = 2
while(currentInd <= n){
  times = rpois(1, lambda = Lambda[currentInd + 1, s[t]]) + 1
  C[currentInd + 1:times] = rep(s[t], times)
  x[currentInd + 1:times] = rnorm(times, mu[s[t]], sigma[s[t]])
  
  currentInd = currentInd + times
  t = t+1
}
x = x[1:n]
z = z[1:n]
C = C[1:n]

emp_pmf = empirical_pmf(C, agsizes = c(10, 50))

pdf("./figures/simulation/part4_empirical_pmf.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(1:10, as.numeric(emp_pmf[[1]]), type = "h", bty = "n", col = color[1], lwd = 2,
     ylab = "probabilities", xlab = "dwell time", main = "state 1", xlim = c(0,10))
plot(1:25, as.numeric(emp_pmf[[2]][1:25]), type = "h", bty = "n", col = color[2], lwd = 2,
     ylab = "probabilities", xlab = "dwell time", main = "state 2", xlim = c(0,25))
dev.off()


# Writing likelihood functions --------------------------------------------

nllHMM = function(par){
  getAll(par, dat)
  # parameter transformations
  sigma = exp(logsigma); REPORT(sigma)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  # evaluating state-dependent densities
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # forward algorithm
  -forward(delta, Gamma, allprobs, ad = TRUE)
}

nlliHSMM = function(par){
  getAll(par, dat)
  # parameter transformations
  sigma = exp(logsigma); REPORT(sigma)
  # calculating mean dwell times
  Lambda = exp(Z %*% t(beta)); REPORT(Lambda) # dwell time means
  # dwell-time pmfs and embedded tpm
  dm = lapply(1:N, function(i) sapply(1:agsizes[i]-1, dpois, lambda = Lambda[,i])); REPORT(dm)
  omega = if (N == 2) tpm_emb() else tpm_emb(logitomega) # case distinction: 2 states -> no parameters to vary
  # evaluating state-dependent densities
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # forward algorithm for inhomogeneous HSMMs
  -forward_ihsmm(dm, omega, allprobs, eps = 1e-20, 
                 startInd = startInd) # setting startInd to make models with different aggregate sizes comparible.
}



# Fitting homogeneous HMM -------------------------------------------------

par1 = list(
  mu = mu,
  logsigma = log(sigma),
  eta = rep(-2, 2)
)

dat = list(x = x, N = 2)

obj1 = MakeADFun(nllHMM, par1)
opt1 = nlminb(obj1$par, obj1$fn, obj1$gr)
mod1 = obj1$report()

# model-implied dwell-time distributions
par(mfrow = c(1,2))
plot(1:20, dgeom(1:20, 1-diag(mod1$Gamma)[1]), type = "h", bty = "n", xlab = "dwell time", ylab = "probabilities")
plot(1:40, dgeom(1:40, 1-diag(mod1$Gamma)[2]), type = "h", bty = "n", xlab = "dwell time", ylab = "probabilities")

# model-implied dwell-time distribution's 99%-quantiles
qgeom(0.99, prob = 1-diag(mod1$Gamma))

# decoding the states
mod1$states = viterbi(mod1$delta, mod1$Gamma, mod1$allprobs)

# computing the empircal dwell-time distribution of decoded states
emp_pmf1 = empirical_pmf(mod1$states, c(20, 40))

# largest non-zero entries
sapply(emp_pmf1, function(x) max(which(x > 0)))


# Fitting inhomogeneous HSMM ----------------------------------------------

# aggregate size grid for state 2

time = llks = pars = list()

# initial HSMM parameter
par2 = list(
  mu = mu, 
  logsigma = log(sigma),
  beta = matrix(c(log(3), log(5), 0, 2), nrow = 2)
)

# I will now fix the first aggregate size at 10 and vary the second one
# from 15 to 55
dat = list(x = x, N = 2, Z = Z, agsizes = c(10, 15), 
           startInd = 60) # fixing startInd at something larger than the largest aggregate size to make models comparible

for(i in 0:40){
  dat$agsizes[2] = 15 + i
  
  obj = MakeADFun(nlliHSMM, par2)
  
  time[[i+1]] = system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(rel.tol = 1e-15))
  )[3]
  
  print(dat$agsizes)
  
  llks[[i+1]] = -obj$fn(opt$par)
  print(-llks[[i+1]])
  pars[[i+1]] = opt$par
}

time = unlist(time)
llks = unlist(llks)
pars = do.call(rbind, pars)


# pdf("./figures/simulation/agsize_est_time.pdf", width = 8, height = 3.5)
par(mfrow = c(1,2), mar = c(4.5,4,3,1) + 0.1)
plot(15:55, time, type = "l", lwd = 2, xlab = "aggregate size in state 2", 
     ylab = "estimation time (seconds)", bty = "n", xaxt = "n")
axis(1, at = seq(15, 55, by = 10), labels = seq(15, 55, by = 10))
abline(lm(time ~ x, data = data.frame(time = time, x = 15:55)), col = "orange", lty = 2)
mod = lm(time ~ x + I(x^2), data = data.frame(time = time, x = 15:55))
lines(seq(15,55, length = 200), predict(mod, data.frame(x = seq(15,55, length = 200))), col = "deepskyblue", lty = 2)
legend("top", legend = c("linear", "quadratic"), col = color, lty = 2, bty = "n")

plot(15:55, llks, type = "l", lwd = 2, xlab = "aggregate size in state 2", 
     ylab = "log-likelihood", bty = "n", xaxt = "n")
axis(1, at = seq(15, 55, by = 10), labels = seq(15, 55, by = 10))
# dev.off()



# Visualising the dwell-time distribution ---------------------------------

# refit the model with agsizes c(10, 40)

agsizes = c(10, 40)

dat = list(x = x, 
           N = 2, 
           Z = Z, 
           agsizes = agsizes, 
           startInd = NULL)

# creating the objective function
obj = MakeADFun(nlliHSMM, par2)
opt = nlminb(obj$par, obj$fn, obj$gr)

# estimated parameter as list
par = as.list(sdreport(obj), "Estimate")

# reporting
mod = obj$report()

# make sparse objects dense for viterbi
delta = as.vector(mod$delta_sparse)
Gamma = simplify2array(lapply(mod$Gamma_sparse, as.matrix))

# decode states
mod$states = viterbi_g(delta, Gamma[,,-dim(Gamma)[3]], mod$allprobs[-(1:39), rep(1:2, times = agsizes)])

# discarding first covariate values
z2 = z[-(1:39)]

cov = data.frame(z = z2, state = mod$states)

# computing quantiles of covariate
quants = quantile(
  subset(cov, state %in% (agsizes[1] + 1:agsizes[2]))$z, 
  probs = c(seq(0, 1, length = 6))
  )
names(quants) = c("min", "20 %", "40 %", "60 %", "80 %", "max")


# pdf("./figures/simulation/pmf_quants.pdf", width = 7, height = 4)
par(mfrow = c(2,3))
for(i in 1:length(quants)){
  pmft = dpois(1:40-1, lambda = exp(par$beta[2,1] + par$beta[2,2] * quants[i]))
  plot(1:40, pmft, type = "h", lwd = 1.5, bty = "n", col = color[2], xlim = c(0, 40), 
       main = names(quants)[i], xlab = "dwell time", ylab = "probabilities")
}
#dev.off()
