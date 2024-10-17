
# Data --------------------------------------------------------------------

# install.packages(moveHMM)
library(moveHMM) # prepping data
# install.packages(PHSMM)
library(PHSMM) # data set

# prepping data -> calculating step lengths and turning angles from coordinates
muskox = moveHMM::prepData(muskox, type = "UTM")
muskox$tod = muskox$tday+1 # time of day variable in 1:24
muskox = muskox[,-c(1,7,8)] # excluding unnecessary columns

head(muskox)

(n = nrow(muskox)) # number of observations
sum(is.na(muskox$step)) # number of missing step lengths
sum(muskox$step[!is.na(muskox$step)]==0) # no need for zero-inflated gamma distribution
sum(is.na(muskox$angle)) # number of missing turning angles



# Defining colors ---------------------------------------------------------

color = c("orange", "deepskyblue", "seagreen2")


# EDA ---------------------------------------------------------------------

library(sf) # spatial data
library(ggplot2) # plotting
# install.packages("rnaturalearth")
library(rnaturalearth) # map plotting
# install.packages("patchwork")
library(patchwork) # multiple plots as one


### first plot: track on map ###
# Create an sf object from UTM data
utm_sf = st_as_sf(muskox[,c("x", "y")], 
                   coords = c("x", "y"), 
                   crs = paste0("+proj=utm +zone=", 27, " +datum=WGS84"),
                   na.fail = FALSE)
# Transform the UTM coordinates to latitude and longitude (EPSG: 4326 is WGS84)
latlon_sf = st_transform(utm_sf, crs = 4326)
# Extract the latitude and longitude columns
coords = data.frame(st_coordinates(latlon_sf))
colnames(coords) = c("longitude", "latitude")
coords$longitude[is.nan(coords$longitude)] = NA
coords$latitude[is.nan(coords$latitude)] = NA

# Create a GPS points spatial object
gps_data <- data.frame(lon = coords$longitude, lat = coords$latitude)
gps_sf <- st_as_sf(gps_data, coords = c("lon", "lat"), crs = 4326, na.fail = FALSE)

# Get world map data from rnaturalearth
world <- ne_countries(scale = "large", returnclass = "sf")

xlim = range(coords$longitude, na.rm = T) + 0.5*c(-1, 1)
ylim = range(coords$latitude, na.rm = T) + 0.2*c(-1, 1)

## region plot showing greenland
plot_coarse = ggplot() +
  geom_sf(data = world, fill = "gray95") +
  # geom_sf(data = gps_sf[1,], color = "black", size = 5) +
  geom_rect(aes(xmin = xlim[1], xmax = xlim[2],                # Rectangle based on xlim
                ymin = ylim[1], ymax = ylim[2]), 
            alpha = 0.7, fill = "deepskyblue") +  
  coord_sf(xlim = c(-55, -12), ylim = c(60,78), expand = FALSE) +
  theme_minimal()
## finer plot showing track

plot_fine = ggplot() +
  geom_sf(data = world, fill = "gray95") +
  # geom_sf(data = gps_line, color = "black", size = .1) +  # Add the connecting line
  geom_sf(data = gps_sf, color = "black", size = .3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  # fixed_plot_aspect(ratio = 1) +
  theme_minimal()# +

# Combine plots side by side
combined_plot <- plot_coarse + plot_fine + theme(plot.margin = margin(0, 0, 0, 0))

# ggsave("./figures/maps.pdf", combined_plot, width = 10, height = 5, units = "in")


## step length and turning angle time series

# pdf("./figures/muskox_step.pdf", width = 8, height = 4)
par(mar = c(4.5,4,2,2))
plot(muskox$step[1:1000], type = "h", bty = "n", ylab = "step length (meters)", xlab = "time")
# dev.off()

# plot(muskox$angle, type = "h", bty = "n", ylab = "turning angle (radians)")

# pdf("./figures/muskox_hist.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
hist(muskox$step, prob = T, main = "", bor = "white", xlab = "step length (meters)", ylab = "density", 
     xlim = c(0,500), breaks = 100, ylim = c(0, 0.01))
hist(muskox$angle, prob = T, breaks = 25, main = "", bor = "white", xlab = "turning angle (radians)", ylab = "density", xaxt = "n")
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), 
                expression(0), expression(pi/2), expression(pi)))
# dev.off()

# mean hourly step lengths
meanstep = sapply(1:24, function(t) mean(muskox$step[which(muskox$tod==t)], na.rm=T))

# pdf("./figures/muskox_boxplot.pdf", width = 6, height = 3)
par(mfrow = c(1,1), mar = c(4.5,4,1.5,2)+0.1)
boxplot(muskox$step ~ muskox$tod, pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000010", 
        frame = F, xlab = "time of day", ylab = "step length (meters)", ylim = c(0,400))
points(1:24, meanstep, type = "o", col = color[2], pch = 19)
# lines(1:24, means, type = "b", col = color[2])
legend(x = 18, y = 420, lwd = c(2, 1), pch = c(NA, 19), col = c("black", color[2]),
       legend = c("median", "mean"), bty = "n")
# dev.off()



### Model fitting -----------------------------------------------------------

# Loading functions and packages ------------------------------------------

source("./application/likelihoods.R") # contains the 4 different likelihood functions used
source("./functions/auxiliary_functions.R")

# devtools::install_github("janoleko/LaMa")
library(LaMa) # HMM and HSMM functions
# install.packages("RTMB")
library(RTMB) # automatic differentiation


# Fitting HMMs ------------------------------------------------------------

N = 3 # all models will have 3 states
estTime = list()

### homogeneous ###

## K = 0
# initial parameter values
par1 = list(eta = rep(-2, 6), # state process intercepts
            logmu = log(c(5, 80, 400)), logsigma = log(c(4, 70, 400)), # state-dependent process step
            mu.turn = c(pi, 0, 0), logkappa = log(c(0.4, 0.5, 1.5))) # state-dependent process angle
# data and hyperparameters
dat = list(step = muskox$step, angle = muskox$angle, N = 3)

# map argument to fix mu.turn at initial value c(pi, 0, 0)
map = list(mu.turn = factor(rep(NA, 3)))

# creating the AD objective function
obj1 = MakeADFun(nllHMM, par1, 
                 map = map) # fix turning angle mean
estTime[[1]] = system.time(
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr) # model fitting
)
mod1 = obj1$report() # reporting
(IC1 = AIC_BIC(opt1, n)) # computing information criteria


### inhomogeneous ###

## K = 1
par2 = list(beta = matrix(c(rep(-2, 6), rep(0, 2*6)), nrow = 6), # now matrix (intercepts + slopes)
            logmu = log(c(5, 80, 400)), logsigma = log(c(4, 70, 400)),
            mu.turn = c(pi, 0, 0), logkappa = log(c(0.4, 0.5, 1.5)))
dat = list(step = muskox$step, angle = muskox$angle, N = 3, 
           Z = trigBasisExp(1:24, degree = 1), tod = muskox$tod)

obj2 = MakeADFun(nllpHMM, par2, map = map)
estTime[[2]] = system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
)
mod2 = obj2$report()
(IC2 = AIC_BIC(opt2, n))


## K = 2
par3 = par2
par3$beta = matrix(c(rep(-2, 6), rep(0, 4*6)), nrow = 6)
dat$Z = trigBasisExp(1:24, degree = 2)

obj3 = MakeADFun(nllpHMM, par3, map = map)
estTime[[3]] = system.time(
  opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
)
mod3 = obj3$report()
(IC3 = AIC_BIC(opt3, n))


# Fitting HSMMs -----------------------------------------------------------

### homogeneous ###
## K = c(0,0,0)
# getting an intuition on dwell time starting values
1/(1-diag(mod1$Gamma)) # mean (geometric) dwell times of the homogeneous HMM
# intuition for aggreagte sizes:
qexp(0.95, 1-diag(mod1$Gamma)) # 95 % quantile of geometric distribution

par4 = list(logmu_dwell = log(rep(3,3)), logphi_dwell = log(rep(0.2, 3)), # dwell-time mean and dispersion
            logitomega = rep(0,3),                                        # embedded tpm
            logmu = log(c(4, 80, 400)), logsigma = log(c(4, 70, 400)),
            mu.turn = c(pi, 0, 0), logkappa = log(c(0.4, 0.3, 2)))

dat = list(step = muskox$step, angle = muskox$angle, N = 3, 
           agsizes = rep(20, 3)) # state aggregate sizes used for the approximation

obj4 = MakeADFun(nllHSMM, par4, map = map)
estTime[[4]] = system.time(
  opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr, control = list(iter.max = 1e3))
)
mod4 = obj4$report()
(IC4 = AIC_BIC(opt4, n))

par4 = as.list(sdreport(obj4), "Estimate") # get estimated parameter in list format


### inhomogeneous ###

## K = c(1,0,0)
par5 = list(betaMu = matrix(c(par4$logmu_dwell, rep(0, 2*3)), nrow = 3), # mean matrix
            betaPhi = par4$logphi_dwell, beta_omega = par4$logitomega,
            logmu = log(c(4, 80, 400)), logsigma = log(c(4, 70, 400)),
            mu.turn = c(pi, 0, 0), logkappa = log(c(0.4, 0.3, 2)))

dat = list(step = muskox$step, angle = muskox$angle, tod = muskox$tod, N = 3, agsizes = rep(30, 3), # increased aggregate sizes
           Z_Mu = cbind(1, trigBasisExp(1:24, degree = 1)), # dwell-time mean design matrix
           Z_Phi = rep(1, 24), Z_omega = rep(1, 24)) # dwell-time dispersion and embedded tpm design matrix are constant

obj5 = MakeADFun(nllpHSMM, par5, map = map)
estTime[[5]] = system.time(
  opt5 <- nlminb(obj5$par, obj5$fn, obj5$gr, control = list(iter.max = 1e3))
)
mod5 = obj5$report()
(IC5 = AIC_BIC(opt5, n))

par5 = as.list(sdreport(obj5), "Estimate")


## K = c(2,0,0)
par6 = par5 # initialise with previous values
par6$betaMu = cbind(par5$betaMu, 0, 0) # add zero columns for new sine and cosine frequencies

dat$Z_Mu = cbind(1, trigBasisExp(1:24, degree = 2)) # adding 2 columns with additional sine and cosine

obj6 = MakeADFun(nllpHSMM, par6, map = map)
estTime[[6]] = system.time(
  opt6 <- nlminb(obj6$par, obj6$fn, obj6$gr, 
                control = list(iter.max = 1e3, rel.tol = 1e-15))
)
mod6 = obj6$report()
(IC6 = AIC_BIC(opt6, n))

par6 = as.list(sdreport(obj6), "Estimate")


## K = c(1,1,0): inhomogeneity in dispersion parameter
par7 = par5 # parameters from K = c(1,0,0) model
par7$betaPhi = cbind(par5$betaPhi, 0, 0) # adding zeros for dispersion slopes

dat = list(step = muskox$step, angle = muskox$angle, tod = muskox$tod, N = 3, agsizes = rep(35,3),
           Z_Mu = cbind(1, trigBasisExp(1:24, degree = 1)), # dwell-time mean design matrix
           Z_Phi = cbind(1, trigBasisExp(1:24, degree = 1)), Z_omega = rep(1, 24)) # dwell-time dispersion design matrix

obj7 = MakeADFun(nllpHSMM, par7, map = map)
estTime[[7]] = system.time(
  opt7 <- nlminb(obj7$par, obj7$fn, obj7$gr, 
                control = list(iter.max = 1e4, rel.tol = 1e-15, step.max = 1e-1))
)
# gradient does not fully go to zero -> models with time-varying dispersion are rather unstable
mod7 = obj7$report()
(IC7 = AIC_BIC(opt7, n))

## visualising model 7 a bit

par(mfrow = c(1,2))
# visualizing time-varying mean
plot(mod7$Mu_dwell[,1]+1, type = "b", lwd = 2, pch = 19, col = color[1], bty = "n", 
     ylim = c(2,9), xlim = c(0,24), main = "mean",
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
for(j in 2:3) lines(mod7$Mu_dwell[,j]+1, type = "b", lwd = 2, pch = 19, col = color[j])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))

# visualizing time-varying dispersion
plot(mod7$Phi_dwell[,1], type = "b", lwd = 2, pch = 19, col = color[1], bty = "n", 
     ylim = c(0,2), xlim = c(0,24), main = "dispersion",
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
for(j in 2:3) lines(mod7$Phi_dwell[,j], type = "b", lwd = 2, pch = 19, col = color[j])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))


### overall dwell-time distributino
source("./functions/dwell_functions.R")
# turn list of sparse tpms into dense array
mod7$Gamma = simplify2array(lapply(mod7$Gamma_sparse, as.matrix))

# compute overall dwell-time distribution
pmf7 = lapply(1:3, function(j) ddwell_hsmm(1:30, j, mod7$dm, mod7$Gamma, rep(30,3)))

par(mfrow = c(1,3))
plot(1:15, pmf7[[1]][1:15], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "resting")
plot(1:15, pmf7[[2]][1:15], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "foraging")
plot(1:15, pmf7[[3]][1:15], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "travelling")

# time-varying dwell time mean and standard deviation
par(mfrow = c(1,2))
todseq = seq(0,24,length=200)
Z = trigBasisExp(todseq, 24, 1)
Museq = exp(cbind(1,Z)%*%t(par7$betaMu))+1
plot(todseq, Museq[,1], ylim = c(0,10), type = "l", lwd = 2, col = color[1],
     bty = "n", xaxt = "n", xlab = "time of day", ylab = "mean dwell time")
for(j in 2:3) lines(todseq, Museq[,j], lwd = 2, col = color[j])
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))

Phiseq = exp(cbind(1,Z)%*%t(par7$betaPhi))
Varseq = Museq+(Museq^2)*Phiseq
Sdseq = sqrt(Varseq)
plot(todseq, Sdseq[,1], ylim = c(0,15), type = "l", lwd = 2, col = color[1],
     bty = "n", xaxt = "n", xlab = "time of day", ylab = "dwell time standard deviation")
for(j in 2:3) lines(todseq, Sdseq[,j], lwd = 2, col = color[j])
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))


## K = c(1,0,1): modelling the conditional tpm
par8 = par5
par8$beta_omega = cbind(par5$beta_omega, 0, 0)

dat = list(step = muskox$step, angle = muskox$angle, tod = muskox$tod, N = 3, agsizes = rep(30, 3),
           Z_Mu = cbind(1, trigBasisExp(1:24, degree = 1)), # dwell-time mean design matrix
           Z_Phi = rep(1, 24), # dwell-time dispersion design matrix (constant)
           Z_omega = cbind(1, trigBasisExp(1:24, degree = 1))) # embedded tpm design matrix

obj8 = MakeADFun(nllpHSMM, par8, map = map)
estTime[[8]] = system.time(
  opt8 <- nlminb(obj8$par, obj8$fn, obj8$gr, 
                control = list(iter.max = 1e4))
)
mod8 = obj8$report()
(IC8 = AIC_BIC(opt8, n))

par8 = as.list(sdreport(obj8), "Estimate")

# conditional tpm Omega as function of the time of day
par(mfrow = c(3,3))
for(i in 1:N){ for(j in 1:N){
    plot(mod8$omega[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", 
         ylab = expression(omega), xlab = "time of day")
}}
# no strong effects here


## K = c(1,1,1): time-varying mean, dispersion and embedded tpm
par9 = par8
par9$betaPhi = par7$betaPhi

dat$Z_Phi = cbind(1, trigBasisExp(1:24, degree = 1))
dat$agsizes = rep(35, 3)

obj9 = MakeADFun(nllpHSMM, par9, map = map)
estTime[[9]] = system.time(
  opt9 <- nlminb(obj9$par, obj9$fn, obj9$gr,
                control = list(iter.max = 1e4, step.max = 1e-1))
)
# gradient does not fully go to zero -> models with time-varying dispersion are rather unstable
mod9 = obj9$report()
(IC9 = AIC_BIC(opt9, n))

# conditional tpm Omega as function of the time of day again
par(mfrow = c(3,3))
for(i in 1:N){ for(j in 1:N){
    plot(mod9$omega[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", 
         ylab = expression(omega), xlab = "time of day")
}}
# again, no strong effects


# Information criteria ----------------------------------------------------

ICs = cbind(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
colnames(ICs) = 1:9
tab = rbind(round(c(-opt1$objective, -opt2$objective, -opt3$objective, -opt4$objective,
              -opt5$objective, -opt6$objective, -opt7$objective, -opt8$objective, -opt9$objective), 1),
            round(ICs,1),
            round(sapply(1:9, function(i) estTime[[i]][3]),2))
rownames(tab) = c("llk", "AIC", "BIC", "time")
t(tab)


### Visualising the results from models 2, 4, and 5 -------------------------

mod5$Gamma = simplify2array(lapply(mod5$Gamma_sparse, as.matrix)) # dense tpm array


# State-dependent distributions of model 5 --------------------------------

# aggregating stationary distribution for weighting of component distributions in plot
Delta5_tilde = stationary_p(mod5$Gamma)
agsizes = rep(30,3)
Delta5 = sapply(1:N, function(j) rowSums(Delta5_tilde[,c(0,cumsum(agsizes))[j]+1:agsizes[j]]))
delta5 = colMeans(Delta5) # overall state distribution


# pdf("./figures/marginal.pdf", width = 7, height = 3.5)
m = matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.15, 1, 1))
par(mar = c(0,2,0.5,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,10), ylim = c(1,10))

legend(x = 1.4, y = 9.5, inset = c(0.3,0), lwd = 2, lty = c(rep(1,3), 2),
       col = c(color, "black"),
       legend = c("resting", "foraging", "travelling", "marginal"), bty = "n", 
       horiz = TRUE, text.width = c(1, 1, 1.2, 1))

par(mar = c(5,4,2.2,2)+0.1)

hist(muskox$step, prob = T, breaks = 100, bor = "white", xlim = c(0,600), ylim = c(0,0.01),
     main = "", xlab = "step length", ylab = "density")
for(j in 1:N) curve(delta5[j]*dgamma2(x, mod5$mu[j], mod5$sigma[j]), add = T, lwd = 2, col = color[j], n = 500)
curve(delta5[1]*dgamma2(x, mod5$mu[1], mod5$sigma[1])+
        delta5[2]*dgamma2(x, mod5$mu[2], mod5$sigma[2])+
        delta5[3]*dgamma2(x, mod5$mu[3], mod5$sigma[3]), add = T, lwd = 2, lty = 2)

hist(muskox$angle, prob = T, breaks = 20, bor = "white", xlim = c(-pi,pi), ylim = c(0,0.3),
     main = "", xlab = "turning angle", ylab = "density", xaxt = "n")
for(j in 1:N) curve(delta5[1]*dvm(x, mod5$mu.turn[j], mod5$kappa[j]), add = T, lwd = 2, col = color[j], n = 500)
curve(delta5[1]*dvm(x, mod5$mu.turn[1], mod5$kappa[1])+
        delta5[2]*dvm(x, mod5$mu.turn[2], mod5$kappa[2])+
        delta5[3]*dvm(x, mod5$mu.turn[3], mod5$kappa[3]), add = T, lwd = 2, lty = 2, n = 500)
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), 
                expression(0), expression(pi/2), expression(pi)))
# dev.off()



# Periodically stationary distribution with confidence bands --------------

library(mvtnorm)

# get skeleton of parlist excluding fixpars to convert sampled parvectors back to named lists
skeletonHSMM = as.relistable(par5[names(par5) != "mu.turn"])
skeletonHMM = as.relistable(par2[names(par2) != "mu.turn"])

# compute hessians for pHMM and pHSMM
H2 = obj2$he(opt2$par)
H5 = obj5$he(opt5$par)

# sample parameter vectors from multivariate normal distribution
B = 1000
set.seed(123)
allparsHSMM = rmvnorm(B, opt5$par, solve(H5)) # matrix
allparsHSMM = lapply(1:B, function(b) relist(allparsHSMM[b,], skeletonHSMM))

allparsHMM = rmvnorm(B, opt2$par, solve(H2)) # matrix
allparsHMM = lapply(1:B, function(b) relist(allparsHMM[b,], skeletonHMM))

# function that computes the HSMM's stationary distribution at a given time point based on parameter
get_delta = function(par, t, agsizes = rep(30,3), startInds = c(0, cumsum(agsizes[1:(length(agsizes)-1)]))){
  todseq = (t + 1:24) %% 24
  Mu_dwell = exp(cbind(1, trigBasisExp(todseq, degree=1)) %*% t(par$betaMu))
  Size = 1 / exp(par$betaPhi)
  dm = lapply(1:N, function(i) sapply(1:agsizes[i]-1, dnbinom, Size[i], Size[i]/(Size[i] + Mu_dwell[,i])))
  omega = tpm_emb(par$beta_omega) # case distinction: 2 states -> no parameters to vary
  Gamma = tpm_phsmm(omega, dm)
  Gamma = simplify2array(lapply(Gamma, as.matrix))
  bigDelta = stationary_p(Gamma, t = 1)
  sapply(1:N, function(j) sum(bigDelta[startInds[j]+1:agsizes[j]]))
}

# function that computes the HMM's stationary distribution at a given time point based on parameter
get_deltaHMM = function(par, t){
  todseq = (t + 1:24) %% 24
  Gamma = tpm_p(tod=todseq, L=24, beta=par$beta, degree=1)
  stationary_p(Gamma, t=1)
}

# point estimates
npoints = 150
timeseq = seq(0, 24, length = npoints)
Delta_cont = Delta_cont_hmm = matrix(nrow = npoints, ncol = 3)
for(i in 1:npoints){
  Delta_cont[i,] = get_delta(relist(opt5$par, skeletonHSMM), t = timeseq[i])
  Delta_cont_hmm[i,] = get_deltaHMM(relist(opt2$par, skeletonHMM), t = timeseq[i])
}

# # sampling for confidence bands
# Deltas = Deltas_hmm = array(dim = c(npoints, 3, B))
# for(t in 1:npoints){
#   print(t)
#   for(b in 1:B){
#     Deltas[t,,b] = get_delta(allparsHSMM[[b]], t = timeseq[t])
#     Deltas_hmm[t,,b] = get_deltaHMM(allparsHMM[[b]], t = timeseq[t])
#   }
# }
# saveRDS(Deltas, file = "./application/objects/Deltas.rds")
# saveRDS(Deltas_hmm, file = "./application/objects/Deltas_hmm.rds")

Deltas = readRDS("./application/objects/Deltas.rds")
Deltas_hmm = readRDS("./application/objects/Deltas_hmm.rds")
DeltaCI = apply(Deltas, c(1,2), quantile, probs = c(0.025, 0.975))
DeltaCI_hmm = apply(Deltas_hmm, c(1,2), quantile, probs = c(0.025, 0.975))


# pdf("./figures/stationary_p_muskox.pdf", width = 7, height = 4)

m = matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.15, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,10), ylim = c(1,10))
legend(x = 2.3, y = 9.5, inset = c(0.3,0), lwd = 2,
       col = color,
       legend = c("resting", "foraging", "travelling"), bty = "n", 
       horiz = TRUE, text.width = c(1, 1, 1.2))
par(mar = c(5,4,2,2)+0.1)

# phsmm
plot(timeseq, Delta_cont[,1], bty = "n", type = "l", lwd = 2, col = color[1], xlim = c(0,24), ylim = c(0,0.8),
     ylab = "stationary state probabilities", xlab = "time of day", xaxt = "n")
# lines(timeseq, Delta_cont[,1], lwd = 2, col = color[1])
polygon(c(timeseq, rev(timeseq)), c(DeltaCI[1,,1], rev(DeltaCI[2,,1])), col=scales::alpha(color[1],0.2), border = F)#,border=color[1])
for(j in 2:3){
  lines(timeseq, Delta_cont[,j], lwd = 2, col = color[j])
  polygon(c(timeseq, rev(timeseq)), c(DeltaCI[1,,j], rev(DeltaCI[2,,j])), col=scales::alpha(color[j],0.2), border = F)
}
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
legend("top", bty = "n", legend = "inhomogeneous HSMM")

# phmm
plot(timeseq, Delta_cont_hmm[,1], bty = "n", type = "l", lwd = 2, col = color[1], xlim = c(0,24), ylim = c(0,0.8),
     ylab = "stationary state probabilities", xlab = "time of day", xaxt = "n")
# lines(timeseq, Delta_cont_hmm[,1], lwd = 2, col = color[1])
polygon(c(timeseq, rev(timeseq)), c(DeltaCI_hmm[1,,1], rev(DeltaCI_hmm[2,,1])), col=scales::alpha(color[1],0.2), border = F)#,border=color[1])
for(j in 2:3){
  lines(timeseq, Delta_cont_hmm[,j], lwd = 2, col = color[j])
  polygon(c(timeseq, rev(timeseq)), c(DeltaCI_hmm[1,,j], rev(DeltaCI_hmm[2,,j])), col=scales::alpha(color[j],0.2), border = F)
}
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
legend("top", bty = "n", legend = "inhomogeneous HMM")

# dev.off()



# Dwell-time distributions ------------------------------------------------

## Mean dwell times
par(mfrow = c(1,1))
plot(mod5$Mu_dwell[,1]+1, type = "b", lwd = 2, col = color[1], bty = "n", 
     ylim = c(2,8), xlim = c(0,24),
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
lines(mod5$Mu_dwell[,2]+1, type = "b", lwd = 2, col = color[2])
lines(mod5$Mu_dwell[,3]+1, type = "b", lwd = 2, col = color[3])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))


# Time-varying dwell-time distribution
todseq = seq(0,24,length=200)
Z = trigBasisExp(todseq, 24, 1)
Museq = exp(cbind(1,Z)%*%t(par5$betaMu))+1

statenames = c("resting", "foraging", "travelling")

# pdf("./figures/time_varying_distr_heat5.pdf", width = 7, height = 2.7)

par(mfrow = c(1,3), mar = c(5,4,4,1.5))
for(j in 1:3){
  plot(rep(1, 20), 1:20, ylim = c(0,16), xlim = c(0,24), pch = 16, 
       col = scales::alpha(color[j], mod5$dm[[j]][1,]/max(mod5$dm[[j]][1,])*0.7), 
       bty = "n", xaxt = "n", xlab = "time of day", ylab = "dwell time", main = statenames[j])
  for(t in 2:24){
    points(rep(t,20), 1:20, pch = 16, 
           col = scales::alpha(color[j], mod5$dm[[j]][t,]/max(mod5$dm[[j]][t,])*0.7))
  }
  lines(todseq, Museq[,j], lwd = 4, col = "#ffffff")
  lines(todseq, Museq[,j], lwd = 2, col = color[j])
  # lines(todseq, Sdseq[,j], lwd = 1, lty = 2, col = color[j])
  
  axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
  # axis(4, at = seq(0,15,by=5), labels = seq(0,15,by=5))
  # mtext("standard deviation", side=4, line = 3)
  if(j==3) legend(x=12, y=16, bty = "n", lwd = 2, lty = 1, legend = c("mean"))
}

# dev.off()


## Overall dwell-time distribution

# periodic HSMM (mod5)
pmf5 = lapply(1:3, function(j) ddwell_hsmm(1:30, j, mod5$dm, mod5$Gamma, rep(30,3)))

# periodic HMM (mod2)
pmf2 = lapply(1:3, function(j) ddwell(1:30, j, mod2$Gamma))

## decode states for empirical dwell-time distributions

source("./functions/dwell_functions.R")

# pHSMM
states5_large = viterbi_p(as.vector(mod5$delta_sparse), mod5$Gamma, mod5$allprobs[,rep(1:3,times=agsizes)], muskox$tod)
mod5$states = rep(NA, nrow(muskox))
mod5$states[which(states5_large %in% 1:30)] = 1
mod5$states[which(states5_large %in% 31:60)] = 2
mod5$states[which(states5_large %in% 61:90)] = 3

mod5$empmpf = empirical_pmf(mod5$states, rep(30,3))

# pHMM
mod2$states = viterbi_p(mod2$delta, mod2$Gamma_p, mod2$allprobs, muskox$tod)

mod2$empmpf = empirical_pmf(mod2$states, rep(30,3))

# HSMM
mod4$Gamma = as.matrix(mod4$Gamma_sparse)
states4_large = viterbi(as.vector(mod4$delta_sparse), mod4$Gamma, mod4$allprobs[,rep(1:3,times=rep(20,3))])

mod4$states = rep(NA, nrow(muskox))
mod4$states[which(states4_large %in% 1:20)] = 1
mod4$states[which(states4_large %in% 21:40)] = 2
mod4$states[which(states4_large %in% 41:60)] = 3

mod4$empmpf = empirical_pmf(mod4$states, rep(30,3))

# pdf("./figures/overall_ddwell_muskox.pdf", width = 8, height = 4.5)

par(mfrow = c(3,3), mar = c(4,4,1.5,1)+0.1)
plot(1:12, pmf5[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "", ylab = "probabilities", main = "resting", ylim = c(0,0.25))
points(1:12, mod5$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, pmf5[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "", ylab = "probabilities", main = "foraging", ylim = c(0,0.2))
points(1:14, mod5$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, pmf5[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "", ylab = "probabilities", main = "travelling", ylim = c(0,0.2))
points(1:14, mod5$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "inhomogeneous HSMM")

plot(1:12, pmf2[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.3))
points(1:12, mod2$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, pmf2[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:14, mod2$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, pmf2[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:14, mod2$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "inhomogeneous HMM")

plot(1:12, mod4$dm[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:12, mod4$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, mod4$dm[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.2))
points(1:14, mod4$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, mod4$dm[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.2))
points(1:14, mod4$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "homogeneous HSMM")

# dev.off()


