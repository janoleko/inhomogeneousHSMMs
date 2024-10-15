
# Packages ----------------------------------------------------------------

library(moveHMM)
library(PHSMM)
library(CircStats)
library(optimParallel)
library(LaMa)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(patchwork)


# Data --------------------------------------------------------------------

muskox = moveHMM::prepData(muskox, type = "UTM")
muskox$tod = muskox$tday+1
muskox = muskox[,-c(1,7,8)]

head(muskox)

nrow(muskox)
sum(is.na(muskox$step))
sum(muskox$step[!is.na(muskox$step)]==0) # no need for zero-inflated gamma distribution
sum(is.na(muskox$angle))


# EDA ---------------------------------------------------------------------

# track on map
## first convert UTM to lontitude and latitude

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

# Plot using ggplot2
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

plot_fine = ggplot() +
  geom_sf(data = world, fill = "gray95") +
  # geom_sf(data = gps_line, color = "black", size = .1) +  # Add the connecting line
  geom_sf(data = gps_sf, color = "black", size = .3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  fixed_plot_aspect(ratio = 1) +
  theme_minimal()# +

# Combine plots side by side
combined_plot <- plot_coarse + plot_fine

# Adjust the theme to reduce margins
combined_plot <- combined_plot + theme(plot.margin = margin(0, 0, 0, 0))

ggsave("./figures/maps.pdf", combined_plot, width = 10, height = 5, units = "in")





# raw time series
# plot(muskox$x, muskox$y, xlab = "UTM easting", ylab = "UTM northing", type = "l", bty = "n")

# pdf("./figures/muskox_step.pdf", width = 8, height = 4)
par(mar = c(4.5,4,2,2))
plot(muskox$step[1:1000], type = "h", bty = "n", ylab = "step length (meters)", xlab = "time")
# dev.off()

# plot(muskox$step, type = "h", bty = "n", ylab = "step length (meters)")
plot(muskox$angle, type = "h", bty = "n", ylab = "turning angle (radians)")

# pdf("./figures/muskox_hist.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
hist(muskox$step, prob = T, main = "", bor = "white", xlab = "step length (meters)", ylab = "density", 
     xlim = c(0,500), breaks = 100, ylim = c(0, 0.01))
hist(muskox$angle, prob = T, breaks = 25, main = "", bor = "white", xlab = "turning angle (radians)", ylab = "density", xaxt = "n")
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), 
                expression(0), expression(pi/2), expression(pi)))
# dev.off()

means = rep(NA,24)
for(t in 1:24){
  means[t] = mean(muskox$step[which(muskox$tod==t)], na.rm=T)
}
plot(1:24, means, type = "h", lwd = 4)

color = c("orange", "deepskyblue", "seagreen2")

# pdf("./figures/muskox_boxplot.pdf", width = 6, height = 3)
par(mfrow = c(1,1), mar = c(4.5,4,1.5,2)+0.1)
boxplot(muskox$step ~ muskox$tod, pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000010", 
        frame = F, xlab = "time of day", ylab = "step length (meters)", ylim = c(0,400))
points(1:24, means, type = "o", col = color[2], pch = 19)
# lines(1:24, means, type = "b", col = color[2])
legend(x = 18, y = 420, lwd = c(2, 1), pch = c(NA, 19), col = c("black", color[2]),
       legend = c("median", "mean"), bty = "n")
# dev.off()




# Loading functions and colors --------------------------------------------

color = c("orange", "deepskyblue", "seagreen2")
n = nrow(muskox)
source("./functions/likelihoods.R")
source("./functions/auxiliary_functions.R")

# Fitting HMMs ------------------------------------------------------------

## homogeneous
theta.star = c(rep(-2,6), log(c(5, 40, 500, 5, 40, 500, 0.1, 0.5, 2)))
mod1 = nlm(mllk, theta.star, X=muskox, N=3, ptransform=partransform, 
           iterlim = 1e3, print.level = 2)
(p1_hat = partransform(mod1$estimate))
(IC1 = AIC_BIC(mod1, n))

## inhomogeneous
### K = 1
theta.star = c(rep(-2,6), rep(0, 2*6), log(c(4, 80, 400, 4, 70, 400, 0.4, 0.2, 1.5)))
mod2 = nlm(mllk_p, theta.star, X=muskox, N=3, L=24, K=1, ptransform_p=partransform_p, 
           iterlim = 1e3, print.level = 2, stepmax = 10, hessian = TRUE)
(p2_hat = partransform_p(mod2$estimate, L=24, K=1, t1=muskox$tod[1]))  
(IC2 = AIC_BIC(mod2, n))

Delta2_hat = stationary_p(p2_hat$Gamma)

plot(1:24, Delta2_hat[,1], type = "l", lwd = 2, bty = "n", col = color[1], 
     ylim = c(0,0.8), xlim = c(0,24), xaxt = "n", xlab = "time of day", ylab = "stationary state probabilities")
points(1:24, Delta2_hat[,1], pch = 19, col = color[1])
lines(1:24, Delta2_hat[,2], type = "l", lwd = 2, col = color[2])
points(1:24, Delta2_hat[,2], pch = 19, col = color[2])
lines(1:24, Delta2_hat[,3], type = "l", lwd = 2, col = color[3])
points(1:24, Delta2_hat[,3], pch = 19, col = color[3])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))

### K = 2
theta.star = c(rep(-2,6), rep(0, 4*6), log(c(5, 80, 410, 3.5, 70, 380, 0.4, 0.2, 1.5)))
mod3 = nlm(mllk_p, theta.star, X=muskox, N=3, L=24, K=2, ptransform_p=partransform_p, 
           iterlim = 1e3, print.level = 2, stepmax = 10)
(p3_hat = partransform_p(mod3$estimate, L=24, K=2, t1=muskox$tod[1])) 
(IC3 = AIC_BIC(mod3, n)) # best AIC but worst BIC



# Fitting HSMMs -----------------------------------------------------------

## homogeneous
# getting an intuition on dwell time starting values
# mean (geometric) dwell times:
1/(1-diag(p1_hat$Gamma))

theta.star = c(log(c(rep(4,3), rep(0.1,3))), log(c(4, 80, 400, 5, 70, 400, 0.4, 0.2, 2)), rep(0,3))

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_s")
# mod4.o = optimParallel(theta.star, mllk_s, X = muskox, N=3, agsizes = c(20,20,20),
#                      parallel=list(forward=T),control = list(trace=4, maxit=1e4, ndeps=rep(1e-4,length(theta.star))))
# stopCluster(cl)
# mod4 = nlm(mllk_s, mod4.o$par, X=muskox, N=3, agsizes = c(20,20,20),
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod4, "./application/models/mod4.rds")

mod4 = readRDS("./application/models/mod4.rds")

p4_hat = partransform_s(mod4$estimate, agsizes = c(20,20,20))
(IC4 = AIC_BIC(mod4, n)) # best so far

par(mfrow = c(1,3))
plot(1:20, p4_hat$dm[[1]], type = "h", lwd = 2, bty = "n", col = color[1], main = "resting",
     xlab = "dwell time", ylab = "probabilities", xlim = c(0,20), ylim = c(0,0.25))
plot(1:20, p4_hat$dm[[2]], type = "h", lwd = 2, bty = "n", col = color[2], main = "foraging",
     xlab = "dwell time", ylab = "probabilities", xlim = c(0,20), ylim = c(0,0.25))
plot(1:20, p4_hat$dm[[3]], type = "h", lwd = 2, bty = "n", col = color[3], main = "travelling",
     xlab = "dwell time", ylab = "probabilities", xlim = c(0,20), ylim = c(0,0.25))


## inhomogeneous
# K = 1
theta.star = c(log(rep(4,3)), rep(0,3*2), log(rep(0.1,3)), log(c(4, 80, 400, 5, 70, 400, 0.4, 0.2, 2)), rep(0,3))

# mllk_sp(theta.star, X=muskox, N=3, L=24, K=1, agsizes = rep(30,3), ptransform_sp = partransform_sp)

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_sp")
# mod5_optim = optimParallel(par=theta.star, fn=mllk_sp, X=muskox, N=3, L=24, K=1, agsizes = rep(30,3),
#                        control = list(trace = 4, maxit = 1e4, ndeps=rep(1e-5, length(theta.star))))
# stopCluster(cl)

# mod5 = nlm(mllk_sp, mod5_optim$par, X=muskox, N=3, L=24, K=1, agsizes = rep(30,3),
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod5, "./application/models/mod5.rds")

mod5 = readRDS("./application/models/mod5.rds")
p5_hat = partransform_sp(mod5$estimate, K=1, agsizes = rep(30,3), t1 = muskox$tod[1])
(IC5 = AIC_BIC(mod5, n)) # a little better again



# K = 2
theta.star = c(log(rep(3,3)), rep(0, 3*2*2), log(rep(10,3)), log(c(4, 80, 400, 5, 70, 400, 0.4, 0.2, 2)), rep(0,3))
theta.star = c(mod5$estimate[1:9], rep(0, 6), mod5$estimate[-(1:9)])

# mllk_sp(theta.star, X=muskox, N=3, L=24, K=2, agsizes = rep(35,3))

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_sp")
# mod.o = optimParallel(theta.star, mllk_sp, X=muskox, N=3, L=24, K=2, agsizes = rep(35,3),
#                        parallel=list(forward=T),control=list(trace=4, maxit=1e4, ndeps=rep(1e-5,length(theta.star))))
# stopCluster(cl)
# mod6_optim = mod.o
# 
# mod6 = nlm(mllk_sp, mod6_optim$par, X=muskox, N=3, L=24, K=2, agsizes = rep(35,3),
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod6, "./application/models/mod6.rds")

mod6 = readRDS("./application/models/mod6.rds")

p6_hat = partransform_sp(mod6$estimate, K=2, agsizes = rep(35,3), t1 = muskox$tod[1])
(IC6 = AIC_BIC(mod6, n)) # worse


## inhomogeneity in dispersion parameter
# K = c(1,1)
K = c(1,1)

# initialize with very little overdispersion
theta.star = c(mod5$estimate[1:9], log(p5_hat$phi_dwell), rep(0,3*2), 
               mod5$estimate[-(1:12)])

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_sp2")
# mod7_optim = optimParallel(par=theta.star, fn=mllk_sp2, X=muskox, N=3, L=24, K=c(1,1), agsizes = rep(30,3),
#                            parallel = list(forward=T),control = list(trace=4, maxit=1e4, ndeps=rep(1e-7,length(theta.star))))
# stopCluster(cl)
# mod7 = nlm(mllk_sp2, mod7_optim$par, X=muskox, N=3, L=24, K=c(1,1), agsizes = rep(40,3), 
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod7, "./application/models/mod7.rds")

mod7 = readRDS("./application/models/mod7.rds")

p7_hat = partransform_sp2(mod7$estimate, K=c(1,1), agsizes = rep(40,3), t1 = muskox$tod[1])
(IC7 = AIC_BIC(mod7, n)) # better AIC but worse BIC than model with K = 1

par(mfrow = c(1,2))
# visualizing time-varying mean
plot(p7_hat$Mu_dwell[,1]+1, type = "l", lwd = 2, col = color[1], bty = "n", 
     ylim = c(2,10), xlim = c(0,24),
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
points(p7_hat$Mu_dwell[,1]+1, pch = 19, col = color[1])
lines(p7_hat$Mu_dwell[,2]+1, type = "l", lwd = 2, col = color[2])
points(p7_hat$Mu_dwell[,2]+1, pch = 19, col = color[2])
lines(p7_hat$Mu_dwell[,3]+1, type = "l", lwd = 2, col = color[3])
points(p7_hat$Mu_dwell[,3]+1, pch = 19, col = color[3])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))
# visualizing time-varying dispersion
plot(p7_hat$Phi_dwell[,1], type = "l", lwd = 2, col = color[1], bty = "n", 
     ylim = c(0,2), xlim = c(0,24),
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
points(p7_hat$Phi_dwell[,1], pch = 19, col = color[1])
lines(p7_hat$Phi_dwell[,2], type = "l", lwd = 2, col = color[2])
points(p7_hat$Phi_dwell[,2], pch = 19, col = color[2])
lines(p7_hat$Phi_dwell[,3], type = "l", lwd = 2, col = color[3])
points(p7_hat$Phi_dwell[,3], pch = 19, col = color[3])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))

### overall dwell-time distributino
source("./functions/dwell_functions.R")

pmf7 = list()
N=3
for(j in 1:N){
  pmf7[[j]] = ddwell_hsmm(1:40, j, p7_hat$dm, p7_hat$Gamma, rep(40,3))
}

par(mfrow = c(1,3))
plot(1:15, pmf7[[1]][1:15], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "resting")
plot(1:20, pmf7[[2]][1:20], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "foraging")
plot(1:25, pmf7[[3]][1:25], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "dwell time", ylab = "probabilities", main = "travelling")


# time-varying dwell time mean and standard deviation
par(mfrow = c(1,2))
todseq = seq(0,24,length=200)
Z = trigBasisExp(todseq, 24, 1)
Museq = exp(cbind(1,Z)%*%t(p7_hat$beta_Mu))+1
plot(todseq, Museq[,1], ylim = c(0,10), type = "l", lwd = 2, col = color[1], 
     bty = "n", xaxt = "n", xlab = "time of day", ylab = "mean dwell time")
lines(todseq, Museq[,2], lwd = 2, col = color[2])
lines(todseq, Museq[,3], lwd = 2, col = color[3])
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
Phiseq = exp(cbind(1,Z)%*%t(p7_hat$beta_Phi))
Varseq = Museq+(Museq^2)*Phiseq
Sdseq = sqrt(Varseq)
plot(todseq, Sdseq[,1], ylim = c(0,15), type = "l", lwd = 2, col = color[1], 
     bty = "n", xaxt = "n", xlab = "time of day", ylab = "dwell time standard deviation")
lines(todseq, Sdseq[,2], lwd = 2, col = color[2])
lines(todseq, Sdseq[,3], lwd = 2, col = color[3])
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))

dm = p7_hat$dm
dmq = array(dim = c(24, 3, 3))
for(j in 1:3){
  for(t in 1:24){
    if(dm[[j]][t,1]<0.25){
      qlow = which(cumsum(dm[[j]][t,])<0.25)
      qlow = qlow[length(qlow)]
    } else{ qlow = 1}
    if(dm[[j]][t,1]<0.5){
      med = which(cumsum(dm[[j]][t,])<0.5)
      med = med[length(med)]
    } else{ med = 1}
    qup = which(cumsum(dm[[j]][t,])>0.75)
    qup = qup[1]
    dmq[t,j,] = c(qlow, med, qup) 
  }
}


# todseq = seq(0,24,length=200)
# Z = trigBasisExp(todseq, 24, 1)
# Museq = exp(cbind(1,Z)%*%t(p7_hat$beta_Mu))+1

# plot the time-varying distribution implied by model 7
pdf("./figures/time_varying_distr_heat.pdf", width = 7, height = 2.7)
statenames = c("resting", "foraging", "travelling")
par(mfrow = c(1,3), mar = c(5,4,4,1.5))
for(j in 1:3){
  plot(rep(1, 20), 1:20, ylim = c(0,16), xlim = c(0,24), pch = 16, 
       col = scales::alpha(color[j],p7_hat$dm[[j]][1,]/max(p7_hat$dm[[j]][1,])*0.7), 
       bty = "n", xaxt = "n", xlab = "time of day", ylab = "dwell time", main = statenames[j])
  for(t in 2:24){
    points(rep(t,20), 1:20, pch = 16, 
           col = scales::alpha(color[j],p7_hat$dm[[j]][t,]/max(p7_hat$dm[[j]][t,])*0.7))
  }
  lines(todseq, Museq[,j], lwd = 2, col = color[j])
  lines(todseq, Sdseq[,j], lwd = 1, lty = 2, col = color[j])
  
  axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
  # axis(4, at = seq(0,15,by=5), labels = seq(0,15,by=5))
  # mtext("standard deviation", side=4, line = 3)
  if(j==3) legend(x=12, y=16, bty = "n", lwd = c(2,1), lty = c(1,2), legend = c("mean", "sd"))
}
dev.off()


## modelling the conditional tpm
# K = c(1,1)
K = 1
theta.star = c(mod5$estimate[1:21], rep(0,3*3))
theta.star = c(mod5$estimate[1:21], runif(9,-1,1))

mllk_sp3(theta.star, X=muskox, K=1, agsizes = rep(30,3))

# cl = makeCluster(8); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_sp3")
# mod8_optim = optimParallel(par=theta.star, fn=mllk_sp3, X=muskox, N=3, L=24, K=1, agsizes = rep(30,3),
#                            parallel=list(forward=T), control = list(trace=4, maxit=1e4, ndeps=rep(1e-6, length(theta.star))))
# stopCluster(cl)
# mod8 = nlm(mllk_sp3, mod8_optim$par, X=muskox, N=3, L=24, K=1, agsizes = rep(30,3),
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod8, "./application/models/mod8.rds")

mod8 = readRDS("./application/models/mod8.rds")

p8_hat = partransform_sp3(mod8$estimate, K=1, agsizes = rep(30,3), t1 = muskox$tod[1])
(IC8 = AIC_BIC(mod8, n)) # better AIC but worse BIC than model with K = 1


pmf = list()
for(j in 1:N){
  pmf[[j]] = ddwell_hsmm(1:30, j, p8_hat$dm, p8_hat$Gamma, rep(30,3))
}


# conditional tpm Omega as function of the time of day
par(mfrow = c(3,3))
for(i in 1:N){
  for(j in 1:N){
    plot(p8_hat$omega[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", 
         ylab = expression(omega), xlab = "time of day")
  }
}
# no strong effects here

# K = 1 for mean, dispersion and omega
theta.star = c(mod7$estimate[1:27], mod8$estimate[22:30])

# cl = makeCluster(8); setDefaultCluster(cl=cl)
# clusterExport(cl, "partransform_sp4")
# mod9_optim = optimParallel(par=theta.star, fn=mllk_sp4, X=muskox, N=3, L=24, K=c(1,1,1), agsizes = rep(40,3),
#                            control = list(trace=4, maxit=1e4, ndeps=rep(1e-7,length(theta.star))))
# stopCluster(cl)
# mod9 = nlm(mllk_sp4, mod9_optim$par, X=muskox, N=3, L=24, K=c(1,1,1), agsizes = rep(45,3), 
#            iterlim = 1e3, print.level = 2)
# saveRDS(mod9, "./application/models/mod9.rds")

mod9 = readRDS("./application/models/mod9.rds")

p9_hat = partransform_sp4(mod9$estimate, K=c(1,1,1), agsizes = rep(45,3), t1 = muskox$tod[1])
(IC9 = AIC_BIC(mod9, n))

# conditional tpm Omega as function of the time of day again
par(mfrow = c(3,3))
for(i in 1:N){
  for(j in 1:N){
    plot(p9_hat$omega[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", 
         ylab = expression(omega), xlab = "time of day")
  }
}
# again, no strong effects

# plot the probability that a dwell time is started

# get_weights = function(state, # which state to compute
#                        Gamma, # list of L pre-calculated periodic state-aggregate Gamma matrices
#                        aggr_sizes # vector of state-aggregate sizes
# ){
#   L = dim(Gamma)[3]
#   largeN = sum(aggr_sizes)
#   I_minus = ifelse(state>1, aggr_sizes[state-1]+1, 1) # index of state to transition to
#   # index set for sum
#   if(state > 1){aggr_ind = cumsum(aggr_sizes)[state-1]+1:aggr_sizes[state]
#   } else{aggr_ind = 1:aggr_sizes[1]}
#   I = setdiff(1:sum(aggr_sizes), aggr_ind)
#   delta = LaMa::stationary_p(Gamma)
#   weights = numeric(L)
#   weights[1] = sum(delta[L, I]*Gamma[I,I_minus,L])
#   for (t in 2:L){ 
#     weights[t] = sum(delta[t-1, I]*Gamma[I,I_minus,t-1]) 
#   }
#   weights = weights/sum(weights)
#   return(weights)
# }



# Information criteria ----------------------------------------------------

ICs = cbind(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
colnames(ICs) = paste0("mod", 1:9)
ICs
tab = rbind(ICs,
            c(-mod1$minimum, -mod2$minimum, -mod3$minimum, -mod4$minimum,
              -mod5$minimum, -mod6$minimum, -mod7$minimum, -mod8$minimum, -mod9$minimum))
rownames(tab) = c("AIC", "BIC", "llk")
round(t(tab),1)


# # Visualizing mod7 --------------------------------------------------------
# 
# color = c("orange", "deepskyblue", "seagreen2")
# 
# ## State-dependent distribution
# 
# Delta_hat7_tilde = stationary_p(p7_hat$Gamma)
# Delta_hat7 = matrix(nrow = 24, ncol = 3)
# agsizes = rep(40,3);N=3
# for(j in 1:N){
#   Delta_hat7[,j] = apply(Delta_hat7_tilde[,c(0,cumsum(agsizes))[j]+1:agsizes[j]], 1, sum)
# }
# delta_hat = colMeans(Delta_hat7)
# 
# shape = p7_hat$mu^2/p7_hat$sigma^2; scale = p7_hat$sigma^2/p7_hat$mu
# 
# hist(muskox$step, prob = T, breaks = 100, bor = "white", xlim = c(0,800), main = "")
# curve(delta_hat[1]*dgamma(x, shape = shape[1], scale = scale[1]),
#       add = T, lwd = 2, col = color[1], n = 500)
# curve(delta_hat[2]*dgamma(x, shape = shape[2], scale = scale[2]),
#       add = T, lwd = 2, col = color[2], n = 500)
# curve(delta_hat[3]*dgamma(x, shape = shape[3], scale = scale[3]),
#       add = T, lwd = 2, col = color[3], n = 1000)
# 
# library(CircStats)
# hist(muskox$angle, prob = T, breaks = 30, bor = "white", xlim = c(-pi,pi), main = "")
# curve(delta_hat[1]*dvm(x, p7_hat$mu.turn[1], p7_hat$kappa[1]),
#       add = T, lwd = 2, col = color[1], n = 500)
# curve(delta_hat[2]*dvm(x, p7_hat$mu.turn[2], p7_hat$kappa[2]),
#       add = T, lwd = 2, col = color[2], n = 500)
# curve(delta_hat[3]*dvm(x, p7_hat$mu.turn[3], p7_hat$kappa[3]),
#       add = T, lwd = 2, col = color[3], n = 500)
# curve(delta_hat[1]*dvm(x, p7_hat$mu.turn[1], p7_hat$kappa[1])+
#         delta_hat[2]*dvm(x, p7_hat$mu.turn[2], p7_hat$kappa[2])+
#         delta_hat[3]*dvm(x, p7_hat$mu.turn[3], p7_hat$kappa[3]), add = T, lwd = 2, lty = 2, n = 500)
# 
# 
# 
# # Periodically stationary distribution ------------------------------------
# 
# theta.star = mod7$estimate
# N = 3; K = c(1,1)
# 
# get_delta = function(theta.star, N=3, L=24, agsizes = rep(40,N)){
#   beta_Mu = matrix(theta.star[1:(N*(1+2*K[1]))], nrow = N)
#   beta_Phi = matrix(theta.star[N*(1+2*K[1])+1:(N*(1+2*K[2]))], nrow = N)
#   if(N>2){ # only needed if N>2
#     omega = matrix(0,N,N)
#     omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
#                             exp(theta.star[N*(2+2*(K[1]+K[2]))+3*N+1:(N*(N-2))])),N,N-1)))
#     omega = t(omega)/apply(omega,2,sum)
#   }else{ omega = matrix(c(0,1,1,0),2,2) }
#   Z1 = LaMa::trigBasisExp(tod=1:L, L=L, degree=K[1])
#   Z2 = LaMa::trigBasisExp(tod=1:L, L=L, degree=K[2])
#   Mu_dwell = exp(cbind(1,Z1)%*%t(beta_Mu))
#   Phi_dwell = exp(cbind(1,Z2)%*%t(beta_Phi))
#   dm = list()
#   for(j in 1:N){
#     dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu=Mu_dwell[,j], size=1/Phi_dwell[,j])
#   }
#   Gamma = LaMa::tpm_phsmm(omega, dm)
#   Delta_large = LaMa::stationary_p(Gamma)
#   startInds = c(0, cumsum(agsizes[1:(N-1)]))
#   Delta = matrix(nrow=L, ncol = N)
#   for(j in 1:N){
#     Delta[,j] = apply(Delta_large[,startInds[j]+1:agsizes[j]],1,sum)
#   }
#   return(Delta)
# }
# 
# Delta = get_delta(mod7$estimate)
# 
# get_delta_cont = function(theta.star, t, N=3, L=24, agsizes = rep(40,N)){
#   beta_Mu = matrix(theta.star[1:(N*(1+2*K[1]))], nrow = N)
#   beta_Phi = matrix(theta.star[N*(1+2*K[1])+1:(N*(1+2*K[2]))], nrow = N)
#   if(N>2){ # only needed if N>2
#     omega = matrix(0,N,N)
#     omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
#                                            exp(theta.star[N*(2+2*(K[1]+K[2]))+3*N+1:(N*(N-2))])),N,N-1)))
#     omega = t(omega)/apply(omega,2,sum)
#   }else{ omega = matrix(c(0,1,1,0),2,2) }
#   todseq = (t + 1:L) %% L
#   Z1 = LaMa::trigBasisExp(tod=todseq, L=L, degree=K[1])
#   Z2 = LaMa::trigBasisExp(tod=todseq, L=L, degree=K[2])
#   Mu_dwell = exp(cbind(1,Z1)%*%t(beta_Mu))
#   Phi_dwell = exp(cbind(1,Z2)%*%t(beta_Phi))
#   dm = list()
#   for(j in 1:N){
#     dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu=Mu_dwell[,j], size=1/Phi_dwell[,j])
#   }
#   Gamma = LaMa::tpm_phsmm(omega, dm)
#   delta_large = LaMa::stationary_p(Gamma, t=1)
#   startInds = c(0, cumsum(agsizes[1:(N-1)]))
#   delta = rep(NA,N)
#   for(j in 1:N){
#     delta[j] = sum(delta_large[startInds[j]+1:agsizes[j]])
#   }
#   return(delta)
# }
# 
# 
# timeseq = seq(0,24,length=200)
# Delta_cont = matrix(nrow = 200, ncol = 3)
# for(i in 1:200){
#   Delta_cont[i,] = get_delta_cont(mod7$estimate, t = timeseq[i])
# }
# 
# ## confidence bands via Monte Carlo:
# H7 = numDeriv::hessian(mllk_sp2, mod7$estimate, X=muskox, agsizes = rep(40,3))
# # approximate normal distribution of MLE inverse Hessian as approximate Fisher Info
# saveRDS(H7, "./models/mod7_hessian.rds")
# I = solve(H7)
# thetas = mvtnorm::rmvnorm(1000, mean = mod7$estimate, sigma = I)
# 
# Deltas = array(dim = c(200, 3, 1000))
# for(t in 1:200){
#   print(t)
#   for(i in 1:1000){
#     Deltas[t,,i] = get_delta_cont(thetas[i,], t = timeseq[t])
#   }
# }
# saveRDS(Deltas, "./models/stationary_MCsim.rds")
# 
# DeltaCI = array(dim = c(200, 3, 2))
# for(t in 1:200){
#   for(j in 1:3){
#     DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
#   }
# }
# 
# pdf("./figures/stationary_p_muskox.pdf", width = 7, height = 4)
# 
# plot(0:23, Delta[,1], bty = "n", pch = 19, col = "white", xlim = c(0,24), ylim = c(0,0.8),
#      ylab = "stationary probabilities", xlab = "time of day", xaxt = "n")
# lines(timeseq, Delta_cont[,1], lwd = 2, col = color[1])
# polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,1,1], rev(DeltaCI[,1,2])), col=alpha(color[1],0.2), border = F)#,border=color[1])
# 
# lines(timeseq, Delta_cont[,2], lwd = 2, col = color[2])
# polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,2,1], rev(DeltaCI[,2,2])), col=alpha(color[2],0.2), border = F)#,border=color[2])
# # points(0:23, Delta[,2], pch = 19, col = color[2])
# 
# lines(timeseq, Delta_cont[,3], lwd = 2, col = color[3])
# polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,3,1], rev(DeltaCI[,3,2])), col=alpha(color[3],0.2), border = F)#,border=color[3])
# # points(0:23, Delta[,3], pch = 19, col = color[3])
# 
# axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
# 
# legend("top", bty="n", lwd=2, col=color, horiz=T,
#        legend = c("resting", "foraging", "travelling"))
# 
# dev.off()



# Visualizing model 5 -----------------------------------------------------

color = c("orange", "deepskyblue", "seagreen2")

## State-dependent distribution

# aggregating stationary distribution for weighting of component distributions in plot
Delta_hat5_tilde = stationary_p(p5_hat$Gamma)
Delta_hat5 = matrix(nrow = 24, ncol = 3)
agsizes = rep(30,3);N=3
for(j in 1:N){
  Delta_hat5[,j] = apply(Delta_hat5_tilde[,c(0,cumsum(agsizes))[j]+1:agsizes[j]], 1, sum)
}
delta_hat = colMeans(Delta_hat5)
shape = p5_hat$mu^2/p5_hat$sigma^2; scale = p5_hat$sigma^2/p5_hat$mu

# plot state-dependent distributions
library(CircStats)

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
curve(delta_hat[1]*dgamma(x, shape = shape[1], scale = scale[1]),
      add = T, lwd = 2, col = color[1], n = 500)
curve(delta_hat[2]*dgamma(x, shape = shape[2], scale = scale[2]),
      add = T, lwd = 2, col = color[2], n = 500)
curve(delta_hat[3]*dgamma(x, shape = shape[3], scale = scale[3]),
      add = T, lwd = 2, col = color[3], n = 1000)
curve(delta_hat[1]*dgamma(x, shape = shape[1], scale = scale[1])+
        delta_hat[2]*dgamma(x, shape = shape[2], scale = scale[2])+
        delta_hat[3]*dgamma(x, shape = shape[3], scale = scale[3]), add = T, lwd = 2, lty = 2)

hist(muskox$angle, prob = T, breaks = 20, bor = "white", xlim = c(-pi,pi), ylim = c(0,0.3),
     main = "", xlab = "turning angle", ylab = "density", xaxt = "n")
curve(delta_hat[1]*dvm(x, p5_hat$mu.turn[1], p5_hat$kappa[1]),
      add = T, lwd = 2, col = color[1], n = 500)
curve(delta_hat[2]*dvm(x, p5_hat$mu.turn[2], p5_hat$kappa[2]),
      add = T, lwd = 2, col = color[2], n = 500)
curve(delta_hat[3]*dvm(x, p5_hat$mu.turn[3], p5_hat$kappa[3]),
      add = T, lwd = 2, col = color[3], n = 500)
curve(delta_hat[1]*dvm(x, p5_hat$mu.turn[1], p5_hat$kappa[1])+
        delta_hat[2]*dvm(x, p5_hat$mu.turn[2], p5_hat$kappa[2])+
        delta_hat[3]*dvm(x, p5_hat$mu.turn[3], p5_hat$kappa[3]), add = T, lwd = 2, lty = 2, n = 500)
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), 
                expression(0), expression(pi/2), expression(pi)))
# dev.off()


## periodically stationary distribution

theta.star = mod5$estimate
N = 3;

get_delta = function(theta.star, N=3, L=24, K=1, agsizes = rep(30,N)){
  beta_mu = matrix(theta.star[1:(N*(1+2*K))], nrow = N)
  phi_dwell = exp(theta.star[N*(1+2*K)+1:N])
  # parametr transform: omega
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
                                           exp(theta.star[N*(1+2*K)+4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  # dwell-time distributions and approximating tpm
  Z = LaMa::trigBasisExp(tod=1:L, L=L, degree = K)
  Mu_dwell = exp(cbind(1,Z)%*%t(beta_mu))
  dm = list()
  for(j in 1:N){
    dm[[j]] = matrix(nrow=L, ncol = agsizes[j])
    for(t in 1:L){
      dm[[j]][t,] = dnbinom(1:agsizes[j]-1, mu = Mu_dwell[t,j], size = 1/phi_dwell[j])
    }
  }
  Gamma = LaMa::tpm_phsmm(omega, dm)
  Delta_large = LaMa::stationary_p(Gamma)
  startInds = c(0, cumsum(agsizes[1:(N-1)]))
  Delta = matrix(nrow=L, ncol = N)
  for(j in 1:N){
    Delta[,j] = apply(Delta_large[,startInds[j]+1:agsizes[j]],1,sum)
  }
  return(Delta)
}

Delta = get_delta(mod5$estimate)

get_delta_cont = function(theta.star, t, N=3, L=24, K=1, agsizes = rep(30,N)){
  beta_mu = matrix(theta.star[1:(N*(1+2*K))], nrow = N)
  phi_dwell = exp(theta.star[N*(1+2*K)+1:N])
  # parametr transform: omega
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),
                          exp(theta.star[N*(1+2*K)+4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  }else{ omega = matrix(c(0,1,1,0),2,2) }
  todseq = (t + 1:L) %% L
  Z = LaMa::trigBasisExp(tod=todseq, L=L, degree=K)
  Mu_dwell = exp(cbind(1,Z)%*%t(beta_mu))
  dm = list()
  for(j in 1:N){
    dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = Mu_dwell[,j], size=1/phi_dwell[j])
  }
  Gamma = LaMa::tpm_phsmm(omega, dm)
  delta_large = LaMa::stationary_p(Gamma, t=1)
  startInds = c(0, cumsum(agsizes[1:(N-1)]))
  delta = rep(NA,N)
  for(j in 1:N){
    delta[j] = sum(delta_large[startInds[j]+1:agsizes[j]])
  }
  return(delta)
}

get_delta_cont_hmm = function(theta.star, t, N=3, L=24, K=1){
  beta = matrix(theta.star[1:(N*(N-1)*(1+2*K))], nrow=N*(N-1), ncol=1+2*K)
  todseq = (t + 1:L) %% L
  Gamma = LaMa::tpm_p(tod=todseq, L=L, beta=beta, degree=K)
  delta = LaMa::stationary_p(Gamma, t=1)
  delta
}

npoints = 150
timeseq = seq(0,24,length = npoints)
Delta_cont = Delta_cont_hmm = matrix(nrow = npoints, ncol = 3)
for(i in 1:npoints){
  Delta_cont[i,] = get_delta_cont(mod5$estimate, t = timeseq[i])
  Delta_cont_hmm[i,] = get_delta_cont_hmm(mod2$estimate, t = timeseq[i])
}

## confidence bands via Monte Carlo (takes very long)

# calculate approximate Hessian 
# H5 = numDeriv::hessian(mllk_sp, mod5$estimate, X=muskox, K = 1, agsizes = rep(30,3))

# approximate normal distribution of MLE: inverse Hessian is obseverd Fisher Info
# saveRDS(H5, "./models/mod5_hessian.rds")

H5 = readRDS("./application/models/mod5_hessian.rds")
I = solve(H5)
I_hmm = solve(mod2$hessian)
thetas = mvtnorm::rmvnorm(1000, mean = mod5$estimate, sigma = I)
thetas_hmm = mvtnorm::rmvnorm(1000, mean = mod2$estimate, sigma = I_hmm)

# Deltas = array(dim = c(npoints, 3, 1000))
# for(t in 1:npoints){
#   print(t)
#   for(i in 1:1000){
#     Deltas[t,,i] = get_delta_cont(thetas[i,], t = timeseq[t])
#   }
# }
# saveRDS(Deltas, "./models/stationary_MCsim5.rds")

Deltas = readRDS("./application/models/stationary_MCsim5.rds")
Deltas_hmm = array(dim = c(npoints, 3, 1000))
for(t in 1:npoints){
  print(t)
  for(i in 1:1000){
    Deltas_hmm[t,,i] = get_delta_cont_hmm(thetas_hmm[i,], t = timeseq[t])
  }
}

DeltaCI = DeltaCI_hmm = array(dim = c(npoints, 3, 2))
for(t in 1:npoints){
  for(j in 1:3){
    DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
    DeltaCI_hmm[t,j,] = quantile(Deltas_hmm[t,j,], probs = c(0.025, 0.975))
  }
}

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
plot(0:23, Delta[,1], bty = "n", pch = 19, col = "white", xlim = c(0,24), ylim = c(0,0.8),
     ylab = "stationary state probabilities", xlab = "time of day", xaxt = "n")
lines(timeseq, Delta_cont[,1], lwd = 2, col = color[1])
polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,1,1], rev(DeltaCI[,1,2])), col=scales::alpha(color[1],0.2), border = F)#,border=color[1])
for(j in 2:3){
  lines(timeseq, Delta_cont[,j], lwd = 2, col = color[j])
  polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,j,1], rev(DeltaCI[,j,2])), col=scales::alpha(color[j],0.2), border = F)
}
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
legend("top", bty = "n", legend = "inhomogeneous HSMM")

# phmm
plot(0:23, Delta[,1], bty = "n", pch = 19, col = "white", xlim = c(0,24), ylim = c(0,0.8),
     ylab = "stationary state probabilities", xlab = "time of day", xaxt = "n")
lines(timeseq, Delta_cont_hmm[,1], lwd = 2, col = color[1])
polygon(c(timeseq, rev(timeseq)), c(DeltaCI_hmm[,1,1], rev(DeltaCI_hmm[,1,2])), col=scales::alpha(color[1],0.2), border = F)#,border=color[1])
for(j in 2:3){
  lines(timeseq, Delta_cont_hmm[,j], lwd = 2, col = color[j])
  polygon(c(timeseq, rev(timeseq)), c(DeltaCI_hmm[,j,1], rev(DeltaCI_hmm[,j,2])), col=scales::alpha(color[j],0.2), border = F)
}
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
legend("top", bty = "n", legend = "inhomogeneous HMM")

# dev.off()



## dwell-time distributions

# mean dwell times
par(mfrow = c(1,1))
plot(p5_hat$Mu_dwell[,1]+1, type = "l", lwd = 2, col = color[1], bty = "n", 
     ylim = c(2,8), xlim = c(0,24),
     ylab = "mean dwell time", xlab = "time of day", xaxt = "n")
points(p5_hat$Mu_dwell[,1]+1, pch = 19, col = color[1])
lines(p5_hat$Mu_dwell[,2]+1, type = "l", lwd = 2, col = color[2])
points(p5_hat$Mu_dwell[,2]+1, pch = 19, col = color[2])
lines(p5_hat$Mu_dwell[,3]+1, type = "l", lwd = 2, col = color[3])
points(p5_hat$Mu_dwell[,3]+1, pch = 19, col = color[3])
axis(1, at = seq(0,24,by = 4), labels = seq(0,24,by = 4))

# time varying dwell-time distribution

todseq = seq(0,24,length=200)
Z = trigBasisExp(todseq, 24, 1)
Museq = exp(cbind(1,Z)%*%t(p5_hat$beta_mu))+1

color = c("orange", "deepskyblue", "seagreen2")
statenames = c("resting", "foraging", "travelling")

# pdf("./figures/time_varying_distr_heat5.pdf", width = 7, height = 2.7)

par(mfrow = c(1,3), mar = c(5,4,4,1.5))
for(j in 1:3){
  plot(rep(1, 20), 1:20, ylim = c(0,16), xlim = c(0,24), pch = 16, 
       col = scales::alpha(color[j],p5_hat$dm[[j]][1,]/max(p5_hat$dm[[j]][1,])*0.7), 
       bty = "n", xaxt = "n", xlab = "time of day", ylab = "dwell time", main = statenames[j])
  for(t in 2:24){
    points(rep(t,20), 1:20, pch = 16, 
           col = scales::alpha(color[j],p5_hat$dm[[j]][t,]/max(p5_hat$dm[[j]][t,])*0.7))
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


# overall dwell-time distribution

# periodic HSMM (mod5)
pmf5 = list()
for(j in 1:N){
  pmf5[[j]] = ddwell_hsmm(1:30, j, p5_hat$dm, p5_hat$Gamma, rep(30,3))
}

# periodic HMM (mod2)
pmf2 = list()
for(j in 1:N){
  pmf2[[j]] = ddwell(1:30, j, p2_hat$Gamma)
}


## decode states for empirical dwell-time distributions

source("./functions/dwell_functions.R")

# pHSMM
p5_hat$allprobs = matrix(1, nrow = nrow(muskox), ncol = N)
ind = which(!is.na(muskox$step) & !is.na(muskox$angle))
for(j in 1:N){
  p5_hat$allprobs[ind,j] = dgamma(muskox$step[ind],shape=p5_hat$mu[j]^2/p5_hat$sigma[j]^2,
                           scale=p5_hat$sigma[j]^2/p5_hat$mu[j])*
    CircStats::dvm(muskox$angle[ind], p5_hat$mu.turn[j], p5_hat$kappa[j])
}
states5_large = LaMa::viterbi_p(p5_hat$delta, p5_hat$Gamma, 
                          t(apply(p5_hat$allprobs, 1, rep, times = rep(30,3))), muskox$tod)
p5_hat$states = rep(NA, nrow(muskox))
p5_hat$states[which(states5_large %in% 1:30)] = 1
p5_hat$states[which(states5_large %in% 31:60)] = 2
p5_hat$states[which(states5_large %in% 61:90)] = 3

p5_hat$empmpf = empirical_pmf(p5_hat$states, rep(30,3))

# pHMM
p2_hat$allprobs = matrix(1, nrow = nrow(muskox), ncol = N)
ind = which(!is.na(muskox$step) & !is.na(muskox$angle))
for(j in 1:N){
  p2_hat$allprobs[ind,j] = dgamma(muskox$step[ind],shape=p2_hat$mu[j]^2/p2_hat$sigma[j]^2,
                                  scale=p2_hat$sigma[j]^2/p2_hat$mu[j])*
    CircStats::dvm(muskox$angle[ind], p2_hat$mu.turn[j], p2_hat$kappa[j])
}
p2_hat$states = LaMa::viterbi_p(p2_hat$delta, p2_hat$Gamma, p2_hat$allprobs, muskox$tod)
p2_hat$empmpf = empirical_pmf(p2_hat$states, rep(30,3))

# HSMM
p4_hat$allprobs = matrix(1, nrow = nrow(muskox), ncol = N)
ind = which(!is.na(muskox$step) & !is.na(muskox$angle))
for(j in 1:N){
  p4_hat$allprobs[ind,j] = dgamma(muskox$step[ind],shape=p4_hat$mu[j]^2/p4_hat$sigma[j]^2,
                                  scale=p4_hat$sigma[j]^2/p4_hat$mu[j])*
    CircStats::dvm(muskox$angle[ind], p4_hat$mu.turn[j], p4_hat$kappa[j])
}
states4_large = LaMa::viterbi(p4_hat$delta, p4_hat$Gamma, 
                                t(apply(p4_hat$allprobs, 1, rep, times = rep(20,3))))
p4_hat$states = rep(NA, nrow(muskox))
p4_hat$states[which(states4_large %in% 1:20)] = 1
p4_hat$states[which(states4_large %in% 21:40)] = 2
p4_hat$states[which(states4_large %in% 41:60)] = 3
p4_hat$empmpf = empirical_pmf(p4_hat$states, rep(30,3))

# pdf("./figures/overall_ddwell_muskox.pdf", width = 8, height = 4.5)

par(mfrow = c(3,3), mar = c(4,4,1.5,1)+0.1)
plot(1:12, pmf5[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "", ylab = "probabilities", main = "resting", ylim = c(0,0.25))
points(1:12, p5_hat$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, pmf5[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "", ylab = "probabilities", main = "foraging", ylim = c(0,0.2))
points(1:14, p5_hat$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, pmf5[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "", ylab = "probabilities", main = "travelling", ylim = c(0,0.2))
points(1:14, p5_hat$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "inhomogeneous HSMM")

plot(1:12, pmf2[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.3))
points(1:12, p2_hat$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, pmf2[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:14, p2_hat$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, pmf2[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:14, p2_hat$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "inhomogeneous HMM")

plot(1:12, p4_hat$dm[[1]][1:12], type = "h", lwd = 2, col = color[1], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.25))
points(1:12, p4_hat$empmpf[[1]][1:12], pch = 19, col = "#00000030")
plot(1:14, p4_hat$dm[[2]][1:14], type = "h", lwd = 2, col = color[2], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.2))
points(1:14, p4_hat$empmpf[[2]][1:14], pch = 19, col = "#00000030")
plot(1:14, p4_hat$dm[[3]][1:14], type = "h", lwd = 2, col = color[3], bty = "n", 
     xlab = "dwell time (hours)", ylab = "probabilities", main = "", ylim = c(0,0.2))
points(1:14, p4_hat$empmpf[[3]][1:14], pch = 19, col = "#00000030")
legend("topright", bty = "n", legend = "homogeneous HSMM")

# dev.off()


