#Written report SBL.06003 Classical models in Biology (exercises) 
#Spring semester 2025

#Student 1
#First name: Mathilde
#Family name: JACQUEY 
#Student number: 20-223-061

#Student 2
#First name: Florence
#Family name: KURZ
#Student number: 

#Student 3
#First name:
#Family name:
#Student number:

library(deSolve)
#!only de library "deSolve" is allowed!

########################
### Topic 4: Lotka-Volterra model for 2 competing species
rm(list=ls())

#R-code

#pré-requis: all parameters must be > 0
#r1 = intrinsic growth rate of species 1
#r2 = intrinsic growth rate of species 2
#alpha11 = intraspecific competition
#alpha12 = interspecific competition
#alpha21 = interspecific competition
#alpha22 = intraspecific competition

##1. Stable Coexistence - convergence to a stable equilibrium
#define the function with several parameters: time, population size and parameters (alpha's and r's)
LV_model <-function(t,N,p){
  N1 <- N[1]
  N2 <- N[2]
  dN1 <- N1 * (p$r1 - p$alpha11 * N1 - p$alpha12 * N2) 
  #p$alpha12 * N2 = competition effect from the second species
  dN2 <- N2 * (p$r2 - p$alpha21 * N1 - p$alpha22 * N2) 
  #N2*(p$r2 - p$alpha21 * N1) is the so-called logistic growth
  return(list(c(dN1,dN2)))
}

#pré-requis: (alpha12 / alpha22) < (r1 / r2) < (alpha11 / alpha21) and 
#(alpha21 / alpha11) < (r2 / r1) < (alpha22 / alpha12)
p <- list(r1 = 1.2, 
          r2 = 1, 
          alpha11 = 1, 
          alpha12 = 0.3, 
          alpha21 = 0.2, 
          alpha22 = 1)

#initial value (N when t=0)
N0 <- c(2.5,0.5)
#time steps at which we aim an numerical evaluation of the solution
time_steps <- seq(0,50,by=0.01)

#run the numerical integrator
out <- ode(y = N0, times = time_steps, func = LV_model, parms = p)
par(mfrow=c(1,1))

#to plot in the phase space
plot(out[,2],out[,3],xlab ='N1',ylab ='N2', type='l',xlim=c(0,6),ylim=c(0,6), lwd=2, 
     main="Phase space - Stable coexistence")
points(out[length(time_steps),2],out[length(time_steps),3],col='black',pch=16,cex=1.5)

#add the non-trivial zero growth isocline for the species 1
#intercept is r1/alpha12 and slope is -alpha11/alpha12
abline(a = p$r1 / p$alpha12, b =  -p$alpha11 / p$alpha12, lty=2, col='blue', lwd=2)
#add the non-trivial zero growth isocline for the species 2
#intercept is r2/alpha22 and slope is -alpha21/alpha22
abline(a = p$r2 / p$alpha22, b =  -p$alpha21 / p$alpha22, lty=2, col='green', lwd=2)

#to add the vector field
N1.g <- seq(0.2,6,0.4)
N2.g <- seq(0.2,6,0.4)
for(i in 1:length(N1.g)){
  for(j in 1:length(N2.g)){
    N1 <- N1.g[i]
    N2 <- N2.g[j]
    dN1 <- N1 * (p$r1 - p$alpha11 * N1 - p$alpha12 * N2)
    dN2 <- N2 * (p$r2 - p$alpha21 * N1 - p$alpha22 * N2)
    arrows(N1, N2, N1+0.1*dN1, N2+0.1*dN2, length = 0.04)
  }
}

#2. Competitive Exclusion - exclusion of one of the species
#same function as above (because this is the same Lotka-Voltera equations)
LV_model <-function(t,N,p){
  N1 <- N[1]
  N2 <- N[2]
  dN1 <- N1 * (p$r1 - p$alpha11 * N1 - p$alpha12 * N2)
  dN2 <- N2 * (p$r2 - p$alpha21 * N1 - p$alpha22 * N2)
  return(list(c(dN1,dN2)))
}

#time steps at which we aim an numerical evaluation of the solution (same as above)
time_steps <- seq(0,50,by=0.01)

#pré-requis: alpha 11 > alpha12
p_1_excludes2 <- list(r1 = 2, r2 = 1, alpha11 = 0.5, alpha12 = 0.3, alpha21 = 1, alpha22 = 1)
#initial time (N when t=0)
N0_1 <- c(0.2,3)
#run the numerical integrator for species 1
out_1 <- ode(y = N0_1, times = time_steps, func = LV_model, parms = p_1_excludes2)

#pré-requis: alpha22 > alpha21
p_2_excludes1 <- list(r1 = 1, r2 = 2, alpha11 = 1, alpha12 = 1, alpha21 = 0.3, alpha22 = 0.5)
#initial time (N when t=0)
N0_2 <- c(1.5,0.1)
#run the numerical integrator for species 2
out_2 <- ode(y = N0_2, times = time_steps, func = LV_model, parms = p_2_excludes1)

#to plot both phase space
par(mfrow=c(1,1))
#phase space for species 1
plot(out_1[,2],out_1[,3],xlab ='N1',ylab ='N2', type='l',xlim=c(0,8),ylim=c(0,8), lwd=2, 
     main="Phase space - exclusion of species 2 & species 1 wins")
points(out_1[length(time_steps),2],out_1[length(time_steps),3],col='black',pch=16,cex=1.5)

#to add the non-trivial zero growth isoclines (for species 1 and 2)
abline(a = p_1_excludes2$r1 / p_1_excludes2$alpha12, 
       b =  -p_1_excludes2$alpha11 / p_1_excludes2$alpha12, lty=2, col='blue', lwd=2)
abline(a = p_1_excludes2$r2 / p_1_excludes2$alpha22, 
       b =  -p_1_excludes2$alpha21 / p_1_excludes2$alpha22, lty=2, col='green', lwd=2)

#to add the vector field
N1.g <- seq(0.2,8,0.4)
N2.g <- seq(0.2,8,0.4)
for(i in 1:length(N1.g)){
  for(j in 1:length(N2.g)){
    N1 <- N1.g[i]
    N2 <- N2.g[j]
    dN1 <- N1 * (p_1_excludes2$r1 - p_1_excludes2$alpha11 * N1 - p_1_excludes2$alpha12 * N2)
    dN2 <- N2 * (p_1_excludes2$r2 - p_1_excludes2$alpha21 * N1 - p_1_excludes2$alpha22 * N2)
    arrows(N1, N2, N1+0.01*dN1, N2+0.01*dN2, length = 0.04)
  }
}

#phase space for species 2
plot(out_2[,2],out_2[,3],xlab ='N1',ylab ='N2', type='l',xlim=c(0,8),ylim=c(0,8), lwd=2, 
     main="Phase space - exclusion of species 1 & species 2 wins")
points(out_2[length(time_steps),2],out_2[length(time_steps),3],col='black',pch=16,cex=1.5)

#to add the non-trivial zero growth isoclines (for species 1 and 2)
abline(a = p_2_excludes1$r1 / p_2_excludes1$alpha12, b =  -p_2_excludes1$alpha11 / p_2_excludes1$alpha12, 
       lty=2, col='blue', lwd=2)
abline(a = p_2_excludes1$r2 / p_2_excludes1$alpha22, b =  -p_2_excludes1$alpha21 / p_2_excludes1$alpha22, 
       lty=2, col='green', lwd=2)

#to add the vector field
N1.g <- seq(0.2,8,0.5)
N2.g <- seq(0.2,8,0.5)
for(i in 1:length(N1.g)){
  for(j in 1:length(N2.g)){
    N1 <- N1.g[i]
    N2 <- N2.g[j]
    dN1 <- N1 * (p_2_excludes1$r1 - p_2_excludes1$alpha11 * N1 - p_2_excludes1$alpha12 * N2)
    dN2 <- N2 * (p_2_excludes1$r2 - p_2_excludes1$alpha21 * N1 - p_2_excludes1$alpha22 * N2)
    arrows(N1, N2, N1+0.01*dN1, N2+0.01*dN2, length = 0.04)}}

#3. Priority Effects and Alternative Stable States - one species will win and the other will lose
#pré-requis: alpha12 > alpha11 and alpha21 > alpha22
LV_model <-function(t,N,p){
  N1 <- N[1]
  N2 <- N[2]
  dN1 <- N1 * (p$r1 - p$alpha11 * N1 - p$alpha12 * N2)
  dN2 <- N2 * (p$r2 - p$alpha21 * N1 - p$alpha22 * N2)
  return(list(c(dN1,dN2)))
}
p <- list(r1 = 2, r2 = 2, alpha11 = 0.5, alpha12 = 1.5, alpha21 = 1.5, alpha22 = 0.5)

#initial values for both species (when t=0)
N0_1 <- c(3.5,3)
N0_2 <- c(3,3.5)
#time steps at which we aim an numerical evaluation of the solution (same as above)
time_steps <- seq(0,50,by=0.01)
#run the numerical integrator for species 1 and 2
out_1 <- ode(y = N0_1, times = time_steps, func = LV_model, parms = p)
out_2 <- ode(y = N0_2, times = time_steps, func = LV_model, parms = p)

par(mfrow=c(1,1))
#to plot the phase space
plot(out_1[,2],out_1[,3],xlab ='N1',ylab ='N2', type='l',xlim=c(0,5),ylim=c(0,5), lwd=2, 
     main="Phase space - Priority effect")
points(out_1[length(time_steps),2],out_1[length(time_steps),3],col='black',pch=16,cex=1.5)
lines(out_2[,2],out_2[,3], type='l', lwd=2)
points(out_2[length(time_steps),2],out_2[length(time_steps),3],col='black',pch=16,cex=1.5)

#to add the non-trivial zero growth isocline for both species (1 and 2)
abline(a = p$r1 / p$alpha12, b =  -p$alpha11 / p$alpha12, lty=2, col='blue', lwd=2)
abline(a = p$r2 / p$alpha22, b =  -p$alpha21 / p$alpha22, lty=2, col='green', lwd=2)

#to add the vector field
N1.g <- seq(0.2,5,0.3)
N2.g <- seq(0.2,5,0.3)
for(i in 1:length(N1.g)){
  for(j in 1:length(N2.g)){
    N1 <- N1.g[i]
    N2 <- N2.g[j]
    dN1 <- N1 * (p$r1 - p$alpha11 * N1 - p$alpha12 * N2)
    dN2 <- N2 * (p$r2 - p$alpha21 * N1 - p$alpha22 * N2)
    arrows(N1, N2, N1+0.03*dN1, N2+0.03*dN2, length = 0.04)
  }
}

#Discussion of the results 

#1. Stable Coexistence
#We demonstrated that the system converges toward a stable equilibrium point (black point), indicating 
#coexistence of both species. 
#This occurs under the classical conditions where intra-specific (alpha11, alpha22) competition is stronger 
#than inter-specific (alpha12, alpha21) competition for each species, and the zero-growth iso-clines 
#intersect in such a way 
#that each species limits itself more than it limits the other.
#This is ecologically plausible in systems where species partition resources efficiently. 
#The vector field confirms that regardless of the initial population sizes,both species stabilize over time. 

#2. Competitive Exclusion
#We explored competitive exclusion where one species systematically drives the other to extinction. 
#In these cases, interference or exploitative competition from one species exceeds the self-limiting 
#forces within that species, leading to unstable coexistence and the eventual dominance of one competitor. 
#This outcome is consistent with the competitive exclusion principle, which states that two species 
#competing for identical resources cannot stably coexist. 
#The resulting dynamics strongly depend on the initial conditions and parameter values, especially the 
#strength of cross-species competition relative to intra-specific regulation.

#3. Priority Effects and Alternative Stable States
#We demonstrated priority effects: the final outcome of competition depends on the initial densities of the 
#species. 
#These dynamics are associated with alternative stable states, where each species can dominate depending 
#on starting conditions. 
#Here, inter-specific (alpha12, alpha21) competition exceeds intra-specific (alpha11, alpha22) competition 
#(i.e., both species are better at suppressing the other than limiting themselves), creating bistability. 
#This situation is ecologically relevant in cases such as invasive species, where early arrival or 
#slight numerical advantage can lead to monopolization of resources, even when both species have 
#otherwise similar ecological traits.






########################
### Topic 5: Lotka-Volterra model for 2 mutualistic species
rm(list=ls())

#R-code

#define the function with several parameters: time, population size and parameters (alpha's, gamma's and r's)
LV_mutualistic <- function(t, N, p) {
  N1 <- N[1]
  N2 <- N[2]
  dN1 <- N1 * (p$r1 - p$alpha11 * N1 + p$gamma12 * N2) #equation given with gamma12 and gamma21 
  #representing the mutualism terms 
  dN2 <- N2 * (p$r2 + p$gamma21 * N1 - p$alpha22 * N2)
  return(list(c(dN1, dN2)))
}

#1. For a stable equilibrium point
#pré-requis: gamma12, gamma21 < alpha11, alpha22
params_stable <- list(r1 = 1, 
                      r2 = 0.8, 
                      alpha11 = 1, 
                      alpha22 = 1, 
                      gamma12 = 0.5, 
                      gamma21 = 0.4)
#initial value (N when t=0)
N0_stable <- c(N1 = 0.5, N2 = 0.5)
#time steps at which we aim an numerical evaluation of the solution (same as above)
times_steps <- seq(0, 50, by = 0.1)
#run the numerical integrator
out_stable <- ode(y = N0_stable, times = times_steps, func = LV_mutualistic, parms = params_stable)

#plot the phase space for a stable equilibrium
plot(out_stable[, "N1"], out_stable[, "N2"], type = "l", lwd = 2, xlab = "N1", ylab = "N2",
     xlim = c(0, 6), ylim = c(0, 6), main = "Phase space - Stable equilibrium with mutualism")
points(tail(out_stable[, "N1"], 1), tail(out_stable[, "N2"], 1), pch = 16, col = "black", cex = 1.5)

#to plot the non-trivial zero growth isocline (dN1/dt = 0): N2 = (alpha11*N1 - r1)/gamma12
curve((params_stable$alpha11 * x - params_stable$r1) / params_stable$gamma12,
      from = 0, to = 8, add = TRUE, col = "green", lty = 2, lwd = 1.5)

#to plot the non-trivial zero growth isocline (dN2/dt = 0): N2 = (gamma21*N1 + r2)/alpha22
curve((params_stable$gamma21 * x + params_stable$r2) / params_stable$alpha22,
      from = 0, to = 8, add = TRUE, col = "blue", lty = 2, lwd = 1.5)

#to add a vector field
N1_seq <- seq(0.1, 6, by = 0.4)
N2_seq <- seq(0.1, 6, by = 0.4)
for (N1 in N1_seq) {
  for (N2 in N2_seq) {
    dN1 <- N1 * (params_stable$r1 - params_stable$alpha11 * N1 + params_stable$gamma12 * N2)
    dN2 <- N2 * (params_stable$r2 + params_stable$gamma21 * N1 - params_stable$alpha22 * N2)
    norm <- sqrt(dN1^2 + dN2^2)
    arrows(N1, N2,
           N1 + 0.2 * dN1 / norm,
           N2 + 0.2 * dN2 / norm,
           length = 0.05)
  }
}

#2. For an unstable equilibrium point
#pré-requis: gamma12, gamma21 > alpha11, alpha22
params_unstable <- list(r1 = 0.5, 
                        r2 = 0.5, 
                        alpha11 = 0.5, 
                        alpha22 = 0.5, 
                        gamma12 = 0.8, 
                        gamma21 = 0.8)
#initial value (when t=0)
N0_unstable <- c(N1 = 0.5, N2 = 0.5)
#time at which we aim an numerical evaluation of the solution (not the same as above)
times <- seq(0, 2, by = 0.01)
#run the numerical integrator
out_unstable <- ode(y = N0_unstable, times = times, func = LV_mutualistic, parms = params_unstable, 
                    maxsteps = 1e5)

#to plot the phase space
plot(out_unstable[, "N1"], out_unstable[, "N2"], type = "l", lwd = 2, xlab = "N1", ylab = "N2",
     xlim = c(0, 6), ylim = c(0, 6), main = "Phase space - Unstable equilibrium with mutualism")
points(tail(out_unstable[, "N1"], 1), tail(out_unstable[, "N2"], 1), pch = 16, col = "black", cex = 1.5)

#to add the non-trivial zero growth isoclines (for both species)
abline(a = params_unstable$r1 / params_unstable$gamma12,
       b = params_unstable$alpha11 / params_unstable$gamma12, col = "green", lty = 2)
abline(a = params_unstable$r2 / params_unstable$gamma21,
       b = params_unstable$gamma21 / params_unstable$alpha22, col = "blue", lty = 2)

#add a vector field
N1_seq <- seq(0.1, 6, by = 0.5)
N2_seq <- seq(0.1, 6, by = 0.5)
for (N1 in N1_seq) {
  for (N2 in N2_seq) {
    dN1 <- N1 * (params_unstable$r1 - params_unstable$alpha11 * N1 + params_unstable$gamma12 * N2)
    dN2 <- N2 * (params_unstable$r2 + params_unstable$gamma21 * N1 - params_unstable$alpha22 * N2)
    arrows(N1, N2, N1 + 0.05 * dN1, N2 + 0.05 * dN2, length = 0.05)
  }
}

#Discussion of the results

#1. Case of a stable equilibrium 
#The parameters were chosen so that the mutualistic interaction terms (gamma12, gamma21) are smaller than 
#the intra-specific competition coefficients (alpha11, alpha22). 
#This means that while each species benefits from the presence of the other, the self-limiting 
#effects dominate, preventing runaway population growth.
#The results show convergence toward a stable equilibrium point in the phase space, indicating a 
#situation of stable coexistence between the two species despite their positive interaction. 
#The vector field shows trajectories pointing inward toward the equilibrium.
#The non-trivial zero-growth iso-clines intersect in the positive quadrant, confirming the existence 
#of a feasible equilibrium.

#2. Case of an unstable equilibrium 
#In contrast, here it uses mutualism parameters (gamma12, gamma21) that exceed the 
#intra-specific competition terms (alpha11, alpha22). In other words, the positive effects of mutualism are 
#too strong compared to self-regulation (intra-specific competition). 
#This leads to explosive population growth, as seen in the phase trajectories, which rapidly diverge 
#from the origin.
#The vector field indicates divergence away from the origin and no equilibrium point, clearly 
#demonstrating that the system is dynamically unstable. 
#This instability reflects a poorly regulated obligate mutualism, where the mutual benefits cause 
#unchecked population increases—biologically unrealistic in systems with finite resources.
