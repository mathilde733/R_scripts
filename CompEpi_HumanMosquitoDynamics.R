#-----------------------
# Libraries
#-----------------------
library(deSolve) #ODE solver
library(reshape2) #For data reshaping
library(ggplot2) #For plotting

#-----------------------
# SIR model (host, human) + SI model (vector, mosquito)
#-----------------------
sir <- function(time, state, parameters) {
  # Unpack state variables and parameters
  with(as.list(c(state, parameters)), {
    # Total populations 
    Nh <- (Sh + Ih + R) # humans
    Nv <- (Sv + Iv) # mosquitoes
    # Host dynamics 
    dSh <- muh * Nh - (beta1*(1-fm*em)*b * Sh * Iv)/Nh - p*Sh - muh * Sh
    dIh <-  (beta1*(1-fm*em)*b * Sh * Iv)/Nh - (gamma + muh )* Ih 
    dRh <-                 gamma * Ih - muh * R + p*Sh
    # Vector dynamics
    dSv <- (1-c)*muv *Nv  - (beta2 *(1-fm*em)* b* Sv * Ih)/Nh -  (1 + k)*muv*Sv
    dIv <- (beta2 * (1-fm*em)*b* Sv*Ih)/Nh - (1 + k)*muv*Iv 
    # Return derivatives
    return(list(c(dSh, dIh, dRh, dSv,dIv)))
  })
}

#-----------------------
# Initial conditions (proportions)
#-----------------------
init <- c(Sh=1-1e-6, 
           Ih=1e-6, 
           R=0, 
           Sv=1-1e-6, 
           Iv=1e-6)

#-----------------------
# Parameters

# beta1 = transmission probability after biting between infected vector to susceptible human
# beta2 = transmission probability after biting between infected human to susceptible vector
# gamma = recovery rate
# muh = birth and death rate for host
# muv = birth and death rate for the vector

# Protection measures parameters
# p = proportion of vaccinated host per unit of time 
# em = efficiency of the mosquito net from 0 to 1 
# fm = fraction population with mosquito net from 0 to 1
# k = percentage of increase of the death rate due to the insecticide
# c = percentage of decrease of the birth rate due to the birth control, from 0 to 1
#-----------------------
parameters <- c(beta1 = 0.8, 
                beta2 = 0.8, 
                gamma = 0.2, 
                muh= 0.1,
                muv = 0.5,
                b= 2, 
                p = 0,
                em=0,
                fm=0,
                k=0,
                c=0)
# Time vector
times <- seq(0, 100, by = 0.01)

#-----------------------
# Model simulations: baseline
#-----------------------
out <- ode(y = init, times = times, func = sir, parms = parameters)
out <- as.data.frame(out)

# Reshape for plotting
dat <- reshape2::melt(out, id.vars = c("time"))
names(dat) <- c("time", "Group", "value")

# Add population type (human vs vector)
dat$Population[dat$Group %in% c("Sh", "Ih", "R")] <- "Human"
dat$Population[dat$Group %in% c("Sv", "Iv")] <- "Vector"

# Plot: Faceted by population
ggplot(data = dat, aes(x = time, y = value)) +
  geom_line(aes(color = Group, group = Group), size = 1.2) +
  facet_wrap(~Population, ncol = 1) +
  theme_classic(base_size = 18) +
  labs(title = "Baseline model",
       x = "Time",
       y = "Proportion of Population") +
  theme(legend.position = "bottom")

#-----------------------
# Model simulations: intervention measures
#-----------------------
# To visualize the time-series for each intervention model, choose the appropriate parameters values:
# Vaccination:
# parameters <- c(muh=0.1, muv=0.5, beta1 = 0.8, beta2 = 0.8, b=2, gamma=0.2, p=0.1, fm=0, em=0, k=0, c=0)
# Mosquito nets:
# parameters <- c(muh=0.1, muv=0.5, beta1 = 0.8, beta2 = 0.8, b=2, gamma=0.2, p=0, fm=0.7, em=0.95, k=0, c=0)
# Insecticide / Death rate control:
parameters <- c(muh=0.1, muv=0.5, beta1 = 0.8, beta2 = 0.8, b=2, gamma=0.2, p=0, fm=0, em=0, k=0.25, c=0)
# Vector Population birth rate control:
# parameters <- c(muh=0.1, muv=0.5, beta1 = 0.8, beta2 = 0.8, b=2, gamma=0.2, p=0, fm=0, em=0, k=0, c=0.25)

out <- ode(y = init, times = times, func = sir, parms = parameters)
out <- as.data.frame(out)

# Reshape for plotting
dat <- reshape2::melt(out, id.vars = c("time"))
names(dat) <- c("time", "Group", "value")

# Add population type (human vs vector)
dat$Population[dat$Group %in% c("Sh", "Ih", "R")] <- "Human"
dat$Population[dat$Group %in% c("Sv", "Iv")] <- "Vector"

#Creating title plot with the current parameters values for the intervention measures
plot_title <- paste0("Intervention measures model: p = ", parameters["p"], 
       ", em = ", parameters["em"],
       ", fm = ", parameters["fm"],
       ", k = ", parameters["k"], 
       ", c = ", parameters["c"]
       )

# Plot: Faceted by population
ggplot(data = dat, aes(x = time, y = value)) +
  geom_line(aes(color = Group), size = 1.2) +
  facet_wrap(~Population, ncol = 1) +
  theme_classic(base_size = 18) +
  labs(title = plot_title,
       x = "Time",
       y = "Proportion of Population") +
  theme(legend.position = "bottom")


