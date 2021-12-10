## IPM for Gizzard Shad (Fall 2021)
# Remember to set the working directory:
# Session -> Set Working Directory -> To Source Location
library(tidyverse)  #used to make pretty graphs
library(lubridate) # used to change date format of data
library(reshape2)

##############################################################
## Section 1 - Define the demographic functions and parameters
##############################################################
m_par <- tibble(
  ## growth from Michaletz (2017) paper - Table 1
  grow_rate = 0.26, # growth rate
  grow_max  =   394.3, # maximum length in mm: L_inf = 394.30
  grow_sd   =   25,  # growth sd (Max TL - L_inf)
  surv_min  =  0.002, # min survival - Bodola (1955)
  # Then et al (2015): surv_max computed from 
  # 1-natural mortality = 1 - 8.872*K^.73 L^-.33
  surv_max = 1 - 8.872*grow_rate^.73*grow_max^(-.33), 
  # inflection point: will be temp dependent
  # computed for La Grange Reach
  surv_alpha = 103.5261, 
  surv_beta = -943.889, # slope
  ## New recruit from Michaletz (2017)
  recruit_mean = 105,
  recruit_sd = 25, # same as grow_sd
  ## From Bodola (1955):
  egg_viable = 0.002,
  ## Estimated from Jons and Miranda (1997)
  egg_slope = -4.361915,
  egg_max = 41540.608025,
  egg_infl = 935.528239,
  ## Spawning Probability - Estimated from Michaletz (2009)
  prob_spawn = 0.90,
  surv0_int = 0.2686,
  surv0_decay = 0.0030
  )

##########################
## Section 2: Model Set-up
##########################

## Growth function
# Given you are size z now returns the pdf of size z1 next time
#  - Computed from von Bertanaffy equation z(t) = L_inf(1-e^K(t-t0))
# -  To find z(t+1) = L_inf*(1-e^(-K)) + e^(-K)*z(t)

g_z1z <- function(z1, z, m_par) {
  mu <- m_par$grow_max * (1 - exp(- m_par$grow_rate)) +
    exp(-m_par$grow_rate) * z           # mean size next year
  sig <- m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

## Adult Survival function, 4-parameter logistic
s_z <- function(z, m_par) {
  m_par$surv_min + (m_par$surv_max - m_par$surv_min) /
    (1 + exp(m_par$surv_beta * (log(z) - log(m_par$surv_alpha))))
}

## Reproduction, 3-parameter logistic
eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  return(1000 * m_par$egg_max / (1
                                 + exp(m_par$egg_slope*
                                         (log(z)-log(m_par$egg_infl)))))
}

## Recruit size pdf
c_1z1 <- function(z1, m_par) {
  mu <- m_par$recruit_mean
  sig <- m_par$recruit_sd
  p_den_recruit <- dnorm(z1, mean = mu, sd = sig)
  p_den_recruit <- p_den_recruit / (sum(p_den_recruit) * delta_z)
  return(p_den_recruit)
}

#####################################################
## Section 3 - Functions to build IPM kernels F and P
#####################################################

## Fecundity Kernel
# Density-depoendent age-0 survival
surv_density <- function(d, m_par) { # probability of survival of age-0 fish dependent
  # on density
  return(m_par$surv0_int * exp(- m_par$surv0_decay * d))
}

surv_age0 <- function(n, z, m_par) { # survival age-0
  # distribution of VIABLE age0 from current population n
  age0_dist <- m_par$egg_viable * eggs_z(z, m_par) *
    m_par$prob_spawn * n
  age0_density <- 10 ** (-3) * (sum(age0_dist * delta_z))
  return(surv_density(age0_density, m_par))
}

f_z1z <- function(z1, z, n, m_par) {
  age1_dist <- m_par$prob_spawn * eggs_z(z, m_par) *
    m_par$egg_viable * surv_age0(n, z, m_par)
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(c_1z1(z1, m_par), age1_dist))
}

## Growth and Survival Kernel
p_z1z <- function(z1, z, m_par) {
  g_matrix <- matrix(0, N, N)
  for (x in 1:N) {
    g_matrix[, x] <- g_z1z(z, rep(z[x], times = N), m_par)
    g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * delta_z)
  }
  return(g_matrix %*% diag(s_z(z, m_par)))
}