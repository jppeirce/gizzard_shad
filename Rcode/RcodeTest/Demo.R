############################################
#Estimating Parameters from Available Data #
############################################

source("gizshadmodel.R")


# I.EGG PRODUCTION
# Using the batch fecundity vs length (Figure 1a) from Jons & Miranda (1997)
egg_size_data <- read.csv('./EggSizeData.csv',header=TRUE,sep=",")
# assume zero eggs below lowest length in data
zero_eggs <- data.frame(x = seq(floor(min(egg_size_data$x))), 
                        EggLengthData=rep(0,floor(min(egg_size_data$x))))
# assume max eggs for length greater than recorded in data
max_eggs <-  data.frame(x = seq(ceiling(max(egg_size_data$x)), 400), 
                        EggLengthData=max(egg_size_data$EggLengthData))
egg_size_data_ext <- rbind(zero_eggs, egg_size_data, max_eggs)

plot(egg_size_data)
plot(egg_size_data_ext)

### We fit a 3 parameter (max, slope, inflection) logit model to the eggs produced data 
# since biological observations suggest that eggs are not produced by females below size 140mm.
egg_extended_nls <- nls(EggLengthData~egg_max/(1+exp(egg_slope*(log(x)-log(egg_inf)))),
                        data = egg_size_data_ext,
                        start = list(
                          egg_max = max(egg_size_data$x),
                          egg_slope = -7,
                          egg_inf = 313
                        ))

m_par$egg_slope = coef(egg_extended_nls)[2] 
m_par$egg_max =  coef(egg_extended_nls)[1]
m_par$egg_infl = coef(egg_extended_nls)[3] 
eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  return(1000 * m_par$egg_max / (1
                                 + exp(m_par$egg_slope*
                                         (log(z)-log(m_par$egg_infl)))))
}

plot_df <- data.frame(z = zmesh, eggs =  eggs_z(zmesh, m_par)/1000 )
ggplot(data = plot_df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced") +
  #     subtitle = "Gizzard Shad")  
  scale_x_continuous(limits = c(0,400), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,700), breaks = seq(0,700,100))+
  geom_point(data = egg_size_data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

plot_df <- data.frame(z = zmesh, eggs =  eggs_z(zmesh, m_par)/1000 )
ggplot(data = plot_df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced") +
  #     subtitle = "Gizzard Shad")  
  scale_x_continuous(limits = c(0,500), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,1200), breaks = seq(0,1200,200))+
  geom_point(data = egg_size_data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

# I.EGG PRODUCTION
# Using the batch fecundity vs length (Figure 1a) from Jons & Miranda (1997)
egg_size_data <- read.csv('./EggSizeData.csv',header=TRUE,sep=",")
# assume zero eggs below lowest length in data
zero_eggs <- data.frame(x = seq(floor(min(egg_size_data$x))), 
                        EggLengthData=rep(0,floor(min(egg_size_data$x))))
# assume max eggs for length greater than recorded in data
max_eggs <-  data.frame(x = seq(ceiling(max(egg_size_data$x)), 500), 
                        EggLengthData=max(egg_size_data$EggLengthData))
egg_size_data_ext <- rbind(zero_eggs, egg_size_data, max_eggs)

plot(egg_size_data_ext)

### We fit a 3 parameter (max, slope, inflection) logit model to the eggs produced data 
# since biological observations suggest that eggs are not produced by females below size 140mm.
egg_extended_nls <- nls(EggLengthData~egg_max/(1+exp(egg_slope*(log(x)-log(egg_inf)))),
                        data = egg_size_data_ext,
                        start = list(
                          egg_max = max(egg_size_data$x),
                          egg_slope = -7,
                          egg_inf = 313
                        ))

m_par$egg_slope = coef(egg_extended_nls)[2] 
m_par$egg_max =  coef(egg_extended_nls)[1]
m_par$egg_infl = coef(egg_extended_nls)[3] 
eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  return(1000 * m_par$egg_max / (1
                                 + exp(m_par$egg_slope*
                                         (log(z)-log(m_par$egg_infl)))))
}

plot_df <- data.frame(z = zmesh, eggs =  eggs_z(zmesh, m_par)/1000 )
ggplot(data = plot_df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced") +
  #     subtitle = "Gizzard Shad")  
  scale_x_continuous(limits = c(0,500), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,700), breaks = seq(0,1200,200))+
  geom_point(data = egg_size_data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

# Stable equilibrium value
m_par_eq <- tibble(
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
  surv_alpha = 103.5187, 
  surv_beta = -1006.1837, # slope
  ## New recruit from Michaletz (2017)
  recruit_mean = 105,
  recruit_sd = 25, # same as grow_sd
  ## From Bodola (1955):
  egg_viable = 0.002,
  ## Estimated from Jons and Miranda (1997)
  egg_slope = -4.723514, # -4.361915
  egg_max =  1658.561222, # 41540.608025
  egg_infl = 417.801634, # 935.528239
  ## Spawning Probability - Estimated from Michaletz (2009)
  prob_spawn = 0.90,
  surv0_int = coef(surv_den_exp)[1], # 0.2686
  surv0_decay = coef(surv_den_exp)[2], # 0.0030
)

m_par_periodic <- tibble(
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
  surv_alpha = 107.3992, 
  surv_beta = -12, # slope
  ## New recruit from Michaletz (2017)
  recruit_mean = 105,
  recruit_sd = 25, # same as grow_sd
  ## From Bodola (1955):
  egg_viable = 0.002,
  ## Estimated from Jons and Miranda (1997)
  egg_slope = -7.181562,
  egg_max =  737.512104, # 41540.608025
  egg_infl = 313.742703, # 935.528239
  ## Spawning Probability - Estimated from Michaletz (2009)
  prob_spawn = 0.90,
  surv0_int = coef(surv_den_exp)[1], # 0.2686
  surv0_decay = coef(surv_den_exp)[2], # 0.0030
)

m_par <- m_par_eq

eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  return(1000 * m_par$egg_max / (1
                                 + exp(m_par$egg_slope*
                                         (log(z)-log(m_par$egg_infl)))))
}

n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5*m_par$grow_max, sd = 30)
#normal like LTRM 1994
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z

for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

# Population size vs time
n_total <- rep(0, tf)
for (i in 1:tf) {
  n_total[i] <- sum(n[, i]) * delta_z
}

plot_df <- data.frame(year = 1:tf, total = n_total)
ggplot(plot_df,
       aes(x = year, y = total)) + 
  geom_line() +
  labs(x = "years",
       y = "total") +
  labs(x = "time (in years)",
       y = "total density",
       title = "Total n(z,t)",
       color = "Legend") +
  theme_bw() +  
  theme(legend.position = c(0.8, 0.4))+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
