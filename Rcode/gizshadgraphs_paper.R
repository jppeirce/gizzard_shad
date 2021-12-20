### Graph And Analysis of Gizzard Shad Model
# Remember to set the working directory:
# Session -> Set Working Directory -> To Source Location
source("gizshadmodel.R")

#######################
## Graphing Fun Time!
#######################

## Parameter Graphs
plot_df <- data.frame(z = zmesh, eggs =  eggs_z(zmesh, m_par)/1000 )
ggplot(data = plot_df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced") +
  #     subtitle = "Gizzard Shad")  
  scale_x_continuous(limits = c(0,as.numeric(m_par["grow_max"])), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,700), breaks = seq(0,700,100))+
  geom_point(data = egg_size_data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/eggs.png")

dmesh <- seq( from = min(Michaletz_data$density), 
              to =max(Michaletz_data$density), 
              length.out = length(Michaletz_data$density) )

plot_df <- data.frame(x = dmesh , prob = surv_density(dmesh, m_par) )
ggplot(data = plot_df,
       aes( x = x, y = prob))+
  geom_line(color = "blue", size = 1)+
  labs(x = "density (age-0 per 1000 m^3)",
       y = "probability of survival",
       title = "Survival Probability for Age-0")+
  scale_x_continuous(limits = c(0,max(Michaletz_data$density)), 
                     breaks = seq(0,max(Michaletz_data$density),200))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.1))+
  geom_point(data = Michaletz_data, 
             aes(x = density, y = survival/100))+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/age0surv.png")


## Model Simulation Setup ##
# Normal Distribution 
N <- 50 # number of size classes
l_shad <- 0.00   # lower size limit in mm
u_shad <- 450.0    # upper size limit in mm - we want this to be
                   # larger than L-infty
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
tf <- 500 # number of years

# Initial length distribution
n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5*m_par$grow_max, sd = 30)
    #normal like LTRM 1994
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Note: sum(n[,1])*delta_z = n0_total

# Dynamical System
# update survival alpha and beta parameters if needed
# use if gizshadmode_surv_a.R was recently run
m_par$surv_alpha <- opt$par[1]
m_par$surv_beta <- opt$par[2]

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
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/ntotal.png")

#####################################################
## Length FREQUENCY distributions for first few years
#####################################################
show_years <- 1:4
#show_years <- c(160, 170, 182)
#show_years <- 41:51
n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
plot_df <- data.frame(z = zmesh, year = n_freq[,show_years])
plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'year')
ggplot(plot_df,
       aes(z, value)) + 
  geom_line(aes(color = year)) + 
  labs(x = "length (in mm)",
       y = "relative frequency",
       title = "n(z,t)/total",
       color = "Legend") +
  theme_bw() +  
  theme(legend.position = c(0.8, 0.4))+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/sim.png")

### Age-0 survival vs time
surv_t <- rep(0, times = tf)
for (i in 1:tf) {
  surv_t[i] <- surv_age0(n = n[, i], z = zmesh, m_par)
}

plot_df <- tibble(time_years = 1:tf, prob = surv_t[1:tf])
ggplot(data = plot_df,
       aes(x = time_years, y = prob)) +
  geom_line()+
  labs(x = "time (in years)",
       y = "probability of survival",
       title = "Survival Probability of Age-0")+
  #       subtitle = paste("Inital Total Density =", n0_total / 1000,
  #                        "per m^3")) +
  scale_x_continuous(limits = c(5, tf), breaks = seq(0, tf, 50)) +
  scale_y_continuous(limits = c(0, .04), breaks = seq(0, 0.4, 0.005)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/Figure2a.pdf")

#############################
########
ltrm_gzsd <- read_csv("ltrm_fish_data.csv")
# Remove length 0 and NA
ltrm_gzsd <- ltrm_gzsd[!is.na(ltrm_gzsd$length) & (ltrm_gzsd$length > 0) &
                         !is.na(ltrm_gzsd$fdate), ]
# Convert date into new format
# Then pull year and add it as a new column
ltrm_gzsd$fdate <- as.Date(ltrm_gzsd$fdate, "%m/%d/%Y")
ltrm_gzsd <- ltrm_gzsd %>% mutate(year = year(fdate))
ltrm_gzsd <- ltrm_gzsd %>%
  filter(year != 2107)

## Model equilibrium n(z,100) computed with LTRMP data from La Grange
#### Now Check with La Grange
n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5 * m_par$grow_max, sd = 30)
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Dynamical System
for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

plot_df <- data.frame(z = zmesh, freq = n[, tf]/sum(n[, tf]))
ltrm_gzsd_lg <- ltrm_gzsd %>% 
  mutate(length_round = round(ltrm_gzsd$length, -1))  %>% 
  filter(pool == "LG")
ggplot(data = ltrm_gzsd_lg, aes(x = length_round)) +
  geom_histogram(aes(y = ..density..), bins = 49)+
  #  geom_density(aes(x = length)) +
  geom_line(data = plot_df, aes(x = z, y = freq/delta_z))+
  labs(x = "length (in mm)",
       y = "frequency",
       color = "Legend") +
  scale_color_manual(breaks=c("a","b"))+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/lagrange.pdf")

