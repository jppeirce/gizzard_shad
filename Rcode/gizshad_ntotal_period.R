# compute period of n_total
## Model Simulation Setup ##
# Normal Distribution 
N <- 50 # number of size classes
l_shad <- 0.00   # lower size limit in mm
u_shad <- 500.0    # upper size limit in mm - we want this to be
# larger than L-infty
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
tf <- 1000 # number of years

# Initial length distribution
n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5*m_par$grow_max, sd = 30)
#normal like LTRM 1994
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Note: sum(n[,1])*delta_z = n0_total

# Dynamical System
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

# Make tibble of many periods at/near stable state
many_period <- tibble(
  t = 800:900,
  n_tot = n_total[800:900]
)
# fit sin curve to many_period
n_total_period <- nls(
  n_tot~(mean(many_period$n_tot)+(max(many_period$n_tot)-mean(many_period$n_tot))*sin((2*pi/per)*(t-shift))),
                      data = many_period,
                      start = list(
                        per = 9,
                        shift = 5))
n_total_model <- function(t){
  per <- coef(n_total_period)[1]
  shift <- coef(n_total_period)[2]
  return(mean(many_period$n_tot)+(max(many_period$n_tot)-mean(many_period$n_tot))*sin((2*pi/per)*(t-shift)))
}

plot_df <- data.frame(z = many_period$t, n_tot =  n_total_model(many_period$t))
ggplot(data = plot_df,
       aes( x = z, y = n_tot)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "time (in years)",
       y = "density")+
#  scale_x_continuous(limits = c(0,as.numeric(m_par["grow_max"])), breaks = seq(0,500,100))+
#  scale_y_continuous(limits = c(0,700), breaks = seq(0,700,100))+
  geom_point(data = many_period,
             aes(x = t, y = n_tot),
             color = "black")+
  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
coef(n_total_period)

year_start <- 90
year_end <- year_start + floor(coef(n_total_period)[1]) -1
show_years <- year_start:year_end
n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
plot_df <- data.frame(z = zmesh, year = n_freq[,show_years])
plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'year')
plot_average <- tibble(z = zmesh, 
                       year = rep("mean",N),
                       n_freq = (1/floor(coef(n_total_period)[1]))*rowSums(n_freq[,year_start:year_end]) )
#plot_df <- bind_rows(plot_df, plot_average)
ggplot(plot_df,
       aes(z, value)) + 
  geom_line(aes(color = year)) + 
  geom_line(plot_average,
            mapping = aes(x=z,y=n_freq))+
                          #color = "mean"))+
  #          color = "black", size = 1)+
  labs(x = "length (in mm)",
       y = "relative frequency",
       title = "n(z,t)/total",
       color = "Legend") +
  scale_colour_manual(values = rainbow(8))+
  theme_bw() +  
  theme(legend.position = c(0.8, 0.4))+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

ltrm_gzsd_lg <- ltrm_gzsd %>% 
  mutate(length_round = round(ltrm_gzsd$length, -1))  %>% 
  filter(pool == "LG")
ggplot(data = ltrm_gzsd_lg, aes(x = length_round)) +
  geom_histogram(aes(y = ..density..), bins = 50)+
  #  geom_density(aes(x = length)) +
  geom_line(data = plot_average, aes(x = z, y = n_freq/delta_z))+
  labs(x = "length (in mm)",
       y = "frequency",
       color = "Legend") +
  scale_color_manual(breaks=c("a","b"))+
  #  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)
