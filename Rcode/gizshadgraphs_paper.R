### Graph And Analysis of Gizzard Shad Model
source("gizshadmodel_paperOLD.R")
#######################
## Simulation Fun Time!
#######################

## Normal Distribution ##
N <- 50 # number of size classes
l_shad <- 0.00   # lower size limit in mm
u_shad <- 400.0    # upper size limit in mm - we want this to be
                   # larger than L-infty
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
tf <- 200 # number of years

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

plot_df <- data.frame(year = 1:tf, total = n_total)
ggplot(plot_df,
       aes(x = year, y = total)) + 
  geom_line() +
  labs(x = "years",
       y = "total") +
  theme_classic()

#####################################################
## Length FREQUENCY distributions for first few years
#####################################################
show_years <- 1:4
#show_years <- c(5, 15, 25)
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
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/presentation/figures/sim.png")


# plot.df <- data.frame(z = zmesh, density = n[,1:5]/sum(n[,1:5]) )
# ggplot(data = plot.df,
#        aes( x = z, y = density.1, color = "t=0"))+
#   geom_line(size =1)+
#   geom_line(data = plot.df,
#             aes( x = z, y = density.2, color="t=1"))+
#   geom_line(data = plot.df,
#             aes( x = z, y = density.3, color="t=2"))+
#   geom_line(data = plot.df,
#             aes( x = z, y = density.4, color="t=3"))+
#   #  geom_line(data = plot.df,
#   #            aes( x = z, y = density.5, color="purple"))+
#   labs(color = "Year",
#        x = "length (in mm)",
#        y = "length frequency",
#        title = "Length Distribtion")+ 
#   scale_x_continuous(limits = c(0,u_shad), breaks = seq(0,500,100), 
#                      expand = c(0,0))+
#   scale_y_continuous(limits = c(0,max(plot.df[,2:6])), breaks = seq(0,max(plot.df[,2:6]),.005), 
#                      expand = c(0,0))+
#   scale_color_manual(breaks = c("t=0", "t=1", "t=2", "t=3"),
#                      values = c("t=0" = "blue", "t=1" = "green", 
#                                 "t=2" = "orange", "t=3" = "red")) +
#   theme_classic()+
#   theme(text = element_text(size=16),
#         aspect.ratio = .7)


### surv vs t
# survival0 vs time
surv_t <- rep(0, times = tf)
for (i in 1:(tf - 1)) {
  surv_t[i] <- surv_age0(n = n[, i], z = zmesh, m_par)
}

plot_df <- tibble(time_years = 1:(tf+1), prob = surv_t[1:(tf+1)])
ggplot(data = plot_df,
       aes(x = time_years, y = prob)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "time (in years)",
       y = "probability of survival",
       title = "Survival Probability of Age-0")+
  #       subtitle = paste("Inital Total Density =", n0_total / 1000,
  #                        "per m^3")) +
  scale_x_continuous(limits = c(10, tf-10), breaks = seq(0, tf-10, 5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.002),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/Figure2a.pdf")


# lambda vs time
lambda_t <- rep(0, tf)
for (i in 1:(tf - 1)) {
  lambda_t[i] <- n_total[i + 1] / n_total[i]
}

plot_df <- tibble(time_years = 1:50, lambda = lambda_t[1:50])
ggplot(data = plot_df,
       aes(x = time_years, y = lambda)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 0, intercept =1) +
  labs(x = "time (in years)",
       y = "lambda")+
  scale_x_continuous(limits = c(10, 50), breaks = seq(0, 50, 5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.75, 1.25), expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        aspect.ratio = .7)
ggsave("~/OneDrive - University of Wisconsin-La Crosse/GizzardShad/paper/figures/Figure2b.pdf")



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
  #  theme_bw()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

# Determines the relative frequency of each length group 
ltrm_gzsd_data <- ltrm_gzsd_lg %>% 
  group_by(length_round) %>%
  summarize(count = n(), .groups = 'drop') %>%
  mutate(rfreq = count/sum(count)) %>%
  arrange(length_round)

ggplot(data = ltrm_gzsd_data, aes(x = length_round, y=count)) +
  geom_point()

ggplot(data = ltrm_gzsd_lg, aes(x = length_round)) +
  geom_histogram(aes(y = ..count..), 
                 bins = 50,
                 position="dodge",
                 color="black",
                 fill = "#800000")+
  labs(x = "length (in mm)",
       y = "count") +
  theme(text = element_text(size=16),
        aspect.ratio = .7)
