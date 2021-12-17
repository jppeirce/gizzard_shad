############################################
#Estimating Parameters from Available Data #
############################################
library(drc) # used for logit model fitting

# I.EGG PRODUCTION
# Using the batch fecundity vs length (Figure 1a) from Jons & Miranda (1997)
egg_size_data <- read.csv('./EggSizeData.csv',header=TRUE,sep=",")
# assume zero eggs below lowest length in data
zero_eggs <- data.frame(x = seq(floor(min(egg_size_data$x))), 
                        EggLengthData=rep(0,floor(min(egg_size_data$x))))
# assume max eggs for length greater than recorded in data
max_eggs <-  data.frame(x = seq(ceiling(max(egg_size_data$x)), floor(max(zmesh))), 
                        EggLengthData=max(egg_size_data$EggLengthData))
egg_size_data.ext <- rbind(zero_eggs, egg_size_data, max_eggs)

### USE LL3. We fit a 3 parameter (max, slope, inflection) logit model to the eggs produced data 
# since biological observations suggest that eggs are not produced by females below size 140mm.
egg_extended_m3 <- drm(formula=EggLengthData~x, data=egg_size_data.ext,
                  fct = LL.3(names=c("slope", "max","infection")))

# II. Survival age-0
# import data from Michaletz paper
Michaletz_data <- read.csv('./Michaletz2009.csv',header=TRUE,sep=",")
names(Michaletz_data) <- c('density', 'survival')
surv_den_exp <- nls(survival/100~s0*exp(-b*density), 
                    data = Michaletz_data,
                    start = list(s0=0.3, b=0.01))