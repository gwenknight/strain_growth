### Statistical model for dehydration experiments###

### Background ####
## RQ: Do successful strains survive better after dehydration than unsuccessful strains? Ie is there a difference in log reduction?
## H0: there is no difference in log reduction between successful/unsuccessful strains
## H1: there is a difference in log reduction -> is log reduction lower in successful strains?

##Libraries and working directory ####
# install.packages("ggplot2")
# install.packages("data.table")
# install.packages("lmerTest")
# install.packages("lme4")

##Load libraries
library(ggplot2)
library(data.table)
library(lmerTest)
library(lme4)

#Set home
setwd(here::here())

##############################################################
# MIXED LINEAR AND LOGISTIC MODELS

d <- fread("output/mm_final_data.csv")[, -1]
names(d)

d
lapply(d, class)

# Marshalling
d[, id := factor(strain_name)]
d[, country := factor(country)]
d[, lineage := factor(lineage)]
d[, inoc := log10(as.numeric(gsub("\\^", "e", inoc)))]
d[, success := success_bin == 1]

# Distribution check
hist(d$logred)
summary(d$logred)
# OK, vaguely symmetric, go on

# Marginal dependence checks
boxplot(logred ~ inoc, d) # logred decreased with higher inoculum
boxplot(logred ~ success, d) # independent
boxplot(logred ~ success + inoc, d) # independent conditional on inoculum
boxplot(logred ~ lineage, d) # Wide variations
par(mar = c(4, 12, 2, 2))
boxplot(logred ~ success + lineage, d, las = 2, horizontal = TRUE)

# Predict logred 
summary(lmer(logred ~ success + inoc + (1 | lineage) + (1|country), d)) 
confint.merMod(lmer(logred ~ success + inoc + (1 | lineage) + (1|country), d), level = 0.95)
# No signal for success, but effect of inoc
# summary(lmer(logred ~ success + inoc + (1 | lineage), d))
# No signal

# # Predict success
# summary(glmer(success ~ logred + inoc + (1 | country) + (1 | lineage), d, family = binomial))
# # No signal
# 
# # Lineage-dependent association?
#summary(lmer(logred ~ success*lineage + inoc + (1 | country), d))
#summary(glmer(success ~ logred*lineage + inoc + (1 | country), d, family = binomial))
# No signal


stop()

