### Statistical model for dehydration experiments###

### Background ####
## RQ: Do successful strains survive better after dehydration than unsuccessful strains? Ie is there a difference in log reduction?
## H0: there is no difference in log reduction between successful/unsuccessful strains
## H1: there is a difference in log reduction -> is log reduction lower in successful strains?

## For 98 strains, log reduction is determined after 168h of dehydration, in 3 inoculum sizes for 3 replicates
## Data is imported in long format. Experimental unit = strain. Two levels: between strains and within strains
## Main variables are:
# log reduction > outcome variable, continuous
# success > fixed effect, categorical
# lineage > fixed effect, nominal
# country > fixed effect, nominal
# inoc > random effect, nominal
# rep > random effect?, nominal
# strain > level of random effect, nominal
# other possible variables: infection, year, mec, pvl, spa

##Libraries and working directory ####
# install.packages("ggplot2")
# install.packages("nlme")
# install.packages("pastecs")
# install.packages("reshape2")
# install.packages("WRS2")

##Load libraries
 library(ggplot2)
# library(nlme)
# library(pastecs)
# library(reshape2)
# library(WRS2)

#Set home
setwd(here::here())


##############################################################
# MIXED LINEAR AND LOGISTIC MODELS

library(data.table)
library(lmerTest)
library(lme4)

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
# Signal expected within CC8, possibly CC45

# Predict logred 
summary(lmer(logred ~ success + inoc + (1 | lineage) + (1|country), d))
confint.merMod(lmer(logred ~ success + inoc + (1 | lineage) + (1|country), d), level = 0.95)

summary(lmer(logred ~ success + inoc + (1 | lineage), d))
# No signal

# Whatever, predict success
summary(glmer(success ~ logred + inoc + (1 | country) + (1 | lineage), d, family = binomial))
# No signal

# Lineage-dependent association ?
summary(lmer(logred ~ success*lineage + inoc + (1 | country), d))
summary(glmer(success ~ logred*lineage + inoc + (1 | country), d, family = binomial))

# Some (very weak) signal in CC8


stop()
# ###############################################################
# 
# 
# str(ddata)
# #print(levels(ddata$inoc))
# 
# #Graph
# ggplot(ddata, aes(x=success, y=logred)) +
#   geom_point(aes(group=inoc, shape = inoc, col = inoc)) +
#   facet_wrap(~lineage, nrow = 3) +
#   scale_y_continuous("Mean log reduction") +
#   scale_x_discrete("Inoculum")
# ggsave("output/Split_success.pdf")
# 
# ggplot(ddata, aes(x=inoc, y=logred)) +
#   geom_point(aes(group=success, shape = success, col = success)) +
#   facet_wrap(~lineage, nrow = 3) +
#   scale_y_continuous("Mean log reduction") +
#   scale_x_discrete("Inoculum")
# ggsave("output/Split_inoc.pdf")
# 
# #### Check assumptions ####
# 
# 
# #### Multilevel linear model ####
# 
# ## Check need for multilevel model/random intercept ##
# #Fit model with intercept only
# m_int <- gls(logred ~ 1, data = ddata, method = "ML")
# summary(m_int)
# #Fit model allowing intercepts to vary by strain
# m_rint <-lme(logred ~ 1, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint)
# #Compare models to decide if random effects are needed
# anova(m_int, m_rint) #m_rint is significantly better -> random effects needed
# 
# ## Adding fixed effects ##
# m_rint1 <-lme(logred ~ success, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint1)
# anova(m_rint, m_rint1)
# ##>> including success makes model worse 
# m_rint2 <-lme(logred ~ lineage, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint2)
# anova(m_rint, m_rint2)
# ##>> including lineage makes model worse
# m_rint3 <-lme(logred ~ success + lineage, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint3)
# anova(m_rint, m_rint3)
# ##>> including success and lineage also makes model worse
# m_rint4 <-lme(logred ~ country, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint4)
# anova(m_rint, m_rint4)
# ##>> including country also makes model worse
# m_rint5 <-lme(logred ~ inoc, random = ~1|strain_name, data = ddata, method = "ML")
# summary(m_rint5)
# anova(m_rint, m_rint5)
# ##>> inoculum does improve model, so has effect on log reduction
# 
# # # additional fixed effects, interaction terms?
# # m_rint3 <-lme(logred ~ success + lineage + country, random = ~1|strain_name, data = ddata, method = "ML")
# # summary(m_rint3)
# # m_rint4 <-lme(logred ~ success + lineage + country + lineage:country, random = ~1|strain_name, data = ddata, method = "ML")
# # # >> Error in MEEM(object, conLin, control$niterEM) : 
# # #Singularity in backsolve at level 0, block 1
# # # >> lineage and country are confounded??
# # #summary(m_rint4)
# # anova(m_rint, m_rint1, m_rint2, m_rint3)
# 
# ## Adding random slopes > should random effects be included as fixed effects first?
# m_rsl <-lme(logred ~ inoc, random = ~inoc|strain_name, data = ddata, method = "ML")
# summary(m_rsl)
# anova(m_rint, m_rsl) #model including random slope for inoc, better than just random intercept
# anova(m_rint5, m_rsl) #model including random slope for inoc better than 
# #inoc as fixed effect +random intercept 
# ## FINAl MODEL???
# 
# # test final model including random effects with REML instead of ML, follow steps MMDA course
# m_rsl <-lme(logred ~ inoc, random = ~inoc|strain_name, data = ddata, method = "REML")
# summary(m_rsl)
# ## Final model estimates and confidence intervals
# intervals(m_rsl, 0.95)

# option to investigate interaction term effect by subset analysis for single groups of interaction variable, country
