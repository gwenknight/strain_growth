### Multilevel modelling

library(tidyverse)
library(nlme)
library(lme4)

## Data
dd <- read_csv("output/mm_final_data.csv")[,-1]



## Check need for multilevel model/random intercept ##
#Fit model with intercept only
m_int <- gls(success ~ 1, data = dd, method = "ML")
summary(m_int)
#Fit model allowing intercepts to vary by strain
m_rint <-lme(logred ~ 1, random = ~1|strain_name, data = dd, method = "ML")
summary(m_rint)
#Compare models to decide if random effects are needed
anova(m_int, m_rint) #if significant > random effects needed



### what are the clusters? 
# (1) ignore groupings by inoculum
# (2) clusters = within a rep more similar than between replicates
model <- glmer(success_bin ~ logred + (1 + logred | strain_name), data=dd, family=binomial(link="logit"))
summary(model)
