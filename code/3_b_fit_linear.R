### FIT LINEAR MODEL to predict inoculum 

library(tidyverse)
theme_set(theme_bw(base_size = 11))
source("function_linear_model.R")
dir.create(paste0("plots/fit"), showWarnings = FALSE)

##### READ IN DATA
### CHOOSE IN LINE WITH DATA CLEANING IN 2_analysis.R
#name_code <- "1_2_6_7_"
#name_code <- "all_(1_13)_"

param <- read.csv("output/param_labelled_repst.csv", stringsAsFactors = FALSE)[,-1]


##### Check exp growth OK across inocula of the same strain at set drying times
perc <- 0.3
param_exp_gr_lab <- param %>%
  group_by(strain_name, rep) %>% 
  mutate(mean_peak_exp_gr = mean(cut_exp),
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
         outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0)) %>%
  ungroup() %>%
  group_by(strain_name,rep) %>% # group by strain and replicate
  dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate

#### check for exp_gr - remove if only one is "wrong" 
param_expok <- param_exp_gr_lab %>%
  filter(sum_outside_s_r < 2) # some have only one in the strain and replicate

length(unique(param_expok$strain_name)) 
length(unique(param$strain_name)) 
# Only 13 strains removed 

strains_typical = unique(param_expok$strain_name) # PERFECT strains

param_exp_gr_lab <- param_exp_gr_lab %>% mutate(typ = ifelse(sum_outside_s_r < 2,0,1))

ggplot(param_exp_gr_lab, aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(typ)),lty = "dashed") + 
  ggtitle(paste0("All strains. Red = outside of limits (",100*perc,"%)"))
ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)

ggplot(param_exp_gr_lab %>% filter(odd_strains == 0), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(typ)),lty = "dashed") + 
  ggtitle(paste0("Typical strains. Red = outside of limits (",100*perc,"%)"))
ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth_typical.pdf"),width = 20, height = 20)

ggplot(param_exp_gr_lab %>% filter(odd_strains == 1), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(typ))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(typ)),lty = "dashed") + 
  ggtitle(paste0("Odd strains. Red = outside of limits (",100*perc,"%)"))
ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth_oddstrains.pdf"),width = 20, height = 20)


######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** LINEAR MODEL 
######****************************######################################################################

# Parameters where the exponential growth is within the agreed range for linear model 
param_expok <- param_exp_gr_lab %>%
  filter(sum_outside_s_r < 2) # some have only one in the strain and replicate


######****** Predicting **********######################################################################
ggplot(param_expok, aes(y=inocl,x = t_m_h_flow, group = strain_name, colour=factor(rep_st))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_x_continuous("Time to peak") + scale_y_continuous("Inoculum size (10^x)") + 
  scale_color_discrete("Experiment") 
ggsave(paste0("plots/fit/time_to_peak_all.pdf"), width = 16, height =10 )

### PLOT variation in other indicators
#ggplot(param_typical, aes(y=inocl,x = v_m_h_flow, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = lag, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = exp_gr, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = auc, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)

######****** Linear model **********######################################################################
# Do for each rep, for each strain_name
reps <- unique(param_expok$rep)
strains <- unique(param_expok$strain_name)

ggplot(param_expok, aes(x=t_m_h_flow,y = log10(scalar), group = strain_name, colour=factor(drytime))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","168hr drying")) 
ggsave(paste0("plots/fit/time_to_peak_all_as_linear_model.pdf"), width = 16, height =10 )

#### Fit linear model 
# Remove 10^2 and 10^6: not for reduction analysis
w26<-c(which(param_expok$inocl == 2), which(param_expok$inocl == 6))
if(length(w26) > 0){param_expok <- param_expok[-w26,]}

reductions_fit <- fit_line_model(reps, strains, param_expok, "t_m_h_flow","Time to max heat flow")

# in reductions. Meas column key: 
# meas = 1 = log reduction 
# meas = 2 = percentage reduction 
# meas = 3 = inoculum
# meas = 4 = predicted inoculum

ggplot(subset(reductions_fit$reductions, meas == 1), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_wrap(~strain_name, nrow = 10) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = pmax(0,mean - sd)),position = "dodge")
ggsave(paste0("plots/fit/log_reductions_byrep.pdf"), width = 15, height = 10)

#### ISSUE IS WHAT IS THE ERRORBAR? mean of sd or sd of mean? 
## Similar to reductions_fit$av_for_strain except for sd calc
reductions_fit$reductions %>% subset(meas == 1) %>% group_by(strain_name, dry) %>% 
  summarise(mean_over_rep = mean(mean), sd_over_rep = sd(mean)) %>%
  ggplot(aes(x=strain_name, y = mean_over_rep, fill = strain_name)) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(~dry) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Strain") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_discrete("Strain") + 
  geom_errorbar(aes(ymax = mean_over_rep + sd_over_rep, ymin = mean_over_rep - sd_over_rep),position = "dodge")
ggsave(paste0("plots/fit/log_reductions_mean_by_drytime.pdf"), width = 20, height = 10)

ggplot(subset(reductions, meas == 2), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Percentage reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave(paste0("plots/fit/",name_code,"perc_reductions.pdf"))


# By inoculum 
ggplot(reductions_fit$av_for_inoculum, aes(x = name, y = mean_inoc)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  facet_wrap(~strain_name, scales = "free") + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 3, lty = "dashed") + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Log reduction by inoculum")
ggsave(paste0("plots/fit/log_reduction_by_inoc.pdf"))

ggplot(reductions_fit$av_for_inoculum, aes(x = name, y = mean_inoc)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  facet_wrap(~strain_name) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 3, lty = "dashed") + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Log reduction by inoculum")
ggsave(paste0("plots/fit/log_reduction_by_inoc_fixed_scales.pdf"))

# Store tables
write.csv(reductions_fit$reductions,"output/fit_reductions_all_data.csv")
write.csv(reductions_fit$av_for_inoculum, "output/fit_av_by_inoc.csv")
write.csv(reductions_fit$av_for_strain, "output/fit_av_by_strain.csv")


