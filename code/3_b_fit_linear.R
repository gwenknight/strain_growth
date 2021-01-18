### FIT LINEAR MODEL to predict inoculum 

library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 11))
source("code/function_linear_model.R")
dir.create(paste0("plots/fit"), showWarnings = FALSE)
dir.create(paste0("plots/output_fit"), showWarnings = FALSE)

##### READ IN DATA
### CHOOSE IN LINE WITH DATA CLEANING IN 2_analysis.R
#name_code <- "1_2_6_7_"
#name_code <- "all_(1_13)_"

param <- read.csv("output/param_labelled_repst.csv", stringsAsFactors = FALSE)[,-1]


##### Check exp growth OK across inocula of the same strain at set drying times
perc <- 0.34

pp_strain_names <- param %>%
  group_by(strain_name, rep) %>% 
  mutate(mean_peak_exp_gr = mean(cut_exp),
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
         outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0)) %>%
  ungroup() %>% 
  group_by(strain_name, rep,drytime) %>% 
  mutate(total_outside_inrep = sum(outside), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep) %>%
  ungroup() %>% 
  group_by(strain_name, rep) %>% 
  mutate(rep_dt_remove = ifelse(perc_outside > 34, 1, ifelse(total_inrep < 2,1,0)), total_rep_rem = sum(rep_dt_remove)) %>%
  ungroup() %>%
  mutate(keep_rep =ifelse((total_rep_rem == 0),rep,0)) %>%
  group_by(strain) %>%
  mutate(remove_strain = ifelse(n_distinct(keep_rep) > 2,0,1))

## e.g.
eg <- pp_strain_names %>% 
  select(strain, rep, inoc, drytime, cut_exp, mean_peak_exp_gr, outside, perc_outside, rep_dt_remove, total_rep_rem, keep_rep, remove_strain, t_m_h_flow) %>% 
  filter(strain == 11016)

ggplot(eg, aes(x=inoc, y = t_m_h_flow)) + geom_point(aes(col = factor(drytime))) + facet_wrap(~rep)

# param_exp_gr_lab <- param %>%
#   group_by(strain_name, rep) %>% 
#   mutate(mean_peak_exp_gr = mean(cut_exp),
#          mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
#          mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
#          outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0)) %>%
#   ungroup() %>%
#   group_by(strain_name,rep) %>% # group by strain and replicate
#   dplyr::mutate(sum_outside_s_r = sum(outside), total_vals = n()) %>% # number of odd datasets in this strain and replicate
#   dplyr::mutate(typ = ifelse(sum_outside_s_r < 2,0,1)) %>%
#   group_by(strain) %>%
#   dplyr::mutate(tot_typ = sum(typ), remove_strain = ifelse(tot_typ > 1,1,0))
# 
# param_exp_gr_lab %>% select(strain, rep, inoc, cut_exp, mean_peak_exp_gr, outside, sum_outside_s_r, total_vals, typ, tot_typ, remove_strain) %>% filter(strain == 11009)

#### check for exp_gr - remove if only one is "wrong" 
param_expok <- pp_strain_names %>%
  filter(remove_strain == 0) %>% # remove those strains with more than 2 wrong reps
  filter(total_rep_rem == 0) # remove those reps with more than 2 outside
  
  #param_exp_gr_lab %>%
  #filter(remove_strain < 1)
  #filter(sum_outside_s_r < 2) # some have only one in the strain and replicate

length(unique(param_expok$strain_name)) 
length(unique(param$strain_name)) 
# Only 11 strains removed: we can't use these are their exponential growth is too variable

strains_typical = unique(param_expok$strain_name) # PERFECT strains

param_exp_gr_lab <- pp_strain_names 

ggplot(param_exp_gr_lab, aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(remove_strain)),lty = "dashed") + 
  ggtitle(paste0("All strains. Red = outside of limits (",100*perc,"%)"))
ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)

# ggplot(param_exp_gr_lab %>% filter(odd_strains == 0), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
#   scale_color_manual("In the limits?", values = c("black","red")) +
#   scale_shape_discrete("Drytime") + 
#   scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
#   scale_y_continuous("Exponential growth") + 
#   facet_wrap(strain_name ~ rep, ncol = 30) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(typ))) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(typ))) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(typ)),lty = "dashed") + 
#   ggtitle(paste0("Typical strains. Red = outside of limits (",100*perc,"%)"))
# ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth_typical.pdf"),width = 20, height = 20)
# 
# ggplot(param_exp_gr_lab %>% filter(odd_strains == 1), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
#   scale_color_manual("In the limits?", values = c("black","red")) +
#   scale_shape_discrete("Drytime") + 
#   scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
#   scale_y_continuous("Exponential growth") + 
#   facet_wrap(strain_name ~ rep, ncol = 30) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(typ))) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(typ))) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(typ)),lty = "dashed") + 
#   ggtitle(paste0("Odd strains. Red = outside of limits (",100*perc,"%)"))
# ggsave(paste0("plots/exp_growth/perc_from_mean_exponential_growth_oddstrains.pdf"),width = 20, height = 20)


######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** LINEAR MODEL 
######****************************######################################################################

# Parameters where the exponential growth is within the agreed range for linear model 
param_expok <- pp_strain_names %>%
  filter(remove_strain == 0) %>% # remove those strains with more than 2 wrong reps
  filter(total_rep_rem == 0) # remove those reps with more than 2 outside
# removes 80 datasets: 6%
  

######****** Predicting **********######################################################################
ggplot(param_expok, aes(y=inocl,x = t_m_h_flow, group = strain_name, colour=factor(rep_st))) + 
  geom_point(size = 2) + facet_wrap(~strain_name) + 
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

ggplot(param_expok, aes(x=cut_timepeak,y = log10(scalar), group = strain_name, colour=factor(drytime))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","168hr drying")) 
ggsave(paste0("plots/fit/time_to_peak_all_as_linear_model.pdf"), width = 16, height =10 )

#### Fit linear model 
# Remove 10^2 and 10^6: not for reduction analysis
w26<-c(which(param_expok$inocl == 2), which(param_expok$inocl == 6))
if(length(w26) > 0){param_expok <- param_expok[-w26,]}

# If r^2 over 8 then good
r2_threshold = 0.79

reductions_fit <- fit_line_model(reps, strains, param_expok, "cut_timepeak","Time to max heat flow", R_cut = r2_threshold, plot = 1) ## plot = 1 will give the underling fit curves

# ### Check fits
ggplot(reductions_fit$fit, aes(x=strain, y = R2)) +
  geom_point() +
  ggtitle("R2 value for fit of linear model to inoculum size = a*time_to_peak + b")

ggplot(reductions_fit$fit, aes(x=R2)) + geom_histogram(binwidth = 0.02) +
  geom_vline(xintercept = 0.9) + 
  geom_vline(xintercept = 0.8, lty = "dashed")
ggsave("plots/fit/cutoff_for_r2.pdf") # only get this if set R_cut = 0 in above?

fitted_strains <- reductions_fit$fit %>% filter(R2 > r2_threshold) %>%
  select(strain) %>% unlist() %>% as.character() %>% unique()# but not all the replicates for these strains

length(unique(param_expok$strain_name)) # 87 into the function
length(unique(reductions_fit$fit$strain)) # 87 Not filtered on R2
length(fitted_strains) # 86 Filtered on R2

setdiff(unique(param_expok$strain_name),unique(reductions_fit$fit$strain))
setdiff(unique(reductions_fit$fit$strain),fitted_strains)

### Only those that have good R^2
reductions_fit$reductions$rep <- reps[reductions_fit$reductions$rep]

# in reductions. Meas column key: 
# meas = 1 = log reduction 
# meas = 2 = percentage reduction 
# meas = 3 = inoculum
# meas = 4 = predicted inoculum

ggplot(reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_wrap(~strain_name, nrow = 10) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = pmax(0,mean - sd)),position = "dodge")
ggsave(paste0("plots/fit/log_reductions_byrep.pdf"), width = 15, height = 10)


subset(reductions_fit$reductions, meas == 1) %>% filter(strain_name == "11016") # had negative if use old time peak but still here negative at 105


#### ISSUE IS WHAT IS THE ERRORBAR? mean of sd or sd of mean? 
## Similar to reductions_fit$av_for_strain except for sd calc
reductions_fit$reductions %>% filter(meas == 1) %>% group_by(strain_name, dry) %>% 
  summarise(mean_over_rep = mean(mean), sd_over_rep = sd(mean)) %>%
  ggplot(aes(x=strain_name, y = mean_over_rep, fill = strain_name)) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(~dry) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Strain") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_discrete("Strain") + 
  geom_errorbar(aes(ymax = mean_over_rep + sd_over_rep, ymin = mean_over_rep - sd_over_rep),position = "dodge")
ggsave(paste0("plots/fit/log_reductions_mean_by_drytime.pdf"), width = 20, height = 10)

ggplot(subset(reductions_fit$reductions, meas == 2), aes(x=ticker, y = mean, fill = factor(ticker))) + 
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


# Average over all strains by inoculum
av_all <- reductions_fit$reductions %>%
  filter(meas == 1, strain_name %in% fitted_strains) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(name) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), se = sd(value, na.rm = TRUE)/n())

ggplot(av_all, aes(x=name, y = mean)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))

## Success? 
succ <- read.csv("data/success_yn.csv")
succ <- read.csv("data/MACOTRA 100collection success_20210107.csv")
succ$strain_name <- as.character(succ$strain)

succ_go <- left_join(reductions_fit$reductions %>% filter(strain_name %in% fitted_strains),succ, by = "strain_name")

av_all_bys <- succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(name, success) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), se = sd(value, na.rm = TRUE)/sqrt(n))

ggplot(av_all_bys, aes(x=name, y = mean)) + geom_bar(stat = "identity", aes(fill = success)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) + 
  facet_wrap(~success)

ggplot(av_all_bys, aes(x=name, y = mean,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - sd, ymax = mean + sd)) + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/fit/succ_unsucc_sd.pdf")

ggplot(av_all_bys, aes(x=name, y = mean,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/fit/succ_unsucc_se.pdf")

### By success only not inoculum
av_all_byso <- succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(success) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), se = sd(value, na.rm = TRUE)/sqrt(n))

ggplot(av_all_byso, aes(x=success, y = mean,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  scale_x_discrete("") + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("")
ggsave("plots/fit/succ_unsucc_se_only_succ.pdf")

### LINEAGE 

ggplot(succ_go %>% filter(meas == 1) %>%
         pivot_longer(cols = c('10^3':'10^5')), aes(x=name, y = value, group = success)) + 
  geom_point(aes(col = success), position = position_dodge(0.8)) + 
  facet_wrap(~lineage) + 
  scale_y_continuous("Log reduction") + 
  scale_x_discrete("Inoculum")
ggsave("plots/fit/lineage_all_points.pdf")

### Average by lineage and success
av_all_byl <- succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(name, lineage, success) %>% 
  filter(value < 4) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), se = sd(value, na.rm = TRUE)/sqrt(n))

ggplot(av_all_byl, aes(x=name, y = mean, group = success)) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  facet_wrap(~lineage) + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/fit/lineage_succ_unsucc.pdf")

ggplot(av_all_byl, aes(x=name, y = mean, group = success)) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  facet_wrap(~lineage) + 
  geom_text(data = av_all_byl, aes(x = name, y = mean + 1, label = n, group = success), position = position_dodge(0.8)) + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/fit/lineage_succ_unsucc_labelled.pdf")

succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>% 
  filter(value < 0) %>%#, lineage == "CC5") %>%
  dplyr::select(strain_name, rep, name, value, meas)
  
  
  
  dplyr::select(strain_name) %>%
  unique()


av_all_byl4 <- succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(name, lineage, success) %>% 
  filter(value < 2, value > 0) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), se = sd(value, na.rm = TRUE)/sqrt(n))

ggplot(av_all_byl4, aes(x=name, y = mean, group = success)) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  facet_wrap(~lineage) 
ggsave("plots/fit/lineage_succ_unsucc_remove_top.pdf")


### Average by country and success
av_all_byc <- succ_go %>%
  filter(meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(name, country, success) %>% 
  filter(value < 4) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), se = sd(value, na.rm = TRUE)/sqrt(n))

ggplot(av_all_byc, aes(x=name, y = mean, group = success)) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean - se, ymax = mean + se)) + 
  facet_wrap(~country) 
ggsave("plots/fit/country_succ_unsucc.pdf")


# Store tables
write.csv(reductions_fit$reductions,"output/fit_reductions_all_data.csv")
write.csv(reductions_fit$av_for_inoculum, "output/fit_av_by_inoc.csv")
write.csv(reductions_fit$av_for_strain, "output/fit_av_by_strain.csv")




