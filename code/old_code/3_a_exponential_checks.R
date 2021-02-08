### TYPICAL STRAIN ANALYSIS
### From 3_next_steps.R have typical param set
library(tidyverse)
theme_set(theme_bw(base_size = 11))
source("function_linear_model.R")

##### READ IN DATA
### CHOOSE IN LINE WITH DATA CLEANING IN 2_analysis.R
#name_code <- "1_2_6_7_"
name_code <- "all_(1_13)_"

param_typical_notexp_gr <- read.csv(paste0("output/",name_code,"param_typical.csv"), stringsAsFactors = FALSE)[,-1]
param_typical_notexp_gr <- as.data.frame(param_typical_notexp_gr)

##### FOR THOSE WITH ONLY ONE ODD BEHAVIOUR CHECK EXPONENTIAL GROWTH
##### Check exp growth OK across inocula of the same strain at set drying times
perc <- 0.2
param_exp_gr_lab <- param_typical_notexp_gr %>%
  group_by(strain_name, rep, drytime) %>% 
  mutate(mean_peak_exp_gr = mean(exp_gr),
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
         outside = ifelse(exp_gr < mean_peak_exp_gr_m10 | exp_gr > mean_peak_exp_gr_p10, 1, 0)) %>%
  ungroup() %>%
  group_by(strain_name,rep) %>% # group by strain and replicate
  dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate

#### check for exp_gr - remove if only one is "wrong" 
param_typical <- param_exp_gr_lab %>%
  filter(sum_outside_s_r < 2) # some have only one in the strain and replicate

length(unique(param_typical$strain_name)) 
# For 1_2_6_ 9 strains still have typical behaviour

strains_typical = unique(param_typical$strain_name) # PERFECT strains

dim(param_typical_notexp_gr)[1] - dim(param_typical)[1] # 79 datasets removed because of exp_gr - one of those from the one strain that is removed ("RN6390B")

ggplot(param_exp_gr_lab, aes(x=inocl, y = exp_gr)) + geom_point(aes(colour = factor(outside))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  facet_wrap(strain_name ~ rep + drytime, scales = "free") +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10)) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10)) +
  geom_hline(aes(yintercept = mean_peak_exp_gr),lty = "dashed") +
  ggtitle(paste0(100*perc," percent from mean exponential growth"))
ggsave(paste0("plots/",name_code,100*perc,"perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)

### OUTPUT WHAT HAPPENS TO STRAINS
infor <- read.csv(paste0("output/",name_code,"information_all_model_fit_params.csv"), stringsAsFactors = FALSE)
infor <- as.data.frame(infor)
infor$rep <- round(infor$rep,1)
## Label all the datasets removed due to exp_gr issues
infor$removed_typ_but_exp_gr <- 0 
w<-which(infor$typical == 1)
w2<-intersect(w, union(which(infor$removed_dataset == 0),which(infor$removed_rep == 0)))
infor[w2,"removed_typ_but_exp_gr"] <- 1 # Label all typical as removed
for(i in 1:dim(param_typical)[1]){
  wst <- which(infor$strain_name == as.character(param_typical[i,"strain_name"]))
  wstr <- intersect(wst, which(infor$rep == as.numeric(param_typical[i,"rep"])))
  wstrt <- intersect(wstr, which(infor$drytime == as.numeric(param_typical[i,"drytime"])))
  wstrti <- intersect(wstrt, which(infor$inocl == as.numeric(param_typical[i,"inocl"])))
  print(c(i,wstrti))
  infor[wstrti,"removed_typ_but_exp_gr"] <- 0 # say not removed if in typical
}
## CHECKS
# w<-which(infor$removed_typ_but_exp_gr == 1)
# unique(infor[w,"strain_name"]) # which strains have some datasets removed because of exp_gr
# 
# infor %>% filter(removed_dataset == 0, removed_rep == 0, typical == 1, removed_typ_but_exp_gr == 0) %>%
#   summarise(unique(strain_name))
# 
# unique(param_exp_gr_lab$strain_name)
# unique(param_typical$strain_name)

write.csv(infor,paste0("output/",name_code,"information_all_model_fit_params.csv"))

######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** TYPICAL STRAIN ANALYSIS
######****************************######################################################################

param_typical$inocl <- as.numeric(param_typical$inocl) # 7 strains! 

######****** Predicting **********######################################################################
ggplot(param_typical, aes(y=inocl,x = t_m_h_flow, group = strain_name, colour=factor(rep))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_x_continuous("Time to peak") + scale_y_continuous("Inoculum size (10^x)") + 
  scale_color_discrete("Experiment") 
ggsave(paste0("plots/",name_code,"time_to_peak_all.pdf"), width = 16, height =10 )

### PLOT variation in other indicators
#ggplot(param_typical, aes(y=inocl,x = v_m_h_flow, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = lag, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = exp_gr, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = auc, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)

######****** Linear model **********######################################################################
# Do for each rep, for each strain_name
reps <- unique(param_typical$rep)
strains <- unique(param_typical$strain_name)

ggplot(param_typical, aes(x=t_m_h_flow,y = log10(scalar), group = strain_name, colour=factor(drytime))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","24hr drying","168hr drying")) 
ggsave(paste0("plots/",name_code,"time_to_peak_all_as_linear_model.pdf"), width = 16, height =10 )

#### Fit linear model 
# Remove 10^2 and 10^6: not for reduction analysis
w26<-c(which(param_typical$inocl == 2), which(param_typical$inocl == 6))
if(length(w26) > 0){param_typical <- param_typical[-w26,]}

reductions_fit <- fit_line_model(reps, strains, param_typical, "t_m_h_flow","Time to max heat flow")

# in reductions. Meas column key: 
# meas = 1 = log reduction 
# meas = 2 = percentage reduction 
# meas = 3 = inoculum
# meas = 4 = predicted inoculum

ggplot(subset(reductions_fit$reductions, meas == 1), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_wrap(~strain_name, nrow = 2) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = pmax(0,mean - sd)),position = "dodge")
ggsave(paste0("plots/",name_code,"log_reductions_byrep.pdf"), width = 15, height = 10)

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
ggsave(paste0("plots/",name_code,"log_reductions_mean_by_drytime.pdf"), width = 20, height = 10)

ggplot(subset(reductions, meas == 2), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Percentage reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave(paste0("plots/",name_code,"perc_reductions.pdf"))


# By inoculum 
ggplot(reductions_fit$av_for_inoculum, aes(x = name, y = mean_inoc)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  facet_wrap(~strain_name, scales = "free") + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 3, lty = "dashed") + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Log reduction by inoculum")
ggsave(paste0("plots/",name_code,"log_reduction_by_inoc.pdf"))

ggplot(reductions_fit$av_for_inoculum, aes(x = name, y = mean_inoc)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  facet_wrap(~strain_name) + 
  geom_hline(yintercept = 1, lty = "dashed") + 
  geom_hline(yintercept = 3, lty = "dashed") + 
  scale_x_discrete("Inoculum") + 
  scale_y_continuous("Log reduction by inoculum")
ggsave(paste0("plots/",name_code,"log_reduction_by_inoc_fixed_scales.pdf"))


######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** NON-TYPICAL STRAIN ANALYSIS
######****************************######################################################################
strains_typical = unique(param_typical$strain_name) # PERFECT strains

# STRAINS GOOD BUT NOT FILTERD FOR HEIGHT
strains_typical_but_height = setdiff(unique(param_typical_notheight$strain_name),strains_typical) 

param_odd <- param %>%
  filter(!strain_name %in% strains_typical) %>%
  filter(!strain_name %in% strains_typical_but_height)

param_odd_height <- param %>%
  filter(strain_name %in% strains_typical_but_height)

length(unique(param_odd$strain_name))
length(unique(param$strain_name))
length(unique(param_odd_height$strain_name))

ggplot(param_odd_height, aes(x=t_m_h_flow, y=inocl)) + geom_point(aes(col=rep)) + facet_grid(drytime ~ strain_name)
ggplot(param_odd_height, aes(x=v_m_h_flow, y=inocl)) + geom_point(aes(col=rep)) + facet_grid(drytime ~ strain_name)

reps <- unique(param_odd_height$rep)
strains <- unique(param_odd_height$strain_name)

#### Fit linear model 
reductions <- fit_line_model(reps, strains, param_odd_height, plot = 1)

ggplot(subset(reductions, meas == 1), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave("plots/nont_odd_height_log_reductions.pdf")

ggplot(subset(reductions, meas == 2), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Percentage reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave("plots/nont_odd_height_perc_reductions.pdf")

