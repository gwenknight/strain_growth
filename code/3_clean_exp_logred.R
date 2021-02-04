#### Next steps
## Where do we go from the parameter set? 
# (1) Get the parameter sets from 2_analysis.R 

# (2) Clean datasets: 
### How many odd datasets in a drytime for a single rep? 
##### If < 50% of datasets in a rep at a drytime are odd then remove just this dataset
##### If > 50% of datasets are odd in one drytime then label the rep as odd
### How many odd replicates are there for a strain?
##### If more than one replicate is odd then strain is odd
##### If only one then remove the replicate from the strain
### > NOW just label for "removal" but as only use up to peak OK to keep in

# (3) Check exponential growth 
### Should be same across replicate

#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) 
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here:here())

#####*************************** READ IN DATA *******************###############
ddm <- read.csv("output/cut_all_time_series_fit_params.csv")[,-1]
param <- read.csv("output/cut_all_model_fit_params.csv")[,-1]

# Change from inoc to scalar
# scalar = 1 at 10^6, 1/10th of each next dilution
param$scalar <- 0

w<- which(param$inocl == "2")
param[w,"scalar"] = 100
w<- which(param$inocl == "3")
param[w,"scalar"] = 1000
w<- which(param$inocl == "4")
param[w,"scalar"] = 10000
w<- which(param$inocl == "5")
param[w,"scalar"] = 100000
w<- which(param$inocl == "6")
param[w,"scalar"] = 1000000

param$scalar2 <- (param$scalar)^2
strains <- unique(param$strain_name)
reps <- unique(param$rep)
drytimes <- unique(param$drytime)

po <- param %>% group_by(strain_name) %>% dplyr::mutate(maxx = max(rep), minn = min(rep), 
                                                       ones = ifelse(rep == minn, 1, 0), threes = ifelse(rep == maxx,1,0), twos = ifelse(ones == 0, ifelse(threes == 0,1,0),0)) %>%
  mutate(rep_st = case_when((ones == 1) ~ 1,
                            (threes == 1)  ~ 3,
                            (twos == 1) ~ 2)) # tries pmax etc didn't work # Labels reps as 1 2 3

write.csv(po, "output/param_labelled_repst.csv")







##### Diagnostic / exploring plots
ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) +
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_replicates.pdf")

ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_bp.pdf")

ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_bpgrp.pdf")

ggplot(po, aes(x=inocl, y = cut_exp,aes(group = drytime))) + geom_point() +  
  geom_smooth(method = "loess") +  #, #, formula = y ~ a * x + b,method.args = list(start = list(a = 0.1, b = 0.1))) + 
  facet_wrap(~drytime) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_pointsline.pdf")

ggplot(po, aes(x=rep_st, y = t_m_h_flow, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
  scale_y_continuous("Time to peak") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/time_to_peak_across_replicates.pdf")

ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
  scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/across_replicates_zoom.pdf")

ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Replicate/Strain") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep.pdf",width = 30, height = 30)

ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Replicate/Strain") +   
  scale_shape_discrete("Dry time") + 
  scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inocl.pdf",width = 30, height = 30)

ggplot(po, aes(x=interaction(inocl, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Inoculum/Strain") +   
  scale_shape_discrete("Dry time") + 
  scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inoclgrp.pdf",width = 30, height = 30)


ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
  scale_y_continuous("exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/across_replicates_all.pdf")

ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/group_replicates.pdf")

ggplot(po %>% filter(drytime == 0), aes(x=inocl, y = cut_exp, group = interaction(rep,strain))) + 
  geom_line(aes(col = strain))

#### Normalised start
ddm <-ddm %>% group_by(strain,rep, inoc, drytime) %>% dplyr::mutate(initial = first(value_J), value_J_norm = value_J - initial)
ggplot(ddm, aes(x = Time, y = value_J_norm, group= interaction(rep, inoc, drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain)
ggsave("plots/exp_growth/normalised_start.pdf", width = 20, height = 20)


#### Example 
# 11272
d2 <- ddm %>% filter(strain == "11277", rep == 1.1) %>% ungroup() %>% select(Time,value_J)
s <- summary(gcFitSpline(d2$Time, d2$value_J))

# ##### Check exp growth OK across inocula of the same strain at set drying times
# perc <- 0.2
# param_exp_gr_lab <- param %>%
#   group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
#   mutate(mean_peak_exp_gr = mean(exp_gr),
#          mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
#          mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
#          outside = ifelse(exp_gr < mean_peak_exp_gr_m10 | exp_gr > mean_peak_exp_gr_p10, 1, 0)) %>%
#   dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate
# 
# #### check for exp_gr - remove if only one is "wrong"
# param_typical <- param_exp_gr_lab %>%
#   filter(sum_outside_s_r < 2) # some have only one in the strain and replicate
# 
# length(unique(param_typical$strain_name))
# # For 1_2_6_ 9 strains still have typical behaviour
# 
# strains_typical = unique(param_typical$strain_name) # PERFECT strains
# 
# dim(param_typical_notexp_gr)[1] - dim(param_typical)[1] # 79 datasets removed because of exp_gr - one of those from the one strain that is removed ("RN6390B")
# 
# ggplot(param_exp_gr_lab[1:100,], aes(x=inocl, y = exp_gr)) + geom_point(aes(colour = factor(outside))) +
#   scale_color_manual("In the limits?", values = c("black","red")) +
#   facet_wrap(strain_name ~ rep + drytime) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_m10)) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_p10)) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr),lty = "dashed") +
#   ggtitle(paste0(100*perc," percent from mean exponential growth"))
# ggsave(paste0("plots/exp_growth/",name_code,100*perc,"perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)
# 
# param %>% filter(strain_name == "11277", rep == 1.1, drytime == 0)
