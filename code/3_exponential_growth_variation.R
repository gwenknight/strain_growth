#### Exponential growth variation 
## How to choose the threshold for exponential growth? 

library(tidyverse) # library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(patchwork) # for plot combinations
library(here)
theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

#####*************************** READ IN DATA *******************###############
ddm_orig <- read.csv("output/cut_all_time_series_fit_params.csv")[,-1]
ddm <- ddm_orig %>% filter(source == "Macotra") # only look at variaton in Macotra strains

param_orig <- read.csv("output/cut_all_model_fit_params.csv")[,-1]
param <- param_orig %>% filter(strain_name %in% unique(ddm$strain))

#################**************** CHECK EXPONENTIAL GROWTH *******************###############
dir.create(file.path(here(), "plots/exp_growth/"),showWarnings = FALSE) # Create a file for this output

# Look at the distribution of exponential growth for typical strains
for(i in unique(param$strain_name)){
  pa <- param %>% filter(strain_name == i)
  g1 <- ggplot(pa, aes(x=rep, y = cut_exp, group = rep)) + geom_boxplot() + ggtitle(paste0(i," all"))+ scale_x_continuous("Replicate") + scale_y_continuous("Exponential growth") 
  
  g2 <- ggplot(pa, aes(x=rep, y = cut_exp, group = interaction(rep,drytime))) + geom_boxplot(aes(col = factor(drytime))) + ggtitle(paste0("Split by drying time")) + 
    scale_color_discrete("Dry time") + scale_x_continuous("Replicate") + scale_y_continuous("Exponential growth") + 
    theme(legend.position="bottom")
  
  g1 + g2 + plot_layout(widths = c(1, 2)) 
  ggsave(paste0("plots/exp_growth/",i,".pdf"), width = 15)
}




#################**************** How far from the mean is the exponential growth? *******************###############

param_exp_gr_lab <- param %>%
  group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
  mutate(mean_peak_exp_gr = mean(cut_exp)) %>%
  ungroup() %>%
  mutate(outside05 = ifelse((cut_exp > (mean_peak_exp_gr + 0.05*mean_peak_exp_gr)) | (cut_exp < (mean_peak_exp_gr - 0.05*mean_peak_exp_gr)), 1, 0),
         outside10 = ifelse(cut_exp > (mean_peak_exp_gr + 0.1*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.1*mean_peak_exp_gr), 1, 0),
         outside15 = ifelse(cut_exp > (mean_peak_exp_gr + 0.15*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.15*mean_peak_exp_gr), 1, 0),
         outside20 = ifelse(cut_exp > (mean_peak_exp_gr + 0.2*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.2*mean_peak_exp_gr), 1, 0),
         outside22 = ifelse(cut_exp > (mean_peak_exp_gr + 0.22*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.22*mean_peak_exp_gr), 1, 0),
         outside24 = ifelse(cut_exp > (mean_peak_exp_gr + 0.24*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.24*mean_peak_exp_gr), 1, 0),
         outside26 = ifelse(cut_exp > (mean_peak_exp_gr + 0.26*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.26*mean_peak_exp_gr), 1, 0),
         outside28 = ifelse(cut_exp > (mean_peak_exp_gr + 0.28*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.28*mean_peak_exp_gr), 1, 0),
         outside30 = ifelse(cut_exp > (mean_peak_exp_gr + 0.3*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.3*mean_peak_exp_gr), 1, 0),
         outside31 = ifelse(cut_exp > (mean_peak_exp_gr + 0.31*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.31*mean_peak_exp_gr), 1, 0),
         outside32 = ifelse(cut_exp > (mean_peak_exp_gr + 0.32*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.32*mean_peak_exp_gr), 1, 0),
         outside33 = ifelse(cut_exp > (mean_peak_exp_gr + 0.33*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.33*mean_peak_exp_gr), 1, 0),
         outside34 = ifelse(cut_exp > (mean_peak_exp_gr + 0.34*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.34*mean_peak_exp_gr), 1, 0),
         outside35 = ifelse(cut_exp > (mean_peak_exp_gr + 0.35*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.35*mean_peak_exp_gr), 1, 0),
         outside36 = ifelse(cut_exp > (mean_peak_exp_gr + 0.36*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.36*mean_peak_exp_gr), 1, 0),
         outside37 = ifelse(cut_exp > (mean_peak_exp_gr + 0.37*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.37*mean_peak_exp_gr), 1, 0),
         outside38 = ifelse(cut_exp > (mean_peak_exp_gr + 0.38*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.38*mean_peak_exp_gr), 1, 0),
         outside39 = ifelse(cut_exp > (mean_peak_exp_gr + 0.39*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.39*mean_peak_exp_gr), 1, 0),
         outside40 = ifelse(cut_exp > (mean_peak_exp_gr + 0.4*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.4*mean_peak_exp_gr), 1, 0),
         outside50 = ifelse(cut_exp > (mean_peak_exp_gr + 0.5*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.5*mean_peak_exp_gr), 1, 0)) %>% 
  pivot_longer(col = c(outside05:outside50)) %>%
  mutate(thresh = as.numeric(gsub( "outside", "", name))) %>% 
  ungroup()

param_exp_gr_lab[,c("strain_name","rep","cut_exp","mean_peak_exp_gr","name","value","thresh","drytime","inocl")] # to have a look 


#### Group by strains & reps and look at how many datasets out 
# have a look at top 30 rows to see the structure. value = 1 says that the value (cut_exp) is > than the +/- percentage value of the mean
param_exp_gr_lab %>% dplyr::select(strain_name, rep, inocl, drytime, cut_exp, mean_peak_exp_gr, thresh, value) %>% print(n=30)

pp <- param_exp_gr_lab %>% 
  ungroup() %>% 
  group_by(strain_name, rep, thresh) %>% # For each strain and replicate, how many datasets
  dplyr::mutate(total_outside = sum(value), total = n()) %>% # are outside this threshold and how many datasets in total (some have < 6)
  dplyr::select(strain_name, rep, inocl, drytime, total_outside, thresh, total, value) %>%
  ungroup() %>% 
  group_by(thresh, total_outside)

pp_all <- pp %>% group_by(thresh, total_outside) %>% # For each threshold and number outside the range, how many 
  summarise(ns = n_distinct(strain_name, rep)) %>% # distinct strain and replicates are there
  group_by(thresh) %>%
  mutate(tot = sum(ns), perc = 100*ns / tot) # and what percentage of the total is this? 

## Want to look at number of datasets excluded by drytime and rep
pp_strain_dryt <- pp %>%
  group_by(strain_name, rep, thresh, drytime) %>% # At each drytime for this strain and rep
  summarise(total_outside_inrep = sum(value), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep) %>% 
  group_by(thresh, strain_name, rep) %>% # in a replicate how many outside? 
  mutate(rep_dt_remove = ifelse(perc_outside > 34, 1, 0)) %>% # see how many outside: flag to remove a rep if more than a third are outside 
  group_by(thresh, strain_name) %>% # in a strain 
  summarise(rep_removes = sum(rep_dt_remove)) %>% # sum how many reps are outside? 
  group_by(thresh, rep_removes) %>% # by a threshold and replicate removal number
  summarise(ns = n_distinct(strain_name)) %>% # how many strains are in each? 
  ungroup() %>% 
  group_by(thresh) %>% # within at threshold how many strains are there? and what percentage are in each? 
  mutate(tot = sum(ns), perc = 100*ns / tot) %>% 
  ungroup()

## Want to look at number of strains excluded 
pp_strain_names <- pp %>%
  group_by(strain_name, rep, thresh, drytime) %>% # in a rep and drytime 
  summarise(total_outside_inrep = sum(value), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep) %>% # how many outside? 
  group_by(thresh, strain_name, rep) %>% # in a replicate how many outside
  mutate(rep_dt_remove = ifelse(perc_outside > 34, 1, 0)) %>% # see how many outside: flag to remove a rep if more than a third are outside 
  group_by(thresh, strain_name, rep) %>% # in a replicate
  summarise(s = sum(rep_dt_remove), rep_remove = ifelse(s>0,1,0)) %>% # if a replicate to remove
  group_by(thresh, strain_name) %>% # within a strain 
  summarise(rep_removes = sum(rep_remove)) %>% # calculate how many replicates to remove
  group_by(thresh, rep_removes) 

pp_strain <- pp_strain_names%>% # building on the above
  summarise(ns = n_distinct(strain_name)) %>% # how many strains within each rep_removal category? e.g. how many strains have 2 reps removed? 
  ungroup() %>% 
  group_by(thresh) %>% # within a threshold: what percentage of strains? 
  mutate(tot = sum(ns), perc = 100*ns / tot) %>%
  ungroup()


ggplot(pp_all, aes(x=thresh, y = ns)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Number of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded") +
  ggsave("plots/exp_growth/exp_explore_n_datasets.pdf")

g1 <- ggplot(pp_all, aes(x=thresh, y = perc)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("% of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded") + 
  geom_hline(yintercept = 10, lty = "dashed") + 
  geom_hline(yintercept = 5, lty = "dashed") +
  geom_hline(yintercept = 30, lty = "dashed") + 
  geom_vline(xintercept = 36, lty = "dashed")
ggsave("plots/exp_growth/exp_explore_perc_datasets.pdf")

# By strain & drytime
ggplot(pp_strain_dryt, aes(x=thresh, y = ns)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(rep_removes))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Number of strains") + 
  scale_fill_discrete("Number of\nrep x drytimes\nexcluded")
ggsave("plots/exp_growth/exp_explore_n_datasets_strains_dryt.pdf")

g2 <- ggplot(pp_strain_dryt, aes(x=thresh, y = perc)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(rep_removes))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Percentage of reps x drytimes") + 
  geom_hline(yintercept = 10, lty = "dashed") + 
  geom_hline(yintercept = 5, lty = "dashed") +
  geom_hline(yintercept = 30, lty = "dashed") +
  scale_fill_discrete("Number of\nrep x drytimes\nexcluded") + 
  geom_vline(xintercept = 36, lty = "dashed")

# By strain
ggplot(pp_strain, aes(x=thresh, y = ns)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(rep_removes))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Number of strains") + 
  scale_fill_discrete("Number of\nreps\nexcluded")
ggsave("plots/exp_growth/exp_explore_n_datasets_strains.pdf")

g3 <- ggplot(pp_strain, aes(x=thresh, y = perc)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(rep_removes))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Percentage of strains") + 
  scale_fill_discrete("Number of\nreps\nexcluded") + 
  geom_hline(yintercept = 10, lty = "dashed") +
  geom_hline(yintercept = 5, lty = "dashed") +
  geom_hline(yintercept = 30, lty = "dashed") + 
  geom_vline(xintercept = 36, lty = "dashed")

(g1 + g2) / g3 + plot_annotation(tag_levels = 'A')
ggsave("plots/exp_growth/exp_explore_combine.pdf", width = 25)

#################****************  PLOT THE INDIVIDUAL POINT VARIATION *******************###############

cutoff <- 0.34 ## Determined as the cutoff that removes the top 10% of variable strains

pp_strain_names <- param %>%
  group_by(strain_name, rep) %>% # Over the replicate within a strain
  mutate(mean_peak_exp_gr = mean(cut_exp), # what is the mean? 
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + cutoff*mean_peak_exp_gr, # what is the upper limit?
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - cutoff*mean_peak_exp_gr, # what is the lower limit?
         outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0)) %>% # is it outside? 
  ungroup() %>% 
  group_by(strain_name, rep,drytime) %>% # in a strain, rep and drytime
  mutate(total_outside_inrep = sum(outside), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep) %>% # how many / percentage outside? 
  ungroup() %>% 
  group_by(strain_name, rep) %>% # in a replicate
  mutate(rep_dt_remove = ifelse(perc_outside > 34, 1, ifelse(total_inrep < 2,1,0)), total_rep_rem = sum(rep_dt_remove)) %>% # remove if > 34 = third outside
  ungroup() %>%
  mutate(keep_rep =ifelse((total_rep_rem == 0),rep,0)) %>% # if none to remove then put in rep name
  group_by(strain_name) %>%
  mutate(remove_strain = ifelse((n_distinct(keep_rep) - any(keep_rep == 0)) >=2,0,1)) # n_distinct counts 0 so need to remove

theme_set(theme_bw(base_size=6))
ggplot(pp_strain_names, aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(remove_strain)),lty = "dashed") + 
  ggtitle(paste0("All strains. Red = outside of limits (",100*cutoff,"%)"))
ggsave(paste0("plots/exp_growth/cutoff_from_mean_exponential_growth.pdf"),width = 30, height = 30)


ggplot(pp_strain_names %>% filter(strain_name == "11288"), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime)), size = 3) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth", lim = c(0,0.03)) + 
  facet_wrap(strain_name ~ rep, ncol = 30) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(remove_strain)),lty = "dashed") + 
  ggtitle(paste0("11288. Red = outside of limits (",100*cutoff,"%)")) + theme_bw(base_size = 20)
ggsave(paste0("plots/exp_growth/eg_cutoff_from_mean_exponential_growth.pdf"))



ggplot(pp_strain_names %>% filter(remove_strain == 1), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, nrow = 5) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(remove_strain)),lty = "dashed") + 
  ggtitle(paste0("All strains. Red = outside of limits (",100*cutoff,"%)"))
ggsave(paste0("plots/exp_growth/cutoff_from_mean_exponential_growth_outside.pdf"),width = 10, height = 10)

ggplot(pp_strain_names %>% filter(remove_strain == 0), aes(x=inocl, y = cut_exp)) + geom_point(aes(colour = factor(outside),pch = factor(drytime))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  scale_shape_discrete("Drytime") + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5)) + 
  scale_y_continuous("Exponential growth") + 
  facet_wrap(strain_name ~ rep, nrow = 5) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10, col = factor(remove_strain))) +
  geom_hline(aes(yintercept = mean_peak_exp_gr, col = factor(remove_strain)),lty = "dashed") + 
  ggtitle(paste0("All strains. Red = outside of limits (",100*cutoff,"%)"))
ggsave(paste0("plots/exp_growth/cutoff_from_mean_exponential_growth_inside.pdf"),width = 10, height = 10)


# Look at example
# This one is to be removed
pp_strain_names %>% filter(strain_name == "11070") %>% 
  dplyr::select(rep, drytime, inocl, mean_peak_exp_gr_m10, cut_exp,  mean_peak_exp_gr_p10, outside, rep_dt_remove, total_outside_inrep, total_inrep, perc_outside, total_rep_rem, remove_strain)

### Look at normalised start
ddm <-ddm %>% group_by(strain,rep, inoc, drytime) %>% dplyr::mutate(initial = first(value_J), value_J_norm = value_J - initial)
ggplot(ddm, aes(x = Time, y = value_J_norm, group= interaction(rep, inoc, drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain)
ggsave("plots/exp_growth/normalised_start.pdf", width = 20, height = 20)


######### Diagnostic / exploring plots - not necessary 
# # Change replicates from within experiment given names to replicate 1 to 3 for all 
# po <- param %>% group_by(strain_name) %>% dplyr::mutate(maxx = max(rep), minn = min(rep), 
#                                                         ones = ifelse(rep == minn, 1, 0), threes = ifelse(rep == maxx,1,0), twos = ifelse(ones == 0, ifelse(threes == 0,1,0),0)) %>%
#   mutate(rep_st = case_when((ones == 1) ~ 1,
#                             (threes == 1)  ~ 3,
#                             (twos == 1) ~ 2)) # tries pmax etc didn't work # Labels reps as 1 2 3
# 
# 
# ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) +
#   scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
#   scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/exp_growth_across_replicates.pdf")
# 
# ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
#   scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
#   scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/exp_growth_across_inocl_bp.pdf")
# 
# ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
#   scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
#   scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/exp_growth_across_inocl_bpgrp.pdf")
# 
# ggplot(po, aes(x=inocl, y = cut_exp,aes(group = drytime))) + geom_point() +  
#   geom_smooth(method = "loess") +  #, #, formula = y ~ a * x + b,method.args = list(start = list(a = 0.1, b = 0.1))) + 
#   facet_wrap(~drytime) + 
#   scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
#   scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/exp_growth_across_inocl_pointsline.pdf")
# 
# ggplot(po, aes(x=rep_st, y = t_m_h_flow, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
#   scale_y_continuous("Time to peak") + 
#   scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/time_to_peak_across_replicates.pdf")
# 
# ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
#   scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
#   scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/across_replicates_zoom.pdf")
# 
# ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
#   geom_jitter(width = 0.1,aes(col = factor(drytime))) + 
#   facet_wrap(~strain_name, scales = "free_x") + 
#   scale_y_continuous("Exponential growth") + 
#   scale_x_discrete("Replicate/Strain") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
# ggsave("plots/exp_growth/across_replicates_zoom_strain_rep.pdf",width = 30, height = 30)
# 
# ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
#   geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
#   facet_wrap(~strain_name, scales = "free_x") + 
#   scale_y_continuous("Exponential growth") + 
#   scale_x_discrete("Replicate/Strain") +   
#   scale_shape_discrete("Dry time") + 
#   scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
#   theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
# ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inocl.pdf",width = 30, height = 30)
# 
# ggplot(po, aes(x=interaction(inocl, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
#   geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
#   facet_wrap(~strain_name, scales = "free_x") + 
#   scale_y_continuous("Exponential growth") + 
#   scale_x_discrete("Inoculum/Strain") +   
#   scale_shape_discrete("Dry time") + 
#   scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
#   theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
# ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inoclgrp.pdf",width = 30, height = 30)
# 
# 
# ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
#   scale_y_continuous("exponential growth") + 
#   scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/across_replicates_all.pdf")
# 
# ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
#   scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
#   scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
#   theme(legend.position="bottom")
# ggsave("plots/exp_growth/group_replicates.pdf")
# 
# ggplot(po %>% filter(drytime == 0), aes(x=inocl, y = cut_exp, group = interaction(rep,strain))) + 
#   geom_line(aes(col = strain))