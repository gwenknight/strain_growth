# EXPLORE exponential growth variation 
library(tidyverse) # good library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(patchwork) # for plot combinations
theme_set(theme_bw(base_size=6)) # theme setting for plots: black and white (bw) and font size (24)

##### READ IN DATA
param <- read.csv("output/param_labelled.csv")
param_clean <- param #%>% filter(removed_rep == 0, removed_dataset == 0) ## NO longer remove as cutting at peak time 

###################************** CHECK EXPONENTIAL GROWTH OK ##############
## Imipact of original vs cut exp
pdf("plots/exp_growth/impact_of_cut.pdf")
plot(param_clean[w,"exp_gr"],param_clean[w,"cut_exp"], xlab = "Pre cut exponential growth", ylab = "After cut exponential growth")
lines(seq(0,0.04,0.001),seq(0,0.04,0.001))
text(cut_exp ~exp_gr, labels=param_clean$strain,data=param_clean, cex=0.9, font=2)
dev.off()


##### Check exp growth OK across inocula of the same strain at set drying times
param_exp_gr_lab <- param_clean %>%
  group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
  mutate(mean_peak_exp_gr = mean(cut_exp)) %>%
  ungroup() %>%
  mutate(outside05 = ifelse((cut_exp > (mean_peak_exp_gr + 0.05*mean_peak_exp_gr)) | (cut_exp < (mean_peak_exp_gr - 0.05*mean_peak_exp_gr)), 1, 0),
         outside10 = ifelse(cut_exp > (mean_peak_exp_gr + 0.1*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.1*mean_peak_exp_gr), 1, 0),
         outside15 = ifelse(cut_exp > (mean_peak_exp_gr + 0.15*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.15*mean_peak_exp_gr), 1, 0),
         outside20 = ifelse(cut_exp > (mean_peak_exp_gr + 0.2*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.2*mean_peak_exp_gr), 1, 0),
         outside25 = ifelse(cut_exp > (mean_peak_exp_gr + 0.25*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.25*mean_peak_exp_gr), 1, 0),
         outside30 = ifelse(cut_exp > (mean_peak_exp_gr + 0.3*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.3*mean_peak_exp_gr), 1, 0),
         outside35 = ifelse(cut_exp > (mean_peak_exp_gr + 0.35*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.35*mean_peak_exp_gr), 1, 0),
         outside40 = ifelse(cut_exp > (mean_peak_exp_gr + 0.4*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.4*mean_peak_exp_gr), 1, 0),
         outside50 = ifelse(cut_exp > (mean_peak_exp_gr + 0.5*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.5*mean_peak_exp_gr), 1, 0)) %>% 
  pivot_longer(col = c(outside05:outside50)) %>%
  mutate(thresh = as.numeric(gsub( "outside", "", name))) %>% 
  ungroup()

param_exp_gr_lab[,c("strain_name","rep","cut_exp","mean_peak_exp_gr","name","value","thresh")] # to have a look 

# Group together wtih percentages: sum over the respective columns (1 if outside so sum is total number)

total_split <- param_exp_gr_lab %>% group_by(thresh, odd_strains) %>% summarise(total = sum(value), perc = 100*total/dim(param)[1])
total <- param_exp_gr_lab %>% group_by(thresh) %>% summarise(total = sum(value), perc = 100*total/dim(param)[1])
total$odd_strains <- 2
tot <- full_join(total, total_split)

# OLD pre pivot longer
#total = cbind(c(seq(5,35,5),40,50), 100*colSums(param_exp_gr_lab[,c("outside05","outside10","outside15","outside20",
#                                                           "outside25","outside30","outside35","outside40","outside50")]) / dim(param)[1])


ggplot(tot, aes(x= thresh, y = perc, group = odd_strains)) + geom_point() + 
  geom_line(aes(col = factor(odd_strains))) + scale_x_continuous("Cutoff (within +/- x% of mean)") + 
  scale_y_continuous("Percentage of datasets outside this range", lim = c(0,100)) + 
  geom_hline(yintercept = 30, lty = "dashed") +
  scale_color_discrete("Odd strains?", breaks = c(0,1,2), labels = c("None","Only","All strains"))
ggsave("plots/exp_growth/exp_explore_perc_outside.pdf")

#### Group by strains & reps and look at how many datasets out 
pp <- param_exp_gr_lab %>% 
  group_by(strain_name,rep, thresh) %>% 
  dplyr::summarise(total_outside = sum(value), total = n()) %>%
  ungroup() %>%
  group_by(thresh, total_outside) %>% 
  summarise(ns = n_distinct(strain_name, rep)) %>%
  group_by(thresh) %>%
  mutate(tot = sum(ns), 
            perc = 100*ns / tot)
  
  

ggplot(pp, aes(x=thresh, y = ns)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Number of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded")
ggsave("plots/exp_growth/exp_explore_n_datasets.pdf")

ggplot(pp, aes(x=thresh, y = perc)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("% of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded") + 
  geom_hline(yintercept = 30, lty = "dashed") 
ggsave("plots/exp_growth/exp_explore_perc_datasets.pdf")



##### strains only? not sure what the colours are here... 
pp <- param_exp_gr_lab %>% 
  group_by(strain_name, thresh) %>% 
  dplyr::summarise(total_outside = sum(value), total = n()) %>%
  ungroup() %>%
  group_by(thresh, total_outside) %>% 
  summarise(ns = n_distinct(strain_name)) %>%
  group_by(thresh) %>%
  mutate(tot = sum(ns), 
         perc = 100*ns / tot)


ggplot(pp, aes(x=thresh, y = ns)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("Number of strains") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded")
#ggsave("plots/exp_growth/exp_explore_n_datasets.pdf")

ggplot(pp, aes(x=thresh, y = perc)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(total_outside))) + 
  scale_x_continuous("Exclude if exponential growth is +/- \nthis percentage away from mean") + 
  scale_y_continuous("% of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded") + 
  geom_hline(yintercept = 30, lty = "dashed") 
#ggsave("plots/exp_growth/exp_explore_perc_datasets.pdf")

##### Exponential agains height
p1 <- ggplot(param_exp_gr_lab, aes(x=cut_exp, y = v_m_h_flow)) + 
  geom_point(aes(col = factor(drytime))) + 
  geom_smooth(method=lm) + 
  geom_abline(intercept = 0.03, slope = 1) + 
  scale_y_continuous("Peak heigh value") + 
  scale_x_continuous("Exponential growth") + 
  scale_color_discrete("Drytime") +
  ggtitle("Linear model fit all to data vs. straight correlation (black)")

p2 <- ggplot(param_exp_gr_lab, aes(x=cut_exp, y = v_m_h_flow)) + 
  geom_point(aes(col = factor(drytime))) + 
  geom_smooth() + 
  geom_abline(intercept = 0.03, slope = 1) + 
  scale_y_continuous("Peak heigh value") + 
  scale_x_continuous("Exponential growth") + 
  scale_color_discrete("Drytime") + 
  ggtitle("Loess model fit to all data vs. straight correlation (black)")

p3 <- ggplot(param_exp_gr_lab, aes(x=cut_exp, y = v_m_h_flow,col = factor(drytime))) + 
  geom_point(aes()) + 
  geom_smooth(fullrange=TRUE) + 
  geom_abline(intercept = 0.03, slope = 1) + 
  scale_y_continuous("Peak heigh value") + 
  scale_x_continuous("Exponential growth") + 
  scale_color_discrete("Drytime") + 
  ggtitle("Loess model fit to data by drytime vs. straight correlation (black)")

p4 <- ggplot(param_exp_gr_lab, aes(x=cut_exp, y = v_m_h_flow,col = factor(odd_strains))) + 
  geom_point(aes()) + 
  geom_smooth(fullrange=TRUE) + 
  geom_abline(intercept = 0.03, slope = 1) + 
  scale_y_continuous("Peak heigh value") + 
  scale_x_continuous("Exponential growth") + 
  scale_color_manual("Odd strain?", values = c("green","pink")) + 
  ggtitle("Loess model fit to data by odd strain type vs. straight correlation (black)")

(p1 | p2) / ( p3 | p4 )
ggsave("plots/exp_growth/height_vs_exp.pdf")




