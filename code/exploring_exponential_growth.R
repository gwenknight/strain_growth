# EXPLORE exponential growth variation 
library(tidyverse) # library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(patchwork) # for plot combinations
theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)

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
         outside22 = ifelse(cut_exp > (mean_peak_exp_gr + 0.22*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.22*mean_peak_exp_gr), 1, 0),
         outside24 = ifelse(cut_exp > (mean_peak_exp_gr + 0.24*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.24*mean_peak_exp_gr), 1, 0),
         outside26 = ifelse(cut_exp > (mean_peak_exp_gr + 0.26*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.26*mean_peak_exp_gr), 1, 0),
         outside28 = ifelse(cut_exp > (mean_peak_exp_gr + 0.28*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.28*mean_peak_exp_gr), 1, 0),
         outside30 = ifelse(cut_exp > (mean_peak_exp_gr + 0.3*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.3*mean_peak_exp_gr), 1, 0),
         outside32 = ifelse(cut_exp > (mean_peak_exp_gr + 0.32*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.32*mean_peak_exp_gr), 1, 0),
         outside34 = ifelse(cut_exp > (mean_peak_exp_gr + 0.34*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.34*mean_peak_exp_gr), 1, 0),
         outside36 = ifelse(cut_exp > (mean_peak_exp_gr + 0.36*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.36*mean_peak_exp_gr), 1, 0),
         outside38 = ifelse(cut_exp > (mean_peak_exp_gr + 0.38*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.38*mean_peak_exp_gr), 1, 0),
         outside40 = ifelse(cut_exp > (mean_peak_exp_gr + 0.4*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.4*mean_peak_exp_gr), 1, 0),
         outside50 = ifelse(cut_exp > (mean_peak_exp_gr + 0.5*mean_peak_exp_gr) | cut_exp < (mean_peak_exp_gr - 0.5*mean_peak_exp_gr), 1, 0)) %>% 
  pivot_longer(col = c(outside05:outside50)) %>%
  mutate(thresh = as.numeric(gsub( "outside", "", name))) %>% 
  ungroup()

param_exp_gr_lab[,c("strain_name","rep","cut_exp","mean_peak_exp_gr","name","value","thresh","drytime","inoc")] # to have a look 

# Group together wtih percentages: sum over the respective columns (1 if outside so sum is total number)
n_odd = sum(param_exp_gr_lab)
total_split <- param_exp_gr_lab %>% group_by(thresh, odd_strains) %>% summarise(total = sum(value), number = n(), perc = 100*total/number)
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
  scale_color_discrete("Strain type", breaks = c(0,1,2), labels = c("Typical","Odd","All strains"))
ggsave("plots/exp_growth/exp_explore_perc_outside.pdf")

#### Group by strains & reps and look at how many datasets out 
pp <- param_exp_gr_lab %>% 
  group_by(strain_name, rep, thresh) %>% 
  mutate(total_outside = sum(value), total = n()) %>%
  select(strain_name, rep, inoc, drytime, total_outside, total, value) %>%
  ungroup() 

pp_all <- pp %>% group_by(thresh, total_outside) %>% 
  summarise(ns = n_distinct(strain_name, rep)) %>%
  group_by(thresh) %>%
  mutate(tot = sum(ns), perc = 100*ns / tot)

pp_strain_dryt <- pp %>%
  group_by(strain_name, rep, thresh, drytime) %>% 
  summarise(total_outside_inrep = sum(value), total_inrep = n()) %>%
  group_by(thresh, strain_name, rep) %>% 
  mutate(rep_dt_remove = ifelse(total_outside_inrep > 1, 1, 0)) %>%
  group_by(thresh, strain_name) %>%
  summarise(rep_removes = sum(rep_dt_remove)) %>% 
  group_by(thresh, rep_removes) %>%
  summarise(ns = n_distinct(strain_name)) %>%
  ungroup() %>% 
  group_by(thresh) %>%
  mutate(tot = sum(ns), perc = 100*ns / tot)

pp_strain <- pp %>%
  group_by(strain_name, rep, thresh, drytime) %>% 
  summarise(total_outside_inrep = sum(value), total_inrep = n()) %>%
  group_by(thresh, strain_name, rep) %>% 
  mutate(rep_dt_remove = ifelse(total_outside_inrep > 1, 1, 0)) %>%
  group_by(thresh, strain_name, rep) %>%
  mutate(s=sum(rep_dt_remove), rep_remove = ifelse(s > 0, 1, 0)) %>%
  group_by(thresh, strain_name) %>%
  summarise(rep_removes = sum(rep_remove)/2) %>% 
  group_by(thresh, rep_removes) %>%
  summarise(ns = n_distinct(strain_name)) %>%
  ungroup() %>% 
  group_by(thresh) %>%
  mutate(tot = sum(ns), perc = 100*ns / tot)


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
  geom_hline(yintercept = 30, lty = "dashed") 
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
  geom_hline(yintercept = 30, lty = "dashed") +
  scale_fill_discrete("Number of\nrep x drytimes\nexcluded")

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
  geom_hline(yintercept = 30, lty = "dashed") 

g1 + g2 + g3 
ggsave("plots/exp_growth/exp_explore_combine.pdf", width = 20)


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

################ 
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

####### Exponential against odd
param_exp_gr_lab$oddph <- param_exp_gr_lab$odd_peaks + param_exp_gr_lab$odd_width
ggplot(param_exp_gr_lab, aes(x= inocl, y = cut_exp, group = interaction(inocl,oddph))) + geom_boxplot(aes(col = factor(oddph))) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("New exponential") + 
  scale_color_discrete("Odd?", breaks = c(0,1,2), labels = c("No","Odd peak or width", "Odd peak & width"))
ggsave("plots/exp_growth/odd_strains_diff_exp.pdf")



