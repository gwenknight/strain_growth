# EXPLORE exponential growth variation 
library(tidyverse) # good library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
theme_set(theme_bw(base_size=6)) # theme setting for plots: black and white (bw) and font size (24)

##### READ IN DATA
param <- read.csv("output/param_labelled.csv")
param_clean <- param %>% filter(removed_rep == 0, removed_dataset == 0)

###################************** CHECK EXPONENTIAL GROWTH OK ##############
param_clean$use_exp <- param_clean$exp_gr
w <- which(!is.na(param_clean$cut_exp))
param_clean[w,"use_exp"] <- param_clean[w,"cut_exp"]

pdf("plots/exp_growth/impact_of_cut.pdf")
plot(param_clean[w,"exp_gr"],param_clean[w,"cut_exp"], xlab = "Pre cut exponential growth", ylab = "After cut exponential growth")
lines(seq(0,0.04,0.001),seq(0,0.04,0.001))
dev.off()


##### Check exp growth OK across inocula of the same strain at set drying times
param_exp_gr_lab <- param_clean %>%
  filter(odd_strains == 0) %>%
  group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
  mutate(mean_peak_exp_gr = mean(use_exp)) %>%
  ungroup() %>%
  mutate(outside05 = ifelse((exp_gr > (mean_peak_exp_gr + 0.05*mean_peak_exp_gr)) | (exp_gr < (mean_peak_exp_gr - 0.05*mean_peak_exp_gr)), 1, 0),
         outside10 = ifelse(exp_gr > (mean_peak_exp_gr + 0.1*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.1*mean_peak_exp_gr), 1, 0),
         outside15 = ifelse(exp_gr > (mean_peak_exp_gr + 0.15*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.15*mean_peak_exp_gr), 1, 0),
         outside20 = ifelse(exp_gr > (mean_peak_exp_gr + 0.2*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.2*mean_peak_exp_gr), 1, 0),
         outside25 = ifelse(exp_gr > (mean_peak_exp_gr + 0.25*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.25*mean_peak_exp_gr), 1, 0),
         outside30 = ifelse(exp_gr > (mean_peak_exp_gr + 0.3*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.3*mean_peak_exp_gr), 1, 0),
         outside35 = ifelse(exp_gr > (mean_peak_exp_gr + 0.35*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.35*mean_peak_exp_gr), 1, 0),
         outside40 = ifelse(exp_gr > (mean_peak_exp_gr + 0.4*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.4*mean_peak_exp_gr), 1, 0),
         outside50 = ifelse(exp_gr > (mean_peak_exp_gr + 0.5*mean_peak_exp_gr) | exp_gr < (mean_peak_exp_gr - 0.5*mean_peak_exp_gr), 1, 0)) 

total = cbind(c(seq(5,35,5),40,50), 100*colSums(param_exp_gr_lab[,c("outside05","outside10","outside15","outside20",
                                                           "outside25","outside30","outside35","outside40","outside50")]) / dim(param)[1])

total <- as.data.frame(total)
colnames(total) <- c("range","perc_outside")


ggplot(total, aes(x= range, y = perc_outside)) + geom_point() + 
  geom_line() + scale_x_continuous("Cutoff (within +/- x% of mean)") + 
  scale_y_continuous("Percentage outside this range", lim = c(0,90))
ggsave("plots/exp_growth/exp_explore_perc_outside.pdf")

#### Group by strains and look at how many datasets out 
pp <- param_exp_gr_lab %>% 
  group_by(strain_name,rep) %>% 
  dplyr::summarise(sum_outside_s_r05 = sum(outside05),
                   sum_outside_s_r10 = sum(outside10),
                   sum_outside_s_r15 = sum(outside15),
                   sum_outside_s_r20 = sum(outside20),
                   sum_outside_s_r25 = sum(outside25),
                   sum_outside_s_r30 = sum(outside30),
                   sum_outside_s_r35 = sum(outside35),
                   sum_outside_s_r40 = sum(outside40),
                   sum_outside_s_r50 = sum(outside50),.groups = "drop") # number of odd datasets in this strain and replicate

reps_split <- matrix(0,7,9)
reps_split[1+min(pp$sum_outside_s_r05):max(pp$sum_outside_s_r05),1] <- table(pp$sum_outside_s_r05)
reps_split[1+min(pp$sum_outside_s_r10):max(pp$sum_outside_s_r10),2] <- table(pp$sum_outside_s_r10)
reps_split[1+min(pp$sum_outside_s_r15):max(pp$sum_outside_s_r15),3] <- table(pp$sum_outside_s_r15)
reps_split[1+min(pp$sum_outside_s_r20):max(pp$sum_outside_s_r20),4] <- table(pp$sum_outside_s_r20)
reps_split[1+min(pp$sum_outside_s_r25):max(pp$sum_outside_s_r25),5] <- table(pp$sum_outside_s_r25)
reps_split[1+min(pp$sum_outside_s_r30):max(pp$sum_outside_s_r30),6] <- table(pp$sum_outside_s_r30)
reps_split[1+min(pp$sum_outside_s_r35):max(pp$sum_outside_s_r35),7] <- table(pp$sum_outside_s_r35)
reps_split[1+min(pp$sum_outside_s_r40):max(pp$sum_outside_s_r40),8] <- table(pp$sum_outside_s_r40)
reps_split[1+min(pp$sum_outside_s_r50):max(pp$sum_outside_s_r50),9] <- table(pp$sum_outside_s_r50)
reps_split <- as.data.frame(reps_split)
colnames(reps_split) <- c(seq(5,35,5),40,50)
reps_split$no_datasets <- seq(0,6,1)
mreps_split <- reshape2::melt(reps_split, id.vars = "no_datasets")

## With a cutoff of 5 %, 

ggplot(mreps_split, aes(x=variable, y = value)) + 
  geom_bar(position = "stack",stat = "identity", aes(fill = factor(no_datasets))) + 
  scale_x_discrete("Percentage of exponential growth excluded") + 
  scale_y_continuous("Number of strain and replicate combinations") + 
  scale_fill_discrete("Number of\ndatasets\nexcluded")
ggsave("plots/exp_growth/exp_explore_n_datasets.pdf")



