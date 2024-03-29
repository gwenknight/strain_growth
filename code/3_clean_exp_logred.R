#### Next steps
## Where do we go from the parameter set? 
# (1) Get the parameter sets from 2_analysis.R 

# (2) Clean datasets: 
## Remove those with exponential growth outside the set range (see 3_exponential_growth_variation.R)

# (3) Predict inoculum reduction 


#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) 
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
library(RColorBrewer)
library(here)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

# Read in code for linear model fitting
source("code/function_linear_model.R")

#####*************************** READ IN DATA *******************###############
ddm_orig <- read.csv("output/cut_all_time_series_fit_params.csv")[,-1]
ddm <- ddm_orig %>% filter(source == "Macotra")

param_orig <- read.csv("output/cut_all_model_fit_params.csv")[,-1]
param <- param_orig %>% filter(strain_name %in% unique(ddm$strain))

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

# Change replicates from within experiment given names to replicate 1 to 3 for all 
po <- param %>% group_by(strain_name) %>% dplyr::mutate(maxx = max(rep), minn = min(rep), 
                                                        ones = ifelse(rep == minn, 1, 0), threes = ifelse(rep == maxx,1,0), twos = ifelse(ones == 0, ifelse(threes == 0,1,0),0)) %>%
  mutate(rep_st = case_when((ones == 1) ~ 1,
                            (threes == 1)  ~ 3,
                            (twos == 1) ~ 2)) %>% # tries pmax etc didn't work # Labels reps as 1 2 3
  dplyr::select(-c(maxx,minn, ones, threes, twos))

write.csv(po, "output/param_labelled_repst.csv")

#####*************************** REMOVE THOSE WITH EXPONENTIAL GROWTH OUTSIDE OF RANGE *******************###############
cutoff <- 0.36 ### Determined by analysis in 3_exponential_growth_variation.R

## Iterative filter 
pp_strain_names <- po %>%
  group_by(strain_name, rep_st) %>% 
  dplyr::mutate(mean_peak_exp_gr = mean(cut_exp), # mean growth of first part of curve
                mean_peak_exp_gr_p10 = mean_peak_exp_gr + cutoff*mean_peak_exp_gr, # +/- range
                mean_peak_exp_gr_m10 = mean_peak_exp_gr - cutoff*mean_peak_exp_gr,
                outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0), # is the growth for this dataset outside the range?
                diff = ifelse(outside == 1, abs(cut_exp - mean_peak_exp_gr), 0), # what is the difference between the growth for this dataset vs the mean? 
                max = max(diff), # what is the maximum difference across this 
                remove_dataset = ifelse(diff == max, ifelse(max > 0,1,0), 0)) %>%
  ungroup() %>% 
  group_by(strain_name, rep_st, remove_dataset) %>% # group by remove_dataset now as want to consider those not removed differently
  dplyr::mutate(mean_peak_exp_gr2 = mean(cut_exp),
                mean_peak_exp_gr_p102 = mean_peak_exp_gr2 + cutoff*mean_peak_exp_gr2,
                mean_peak_exp_gr_m102 = mean_peak_exp_gr2 - cutoff*mean_peak_exp_gr2,
                outside2 = ifelse(cut_exp < mean_peak_exp_gr_m102 | cut_exp > mean_peak_exp_gr_p102, 1, 0),
                diff2 = ifelse(outside2 == 1, abs(cut_exp - mean_peak_exp_gr2), 0),
                max2 = max(diff2), 
                remove_dataset2 = ifelse(diff2 == max2, ifelse(max2 > 0,1,0), 0)) %>% 
  group_by(strain_name, rep_st, remove_dataset, remove_dataset2) %>%
  dplyr::mutate(mean_peak_exp_gr3 = mean(cut_exp),
                mean_peak_exp_gr_p103 = mean_peak_exp_gr3 + cutoff*mean_peak_exp_gr3,
                mean_peak_exp_gr_m103 = mean_peak_exp_gr3 - cutoff*mean_peak_exp_gr3,
                outside3 = ifelse(cut_exp < mean_peak_exp_gr_m103 | cut_exp > mean_peak_exp_gr_p103, 1, 0),
                diff3 = ifelse(outside3 == 1, abs(cut_exp - mean_peak_exp_gr3), 0),
                max3 = max(diff3), 
                remove_dataset3 = ifelse(diff3 == max3, ifelse(max3 > 0,1,0), 0))  %>% 
  group_by(strain_name, rep_st, remove_dataset, remove_dataset2, remove_dataset3) %>%
  dplyr::mutate(mean_peak_exp_gr4 = mean(cut_exp),
                mean_peak_exp_gr_p104 = mean_peak_exp_gr4 + cutoff*mean_peak_exp_gr4,
                mean_peak_exp_gr_m104 = mean_peak_exp_gr4 - cutoff*mean_peak_exp_gr4,
                outside4 = ifelse(cut_exp < mean_peak_exp_gr_m104 | cut_exp > mean_peak_exp_gr_p104, 1, 0),
                diff4 = ifelse(outside4 == 1, abs(cut_exp - mean_peak_exp_gr4), 0),
                max4 = max(diff4), 
                remove_dataset4 = ifelse(diff4 == max4, ifelse(max4 > 0,1,0), 0)) %>%
  ungroup() %>%
  dplyr::mutate(remove_dataset_exp_iter = remove_dataset + remove_dataset2 + remove_dataset3 + remove_dataset4) %>% # sum over all: remove if any are > 0
  group_by(strain_name, rep, drytime) %>%  # now group by drytime: matters if have few datasets pre / post drying 
  dplyr::mutate(total_outside_inrep = sum(remove_dataset_exp_iter), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep,# how many of the datasets are removed in each drytime of the rep? 
                rep_remove = ifelse(drytime == 0,ifelse(perc_outside > 34,1,0), ifelse(perc_outside == 100, 1, 0))) %>% # If drytime == 0 then need 3 left for linear model. If drytime == 168 then remove if all gone
  ungroup() %>% 
  group_by(strain_name, rep) %>% 
  dplyr::mutate(total_rep_rem = sum(rep_remove)) %>%
  ungroup() %>%
  dplyr::mutate(keep_rep = ifelse((total_rep_rem == 0),1,0)) 

# CHECK look at those reps removed
pp_strain_names %>% filter(keep_rep == 0) %>% dplyr::select(strain_name, rep, drytime, inocl, cut_exp, remove_dataset_exp_iter, total_outside_inrep, perc_outside, rep_remove, keep_rep, total_rep_rem) %>% print(n=Inf)
# CHECK look at those datasets removed
pp_strain_names %>% filter(remove_dataset_exp_iter == 1) %>% dplyr::select(strain_name, rep, drytime, inocl, cut_exp, remove_dataset_exp_iter, total_outside_inrep, perc_outside, rep_remove, keep_rep, total_rep_rem) %>% print(n=Inf)

#### Remove those with exponential values outside the above range
param_expok <- pp_strain_names %>%
  #filter(remove_dataset == 0, remove_dataset2 == 0, remove_dataset3 == 0, remove_dataset4 == 0) # remove those in the interative calculation of the mean: single datasets
  filter(remove_dataset_exp_iter == 0) %>% # same as above line: this is the sum of each of these indices
  filter(keep_rep == 1) # Keep those reps that have enough datasets
  

## Which are removed? 
dim(pp_strain_names) # 1739 datasets initially
length(which(pp_strain_names$remove_dataset_exp_iter == 1)) # 171 single datasets removed from strains as outside range
length(which(pp_strain_names$keep_rep == 0)) # 216 

length(which(pp_strain_names$remove_dataset == 1)) # 119 Each of these is remove one dataset and then recalculate the mean within a replicate
length(which(pp_strain_names$remove_dataset2 == 1)) # 43
length(which(pp_strain_names$remove_dataset3 == 1)) # 9
length(which(pp_strain_names$remove_dataset4 == 1)) # 0: so we know that we don't need to go any further with the iteration (as it is zero)

dim(pp_strain_names)[1] - dim(param_expok)[1] # total number of datasets removed this is the 171 single datasets + extra removed as not enough with the replicate

outside_datasets <- pp_strain_names %>% filter(remove_dataset_exp_iter == 1) %>% 
  dplyr::select(strain_name, rep, drytime, inocl,remove_dataset_exp_iter, outside)

# No strains are removed
length(unique(param_expok$strain_name)) # New with 5% of strain removed
length(unique(param$strain_name)) # Original total

setdiff(unique(param$strain_name),unique(param_expok$strain_name))
# 0 strains removed: all have at least one replicate

strains_typical = unique(param_expok$strain_name) # PERFECT strains

#####*************************** FILTERED plot - only the clean data *******************###############
all_strains = unique(param$strain)

## Add in label for odd exponential growth

#### CHANGE THIS TO MATCH THE ABOVE pp_strains selection... 
param <- left_join(param, pp_strain_names %>% dplyr::select("strain_name", "rep","drytime","inocl","remove_dataset_exp_iter","total_rep_rem"), by = c("strain_name", "rep","drytime","inocl"))
w1<-which(param$remove_dataset_exp_iter == 1) # single datasets to remove
w2<-which(param$total_rep_rem > 0) # total reps to remove
w <- union(w1,w2)
param[w,"odd_type_db"] <- paste0(param[w,"odd_type_db"],"5")

write_csv(param, "output/param_all_odd_labelled.csv") # save so can read in later 

pp_strain_names$strain <- pp_strain_names$strain_name
pp_strain_names$inoc <- pp_strain_names$inocl

ddm<- left_join(ddm, pp_strain_names %>% dplyr::select("strain", "rep","drytime","inoc","remove_dataset_exp_iter","total_rep_rem"), 
                 by = c("strain", "rep","drytime","inoc"))
w1<-which(ddm$remove_dataset_exp_iter == 1)
w2<-which(ddm$total_rep_rem > 0)
w <- union(w1,w2)
ddm[w,"odd_type_db"] <- paste0(ddm[w,"odd_type_db"],"5")

write.csv(ddm, "output/ddm_all_odd_labelled.csv")

cols = c(1,brewer.pal(n = 11, name = "Set3"),"#EDF8E9")
dir.create(file.path(here(), "plots/final_data_split_highlighted/"),showWarnings = FALSE)

for(jj in 1:length(all_strains)){ # for each strain
  
  # Want to keep the unique combination of replicate and dataset that are clean
  clean = param %>% filter(strain_name == all_strains[jj])
  
  ddm_strain <- ddm %>% filter(strain == all_strains[jj])
  ddm_orig_s <- ddm_orig %>% filter(strain == all_strains[jj])
  
  wc <- c()
  for(i in 1:length(clean[,1])){
    w1 <- intersect(which(ddm_strain$rep == clean[i,"rep"]),which(ddm_strain$inoc == clean[i,"inocl"]))
    wc<-c(wc,intersect(w1,which(ddm_strain$drytime == clean[i,"drytime"])))
  }
  
  dd <- ddm_strain[wc,] %>% group_by(strain, inoc, rep) %>% filter(Time <= shoulder_point_t, Time > 3.5)
  dd$odd_type <- as.character(dd$odd_type)
  ddm_orig_s$odd_type <- as.character(ddm_orig_s$odd_type)
  dd$odd_type_db <- as.character(dd$odd_type_db)
  ddm_orig_s$odd_type_db <- as.character(ddm_orig_s$odd_type_db)
  
  ggplot(dd, aes(x=Time, y = value_J)) + 
    geom_line(aes(group = inoc, col = odd_type_db, linetype = factor(inoc)), lwd = 2) + 
    facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
    scale_color_manual("Odd_type", 
                       breaks = c("0","14","24","34","124","134","05",
                                  "145","245","1245",
                                  "234","1234","345","1345",
                                  "2345", "12345"),
                       labels = c("None","Peak&Double",
                                  "Width&Double","Shoulder&Double","Peak Width&Double","Peak Shoulder&Double","ExpGr",
                                  "Peak Shoulder&ExpGr","Width Shoulder&ExpGr","Peak Width Shoulder&ExpGr",
                                  "Width Shoulder&Double", "Peak Width Shoulder&Double", "Shoulder Double&ExpGr","Peak Shoulder Double&ExpGr",
                                  "Width Shoulder Double&ExpGr","All"),
                       values = cols, drop = FALSE) + 
    scale_linetype_discrete("Inoculum size 10^x") + 
    geom_line(data =  ddm_orig_s, aes(group = inoc, col = odd_type_db, linetype = factor(inoc)), alpha = 0.2, lwd = 2) + 
    geom_point(data = dd, aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    geom_point(data = dd, aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    labs(y = "Heat flow (mW)") +
    ggtitle(all_strains[jj]) + 
    theme(legend.key.width=unit(2,"cm"))

  
  ggsave(paste0("plots/final_data_split_highlighted/",all_strains[jj],"_filtered.pdf"), height = 10, width = 15) # if any to highlight it is shown here
  ggsave(paste0("plots/final_data_split_highlighted/",all_strains[jj],"_filtered.png"), height = 10, width = 15) # if any to highlight it is shown here
}



######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** LINEAR MODEL 
######****************************######################################################################
dir.create(file.path(here(), "plots/linear_fit/"),showWarnings = FALSE)
dir.create(file.path(here(), "plots/output_fit/"),showWarnings = FALSE)
theme_set(theme_bw(base_size=20))

### Need to remove those with only 2 datapoints at time 0: perfect line
#how many? 
#table(table(param_expok %>% filter(drytime == 0) %>% dplyr::select(strain_name, inocl))) 
param_expok_filt <- param_expok %>% 
  group_by(strain_name, drytime, rep) %>% 
  mutate(ndata = n() + ifelse(drytime == 168, 3,0)) %>% # add 3 as only care about little data at time zero
  filter(ndata > 2) # if remove 0 timepoint then won't evaluate in linear fit so only need to remove these

# What strains have all / some data?
table(table(pp_strain_names$strain_name)) # Original data 
table(table(param_expok$strain_name)) # After exponential filtering
table(table(param_expok_filt$strain_name)) # After filtering for more than 2 datasets at baseline

######****** Predicting **********######################################################################
ggplot(po, aes(y=inocl,x = timepeak, group = strain_name, colour=factor(rep_st))) + 
  geom_point(size = 2) + facet_wrap(~strain_name) + 
  scale_x_continuous("Time to peak") + scale_y_continuous("Inoculum size (10^x)") + 
  scale_color_discrete("Experiment") 
ggsave(paste0("plots/linear_fit/time_to_peak_all.pdf"), width = 16, height =10 )

### PLOT variation in other indicators
#ggplot(param_typical, aes(y=inocl,x = v_m_h_flow, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = lag, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = exp_gr, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_typical, aes(y=inocl,x = auc, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)

######****** Linear model **********######################################################################
# Do for each rep, for each strain_name
reps <- unique(param_expok$rep)
strains <- unique(param_expok$strain_name)

ggplot(param_expok, aes(x=timepeak,y = log10(scalar), group = strain_name, colour=factor(drytime))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","168hr drying")) 
ggsave(paste0("plots/linear_fit/time_to_peak_all_as_linear_model.pdf"), width = 20, height = 20)
ggsave(paste0("plots/linear_fit/time_to_peak_all_as_linear_model.png"), width = 15, height = 15)

#### Fit linear model 
# Remove 10^2 and 10^6: not for reduction analysis
w26<-c(which(param_expok$inocl == 2), which(param_expok$inocl == 6))
if(length(w26) > 0){param_expok <- param_expok[-w26,]}

# If r^2 over 8 then good
r2_threshold = 0.75

###########*********** MODEL FIT ********************#########################

### RUN first time
reductions_fit <- fit_line_model(reps, strains, param_expok_filt, "timepeak","Time to max heat flow", R_cut = r2_threshold, plot = 1) ## plot = 1 will give the underling fit curves
# 
write_csv(reductions_fit$reductions, "plots/linear_fit/reductions_fit_reductions.csv") # save so can read in later
write_csv(reductions_fit$fit, "plots/linear_fit/reductions_fit_fit.csv") # save so can read in later

## READ in on subsequent times
#reductions_fit <- c() # initialise
#reductions_fit$reductions <- read_csv("plots/linear_fit/reductions_fit_reductions.csv") 
#reductions_fit$fit <- read_csv("plots/linear_fit/reductions_fit_fit.csv")

###########*###########***********###########***********###########***********
###########*
# ### Check fits
reductions_fit$fit$R2 <- as.numeric(reductions_fit$fit$R2)
ggplot(reductions_fit$fit, aes(x=strain, y = R2)) +
  geom_point() +
  ggtitle("R2 value for fit of linear model to inoculum size = a*time_to_peak + b")

ggplot(reductions_fit$fit, aes(x=R2)) + geom_histogram(binwidth = 0.02) +
  #geom_vline(xintercept = 0.9) + 
  geom_vline(xintercept = 0.75, lty = "dashed") 
ggsave("plots/linear_fit/cutoff_for_r2.pdf", width = 10, height = 10) # only get this if set R_cut = 0 in above?
ggsave("plots/linear_fit/cutoff_for_r2.png", width = 10, height = 10) # only get this if set R_cut = 0 in above?

fitted_strains <- reductions_fit$fit %>% filter(R2 > r2_threshold) %>%
  dplyr::select(strain) %>% unlist() %>% as.character() %>% unique()# but not all the replicates for these strains

# How many remove by filtering on R2? 
dim(reductions_fit$fit)
dim(reductions_fit$fit %>% filter(R2 > r2_threshold))
dim(reductions_fit$fit) - dim(reductions_fit$fit %>% filter(R2 > r2_threshold)) # 7 replicates
reductions_fit$fit %>% filter(R2 <= r2_threshold)

length(unique(param_expok_filt$strain_name)) # 98 into the function
length(unique(reductions_fit$fit$strain)) # 90 have values
length(fitted_strains) # 90 after filtered on R2

setdiff(unique(param_expok$strain_name),unique(reductions_fit$fit$strain))
setdiff(unique(param_expok_filt$strain_name),unique(reductions_fit$fit$strain))
setdiff(unique(reductions_fit$fit$strain),fitted_strains)


###############################******************** RESULTS *****************************###################################
######## Over all strains
# in reductions. Meas column key: 
# meas = 1 = log reduction 
# meas = 2 = percentage reduction 
# meas = 3 = inoculum
# meas = 4 = predicted inoculum
rw <- reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1) %>% pivot_longer(`10^2`:`10^6`) 
rw$inocl <- as.numeric(substr(rw$name,4,4))

theme_set(theme_bw(base_size=14))
#### Plot replicate variation
ggplot(reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1), aes(x=ticker, y = mean)) + 
  geom_bar(stat = "identity", aes(fill = factor(ticker))) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) + 
  #geom_point() + 
  facet_wrap(~strain_name) + 
  scale_fill_discrete("Replicate") + 
  scale_x_continuous("Replicate") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/linear_fit/replicate_variation_all.pdf", width = 15, height = 15)
ggsave("plots/linear_fit/replicate_variation_all.png", width = 15, height = 15)


repv <- reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1) %>% 
  ungroup() %>% 
  group_by(strain_name) %>%
  dplyr::select(strain_name, mean, ticker) %>%
  pivot_wider(names_from = ticker, values_from = mean) %>%
  mutate(diff1 = abs(ifelse(is.na(`2`),NA,ifelse(is.na(`1`),NA,`2`-`1`))),
         diff2 = abs(ifelse(is.na(`3`),NA,ifelse(is.na(`2`),NA,`3`-`2`))),
         diff3 = abs(ifelse(is.na(`3`),NA,ifelse(is.na(`1`),NA,`3`-`1`)))) %>%
  pivot_longer(cols = diff1:diff3) %>% 
  filter(!is.na(value))

pdf("plots/linear_fit/replicate_variable.pdf", width = 10, height = 10)
hist(repv$value, breaks = seq(0,5,0.2)) 
dev.off()

# All replicates
counts <- c()
for(i in seq(0,3,0.01)){
  counts <- c(counts, count(repv$value < i))
}
perc <- as.data.frame(100*counts / length(repv$value))
perc$thresh_diff = seq(0,3,0.01)
colnames(perc) <- c("perct","thresh_diff")

g1 <- ggplot(perc, aes(x=thresh_diff, y = perct)) + geom_point() + 
  geom_smooth(col = "black") + 
  scale_x_continuous("Difference between reduction from each replicate") + 
  scale_y_continuous("Percentage of replicates with this difference") + 
  geom_hline(yintercept =  c(50,90), lty = "dashed")

# Using mean value per strain 
repvm <- repv %>% group_by(strain_name) %>% 
  summarise(mean = mean(value))

counts <- c()
for(i in seq(0,3,0.01)){
  counts <- c(counts, count(repvm$mean < i))
}
perc <- as.data.frame(100*counts / length(repvm$mean))
perc$thresh_diff = seq(0,3,0.01)
colnames(perc) <- c("perct","thresh_diff")

g2 <- ggplot(perc, aes(x=thresh_diff, y = perct)) + geom_point() + 
  geom_smooth(col = "black") + 
  scale_x_continuous("Mean difference between\nreduction from each replicate\nper strain") + 
  scale_y_continuous("Percentage of strains with this difference") + 
  geom_hline(yintercept = c(50,90), lty = "dashed")

g1 + g2 
ggsave("plots/linear_fit/variation_over_reps.pdf")


###### RESULTS: OVERALL Take the average over all values: this is the overall log reduction 
rwo %>% 
  ungroup() %>% 
  group_by(strain_name) %>% 
  dplyr::summarise(mean_strain = mean(mean, na.rm = TRUE), sd_strain = sd(mean, na.rm = TRUE), .groups = "drop") %>% # This is the mean over the replicates for each strain
  ungroup() %>%
  dplyr::summarise(mean_all = mean(mean_strain, na.rm = TRUE), sd_all = sd(mean_strain, na.rm = TRUE), .groups = "drop") # This is the mean over the mean replicates for each strain

###### RESULTS: INOCULUM Take the average by inoculum 
av_inoc_all <- reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1) %>% 
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  dplyr::select(strain_name, name, value) %>% 
  dplyr::group_by(strain_name, name) %>% 
  dplyr::summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE), .groups = "drop")
av_inoc <- av_inoc_all %>%  group_by(name) %>% 
  dplyr::summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE), .groups = "drop")

av_inoc

# ANOVA of impact of inoculum
anova.inoc <- aov(mean_strain ~ name, data = av_inoc_all)
summary(anova.inoc)
TukeyHSD(anova.inoc)

av_inoc$lab = as.numeric(substr(av_inoc$name,4,4))

g1 <- ggplot(av_inoc, aes(x=lab, y = mean_inoc)) + geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  theme(legend.position = "none")

###### RESULTS: SUCCESS
succ <- read_csv("data/MACOTRA_100collection success_20210121.csv")
succ$strain_name <- as.character(succ$strain)

succ_go <- left_join(reductions_fit$reductions %>% filter(strain_name %in% fitted_strains),succ, by = "strain_name")

av_all_bys_all <- succ_go %>%
  filter(r2 > r2_threshold, meas == 1) %>%
  pivot_longer(cols = c('10^3':'10^5')) %>%
  group_by(strain_name, success) %>% 
  dplyr::summarise(mean_strain = mean(mean, na.rm = TRUE), sd_strain = sd(mean, na.rm = TRUE), .groups = "drop") # Mean = mean over inoc for this replicate. Mean_strain = mean over replicates

av_all_bys <- av_all_bys_all %>%  group_by(success) %>% 
  dplyr::summarise(mean_succ = mean(mean_strain, na.rm = TRUE), sd_succ = sd(mean_strain, na.rm = TRUE), .groups = "drop")

av_all_bys 

## T test of difference
t.test(mean_strain ~ success, data = av_all_bys_all)

ggplot(av_all_bys, aes(x=success, y = mean_succ)) + geom_bar(stat = "identity", aes(fill = success)) + 
  geom_errorbar(aes(ymin = mean_succ - sd_succ, ymax = mean_succ + sd_succ)) + 
  scale_y_continuous("Mean log reduction") + 
  scale_x_discrete("") + theme(legend.position = "none")
ggsave("plots/linear_fit/succ_unsucc_sd.pdf", width = 5, height = 5)


###### RESULTS: INOCULUM AND SUCCESS
# Take the average by inoculum and success
av_inoc_succ <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success) %>% 
  dplyr::summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success) %>% 
  dplyr::summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE))

av_inoc_succ$lab = as.numeric(substr(av_inoc_succ$name,4,4))

g2 <- ggplot(av_inoc_succ, aes(x=lab, y = mean_inoc,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")


###### RESULTS: LINEAGE 
av_inoc_succ_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, lineage) %>% 
  dplyr::summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success, lineage) %>% 
  dplyr::summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

av_inoc_succ_lin$lab = as.numeric(substr(av_inoc_succ_lin$name,4,4))

g3 <- ggplot(av_inoc_succ_lin, aes(x=lab, y = mean_inoc,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  facet_wrap(~lineage) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")


###### RESULTS: COUNTRY AND SUCCESS 
av_inoc_succ_country <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, country) %>% 
  dplyr::summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success, country) %>% 
  dplyr::summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

av_inoc_succ_country$lab = as.numeric(substr(av_inoc_succ_country$name,4,4))

g4 <- ggplot(av_inoc_succ_country, aes(x=lab, y = mean_inoc,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  facet_wrap(~country) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")

(g1 + g2) / (g3 + g4) + plot_layout(guides = 'collect', widths = c(1,2)) + plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')
ggsave("plots/final/figure4.pdf", width = 15, height = 15)

###### RESULTS: PLOT INDIVIDUAL DATA
v_inoc_succ_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, lineage, country) %>% 
  dplyr::summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) 

v_inoc_succ_lin$lab = as.numeric(substr(v_inoc_succ_lin$name,4,4))

ggplot(v_inoc_succ_lin, aes(x=lab, y = mean_strain,group = success)) + geom_point(aes(pch = success, col = country),position = position_dodge(0.8), size = 3, alpha = 0.8) + 
  facet_wrap(~lineage) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_colour_discrete("") + 
  scale_shape_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data.pdf", width = 10, height = 8)
ggsave("plots/final/underlying_all_data.png", width = 10, height = 8)

ggplot(v_inoc_succ_lin, aes(x=lab, y = mean_strain,group = success)) + geom_point(aes(pch = success, col = country),position = position_dodge(0.8), size = 3, alpha = 0.8) + 
  facet_wrap(~lineage) + 
  geom_smooth(aes(group = country, col = country, fill = country)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_colour_discrete("") + 
  scale_shape_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_lines.pdf")

ggplot(v_inoc_succ_lin, aes(x=lab, y = mean_strain,group = success)) + geom_point(aes(),position = position_dodge(0.8), size = 3, alpha = 0.8) + 
  facet_wrap(~lineage) + 
  geom_smooth(aes(group = success, col = success, fill = success)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  #scale_colour_discrete("") + 
  scale_shape_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_lines.pdf", width = 10, height = 8)


av_inoc_succ_country_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, country, lineage) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success, country, lineage) %>% 
  summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

av_inoc_succ_country_lin$lab = as.numeric(substr(av_inoc_succ_country_lin$name,4,4))

ggplot(av_inoc_succ_country_lin, aes(x=lab, y = mean_inoc,group = country)) + geom_bar(stat = "identity",position = "dodge", aes(fill = country)) + 
  facet_wrap(lineage~success) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_succ_bar.pdf", width = 10, height = 8)

av_inoc_country_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, country, lineage) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, country, lineage) %>% 
  summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

#av_inoc_country_lin <- av_inoc_country_lin %>% mutate(lab = as.numeric(paste0(substr(name,4,4),ifelse(country == "France",".2",ifelse(country == "Netherlands",".5",".7"))))) # doesn't work
av_inoc_country_lin$lab = as.numeric(substr(av_inoc_country_lin$name,4,4))

g5 <- ggplot(av_inoc_country_lin, aes(x=lab, y = mean_inoc, group = interaction(country, lab, lineage))) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = country)) + 
  facet_wrap(~lineage) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_bar.pdf", width = 10, height = 8)

(g1 + g2) / (g3 + g5) + plot_layout(guides = 'collect', widths = c(1,2)) + plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')
ggsave("plots/final/figure4_alternative.pdf", width = 15, height = 15)
ggsave("plots/final/figure4_alternative.png", width = 15, height = 15)

ggplot(av_inoc_succ_country_lin, aes(x=lab, y = mean_inoc,group = country)) + geom_bar(stat = "identity",position = "dodge", aes(fill = country)) + 
  facet_wrap(lineage~success) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_succ_bar.pdf", width = 10, height = 8)

## For multilevel modelling
mm <- v_inoc_succ_lin %>%
  rename(logred = mean_strain, inoc = name) %>% 
  filter(!is.na(logred)) %>%
  ungroup() %>%
  mutate(success_bin = ifelse(success == "Successful",1,0)) %>% 
  dplyr::select(strain_name, inoc, country, lineage, success, success_bin, logred)

write.csv(mm, "output/mm_final_data.csv")

