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
theme_set(theme_bw(base_size=14)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

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
cutoff <- 0.36

pp_strain_names <- po %>%
  group_by(strain_name, rep_st) %>% 
  dplyr::mutate(mean_peak_exp_gr = mean(cut_exp),
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + cutoff*mean_peak_exp_gr,
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - cutoff*mean_peak_exp_gr,
         outside = ifelse(cut_exp < mean_peak_exp_gr_m10 | cut_exp > mean_peak_exp_gr_p10, 1, 0),
         diff = ifelse(outside == 1, abs(cut_exp - mean_peak_exp_gr), 0),
         max = max(diff), 
         remove_dataset = ifelse(diff == max, ifelse(max > 0,1,0), 0)) %>%
  ungroup() %>% 
  group_by(strain_name, rep,drytime) %>% 
  dplyr::mutate(total_outside_inrep = sum(outside), total_inrep = n(), perc_outside = 100*total_outside_inrep / total_inrep) %>%
  ungroup() %>% 
  group_by(strain_name, rep) %>% 
  dplyr::mutate(rep_dt_remove = ifelse(perc_outside > 34, 1, ifelse(total_inrep < 2,1,0)), total_rep_rem = sum(rep_dt_remove)) %>%
  ungroup() %>%
  dplyr::mutate(keep_rep =ifelse((total_rep_rem == 0),rep,0)) %>%
  group_by(strain_name) %>%
  dplyr::mutate(remove_strain = ifelse((n_distinct(keep_rep) - any(keep_rep == 0)) >=2,0,1)) %>% # n_distinct counts 0 so need to remove
  group_by(strain_name, rep_st, remove_dataset) %>%
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
  dplyr::mutate(remove_dataset_exp_iter = remove_dataset + remove_dataset2 + remove_dataset3 + remove_dataset4)

pp_strain_names %>% dplyr::select(strain_name, rep, drytime, inocl, cut_exp, mean_peak_exp_gr_m10, mean_peak_exp_gr, mean_peak_exp_gr_p10, outside)


pp_strain_names %>% filter(remove_strain == 0, total_rep_rem == 0) %>% 
  dplyr::select(strain_name, rep_st, drytime, inocl, remove_strain, remove_dataset, remove_dataset2,cut_exp,
         #        cut_exp, mean_peak_exp_gr_m10, mean_peak_exp_gr, mean_peak_exp_gr_p10, outside, diff, max, remove_dataset,
         #        mean_peak_exp_gr_m102, mean_peak_exp_gr2, mean_peak_exp_gr_p102, outside2, diff2, max2, remove_dataset2,
         mean_peak_exp_gr_m103, mean_peak_exp_gr3, mean_peak_exp_gr_p103, outside3, diff3, max3, remove_dataset3) %>%
  #       mean_peak_exp_gr_m104, mean_peak_exp_gr4, mean_peak_exp_gr_p104, outside4, diff4, max4, remove_dataset4) %>%
  filter(strain_name == "11271")



#### Remove those with exponential values outside the above range
param_expok <- pp_strain_names %>%
  filter(remove_strain == 0) %>% # remove those strains with more than 2 wrong reps
  filter(total_rep_rem == 0) %>% # remove those reps with more than 2 outside
  #filter(outside == 0) # remove those datasets outside the range: this is with no iterative calculation of the mean
  filter(remove_dataset == 0, remove_dataset2 == 0, remove_dataset3 == 0, remove_dataset4 == 0) # remove those in the interative calculation of the mean

## Which are removed? 
length(which(pp_strain_names$remove_strain == 1)) # 102 datasets removed due to strain being removed
length(which(pp_strain_names$total_rep_rem == 1)) # 4 single reps removed from strains
length(which(pp_strain_names$outside == 1)) # 201 single datasets removed from strains

p2 <- pp_strain_names %>% filter(remove_strain == 0) %>% # remove those strains with more than 2 wrong reps
  filter(total_rep_rem == 0)
length(which(p2$remove_dataset == 1)) # 94 Each of these is remove one dataset and then recalculate the mean within a replicate
length(which(p2$remove_dataset2 == 1)) # 25
length(which(p2$remove_dataset3 == 1)) # 2
length(which(p2$remove_dataset4 == 1)) # 0: so we know that we don't need to go any further with the iteration (as it is zero)

length(which(p2$outside == 1)) # 126 single datasets removed from strains if do not recalculate the mean (after exclude those reps and strains)
length(which(p2$remove_dataset == 1)) + length(which(p2$remove_dataset2 == 1)) + length(which(p2$remove_dataset3 == 1)) + length(which(p2$remove_dataset4 == 1)) # only remove 121

dim(pp_strain_names)[1] - dim(param_expok)[1]

outside_datasets <- pp_strain_names %>% filter(outside == 1) %>% dplyr::select(strain_name, rep_st, drytime, inocl, outside)
outside_datasets2 <- pp_strain_names %>% filter(remove_dataset_exp_iter == 1) %>% dplyr::select(strain_name, rep_st, drytime, inocl)

length(unique(param_expok$strain_name)) # New with 5% of strain removed
length(unique(param$strain_name)) # Original total

setdiff(unique(param$strain_name),unique(param_expok$strain_name))
# Only 6 strains removed: we can't use these are their exponential growth is too variable

strains_typical = unique(param_expok$strain_name) # PERFECT strains

#####*************************** FILTERED plot - only the clean data *******************###############
all_strains = unique(param$strain)

## Add in label for odd exponential growth
param[which(!param$strain_name %in% strains_typical),"odd_type_db"] <- paste0(param[which(!param$strain_name %in% strains_typical),"odd_type_db"],"5")
ddm[which(!ddm$strain %in% strains_typical),"odd_type_db"] <- paste0(ddm[which(!ddm$strain %in% strains_typical),"odd_type_db"],"5")

cols = c(1,brewer.pal(n = 11, name = "Set3"))
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
                                  "145","245",
                                  "234","1234","345"),
                       labels = c("None","Peak&Double",
                                  "Width&Double","Shoulder&Double","Peak Width&Double","Peak Shoulder&Double","ExpGr",
                                  "Peak Shoulder&ExpGr","Width Shoulder&ExpGr",
                                  "Width Shoulder&Double", "Peak Width Shoulder&Double", "Shoulder Double&ExpGr"),
                       #breaks = c("0","1","2","3","4","5",
                       #                        "12","13","23","123","14",
                       #                        "24","34","124","134","234",
                       #                        "1234","05",
                       #                        "15","25","35","45",
                       #                        "145","245","345"), 
                       #                        "1235","235","1235","135",
                       # labels = c("None","Peak","Width","Shoulder","Double","ExpGr",
                       #            "Peak&Width","Peak&Shoulder","Width&Shoulder","Peak Width&Shoulder","Peak&Double",
                       #            "Width&Double","Shoulder&Double","Peak Width&Double","Peak Shoulder&Double","Width Shoulder&Double",
                       #            "Peak Width Shoulder&Double","ExpGr",
                       #            "Peak&ExpGr","Width&ExpGr","Shoulder&ExpGr","Double&ExpGr",
                       #            "Peak Width Shoulder&ExpGr","Width Shoulder&ExpGr","Peak Width Shoulder&ExpGr","Peak Shoulder&ExpGr",
                       #            "Peak Double&ExpGr","Width Double&ExpGr","Shoulder Double ExpGr"),
                       values = cols, drop = FALSE) + 
    scale_linetype_discrete("Inoc.") + 
    geom_line(data =  ddm_orig_s, aes(group = inoc, col = odd_type_db, linetype = factor(inoc)), alpha = 0.2, lwd = 2) + 
    geom_point(data = dd, aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    geom_point(data = dd, aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    ggtitle(all_strains[jj])
  
  # ggplot(dd, aes(x=Time, y = value_J)) + 
  #   geom_line(aes(group = inoc, col = odd_type, linetype = factor(inoc)), lwd = 1) + 
  #   facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
  #   scale_color_manual("Odd_type", breaks = c("0","1","2","3","12","13","23","123"), 
  #                      labels = c("None","Peak","Width","Shoulder","Peak&Width","Peak&Shoulder",
  #                                 "Width&Shoulder","Peak Width&Shoulder"),
  #                      values = cols, drop = FALSE) + 
  #   scale_linetype_discrete("Inoc.") + 
  #   geom_line(data =  ddm_orig_s, aes(group = inoc, col = odd_type, linetype = factor(inoc)), alpha = 0.2, size = 1) + 
  #   geom_point(aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
  #   geom_point(data = dd, aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
  #   ggtitle(all_strains[jj])
  
  ggsave(paste0("plots/final_data_split_highlighted/",all_strains[jj],"_filtered.pdf"), width = 15) # if any to highlight it is shown here
}



######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** LINEAR MODEL 
######****************************######################################################################

### Need to remove those with only 2 datapoints at time 0: perfect line
#how many? 
table(table(param_expok %>% filter(drytime == 0) %>% dplyr::select(strain_name, inocl))) # 62 + 15 = 77... 
param_expok_filt <- param_expok %>% group_by(strain_name, drytime, rep) %>% 
  mutate(ndata = n() + ifelse(drytime == 168, 3,0)) %>% # add 3 as only care about little data at time zero
  filter(ndata > 2)

# What strains have all / some data?
table(table(pp_strain_names$strain_name)) # Original data 
table(table(param_expok$strain_name)) # After exponential filtering
table(table(param_expok_filt$strain_name)) # After exponential filtering

######****** Predicting **********######################################################################
ggplot(po, aes(y=inocl,x = timepeak, group = strain_name, colour=factor(rep_st))) + 
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

ggplot(param_expok, aes(x=timepeak,y = log10(scalar), group = strain_name, colour=factor(drytime))) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","168hr drying")) 
ggsave(paste0("plots/fit/time_to_peak_all_as_linear_model.pdf"), width = 20, height = 20)

#### Fit linear model 
# Remove 10^2 and 10^6: not for reduction analysis
w26<-c(which(param_expok$inocl == 2), which(param_expok$inocl == 6))
if(length(w26) > 0){param_expok <- param_expok[-w26,]}

# If r^2 over 8 then good
r2_threshold = 0.75

###########*********** MODEL FIT ********************#########################


reductions_fit <- fit_line_model(reps, strains, param_expok_filt, "timepeak","Time to max heat flow", R_cut = r2_threshold, plot = 1) ## plot = 1 will give the underling fit curves


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
ggsave("plots/fit/cutoff_for_r2.pdf", width = 10, height = 10) # only get this if set R_cut = 0 in above?

fitted_strains <- reductions_fit$fit %>% filter(R2 > r2_threshold) %>%
  dplyr::select(strain) %>% unlist() %>% as.character() %>% unique()# but not all the replicates for these strains

dim(reductions_fit$fit)
dim(reductions_fit$fit %>% filter(R2 > r2_threshold))
dim(reductions_fit$fit) - dim(reductions_fit$fit %>% filter(R2 > r2_threshold))


length(unique(param_expok_filt$strain_name)) # 92 into the function
length(unique(reductions_fit$fit$strain)) # 83 Not filtered on R2
length(fitted_strains) # 83 Filtered on R2

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


#### Can compare which datasets were removed if use straight outside mean calculation or iterative
# Straight outside calculation
# outside_datasets$ticker <- outside_datasets$rep_st 
# outside_datasets168 <- outside_datasets %>% filter(drytime == 168) %>% dplyr::select(-c(rep_st, drytime))
# outside_datasets168$outside168 <- outside_datasets168$outside
# 
# rwoo <- left_join(rw, outside_datasets %>% filter(drytime == 0) %>% dplyr::select(-c(rep_st, drytime)))
# rwo <- left_join(rwoo, outside_datasets168 %>% dplyr::select(-c(outside))) %>%
#   mutate(outside_f = ifelse(is.na(outside), ifelse(is.na(outside168),"NA",1),1))
# 
# ggplot(rwo, aes(x=inocl, y = value)) + geom_point(aes(col = factor(outside_f))) + 
#   facet_wrap(~strain_name) + 
#   scale_color_discrete("Outside or in\nthe range?", breaks = c(1, "NA"), labels = c("Yes","No")) + 
#   ggtitle("Exclude all initially outside range")
# ggsave("output/reductions_exp_outlier.pdf")

# Iterative outside calculation
outside_datasets2$ticker <- outside_datasets2$rep_st 
outside_datasets2$outside <- 1
outside_datasets2168 <- outside_datasets2 %>% filter(drytime == 168) %>% dplyr::select(-c(rep_st, drytime))
outside_datasets2168$outside168 <- 1

rwoo <- left_join(rw, outside_datasets2 %>% filter(drytime == 0) %>% dplyr::select(-c(rep_st, drytime)))
rwo <- left_join(rwoo, outside_datasets2168 %>% dplyr::select(-c(outside))) %>%
  mutate(outside_f = ifelse(is.na(outside), ifelse(is.na(outside168),"NA",1),1))

ggplot(rwo, aes(x=inocl, y = value)) + geom_point(aes(col = factor(outside_f))) + 
  facet_wrap(~strain_name) + 
  scale_color_discrete("Outside or in\nthe range?", breaks = c(1, "NA"), labels = c("Yes","No")) + 
  ggtitle("Iterative exp exclusion")
ggsave("output/reductions_exp_outlier_iterative.pdf", width = 10, height = 10)


#### Average
# Take the average over all values
rwo %>% 
  ungroup() %>% 
  group_by(strain_name) %>% 
  dplyr::summarise(mean_strain = mean(mean, na.rm = TRUE), sd_strain = sd(mean, na.rm = TRUE), .groups = "drop") %>% # This is the mean over the replicates for each strain
  ungroup() %>%
  dplyr::summarise(mean_all = mean(mean_strain, na.rm = TRUE), sd_all = sd(mean_strain, na.rm = TRUE), .groups = "drop") # This is the mean over the mean replicates for each strain

#### Replicate variation
ggplot(reductions_fit$reductions %>% filter(r2 > r2_threshold, meas == 1), aes(x=ticker, y = mean)) + 
  geom_bar(stat = "identity", aes(fill = factor(ticker))) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) + 
  #geom_point() + 
  facet_wrap(~strain_name) + 
  scale_fill_discrete("Replicate") + 
  scale_x_continuous("Replicate") + 
  scale_y_continuous("Mean log reduction")
ggsave("plots/final/replicate_variation_all.pdf", width = 15, height = 15)


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

pdf("plots/final/replicate_variable.pdf", width = 10, height = 10)
hist(repv$value, breaks = seq(0,5,0.2)) 
dev.off()

ggplot(repv, aes(x=strain_name, y = value)) + geom_point(aes(col = name))

# Take the average by inoculum 
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

## Success? 
succ <- read.csv("data/MACOTRA 100collection success_20210121.csv") #Gwen
#succ <- read.csv2("data/MACOTRA 100collection success_20210121.csv") #Valérie
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
ggsave("plots/fit/succ_unsucc_sd.pdf", width = 5, height = 5)


## Take the average by inoculum and success
av_inoc_succ <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success) %>% 
  summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE))

av_inoc_succ$lab = as.numeric(substr(av_inoc_succ$name,4,4))

g2 <- ggplot(av_inoc_succ, aes(x=lab, y = mean_inoc,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")


### Lineage
av_inoc_succ_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, lineage) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success, lineage) %>% 
  summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

av_inoc_succ_lin$lab = as.numeric(substr(av_inoc_succ_lin$name,4,4))

g3 <- ggplot(av_inoc_succ_lin, aes(x=lab, y = mean_inoc,group = success)) + geom_bar(stat = "identity",position = "dodge", aes(fill = success)) + 
  facet_wrap(~lineage) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")

# individual data
v_inoc_succ_lin <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, lineage, country) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) 

v_inoc_succ_lin$lab = as.numeric(substr(v_inoc_succ_lin$name,4,4))

ggplot(v_inoc_succ_lin, aes(x=lab, y = mean_strain,group = success)) + geom_point(aes(pch = success, col = country),position = position_dodge(0.8), size = 3, alpha = 0.8) + 
  facet_wrap(~lineage) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_colour_discrete("") + 
  scale_shape_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data.pdf", width = 10, height = 8)

ggplot(v_inoc_succ_lin, aes(x=lab, y = mean_strain,group = success)) + geom_point(aes(pch = success, col = country),position = position_dodge(0.8), size = 3, alpha = 0.8) + 
  facet_wrap(~lineage) + 
  geom_smooth(aes(group = country, col = country, fill = country)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  #scale_colour_discrete("") + 
  scale_shape_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_lines.pdf", width = 10, height = 8)

#### Country and success
av_inoc_succ_country <- succ_go %>% 
  filter(r2 > r2_threshold, meas == 1) %>%
  ungroup() %>% 
  pivot_longer(`10^3`:`10^5`) %>%
  group_by(strain_name, name, success, country) %>% 
  summarise(mean_strain = mean(value, na.rm = TRUE), sd_strain = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name, success, country) %>% 
  summarise(mean_inoc = mean(mean_strain, na.rm = TRUE), sd_inoc = sd(mean_strain, na.rm = TRUE)) %>% ungroup()

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
ggsave("plots/final/figure1.pdf", width = 15, height = 15)

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

av_inoc_country_lin$lab = as.numeric(substr(av_inoc_country_lin$name,4,4))

ggplot(av_inoc_country_lin, aes(x=lab, y = mean_inoc,group = country)) + geom_bar(stat = "identity",position = "dodge", aes(fill = country)) + 
  facet_wrap(~lineage) + 
  geom_errorbar(position=position_dodge(),aes(ymin = mean_inoc - sd_inoc, ymax = mean_inoc + sd_inoc)) + 
  scale_x_continuous("Inoculum", breaks = c(3,4,5), labels = function(x) parse(text=paste("10^",x))) + 
  scale_y_continuous("Mean log reduction") + 
  scale_fill_discrete("") + 
  theme(legend.position="bottom")
ggsave("plots/final/underlying_all_data_country_bar.pdf", width = 10, height = 8)

## For multilevel modelling
mm <- v_inoc_succ_lin %>%
  rename(logred = mean_strain, inoc = name) %>% 
  filter(!is.na(logred)) %>%
  ungroup() %>%
  mutate(success_bin = ifelse(success == "Epidemic",1,0)) %>% 
  dplyr::select(strain_name, inoc, country, lineage, success, success_bin, logred)

write.csv(mm, "output/mm_final_data.csv")

