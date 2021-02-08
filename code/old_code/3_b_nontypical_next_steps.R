### TYPICAL STRAIN ANALYSIS
### From 3_next_steps.R have typical param set

source("20_function_linear_model.R")

##### READ IN DATA
### CHOOSE IN LINE WITH DATA CLEANING IN 2_analysis.R
name_code <- "1_2_6_7_"

# read in 
param_nontypical <- read.csv(paste0("output/",name_code,"param_non_typical.csv"), stringsAsFactors = FALSE)[,-1]
param_nontypical <- as.data.frame(param_nontypical)

infor <- read.csv(paste0("output/",name_code,"information_all_model_fit_params.csv"), stringsAsFactors = FALSE)

### (1) Remove those with different behaviour time 0 and time 168
non_typical_strains <- unique(param_nontypical$strain_name)

check_for_db <- c()

for(i in 1:length(non_typical_strains)){
  print(non_typical_strains[i])
  
  pp <- param_nontypical %>% filter(strain_name == non_typical_strains[i])
  reps <- unique(pp$rep)
  for(j in 1:length(reps)){
    ppr0 <- pp %>% filter(rep == reps[j]) %>% filter(drytime == 0)
    ppr168 <- pp %>% filter(rep == reps[j]) %>% filter(drytime == 168)
    
    o0 <- ppr0[which(ppr0$odd_type != 0),"odd_type"]
    o168 <- ppr168[which(ppr168$odd_type != 0),"odd_type"]
    
    print(sum(match(o0, o168, nomatch = 0)))
    
    if(sum(match(o0, o168, nomatch = 0)) < 2){ # if just one match still rubbish (for this rep)
      check_for_db <- rbind(check_for_db, c(reps[j],non_typical_strains[i]))
    }
  }
  
}

tb <- table(check_for_db[,2])


#### double curves 
p_double_curves <- read.csv(paste0("output/",name_code,"p_double_curves.csv"))[,-1]

## TRY DOUBLE curve if odd between reps.
db_fit <- names(which(tb > 1)) # Try a double fit to those with more than one odd rep

# Store information - which are different between 0 and 168 so try a double curve?
infor <- infor %>% mutate(removed_try_db = ifelse(strain_name %in% db_fit,1,0)) 
write.csv(infor,paste0("output/",name_code,"information_all_model_fit_params.csv"))

# ### 
# for(i in 1:length(db_fit)){
#   p <- (p_double_curves %>% filter(strain == db_fit[i]))
#   reps <- unique(p$replicate)
#   for(j in 1:length(reps)){
#     pp <- p %>% filter(replicate == reps[j])
#     x <- pp$time
#     pa <- pp[1,] # paraemters a-f same across rows 
#     y1 <- (pa$a/pa$b)*exp(-(x-pa$c)^2/(2*pa$b^2)) # curve 1
#     y2 <- (pa$d/pa$e)*exp(-(x-pa$f)^2/(2*pa$e^2)) # curve 2
#   }
#   
# }

## What are the strains?
u <- db_fit
## How many replicates? 
r <- unique(param_nontypical$rep) # replicates
# What are the inoculums? 
q <- unique(param_nontypical$inoc)

## Run thru each experiment (columns in original data) for each strain
# Fit separately to each replicate as otherwise don't have enough data for later predictions
# still having up to 3 experimental drytimes despite set2 not having 24hr data
# Where the parameters for each strain are stored
param_n <- matrix(0, length(u)*length(r)*length(q)*3, 4); # number of strains x number of replicates x number of experimental conditions
param <- matrix(0, length(u)*length(r)*length(q)*3, 10); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- c(0,24,168)
param_o <- c()

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to the first
# of the fitted double peaks
for(jj in 1:length(u)){ # for each strain
  for(ii in 1:length(r)){ # for each replicate: fit to all the data, not just each replicate
    for(kk in 1:length(drying_times)){ #each of the three experimental conditions (0, 24, 168)
      for(ll in 1:length(q)){ #each of the inoculums
        print(c(jj,ii,kk,ll))
        
        ## which rows of the data are this strain and replicate?
        wi <- intersect(which(p_double_curves$strain == u[jj]),which(p_double_curves$replicate == r[ii])) # if fit to each replicate
        wj <- intersect(wi, which(p_double_curves$condition == drying_times[kk]))
        w <- intersect(wj, which(p_double_curves$inocl == q[ll]))
        if(length(w)>0){
          print(c(u[jj]))
          dd <- p_double_curves[w,]
          wn1 <- which.max(dd$normal_curve1)
          tn1 <- dd[wn1,"time"]
          wn2 <- which.max(dd$normal_curve2)
          tn2 <- dd[wn2,"time"]
          if(tn1 < tn2){
            gc_fit <- gcFitSpline(dd$time, cumsum(dd$normal_curve1))
            # parameters from this fit
            s <- summary(gc_fit)
            param_o   <- rbind(param_o,
                               c(u[jj], r[ii], drying_times[kk], q[ll],
                                 tn1, dd[wn1,"value_J"], 
                                 s$mu.spline, s$lambda.spline,s$integral.spline))
            dd$used <- dd$normal_curve1
            dd$other <- dd$normal_curve2
          }else{
            gc_fit <- gcFitSpline(dd$time, cumsum(dd$normal_curve2))
            # parameters from this fit
            s <- summary(gc_fit)
            param_o   <- rbind(param_o,
                               c(u[jj], r[ii], drying_times[kk], q[ll],
                                 tn2, dd[wn2,"value_J"], 
                                 s$mu.spline, s$lambda.spline,s$integral.spline))
            dd$used <- dd$normal_curve1
            dd$other <- dd$normal_curve2
          }
          
          ggplot(dd, aes(x=time, y = value_J)) + geom_line() + 
            geom_line(aes(x=time, y = pred),col="blue") + 
            geom_line(aes(x=time, y = used), col = "red") + 
            geom_line(aes(x=time, y = other), col = "pink") + 
            scale_x_continuous("Time") 
          ggsave(paste0("output/", name_code,"dbo_",u[jj],"_",r[ii],"_",drying_times[kk],"_",q[ll],".pdf"))
          
          
        }
      }
    }
  }
} 

param_o <- as.data.frame(param_o)
param_oo <- param_o # backup
colnames(param_o) <- c("strain_name","rep","drytime","inocl","t_m_h_flow","v_m_h_flow","exp_gr","lag","auc")

param_o$v_m_h_flow <- as.numeric(as.character(param_o$v_m_h_flow))
param_o$t_m_h_flow <- as.numeric(as.character(param_o$t_m_h_flow))

##### Check height OK across inocula of the same strain at set drying times
perc <- 0.20
param_height_lab <- param_o %>%
  group_by(strain_name, rep, drytime) %>% 
  mutate(mean_peak_height = mean(v_m_h_flow),
         mean_peak_height_p10 = mean_peak_height + perc*mean_peak_height,
         mean_peak_height_m10 = mean_peak_height - perc*mean_peak_height,
         outside = ifelse(v_m_h_flow < mean_peak_height_m10 | v_m_h_flow > mean_peak_height_p10, 1, 0)) %>%
  ungroup() %>%
  group_by(strain_name,rep) %>% # group by strain and replicate
  dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate

#### check for height - remove if only one is "wrong" 
param_o_typical <- param_height_lab %>%
  filter(sum_outside_s_r < 2) # some have only one in the strain and replicate

length(unique(param_o_typical$strain_name)) 
# "1_2_6_":7 strains with double peaks have typical behaviour: similar heights

strains_o_typical = unique(param_o_typical$strain_name) # Double peak strains





######****************************######################################################################
######****************************######################################################################
######****************************######################################################################
#######****** TYPICAL STRAIN ANALYSIS
######****************************######################################################################

param_o_typical$inocl <- as.numeric(param_o_typical$inocl) 


######****** Predicting **********######################################################################
ggplot(param_o_typical, aes(y=inocl,x = t_m_h_flow, group = strain_name, colour=rep)) + 
  geom_point(size = 3) + facet_grid(drytime~rep + strain_name) + 
  scale_x_continuous("Time to peak") + scale_y_continuous("Inoculum size (10^x)") + 
  scale_color_discrete("Experiment") 
ggsave(paste0("output/",name_code,"dbo_time_to_peak_all.pdf"), width = 16, height =10 )

### PLOT variation in other indicators
#ggplot(param_o_typical, aes(y=inocl,x = v_m_h_flow, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_o_typical, aes(y=inocl,x = lag, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_o_typical, aes(y=inocl,x = exp_gr, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)
#ggplot(param_o_typical, aes(y=inocl,x = auc, group = strain_name,colour=rep)) + geom_point() + facet_wrap(~strain_name)

######****** Linear model **********######################################################################
# Change from inoc to scalar
# scalar = 1 at 10^6, 1/10th of each next dilution
param_o_typical$scalar <- 0

w<- which(param_o_typical$inocl == "2")
param_o_typical[w,"scalar"] = 10
w<- which(param_o_typical$inocl == "3")
param_o_typical[w,"scalar"] = 100
w<- which(param_o_typical$inocl == "4")
param_o_typical[w,"scalar"] = 1000
w<- which(param_o_typical$inocl == "5")
param_o_typical[w,"scalar"] = 10000
w<- which(param_o_typical$inocl == "6")
param_o_typical[w,"scalar"] = 100000

# Do for each rep, for each strain_name
reps <- levels(param_o_typical$rep)
strains <- levels(param_o_typical$strain_name)

ggplot(param_o_typical, aes(x=t_m_h_flow,y = log10(scalar), group = strain_name, colour=drytime)) + 
  geom_point(size = 3) + facet_wrap(~strain_name) + 
  scale_y_continuous("Log(inoculum)") + scale_x_continuous("Time to max heat flow (h)") + 
  scale_color_discrete("Experiment", labels = c("Baseline","24hr drying","168hr drying")) 
ggsave(paste0("output/",name_code,"dbo_time_to_peak_all_as_linear_model.pdf"), width = 16, height =10 )

#### Fit linear model 
reductions <- fit_line_model(reps, strains, param_o_typical)

ggplot(subset(reductions, meas == 1), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Log reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave(paste0("output/",name_code,"dbo_log_reductions.pdf"))

ggplot(subset(reductions, meas == 2), aes(x=ticker, y = mean, fill = factor(ticker))) + 
  geom_bar(stat = "identity", position=position_dodge(2)) + 
  facet_grid(dry~strain_name) + scale_y_continuous("Percentage reduction") + 
  scale_x_discrete("Replicate") + 
  scale_fill_discrete("Replicate") + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),position = "dodge")
ggsave(paste0("output/",name_code,"dbo_perc_reductions.pdf"))



