#### Fit linear model 

fit_line_model <- function(reps, strains, param_here, var, var_name = "Variable name",R_cut = 0.9, plot = 0){
  # reps = how many replicates
  # strains = strain names
  # param_here = the big set of parameters
  # var = the string parameter name 
  # R_cut = upper limit for R: needs to be higher than this
  # plot = 1 if want to plot outputs
  
  # Check name for parameter is a column heading 
  ifelse(var %in% colnames(param_here),"",print("ERROR input parameter"))
  
  var <- as.name(var)
  reductions <- c()
  fit <- c()
  
  for(i in 1:length(reps)){
    for(j in 1:length(strains)){
      
      print(c(reps[i],strains[j]))
      # fit linear model to the data for experiment "a" (no drying)
      w1<-which(param_here$drytime == 0)
      w2<-intersect(w1,which(param_here$strain_name == strains[j]))
      w3<-intersect(w2,which(param_here$rep == reps[i])) 
      
      if(length(w3) > 0){ # If there is data...
        dda <- param_here[w3,] # The data for this experiment for this strain
        
        ## USE variable to predict inoculum
        lm.1 <- lm(log10(dda$scalar) ~ as.numeric(unlist(dda[,var]))) # + dda$scalar2) ## Stuck with linear model as quadratic made no difference
        #lm.2 <- lm(dda$scalar ~ as.numeric(unlist(dda[,var])) + as.numeric(unlist(dda[,var]))^2)
        
        # Check if good fit 
        fit <- rbind(fit, c(reps[i],strains[j],round(summary(lm.1)$r.squared,2)))
        if(round(summary(lm.1)$r.squared,2) < R_cut){print(paste0("Poor fit for ", reps[i],", ", strains[j],". R^2 = ",round(summary(lm.1)$r.squared,2)));}
        
        xx <- c(min(as.numeric(unlist(dda[,var])) ) - 10, as.numeric(unlist(dda[,var])),
                max(as.numeric(unlist(dda[,var])) ) + 10 ) # the variable
        
        xx <- seq(min(as.numeric(unlist(dda[,var])) ) - 10,max(as.numeric(unlist(dda[,var])) ) + 10, length = 100) # the variable
        
        yy <- lm.1$coefficients[1] + lm.1$coefficients[2] * xx # + lm.1$coefficients[3] * xx^2 
        pred <- as.data.frame(cbind(xx,yy)) # prediction
        
        if(plot == 1){
          g1 <- ggplot(dda, aes_string(x = var, y = "scalar")) + geom_point(size = 3) + 
            geom_line(data = pred, aes(x=xx, y = 10^yy), col = "red") + 
            scale_x_continuous(var_name, limits = c(3,15)) + 
            scale_y_log10("Inoculum",limits = c(0.001,120000)) + 
            geom_text(x=7, y=2, label=paste0("R^2 = ", round(summary(lm.1)$r.squared,4))) + 
            ggtitle("Pre-drying: black = data, red = linear model fit")
          
          ### Not a great fit to the raw data
          g2 <- ggplot(dda, aes_string(x = var, y = "scalar")) + geom_point(size = 3) + 
            geom_line(data = pred, aes(x=xx, y = 10^yy), col = "red") + 
            scale_y_continuous("Inoculum", limits = c(0.001,120000)) + scale_x_continuous(var_name,limits = c(3,15)) +
            ggtitle("Pre-drying: black = data, red = linear model fit")
        }
        
        #if have equation of line - can work out reduction due to dehydration
        # 24hrs
        w1<-which(param_here$drytime == 24)
        w2<-intersect(w1,which(param_here$strain_name == strains[j]))
        w3b<-intersect(w2,which(param_here$rep == reps[i]))
        if(length(w3b) > 0 ){ # Not all have 24hr measures
          ddb <- param_here[w3,]
          
          ddb$pred_inoc <- lm.1$coefficients[1] + lm.1$coefficients[2] * as.numeric(unlist(ddb[,var]))
          ddb$pred_i_ten <- 10^ddb$pred_inoc
          
          if(plot == 1){
            g3 <- ggplot(ddb, aes_string(x=var, y = "scalar")) + geom_point(size = 3) + 
              geom_line(data = pred, aes(x=xx, y = 10^yy), col = "red") + 
              geom_point(data = ddb, aes_string(x=var, y = "pred_i_ten"), col="red",size = 3) + 
              scale_x_continuous(var_name, limits = c(3,15)) + 
              scale_y_log10("Inoculum",limits = c(0.001,120000)) 
            
            g4 <- ggplot(ddb, aes_string(x= var, y = "scalar")) + geom_point(size = 3) + 
              geom_line(data = pred, aes(x=xx, y = 10^yy), col = "red") + 
              geom_point(data = ddb, aes_string(x=var, y = "pred_i_ten"), col="red",size = 3) + 
              scale_y_continuous("Inoculum", limits = c(0.001,120000)) + scale_x_continuous(var_name,limits = c(3,15)) 
          }
        }
        # 168hrs
        w1<-which(param_here$drytime == 168)
        w2<-intersect(w1,which(param_here$strain_name == strains[j]))
        w3<-intersect(w2,which(param_here$rep == reps[i]))
        ddc <- param_here[w3,]
        
        ddc$pred_inoc <- lm.1$coefficients[1] + lm.1$coefficients[2] * as.numeric(unlist(ddc[,var]))
        ddc$pred_i_ten <- 10^ddc$pred_inoc
        
        if(plot == 1){
          g5 <- ggplot(ddc, aes_string(x=var, y = "scalar")) + geom_point(size = 3) + 
            geom_line(data = pred, aes(x=xx, y = 10^yy), col = "blue") + 
            geom_point(data = ddc, aes_string(x=var, y = "pred_i_ten"), col="blue",size = 3) + 
            scale_x_continuous(var_name, limits = c(3,15)) + 
            scale_y_log10("Inoculum",limits = c(0.001,120000)) + 
            ggtitle("After 168hr drying:\nblack = data (orig inoc vs new time to max heat)\nblue = linear model prediction of inoc given this time to peak")
          if(length(w3b) > 0 ){ g5 <- g5 + geom_point(data = ddb, aes_string(x=var, y = "pred_i_ten"), col="red",size = 3)}
          
          g6 <- ggplot(ddc, aes_string(x= var, y = "scalar")) + geom_point(size = 3) + 
            geom_line(data = pred, aes(x=xx, y = 10^yy), col = "blue") +
            scale_y_continuous("Inoculum", limits = c(0.001,120000)) + scale_x_continuous(var_name,limits = c(3,15)) +
            ggtitle("After 168hr drying:\nblack = data (orig inoc vs new time to max heat)\nblue = linear model prediction of inoc given this time to peak")
          if(length(w3b) > 0 ){ g6 <- g6 + geom_point(data = ddb, aes(x=var, y = 10^pred_inoc), col="red",size = 3)}
        }
        
        if(plot == 1){
          if(length(w3b) > 0 ){(g1|g2) / (g3|g4) / (g5|g6)
          }else{(g1|g2) /  (g5|g6)}
          ggsave(paste0("plots/output_fit/",strains[j],"_",reps[i],"_pred",".pdf"), width = 18, height = 15)
        }
        
        
        # Log reduction
        log_red_24 <- matrix(NA,3,9)
        perc_red_24 <- matrix(NA,1,5)
        log_red_168 <- matrix(NA,3,9)
        perc_red_168 <- matrix(NA,1,5)
        
        # negative = increase
        if(length(w3b) > 0){
          log_red_24[1,(ddb$inocl-1)] <- ddb$inocl - ddb$pred_inoc
          log_red_24[1,6:9] <-c(i,j,24,1)
          perc_red_24[1,(ddb$inocl-1)] <- 100*(10^ddb$inocl - 10^ddb$pred_inoc)/10^ddb$inocl
          
          log_red_24[2,(ddb$inocl-1)] <- ddb$inocl
          log_red_24[2,6:9] <-c(i,j,24,3)
          log_red_24[3,(ddb$inocl-1)] <- ddb$pred_inoc
          log_red_24[3,6:9] <-c(i,j,24,4)
          
          reductions <- rbind(reductions,log_red_24)
          reductions <- rbind(reductions,
                              c(perc_red_24, i, j, 24, 2))
          
        }
        log_red_168[1,(ddc$inocl-1)] <- ddc$inocl - ddc$pred_inoc
        log_red_168[1,6:9] <-c(i,j,168,1)
        perc_red_168[1,(ddc$inocl-1)] <- 100*(10^ddc$inocl - 10^ddc$pred_inoc)/10^ddc$inocl
        
        log_red_168[2,(ddc$inocl-1)] <- ddc$inocl
        log_red_168[2,6:9] <-c(i,j,168,3)
        log_red_168[3,(ddc$inocl-1)] <- ddc$pred_inoc
        log_red_168[3,6:9] <-c(i,j,168,4)
        
        log_red_168 <- cbind(log_red_168,round(summary(lm.1)$r.squared,3))
        
        reductions <- rbind(reductions,log_red_168)
        reductions <- rbind(reductions,
                            c(perc_red_168, i, j, 168, 2,round(summary(lm.1)$r.squared,3)))
      }else{print(c("no data for:",reps[i], strains[j]))}
      
      # Expect same percentage reduction for each inoculum? 
      # if so then could take the average across these reps to give the average percentage/log
      # reduction for each drying period
    }
  }
  
  reductions <- as.data.frame(reductions)
  colnames(reductions) <- c("10^2","10^3","10^4","10^5","10^6","rep","strain_name","dry","meas","r2") 
  
  reductions$strain_name <- as.character(strains[reductions$strain_name])
  
  reductions$mean <- rowMeans(reductions[c("10^6","10^5","10^4","10^3","10^2")], na.rm=TRUE)
  reductions$sd <- 0
  for(i in 1:length(reductions[,1])){
    reductions[i,"sd"] = sd(reductions[i,c("10^6","10^5","10^4","10^3","10^2")],na.rm = TRUE)}
  
  reductions$number = 1
  reductions <- as.data.frame(reductions) %>%
    group_by(strain_name,dry,meas) %>%
    dplyr::mutate(ticker = dplyr::row_number()) ## Just counts replicates as all had different names
  
  ###### 4 results per strain: (1) log reduction for strain overall (2) each inoculum 
  av_for_strain <- reductions %>% filter(meas == 1) %>% group_by(strain_name, dry)  %>%
    dplyr::summarise(mean_strain = mean(mean, na.rm=TRUE), sd_strain = mean(sd, na.rm = TRUE),n_vals = n(),.groups = 'drop')
  
  av_for_inoculum <- reductions %>% dplyr::select(-c("10^2","10^6")) %>% pivot_longer(cols = c("10^3","10^4","10^5")) %>%
    filter(!is.na(value)) %>%
    group_by(strain_name, dry, name) %>% filter(meas == 1) %>%
    dplyr::summarise(mean_inoc = mean(value, na.rm=TRUE), sd_inoc = sd(value, na.rm = TRUE), 
              n_vals = n(),.groups = 'drop')
  
  ### Fit
  fit <- as.data.frame(fit)
  colnames(fit) <- c("rep", "strain","R2")
  
  return(list(reductions = reductions, av_for_strain = av_for_strain, av_for_inoculum = av_for_inoculum, fit = fit)) 
}
