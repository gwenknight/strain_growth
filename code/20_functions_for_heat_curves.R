#### FUNCTIONS FOR HEAT MAP
## (1) Read in data


######************************************************************************************************************
### (1) Read in data
## loop thru experiment names
## !! Still too set specific so no using very much

read_in_data <- function(name){  
  
  new_data <- c()
  
  da <- paste0("data/Exp",name,"a.txt")
  db <- paste0("data/Exp",name,"b.txt")
  dc <- paste0("data/Exp",name,"c.txt")
  
  if(file.exists(da)){dda <- read.table(da,header = TRUE)} # read in the data 
  if(file.exists(db)){ddb <- read.table(db,header = TRUE)} # read in the data 
  if(file.exists(dc)){ddc <- read.table(dc,header = TRUE)} # read in the data 
  
  # which experiment?
  dda$exp <- "a"
  ddb$exp <- "b"
  ddc$exp <- "c"
  
  # which replicate?
  dda$rep <- name
  ddb$rep <- name
  ddc$rep <- name
  
  new_data <- rbind(rbind(rbind(new_data, dda), ddb), ddc)
}


######************************************************************************************************************
### (2) Fit growth curve to data and return parameters


fit_growth_curve <- function(strain, replicate, condition, inocl, data, plot = 1, 
                             thresh_wide = 80, thresh_peak = 30, thresh_shoulder = 2.2e-06){
  
  ## which rows of the data are this strain and replicate?
  wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
  wj <- intersect(wi, which(data$drytime == condition))
  w <- intersect(wj, which(data$inoc == as.numeric(inocl)))
  
  ## is this strain replicate ODD? 
  odd_peak <- 0 # 0 = one peak
  odd_width <- 0 # 0 = not a wide peak
  odd_shoulder <- 0 # 0 = no shoulder
  max_level <- 0 # height of peak
  
  if(length(w) > 0){ # if this replicate exists for this strain
    data1 <- data[w,]  # just get the data for this experiment (time, value, drying time)
    
    ## GC # (see fig 3 of vv33i07.pdf)
    ## This gives lag time and exponential growth rate
    gc_fit <- gcFitSpline(data1$Time, data1$csum)
    
    ## What is the maximum heat flow and when?
    wmax <- which.max( data1$value_J[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
    time_max_heat_flow <- as.numeric(data1[wmax,"Time"])
    value_max_heat_flow <- as.numeric(data1[wmax,"value_J"])
    
    ## ODD 
    # (1) Is the peak broad? 
    interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
    interval_value <- as.numeric(unlist(data1[c(which(data1$Time == interval_peak[1]), which(data1$Time == interval_peak[2])),"value_J"]))
    
    max_level <- 0;
    max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest
    
    # (2) Are there multiple peaks? 
    # GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
    peaks_index = find_peaks(data1$value_J, m = 2)
    # REMOVE - early ones
    w<-which(data1[peaks_index,"Time"]>3) 
    # REMOVE - late ones
    w<-intersect(w,which(data1[peaks_index,"Time"]< 0.95*max(data1$Time))) 
    peaks_index <- peaks_index[w] # remove the ones at beginning / end
    # Sort by height - want to keep the tallest
    o <- order(data1[peaks_index,"value_J"], decreasing = "TRUE")
    peaks_index <- peaks_index[o]
    
    # When peaks? 
    time_peaks <- as.numeric(unlist(data1[peaks_index,"Time"]))
    time_peaks
    
    time_peaks_diff <- time_peaks - time_peaks[1]
    keep_time_far_apart <- which(abs(time_peaks_diff) > 5) # If multiple far apart then issue: double peaks
    
    if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
    if(any(data1[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low
    
    # If close and same height then odd (peak decline plateau decline OK)
    close_peaks <- which(abs(time_peaks_diff) <= 5)
    close_peaks_i <- 1
    if(length(close_peaks)>1){
      for(i in 2:length(close_peaks)){
        ifelse(data1[peaks_index[close_peaks[i]],"value_J"]/data1[peaks_index[close_peaks[1]],"value_J"]> 0.9,
               close_peaks_i <- c(close_peaks_i,i),"")
      }}
    
    if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks
    
    time_peaks <- time_peaks[c(1,keep_time_far_apart)]
    peaks_index <- peaks_index[c(1,keep_time_far_apart)]
    
    ### Visualise where peaks are if needed
    #plot(data1$Time,data1$value_J, type = 'l')
    #points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'red', pch = 19)
    
    # (3) Is there a shoulder? 
    # GIVES squared distance from line to curve post peak - if far from line then there is a "shoulder"
    ###### SHOULDER
    ### (1) Draw straight line from peak
    time_endline <- max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak
    where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
    value_endline <- data1[where.end, "value_J"] # Find value at this time point
    ## Start of line is where? 
    time_startline <- time_peaks[1]
    value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
    
    # Draw straight line, assuming peak time = time 0. 
    grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline))
    times_line <- as.numeric(unlist(data1[peaks_index[1]:where.end,"Time"])) # The times for the line
    
    pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
    
    # plot against peak point 
    #lines(times_line, pred_points_fit, col= "blue")
    
    #### Check not plateauing before 10 hrs after peak
    dist <- pred_points_fit - data1[peaks_index[1]:where.end,"value_J"]
    
    # Want there to be some wobble on the line - if the predicted line never goes below actual data line then 
    # the distance will always be positive. Want there to be some later negative values or at least a close fit initially
    # Hence the removal of the first few terms - near the peak is always some wobble
    
    # While there are not enough negative points: i.e. not enough overlap between the straight line and the curve
    # And while the time of the endline is not too close to peak
    while(all(round(dist[-c(1:10)],10) >= 0) & time_endline < (time_peaks[1]-3)){ 
      
      time_endline <- as.numeric(time_endline + 1) # 10 hrs from first peak # SUBTRACT ONE EACH TIME
      where.end <- which.min(abs(data1$Time - time_endline))
      value_endline <- as.numeric(data1[where.end, "value_J"]) # Find value nearest this time point
      
      time_startline <- time_peaks[1]
      value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
      
      # straight line 
      grad = (value_endline - value_startline)/(time_endline - time_startline)
      times_line <- as.numeric(unlist(data1[peaks_index[1]:where.end,"Time"]))
      
      pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
      
      dist <- pred_points_fit - as.numeric(unlist(data1[peaks_index[1]:where.end,"value_J"]))
      
    }
    
    squared_dist <- sum(dist^2)
    fpsq <- find_peaks(as.numeric(unlist(dist)), m = 5) # if there are multiple peaks then have the bump
    if(length(fpsq)>1){odd_shoulder <- 1} # otherwise just a smooth curve 
    #if(squared_dist > thresh_shoulder & length(fpsq)>1){odd_shoulder <- 1}
    
    
    
    # max_dist <- max(abs(dist))
    #if(max_dist > 0.0007){odd_shoulder <- 1}
    
    ##### Double curves? 
    
    x <- data1$Time
    y <- data1$value_J
    
    
    startl=list(a=unlist(data1[peaks_index,"value_J"])[1]*1.5, 
                b=diff(interval_peak)/3, 
                c=unlist(data1[peaks_index,"Time"])[1],
                d=unlist(data1[peaks_index,"value_J"])[1]*1.5, 
                e=diff(interval_peak), 
                f=unlist(data1[peaks_index,"Time"])[1])
    
    fit1 <- 0
    plot_p <- 0
    p <- 0
    try(fit1 <- nls(y~(a/b)*exp(-(x-c)^2/(2*b^2))+(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl), silent = TRUE)
    
    if(length(fit1) > 1){
      print(paste("Double curve fit", strain," ", replicate, condition, inocl, sep = " "))
      pred_p <- predict(fit1)
      
      p <- coef(fit1)
      g1p <- (p["a"]/p["b"])*exp(-(x-p["c"])^2/(2*p["b"]^2)) #### PULL OUT THE PARAMETSR
      g2p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
      
      plot_p <- as.data.frame(cbind(x,y,pred_p, g1p, g2p))
      colnames(plot_p) <- c("time","value_J","fit","normal_curve1","normal_curve2")
      plot_p$a <- p["a"]
      plot_p$b <- p["b"]
      plot_p$c <- p["c"]
      plot_p$d <- p["d"]
      plot_p$e <- p["e"]
      plot_p$f <- p["f"]
      
    }
    
    
    
    
    
    ## Plot data1, cumulative and fit
    ## data1:
    if(plot == 1){
      datam <- melt(data1[,c("Time","value_J","csum")], id.vars = "Time")
      gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() + facet_wrap(~strain, scales = "free") 
      ## model fit:
      gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
      colnames(gc_df) <- c("Time","csum")
      gc_df$value_J <- c(gc_df[1,"csum"],diff(gc_df$csum))
      gc_dfm <- melt(gc_df, id.vars = "Time")
      ## add fit to data plot
      gg <- gg + geom_line(data = gc_dfm, aes(x=Time,y=value), col= "red")
      ## save
      ggsave(paste0("output/",strain,"_rep_", replicate,"_",exp_conditions[kk],"_model_fit.pdf"))
    } 
    ## Also save build in plot
    #pdf(paste0("output/",strain,"_rep_", replicate,"_",condition,"_model_fit_builtin.pdf"))
    #plot(gc_fit)
    #dev.off()
    
    # parameters 
    s <- summary(gc_fit)
    
    # ODD? 
    if(max_level >= thresh_wide){odd_width <- 1}
    
    
    ## Required parameters
    param_n_o  <- c(strain,replicate,condition,inocl)
    param_o   <- c(time_max_heat_flow, value_max_heat_flow, 
                   s$mu.spline, s$lambda.spline,s$integral.spline, 
                   odd_peak, odd_width, max_level, odd_shoulder, squared_dist)
    
    return(list(param_n = param_n_o, param = param_o, max_level = max_level, plot_dbl = plot_p, para_dbl = p))
  }
}


### Find peaks function: from https://github.com/stas-g/findPeaks
# if bigger than "m" number of points either side of it
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    # It must be bigger than or equal to all points m to the left and to the right
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

### Modified function - those that are bigger than several to their right fewer on left
find_peaks_right <- function (x, mr = 3, ml = 1){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - ml + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + mr + 1
    w <- ifelse(w < length(x), w, length(x))
    # It must be bigger than or equal to all points m to the left and to the right
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

