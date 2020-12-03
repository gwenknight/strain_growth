#### Function to 
# (1) fit growth curve to data 
# (2) return parameters describing the curve
# (3) return any "odd" behaviour flags
# (4) attempt to fit double curves, flag this and return parameters if can

# Data must have columns: 
#"Time": time of reading
#"rep":  replicate label (e.g. 6.1 etc) 
#"value_J": joules reading level
#"csum": cumulative sum of value_J   
#"strain": strain name   
#"inoc": inoculum level (10^x) 
#"drytime": how long dried for


######************************************************************************************************************
# FUNCTION
# Inputs
# First four identify the data combination
# strain: strain name
# replicate:  experiment name
# condition: number of hours dried
# inocl: inoculum (10^x)

# data: the data from the experiment
# plot: if plot = 1 then an output of this strain's data is given

# For odd behaviours
# thresh_wide: width of the curve peak

### Bug fixing
# strain <- u[jj];
# replicate <- r[ii]
# condition <- drying_times[kk]
# inocl <- q[ll]
# data <- ddm
# fit_growth_curve(strain, replicate, condition, inocl, data)

fit_growth_curve <- function(strain, replicate, condition, inocl, data, 
                             plot = 0, thresh_wide = 80){
  
  ## which rows of the data are this strain and replicate?
  wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
  wj <- intersect(wi, which(data$drytime == condition))
  w <- intersect(wj, which(data$inoc == as.numeric(inocl)))
  
  ## is this strain replicate ODD? Set all ODD indicators to zero initially
  odd_peak <- 0 # 0 = one peak
  odd_width <- 0 # 0 = not a wide peak
  odd_shoulder <- 0 # 0 = no shoulder
  max_level <- 0 # height of peak
  odd_double <- 0 # double curve fits this data
  
  if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
    data1 <- data[w,]  # just get the data for this experiment (strain, time, value, drying time)
    
    ## Growth Curve # (see fig 3 of vv33i07.pdf)
    ## This gives lag time and exponential growth rate cumulative
    gc_fit <- gcFitSpline(data1$Time, data1$value_J)
    # parameters from this fit
    s <- summary(gc_fit)
    
    ## What is the maximum heat flow and when?
    wmax <- which.max( data1$value_J[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
    time_max_heat_flow <- as.numeric(data1[wmax,"Time"])
    value_max_heat_flow <- as.numeric(data1[wmax,"value_J"])
    
    ## ODD 
    # (1) Is the peak broad? 
    interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
    interval_value <- as.numeric(unlist(data1[c( which(round(data1$Time,4) == round(interval_peak[1],4)), which(round(data1$Time,4) == round(interval_peak[2],4))),"value_J"]))
    
    max_level <- 0;
    max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest
    
    # ODD? 
    if(max_level >= thresh_wide){odd_width <- 1}
    
    # (2) Are there multiple peaks? 
    # GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
    peaks_index = find_peaks(data1$value_J, m = 5)
    if(length(peaks_index) > 1){
      # REMOVE - early ones
      we<-which(data1[peaks_index,"Time"]>3) 
      # REMOVE - late ones
      w<-intersect(we,which(data1[peaks_index,"Time"]< 0.90*max(data1$Time))) # remove > 95% of time
      wl<-intersect(w,which(data1[peaks_index,"Time"] > 0.6*max(data1$Time))) # which in the odd 70-95% of the time range
      if(length(wl)>0){ # if a late peak check it is big 
        for(gg in 1:length(wl)){
          if(data1[peaks_index[wl[gg]],"value_J"] < 0.45*max(data1[peaks_index,"value_J"])){ # if not bigger than 40% of peak
            w <-setdiff(w,wl[gg]) }}}   # then remove
      peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
      # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
      o <- order(data1[peaks_index,"value_J"], decreasing = "TRUE")
      peaks_index <- peaks_index[o]
    }
    
    # When are the peaks? 
    time_peaks <- as.numeric(unlist(data1[peaks_index,"Time"]))
    time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
    # If multiple far apart then issue: double peaks
    keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 
    
    if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
    if(any(data1[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low
    
    # If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
    close_peaks <- which(abs(time_peaks_diff) <= 5)
    close_peaks_i <- 1
    if(length(close_peaks)>1){
      for(i in 2:length(close_peaks)){
        ifelse(data1[peaks_index[close_peaks[i]],"value_J"]/data1[peaks_index[close_peaks[1]],"value_J"]> 0.9,
               close_peaks_i <- c(close_peaks_i,i),"")
      }}
    
    if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks
    
    # Only keep those peaks that are far apart
    time_peaks <- time_peaks[c(1,keep_time_far_apart)]
    peaks_index <- peaks_index[c(1,keep_time_far_apart)]
    
    
    ###### SHOULDER
    # (3) Is there a shoulder? 
    # GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
    # First, draw straight line from peak
    # endline = bottom of line near time zero
    time_endline <- -40 #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
    value_endline <- 0
    ## Start of line is where? (top point)
    time_startline <- time_peaks[1] # highest peak
    value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
    
    # Draw straight line, assuming peak time = time 0. 
    grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
    # Only want exponential growth line
    exp_start <- which.min(abs(data1$Time-s$lambda.spline+1))
    times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
    
    
    # Predicted straight line
    pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
    
    ### Visualise where peaks are if needed
    # plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(-40,25))
    #  points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
    ## and can plot line against this if needed too 
    #  lines(times_line, pred_points_fit, col= "blue")
    
    ### Run through lines. 
    odd_shoulder<- 0
    shoulder_point <- time_peaks[1] # shoulder before this
    ### Visualise where peaks are if needed
    # plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(0,25))
    # points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
    
    for(i in -40:5){
      time_endline <- i #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
      value_endline <- 0
      #where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
      #value_endline <- data1[where.end, "value_J"] # Find value at this time point
      ## Start of line is where? (top point)
      time_startline <- time_peaks[1] # highest peak
      value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
      
      # Draw straight line, assuming peak time = time 0. 
      grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
      times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
      
      # Predicted straight line
      pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
      
      
      ## and can plot line against this if needed too 
      # lines(times_line, pred_points_fit, col= "blue")
      
      #### How far is the predicted straight line from the data? 
      #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
      dist <- c(pred_points_fit - data1[peaks_index[1]:exp_start,"value_J"])
      
      # Crossing points
      if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
        if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
          if(max(dist) > 0.001){
            #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
            fpsq <- find_peaks(as.numeric(unlist(dist)), m = 3) # Are there more than 1 peaks?
            if(length(fpsq) != 0){ # of no preaks then not a clear shoulder
              shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
              ws <- which(round(data1$Time,5) == round(shoulder_point1,5)); shoulder_point_v1 <- data1[ws,"value_J"]
              height <- shoulder_point_v1 / data1[peaks_index[1],"value_J"]
              if(height > 0.5 && height < 0.96){ # if the shoulder is greater than half way up but not too close
                shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
                odd_shoulder <- 1;# print(c("new",i)) # Then odd shoulder  
              }
            }
          }
        }
      }
      
    }
    if(odd_shoulder == 0){shoulder_point <- 0; shoulder_point_v <- 0 # no shoulder_point if no shoulder!
    }else{ws <- which(round(data1$Time,5) == round(shoulder_point,5)); shoulder_point_v <- data1[ws,"value_J"]}
    
    #### If no shoulder but multiple peaks, want to grab time of first peak
    if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
      shoulder_point1 = min(time_peaks)
      shoulder_point_v1 = data1[which(round(data1$Time,5) == round(shoulder_point1,5)),"value_J"]
      if((shoulder_point_v1 / data1[peaks_index[1],"value_J"]) > 0.5){ # shoulder needs to be high still! 
        shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
      }
    }
    
    ##### Double curves? 
    x <- data1$Time
    y <- data1$value_J
    
    # Guesses for the parameters in the curve
    startl=list(a=unlist(data1[peaks_index,"value_J"])[1]*1.5, 
                b=diff(interval_peak)/3, 
                c=unlist(data1[peaks_index,"Time"])[1],
                d=unlist(data1[peaks_index,"value_J"])[1]*1.5, 
                e=diff(interval_peak), 
                f=unlist(data1[peaks_index,"Time"])[1])
    # Set indicator/output parameters to zero
    fit1 <- 0
    plot_p <- 0
    p <- 0
    # Try to fit two normal curves using non-linear lear squares methods 
    try(fit1 <- nls(y~(a/b)*exp(-(x-c)^2/(2*b^2))+(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl), silent = TRUE)
    
    # If manage to fit
    if(length(fit1) > 1){
      # print that can 
      print(paste("Double curve fit", strain," ", replicate, condition, inocl, sep = " "))
      odd_double <- 1
      # Save the predictions of the fit 
      pred_p <- predict(fit1)
      # Save the coefficients of the fit
      p <- coef(fit1)
      # Pull out the actual curves using the parameters in p
      g1p <- (p["a"]/p["b"])*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
      g2p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
      # Store the curves
      plot_p <- as.data.frame(cbind(x,y,pred_p, g1p, g2p))
      colnames(plot_p) <- c("time","value_J","fit","normal_curve1","normal_curve2")
      # Store the parameters
      plot_p$a <- p["a"]; plot_p$b <- p["b"]; plot_p$c <- p["c"]
      plot_p$d <- p["d"]; plot_p$e <- p["e"]; plot_p$f <- p["f"]
    }
    
    ## Plot data1, cumulative and fit
    # if functions takes in a command to plot
    if(plot == 1){
      #### FIT TO Csum? 
      # datam <- melt(data1[,c("Time","value_J","csum")], id.vars = "Time")
      # gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() + facet_wrap(~strain, scales = "free") 
      # ## model fit:
      # gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
      # colnames(gc_df) <- c("Time","csum")
      # gc_df$value_J <- c(gc_df[1,"csum"],diff(gc_df$csum))
      # gc_dfm <- melt(gc_df, id.vars = "Time")
      # ## add fit to data plot
      # gg <- gg + geom_line(data = gc_dfm, aes(x=Time,y=value), col= "red")
      # 
      
      ### Fit to value_J : could add cumulative plot in here but not done
      datam <- reshape2::melt(data1[,c("Time","value_J")], id.vars = "Time")
      gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() #+ facet_wrap(~strain, scales = "free") 
      ## model fit:
      gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
      colnames(gc_df) <- c("Time","value")
      #gc_df$value_J <- c(gc_df[1,"csum"],diff(gc_df$csum)) # CHANGE TO CUMULATIVE? 
      gc_dfm <- reshape2::melt(gc_df, id.vars = "Time")
      ## add fit to data plot
      gg <- gg + geom_line(data = gc_df, aes(x=Time,y=value), col= "red")
      
      ## save
      ggsave(paste0("output/",strain,"_rep_", replicate,"_",condition,"_model_fit.pdf"))
    } 
    
    ## Build vectors of required parameters to output
    param_n_o  <- c(strain,replicate,condition,inocl)
    param_o   <- c(time_max_heat_flow, value_max_heat_flow, 
                   s$mu.spline, s$lambda.spline,s$integral.spline, 
                   odd_peak, odd_width, max_level, odd_shoulder, odd_double, shoulder_point, shoulder_point_v)
    
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