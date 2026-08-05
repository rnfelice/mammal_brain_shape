timeslices=.1
interval_rates <- rateMap(
  list(jaws_search),
  summary = "interval",
  control = list(res = 200)
)
obj<-as_tibble(interval_rates$intervals)
# get the rates at timeslices from a processed dataframe
rate.at.time.df.bifrost <- function(timeslices, obj, plot=F, relative.rates=c("F","mean","median")){
  # define the timeslices for extracting trait values, and apply it to your tree
  obj<-obj |> mutate(timestart=depth_start-max(depth_end), timeend=depth_end-max(depth_end))
  ts <- c(seq(from=min(c(obj$timestart,obj$timeend)), to=max(c(obj$timestart,obj$timeend)), by=timeslices), max(c(obj$timestart,obj$timeend)))
  
  # for each timepoint reconstruct the trait value as a function of the distance between the parent and child node
  Yt <- list()
  for (i in 1:(length(ts)-1)) {      
    # which edges overlap timeslice i
    curr.edges <- filter(obj, timestart <= ts[[i]] & timeend >= ts[[i]])
    
    # extract the rates for the branches of interest
    curr.rates <- curr.edges$value
    
    # apply the current slice time
    slice <- rep(ts[i],times = length(curr.rates))
    
    # add all the trait values at timeslice i to the list
    Yt[[i]] <- cbind(curr.rates, slice)
  }
  # remove the fractional period just before the tips (combine the penultimate and ultimate windows)
  ts2 <- ts[-(length(ts)-1)]
  
  # now make a data frame of the trait values at given times
  shuffled.ages <- NULL
  for (k in 1:length(ts2)){
    curr.slice <- unlist(Yt[[k]])
    curr.slice.df <- data.frame(curr.slice, "time" = ts2[k])
    shuffled.ages <- rbind(shuffled.ages, curr.slice.df)
  }
  
  # extract just the trait and timeslices
  rate.time <- shuffled.ages[,c("curr.rates","time")]
  
  # reorder the columns so 'time' comes first
  rate.time <- relocate(rate.time, time)
  colnames(rate.time) <- c("time", "rate")
  
  # rescale rates relative to the mean if requested
  if(relative.rates=="F"){rate.time$rate <- rate.time$rate}
  if(relative.rates=="mean"){rate.time$rate <- rate.time$rate/mean(rate.time$rate)}
  if(relative.rates=="median"){rate.time$rate <- rate.time$rate/median(rate.time$rate)}
  
  # plot the values if you're interested, it should look like the horizontal branches of the tree
  if(plot==T){plot(rate.time$rate ~ rate.time$time, xlab="Time", ylab="Rate", pch=16)}
  
  # give up the object
  return(rate.time)
}

extract.stat <- function(rate.time.obj, stat=c("mean", "median", "scale"), range=c("confidence","quantile"),
                         plot=c("average", "corrected", "sideXside", FALSE)){
  
  # create a data frame of the individual time slices
  time = unique(rate.time.obj$time)
  time.mean <- data.frame(time)
  
  # establish the minimum and maximum values at each slice
  if(stat=="mean"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)}
  if(stat=="median"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), median)}
  if(stat=="scale"){
    rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)
    rate.mean[,2] <- ((rate.mean[,2] - min(rate.mean[,2]))/diff(range(rate.mean[,2])) * 99) + 1
  }
  rate.sd <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), sd)
  
  
  # if you want to estimate quantiles on the rates:
  if(range=="quantile"){
    rate.qt5 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.025)
    rate.qt95 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.975)
    rate.mean <- cbind(rate.mean, rate.sd[,2], rate.qt5[,2], rate.qt95[,2])    
  }
  
  # if you prefer to estimate confidence intervals on the rates:
  if(range=="confidence"){
    rate.all <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN = function(x) Rmisc::CI(x, ci = 0.95))
    rate.mean <- data.frame(time = time.mean[,1], rate = rate.mean[,2], 
                            sd = rate.sd[,2], "5%" = rate.all$x[,3], "95%" = rate.all$x[,1])   
  }
  # regardless, rename the columns
  colnames(rate.mean) <- c("time", "rate", "sd", "5%", "95%")
  
  # establish the number of species living in each time slice
  spp.slice <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), length)
  
  # add the richness information to the data frame
  rate.mean$richness <- spp.slice[,2]
  
  # correct rate by the number of species living at that time
  rate.mean$rate.rich <- rate.mean$rate/rate.mean$richness
  
  # get the name of the current variable
  var.name <- colnames(rate.time.obj)[2]
  
  # plot the results if you'd like
  if(plot=="average"){plot(rate.mean$rate ~ rate.mean$time, 
                           xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="corrected"){plot(rate.mean$rate.rich ~ rate.mean$time, 
                             xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="sideXside"){layout(matrix(1:2,ncol=2))
    plot(rate.mean$rate ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    plot(rate.mean$rate.rich ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    layout(matrix(1))}
  if(plot==F){}
  
  return(rate.mean)
  
}

