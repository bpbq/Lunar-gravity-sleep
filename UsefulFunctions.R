##### Useful function #####

#setwd("/Users/ViridianK/Downloads")

#Open file with gravity data data

##### Function to remove data from a file when more than X consecutive NAs are found #####

#The function removes data from the NAs sequence to the closet "border of the data" for each subject 
#-> Can be improve taking into consideration the quantity of NAs in each subsequence and not only the length

#Parameters: dataFrame  -> is the data frame that you want to modify, it should be ordered such consecutive data points are consecutively ordered
#            ID         -> a vector with Participants IDs (same order as in data)
#            data       -> the column of the data frame where we are going to look for the NAs sequences (same order as in data)
#            n          -> Number of NAs that would cause data removal

removeSubsequentNAs <-function(dataFrame, ID, data, n){
  indexesToKeep=c()
  for (id in unique(ID)){
    index<-which(ID==id)
    NAcount = 0
    #indexCount = 0
    datalength = length(index)
    ilocNAs=c()
    secNAs=c()
    
    for (i in index){
      #indexCount = indexCount + 1
      if(is.na(data[i]) | is.nan(data[i])) {
        NAcount = NAcount+1
      } else{
        if(NAcount >= n){
          ilocNAs<-append(ilocNAs, i-1)
          secNAs= append(secNAs, NAcount)
        }
        NAcount = 0
      }
    }
    ilocNAs<-append(ilocNAs, index[length(index)])
    secNAs= append(secNAs, 0)
    
    current = index[1]
    seqLength = c()
    
    for (i in 1:length(ilocNAs)){
      seqLength = append(seqLength, ((ilocNAs[i]-current)-secNAs[i]))
      current = ilocNAs[i] 
    }
    current = ilocNAs[1] - secNAs[1]
    base = index[1]
    if (length(seqLength) > 1 ){
      if (current != (ilocNAs[which(seqLength==max(seqLength))[1]] - secNAs[1])){
          current = ilocNAs[which(seqLength==max(seqLength))[1]]  - secNAs[which(seqLength==max(seqLength))[1]]
          base = ilocNAs[((which(seqLength==max(seqLength))[1]) - 1)]
      }
    }
    indexesToKeep=c(indexesToKeep, c(base:current))
    

  }
  
  data=dataFrame[indexesToKeep,]
  return (data)
  
}

#### Calculate z-scores (a normalized measure per subject)
#Parameters: data -> data frame
#           variables -> a vector of columnNames

calculate_zscores<-function(data, variables){
  data <- data %>%
    group_by(id) %>%
    mutate_at(variables, zscore) %>%
    ungroup()
  
  return(data)
}


###### Transform onsets and offsets to decimal ######
timedecimal <- function(x) {
  xlist <- strsplit(x,split=":")
  h <- as.numeric(sapply(xlist,"[",1))
  m <- as.numeric(sapply(xlist,"[",2))
  s <- as.numeric(sapply(xlist,"[",3))
  xdec <- h+(m/60)+(s/60/60)
  return(xdec)
}

##### Mimics plyr round_any function #####
round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}


##### Plots Lunar #####

# data -> data to plot
# results -> Estimates of the model
# n -> number of periods included in the model (1 or 2)
# path -> character indicating the path for saving the plots

sinusoid <-function(data, results, n=1, period){
  if (n==1){
    fun = function (x){ results["M"] + (results["A"] * 
                                          sin((2*pi*1/period)* (x - (results["P"]))))}
  }else {
    fun = function (x){ results["M"] + (results["A"] * 
                                          sin((2*pi*1/period)* (x - (results["P"]))) + results["A2"] * 
                                          sin((2*pi*1/(period/2))* (x - (results["P2"])))) 
    }
  }
  return(fun)
}

lunarPlot <- function(data, variable, results, n=1, period){
  
  
  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle), y = variable)) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    stat_function(fun = sinusoid(data, results, n, period),
                  colour = "#D55E00", size = 2) +
    #geom_vline(aes(xintercept = p.ons.zscore.all["P"]+30/4), colour = "#D55E00", size = 1, linetype = 2)+
    stat_summary(fun="mean", colour = "black", pch=16, size=0.8, position = position_dodge(width = .1)) + 
    labs(x = "Night in cycle", y = variable) +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15), 
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text())#+
  
  
  return(plot)
  
}

#### Moving average ####

daily_rollmean <- function(data, var, k, max_gap) {
  # This function calculates the roll mean using by interpolating the NA gaps smaller than
  # max_gap and then applying a moving average with a rectangular window of size k
  #   data: A tibble containing the night and the variable of interest.
  #   var: The variable of interest (for example sleep_time)
  #   k: Size of the rectangular window used in the moving average.
  #   max_gap; Maximum gap in the signal not to be considered as missing data.
  
  head(data)
  # We first complete all the data
  complete_data <- data |>
    #ungroup() |>  # Just in case the data is grouped
    right_join(tibble(nightNumber = seq(min(data$nightNumber), max(data$nightNumber), by = 1)), by = "nightNumber") |>
    arrange(nightNumber) |> # We arrange based on night
    mutate(completed = na.approx({{ var }},  
                                 maxgap = max_gap), # Complete the NA gaps smaller than max_gap size
           "{{ var }}" := rollmean(completed,
                                      k = k,
                                      na.pad = T)) |> # Moving average of k data
    drop_na({{ var }}) |>
    select(-completed)
  
  return(complete_data)
}


####### Moving average: like casiraghi 2020 or more conservative #######

### Expected columns in data: id -> participant id
#                             night -> date of the night of sleep
#                             DayInCycle -> between 0 and 29 or 1 and 30, day of the moon cycle with 0 representing the full moon and 15 the new moon
#                             variable -> the name of the variable under study in data (e.g. Duration)

# type -> 1: fill missing values with linear interpolation and run the moving average. Keep all data points.
#         2: whar daily_rollmean() does

# n -> number of days for the moving average

# k -> max NA window to fill if type == 2

movingAverage <- function(data, type=2, n=7, k = 3){
  
  if(type == 1){
    
    filling=expand.grid(id=unique(data$id), nightNumber=c(min(data$nightNumber):max(data$nightNumber)))
    data<-merge(filling, data, by=c("id", "nightNumber"), all.x=TRUE)  
    data <- data[order(data$id, data$nightNumber),]
    data = data %>%
      group_by (id) %>%
      mutate_at(c("variable", "DayinCycle"), na.approx, na.rm=FALSE)%>%
      ungroup()
    
    data<-na.omit(data)
    
    data = data %>%
      group_by (id) %>%
      mutate( variable = SMA(variable, n)) %>% #picking n days as window
      ungroup() 
  
    data$variable[1:(length(data$variable)-3)]<-data$variable[4:length(data$variable)]
  }else{
    
      data <- data |>
      group_by(id) |>
      nest() |>
      mutate(ma = map(data, ~daily_rollmean(data = ., variable, n, k))) |>
      select(-data) |>
      unnest(c(ma))
    }
  return(data)
}


prepareForMA<-function(data){
  data2<-data
  data2$DayinCycle<-data2$DayinCycle+30
  data<-rbind(data, data2)
  data2<-data[data$DayinCycle>=0 & data$DayinCycle<=2,]
  data2$DayinCycle<-data2$DayinCycle+60
  data<-rbind(data, data2)
  data2<-data[data$DayinCycle>=27 & data$DayinCycle<=29,]
  data2$DayinCycle<-data2$DayinCycle-30
  data<-rbind(data, data2)
  return(data)
}


####### Analysis using nls ######

### Expected columns in data: id -> participant id
#                             night -> date of the night of sleep
#                             DayInCycle -> between 0 and 29 or 1 and 30, day of the moon cycle with 0 representing the full moon and 15 the new moon

# variable -> the name of the variable under study in data (e.g. Duration)

# name -> string used in the plots and for saving them

# ma -> indicates if a moving average is applied or not and which type of ma (0-> not applied, 1/2 -> same ase used in function moving Average)

# n -> number of days for the moving average

# zscores -> TRUE is zscores are calculated FALSE if not (default FALSE)
applyingNLS <- function(data, variable, period = 30, ma = 0, n = 7, k = 3, zscore = FALSE, plot = FALSE, path = getwd()){
  zc=""
  movA=""
  
  
  
  data$night<-as.POSIXct(data$night)
  data$nightNumber <- as.numeric((data$night-min(data$night))/86400)
  
  data<-data.frame(id=data$id,
                   nightNumber = data$nightNumber,
                   DayinCycle = data$DayinCycle,
                   variable=data[,variable])
  
  
  #### Removing outliers
  data<-data[abs(zscore(data$variable))<3,]
  
  if(zscore==TRUE){
    data <- data %>%
      group_by(id) %>%
      mutate(variable = zscore(variable)) %>%
      ungroup()
    zc="zscore"
    }
  
  if(ma != 0){
    data<-movingAverage(data, ma, n, k)
    movA="ma"
  }
  
  data<-na.omit(data)
  
  #### Model with 1 period (period)
  model1 <- nls(variable ~ M+(A * sin((2*pi*1/period)*(DayinCycle - P))), 
                          start = list(A = sd(data$variable), P = (period/2), M=mean(data$variable)),
                          lower = list(A = -Inf, P = 0, M=-Inf),
                          upper = list(A = Inf, P = period, M=Inf), algorithm = "port",
                          data = data, nls.control(maxiter = 20000, warnOnly = TRUE))
  
  results<-as.data.frame(summary(model1)$parameters)
  results$period<-period

  
  # #### Model with 2 periods (period and period/2)
  # model2 <- nls(variable ~ M+(A * sin((2*pi*1/period)*(DayinCycle - P)) +  A2 * sin((2*pi*1/(period/2))*(DayinCycle - P2))), 
  #                          start = list(A = sd(data$variable), P = 10, A2 = sd(data$variable), P2 = 5, M=mean(data$variable)),
  #                          lower = list(A = -Inf, P = 0, P2 = 0, A2=0, M=-Inf),
  #                          upper = list(A = Inf, P = period, P2 = period/2, A2= Inf, M=Inf), algorithm = "port",
  #                          data = data, nls.control(maxiter = 10000, warnOnly = TRUE))
  # 
  # result_model2<-as.data.frame(summary(model2)$parameters)
  # result_model2$period<-paste(period, "and", period/2)
  # 
  # 
  # #### Model with 1 period (period/2)
  # model3 <- nls(variable ~ M+(A * sin((2*pi*1/(period/2))*(DayinCycle - P))), 
  #                          start = list( A = sd(data$variable), P = 7, M=mean(data$variable)),
  #                          lower = list(A = -Inf, P = 0, M=-Inf),
  #                          upper = list(A = Inf, P = period/2, M= Inf), algorithm = "port",
  #                          data = data, nls.control(maxiter = 10000, warnOnly = TRUE))
  # 
  # result_model3<-as.data.frame(summary(model3)$parameters)
  # result_model3$period<-period/2
  # 
  # results<-rbind(result_model1, result_model2, result_model3)
  # 
  # modelComparison<- data.frame(model=(c(period, paste(period, "and", period/2), period/2)),
  #                              stdErr_model= c(summary(model1)$sigma, summary(model2)$sigma, summary(model3)$sigma),
  #                              aic=c(AIC(model1), AIC(model2), AIC(model3)))
  #   
  write.csv(results, paste(path,"_modelResults_",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
  #write.csv(results, paste(path,"_modelComparison_",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
  
  if (plot){
    data2<-data
    data2$DayinCycle<-data2$DayinCycle+30
    data<-rbind(data, data2)
    
    tiff(paste(path,"_",variable,"_",zc,"_",movA,"_",period,".tif", sep=""), width = 2000, height = 1200, res = 300,
         compression = "lzw") # file & size
    print(lunarPlot(data, variable, summary(model1)$parameters[,1], 1, period))
    dev.off() #close
    
    # tiff(paste(path,"_",variable,"_",zc,"_",movA,"_",period,"+",period/2,".tif", sep=""), width = 2000, height = 1200, res = 300,
    #      compression = "lzw") # file & size
    # print(lunarPlot(data, variable, summary(model2)$parameters[,1], 2, period))
    # dev.off() #close
    # 
    # tiff(paste(path,"_",variable,"_",zc,"_",movA,"_",period/2,".tif", sep=""), width = 2000, height = 1200, res = 300,
    #      compression = "lzw") # file & size
    # print(lunarPlot(data, variable, summary(model3)$parameters[,1], 1, period/2))
    # dev.off() #close
  }
  
  
}


##### Cleaning data #####

clean_data <- function(data, variable, period = 30, ma = 0, n = 7, k = 3, zscore = FALSE){
  zc=""
  movA=""
  
  
  
  data<-data.frame(id=data$id,
                   DayinCycle = data$DayinCycle,
                   variable=data[,variable],
                   night = data$night)
  
  #### Removing outliers
  data<-data[abs(zscore(data$variable))<3,]
  
  if(zscore==TRUE){
    data <- data %>%
      group_by(id) %>%
      mutate(variable = zscore(variable)) %>%
      ungroup()
    zc="zscore"
  }
  
  if (ma==1){
  data<-data %>% group_by(id, DayinCycle) %>%
    summarise(variable = mean(variable))
  
  
  data<-prepareForMA(data)
  data$nightNumber<-data$DayinCycle
  }
  
  if(ma==3){
    
    data=data[order(sleep$id, sleep$night),]
    
    data<-data %>%
      group_by(id) %>%
      mutate(night=c(1:length(id)))
    data<-rename(data, nightNumber = night)
  }
  
  if(ma != 0){
    data<-movingAverage(data, ma, n, k)
    movA="ma"
  }
  
  data<-na.omit(data)
  
  return(data)
  }



##### Individual fits using nls

applyingNLS_individuals <- function(data, variable, period = 30, ma = 0, n = 7, k = 3, zscore = FALSE, plot = FALSE, path = getwd()){
  zc=""
  movA=""
  
  
  
  data<-data.frame(id=data$id,
                   DayinCycle = data$DayinCycle,
                   variable=data[,variable])
  
  #### Removing outliers
  data<-data[abs(zscore(data$variable))<3,]
  
  if(zscore==TRUE){
    data <- data %>%
      group_by(id) %>%
      mutate(variable = zscore(variable)) %>%
      ungroup()
    zc="zscore"
  }
  
  data<-data %>% group_by(id, DayinCycle) %>%
    summarise(variable = mean(variable))
  
  
  data<-prepareForMA(data)
  data$nightNumber<-data$DayinCycle
  
  if(ma != 0){
    data<-movingAverage(data, ma, n, k)
    movA="ma"
  }
  
  data<-na.omit(data)
  
  reserve<-data
  
  results_all<-data.frame()
  results_all2<-data.frame()
  modelComparison_all<-data.frame()
  
  for (i in unique(reserve$id)){  
    
      data<-reserve[reserve$id==i,]
      
      #### Model with 1 period (period)
      model1 <- nls(variable ~
                      M+(A * sin((2*pi*1/period)*(DayinCycle - P))), 
                    start = list(A = sd(data$variable), P = (period/2), M=mean(data$variable)),
                    lower = list(A = -Inf, P = 15, M=-Inf),
                    upper = list(A = Inf, P = period, M=Inf), algorithm = "port",
                    data = data, nls.control(maxiter = 20000, warnOnly = TRUE))
      
      results<-as.data.frame(summary(model1)$parameters)
      results$period<-period
      

      #### Model with 2 periods (period and period/2)
      model2 <- nls(variable ~
                    M+(A * sin((2*pi*1/period)*(DayinCycle - P)) +  A2 * sin((2*pi*1/(period/2))*(DayinCycle - P2))),
                    start = list(A = sd(data$variable), P = 7, A2 = sd(data$variable), P2 = 22, M=mean(data$variable)),
                    lower = list(A = -Inf, P = 0, P2 = 0, A2=0, M=-Inf),
                    upper = list(A = Inf, P = period, P2 = period/2, A2= Inf, M=Inf), algorithm = "port",
                    data = data, nls.control(maxiter = 10000, warnOnly = TRUE))

      results2<-as.data.frame(summary(model2)$parameters)
      results2$period<-paste(period, "and", period/2)

      # 
      # 
      # #### Model with 1 period (period/2)
      # model3 <- nls(variable ~
      #                 M+(A * sin((2*pi*1/(period/2))*(DayinCycle - P))), 
      #               start = list( A = sd(data$variable), P = 7, M=mean(data$variable)),
      #               lower = list(A = -Inf, P = 0, M=-Inf),
      #               upper = list(A = Inf, P = period/2, M= Inf), algorithm = "port",
      #               data = data, nls.control(maxiter = 10000, warnOnly = TRUE))
      # 
      # result_model3<-as.data.frame(summary(model3)$parameters)
      # result_model3$period<-period/2
      # 
      # 
      # results<-rbind(result_model1, result_model2, result_model3)
      # results$id<-i
      # 
      # modelComparison<- data.frame(model=(c(period, paste(period, "and", period/2), period/2)),
      #                              stdErr_model= c(summary(model1)$sigma, summary(model2)$sigma, summary(model3)$sigma),
      #                              aic=c(AIC(model1), AIC(model2), AIC(model3)),
      #                              id=i
      #                              )
      #modelComparison_all<-rbind(modelComparison_all, modelComparison)
      
      results_all<-rbind(results_all, results)
      results_all2<-rbind(results_all2, results2)
      
      
      if (plot){
        
        tiff(paste(path,"_",i,"_",variable,"_",zc,"_",movA,"_",period,".tif", sep=""), width = 2000, height = 1200, res = 300,
             compression = "lzw") # file & size
        print(lunarPlot(data, variable, summary(model1)$parameters[,1], 1, period))
        dev.off() #close
        
        tiff(paste(path,"_",i,"_",variable,"_",zc,"_",movA,"_",period,"+",period/2,".tif", sep=""), width = 2000, height = 1200, res = 300,
             compression = "lzw") # file & size
        print(lunarPlot(data, variable, summary(model2)$parameters[,1], 2, period))
        dev.off() #close

        # tiff(paste(path,"_",i,"_",variable,"_",zc,"_",movA,"_",period/2,".tif", sep=""), width = 2000, height = 1200, res = 300,
        #      compression = "lzw") # file & size
        # print(lunarPlot(data, variable, summary(model3)$parameters[,1], 1, (period/2)))
        # dev.off() #close
      }
  }
  
  write.csv(results_all, paste(path,"_modelResults_individual_30days",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
  write.csv(results_all2, paste(path,"_modelResults_individual_30+15days",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
  #write.csv(results, paste(path,"_modelComparison_individual_",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
}

####### Analysis using GLMMcosinor ########

### Expected columns in data: id -> participant id
#                             night -> date of the night of sleep
#                             DayinCycle -> between 0 and 29 or 1 and 30, day of the moon cycle with 0 representing the full moon and 15 the new moon
#                             group -> group if group parameter is set as TRUE
# variable -> the name of the variable under study in data (e.g. Duration)

# name -> string used in the plots and for saving them

# group -> whether we are controlling for a group in the models

# controlVariable -> Variable that is included for control (e.g. typeOfDay -weekend or weekday- or gender) at the mesor level

# ma -> indicates if a moving average is applied or not and which type of ma (0-> not applied, 1/2 -> same ase used in function moving Average)

# n -> number of days for the moving average

# zscores -> TRUE is zscores are calculated FALSE if not (default FALSE)

applyingGLMMCosinor <- function(data, variable, group = FALSE, controlVariable = NaN, components = 1, period = 30, ma = 0, n = 7, k=3, zscore = FALSE, path = getwd()){
    zc=""
    movA=""
    
    data$variable=data[,variable]
    if(!is.nan(controlVariable)){
       data$controlVariable= data[,controlVariable]
    }
    
    #### Removing outliers
    data<-data[abs(zscore(data$variable))<3,]
       
     if(zscore==TRUE){
      data <- data %>%
        group_by(id) %>%
        mutate(variable = zscore(variable)) %>%
        ungroup()
      zc="zscore"
        }
    

    
    data$night<-as.POSIXct(data$night)
    data$nightNumber <- as.numeric((data$night-min(data$night))/86400)
    

    
    if(ma != 0){
      data<-movingAverage(data, ma, n, k)
      movA="ma"
    }
    
    data<-na.omit(data)

    if(group & !is.nan(controlVariable)){
      mixed_mod <- cglmm(
        variable ~ controlVariable + amp_acro(DayinCycle,
                      n_components = 1,
                      period = 30, group = group) + (1 + amp_acro1 | id),
                      data = data
      )
      
  
      
    } else if (group){
      mixed_mod <- cglmm(
        variable ~ amp_acro(DayinCycle,
                      n_components = 1,
                      period = 30, group = group) + (1 + amp_acro1 | id),
                      data = data
      )
      

    } else if (!is.nan(controlVariable)){
       mixed_mod <- cglmm(
        variable ~ controlVariable + amp_acro(DayinCycle,
                     n_components = 1,
                     period = 30) + (1 + amp_acro1 | id),
                     data = data
      )
      
      
    } else{
      mixed_mod <- cglmm(
        variable ~ amp_acro(DayinCycle,
                    n_components = 1,
                    period = 30) + (1 + amp_acro1 | id),
                    data = data
      )
    } 
    
    
    summ<-summary(mixed_mod)
    write.csv(as.data.frame(summary(mixed_mod)$transformed.table), paste(path,"_modelResults_transformedTable_",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
    write.csv(as.data.frame(summary(mixed_mod)$raw.table), paste(path,"_modelResults_rawTable_",variable,"_",zc,"_",movA,"_",period,".csv", sep=""))
    
    tiff(paste(path,"_",variable,"_",zc,"_",movA,"_",period,".tif", sep=""), width = 2000, height = 1200, res = 300,
         compression = "lzw") # file & size
    print(autoplot(mixed_mod))
    dev.off() #close
}
