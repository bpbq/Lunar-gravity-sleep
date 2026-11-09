# Load necessary libraries for data manipulation, time series analysis, and modeling
library(dplyr)
library(TTR) #MA function
library (lubridate) #Date/Time format
library(data.table) #for setDT function
library (tidyr)
library (cgwtools)
library(suncalc)
library(tidyverse) # to work with dataframes
library(mosaic) # to calculate zscores
library(nlstools) # to analyze NLS fits
library(data.table)
library(zoo)
library(stats)
library(AICcmodavg)
library(modelsummary)
library(GLMMcosinor) # Main package for Cosinor analysis
library(tibble)
library(here)

# Clear the R environment
rm(list=ls())

##### Functions #####

# Formats hours-from-midnight as HH:MM clock labels (wraps mod 24 in display only)
hours_to_clock <- function(x) {
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x)
  h   <- x[ok] %% 24
  hh  <- floor(h)
  mm  <- round((h - hh) * 60)
  hh  <- (hh + mm %/% 60) %% 24
  mm  <- mm %% 60
  out[ok] <- sprintf("%02d:%02d", as.integer(hh), as.integer(mm))
  out
}

# Function to normalize/transform the phase of the rhythm
transformPhase <- function(x, period){
  
  x = (((x*180/pi)*period)/360)-(period/4) 
  
  return(x)
}

# Function to generate the fitted sinusoidal curve values based on model results
# Handles both single component (e.g., just 30-day) and dual component (30+15 day) models
sinusoid <-function(data, results, n=1, period){
  if (n==1){
    fun = function (x){ results[,"M"] + (results[,"A"] * sin((2*pi*1/period)* (x - (results[,"P"]))))}
  }else {
    fun = function (x){ results[,"M"] + (results[,"A"] * sin((2*pi*1/period)* (x - (results[,"P"]))) + results[,"A2"] * sin((2*pi*1/(period/2))* (x - (results[,"P2"])))) 
    }
  }
  return(fun)
}

# Function to generate a specific ggplot for lunar cycles
# Plots the raw data points and the fitted curve, marking Full/New moon phases
lunarPlot <- function(data, variable, results, n=1, period,
                      axis_text_size = NULL, axis_title_size = NULL){


  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle), y = variable)) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    stat_function(fun = sinusoid(data, results, n, period),
                  colour = "#D55E00", size = 2) +
    #geom_vline(aes(xintercept = p.ons.zscore.all["P"]+30/4), colour = "#D55E00", size = 1, linetype = 2)+
    stat_summary(fun="mean", colour = "black", pch=16, size=0.8, position = position_dodge(width = .1)) + 
    labs(x = "Night in cycle") +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15),
                       labels = c("Full", "New","Full", "New", "Full"))+
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variable)) sub("hours from midnight", "clock time", variable) else variable,
      labels = if (grepl("hours from midnight", variable)) hours_to_clock else waiver())+
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
          axis.title.x = element_text(size = axis_title_size),
          axis.title.y = element_text(size = axis_title_size),
          axis.text    = element_text(size = axis_text_size))#+
  
  
  return(plot)

}

# Function: makeIndivFitPlots
# Saves one smoothed, double-plotted 30+15 cosinor fit per individual (Fig-1 style)
# into 'outdir'. 'results_30_15' must be the 30+15 results table AFTER
# rownames_to_column (5 rows per id: (Intercept), amp1, amp2, acr1, acr2; the
# estimate is column 2). Skips individuals with an incomplete fit.
makeIndivFitPlots <- function(reserve, results_30_15, ylabel, pop_label, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  n_ok <- 0
  for (idv in unique(results_30_15$id)) {
    rr  <- results_30_15[results_30_15$id == idv, ]
    dat <- reserve[reserve$id == idv, ]
    if (nrow(rr) < 5 || nrow(dat) == 0) next
    rp <- data.frame(M  = rr[1, 2], A = rr[2, 2], A2 = rr[3, 2],
                     P  = transformPhase(rr[4, 2], 30),
                     P2 = transformPhase(rr[5, 2], 15))
    if (any(is.na(unlist(rp)))) next
    safe <- gsub("[^A-Za-z0-9]+", "_", idv)
    tiff(file.path(outdir, paste0(pop_label, "_", safe, "_30+15.tif")),
         width = 1200, height = 900, res = 300, compression = "lzw")
    print(lunarPlot(dat, ylabel, rp, n = 2, 30) +
            ggtitle(paste0(pop_label, ": ", idv)) +
            theme(plot.title = element_text(hjust = 0.5, size = 8)))
    dev.off()
    n_ok <- n_ok + 1
  }
  cat(paste0("  ", pop_label, ": ", n_ok, " individual 30+15 fit plots saved in ", outdir, "/\n"))
}

# Function: saveIndivFitSVG
# Big-font VECTOR (SVG) export of ONE chosen individual, for the Fig. 1 example
# panels. Same logical size as the .tif (1200x900 px @300 dpi -> 4 x 3 in): the
# plot is identical except for larger axis text/titles, NO id title on top, extra
# margin (so descenders like the "g" in "Night" are not clipped), and vector.
# 'results_30_15' must be in the same format makeIndivFitPlots expects.
saveIndivFitSVG <- function(reserve, results_30_15, ylabel, pop_label, idv, outfile,
                            axis_text_size = 8.8 * 1.5, axis_title_size = 11 * 1.3) {
  rr  <- results_30_15[results_30_15$id == idv, ]
  dat <- reserve[reserve$id == idv, ]
  if (nrow(rr) < 5 || nrow(dat) == 0) {
    cat("  saveIndivFitSVG: no complete data for '", idv, "' -- skipped\n", sep = ""); return(invisible())
  }
  rp <- data.frame(M  = rr[1, 2], A = rr[2, 2], A2 = rr[3, 2],
                   P  = transformPhase(rr[4, 2], 30),
                   P2 = transformPhase(rr[5, 2], 15))
  svg(outfile, width = 1200/300, height = 900/300)   # 4 x 3 in, same as the .tif
  print(lunarPlot(dat, ylabel, rp, n = 2, 30,
                  axis_text_size = axis_text_size, axis_title_size = axis_title_size) +
          theme(plot.margin = margin(6, 10, 10, 12)))   # no title; room for axis descenders
  dev.off()
  cat("  saved big-font vector (no title): ", outfile, "\n", sep = "")
}

# Function: saveIndivFitSVG_clock
# Same as saveIndivFitSVG but RE-FITS the 30+15 model on the RAW (non z-scored)
# data of one individual, so the y-axis is shown in CLOCK TIME instead of z-score.
# z-scoring is a linear rescaling, so the fitted-curve SHAPE is identical to the
# z-scored fit; only the y-axis units change. Used for the Fig. 1 Toba example.
saveIndivFitSVG_clock <- function(sleep_df, var, idv, ylabel, outfile,
                                  axis_text_size = 8.8 * 1.5, axis_title_size = 11 * 1.3) {
  reserve_raw <- na.omit(clean_data(sleep_df, var, period = 30, n = 7, k = 3, zscore = FALSE))
  dat <- reserve_raw[reserve_raw$id == idv, ]
  if (nrow(dat) < 30) {
    cat("  saveIndivFitSVG_clock: <30 nights for '", idv, "' -- skipped\n", sep = ""); return(invisible())
  }
  m2 <- cglmm(variable ~ amp_acro(DayinCycle, n_components = 2, period = c(30, 15)), data = dat[1:30, ])
  rr <- rownames_to_column(summary(m2)$transformed.table)
  rp <- data.frame(M  = rr[1, 2], A = rr[2, 2], A2 = rr[3, 2],
                   P  = transformPhase(rr[4, 2], 30),
                   P2 = transformPhase(rr[5, 2], 15))
  svg(outfile, width = 1200/300, height = 900/300)
  print(lunarPlot(dat, ylabel, rp, n = 2, 30,
                  axis_text_size = axis_text_size, axis_title_size = axis_title_size) +
          theme(plot.margin = margin(6, 10, 10, 12)))
  dev.off()
  cat("  saved big-font vector clock-time (no title): ", outfile, "\n", sep = "")
}

# Helper function to pad data at the beginning and end
# This is necessary for circular data (cycles) so the Moving Average calculates correctly at the boundaries
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

# Main Data Cleaning Function
# 1. Filters variable of interest
# 2. Removes outliers (>3 z-score)
# 3. Optional normalization (z-score)
# 4. Interpolates missing nights
# 5. Applies a simple moving average (SMA)
clean_data <- function(data, variable, period = 30, n = 7, k = 3, zscore = FALSE){
  zc=""
  movA=""
  
  
  data<-data.frame(id=data$ID_year,
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
  
  #Average by day in cycle by id
  data<-data %>% group_by(id, DayinCycle) %>%
    summarise(variable = mean(variable))
  
  #Duplicating data previous to moving average
  data<-prepareForMA(data)
  data$nightNumber<-data$DayinCycle
  
  #Filling missing data using linear interpolation
  filling=expand.grid(id=unique(data$id), nightNumber=c(min(data$nightNumber):max(data$nightNumber)))
  data<-merge(filling, data, by=c("id", "nightNumber"), all.x=TRUE)  
  data <- data[order(data$id, data$nightNumber),]
  data = data %>%
    group_by (id) %>%
    mutate_at(c("variable", "DayinCycle"), na.approx, na.rm=FALSE)%>%
    ungroup()
  
  data<-na.omit(data)
  
  #Moving average
  data = data %>%
    group_by (id) %>%
    mutate( variable = SMA(variable, n)) %>% #picking n days as window
    ungroup() 
  
  data$variable[1:(length(data$variable)-3)]<-data$variable[4:length(data$variable)]
  
  
  data<-na.omit(data)
  
  return(data)
}





#### SEATTLE ANALYSIS ####

setwd(here())

sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)

##### Individual analysis with GLMMcosinor #######

#### Variable: Start Time (Sleep Onset) #####

# Clean and smooth data (zscore = FALSE for Seattle)
reserve<-clean_data(sleep, "Start_Time", period = 30, n = 7, k = 3, zscore = FALSE)
reserve<-na.omit(reserve)

# Initialize containers for results
results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")


# Loop through every individual participant
for (i in unique(sleep$ID_year)){  
  data<-reserve[reserve$id==i,]
  
  # --- Model 1: Single component, 30-day period (Lunar) ---
  model1<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 30),
    data = data[1:30,])
  
  # Store results and AIC
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # Save Plot for Model 1
  # tiff(paste("Figures/individual_plots/Seattle",i,"_StartTime_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
  # --- Model 2: Dual component, 30-day + 15-day period (Lunar + Semi-lunar) ---
  model2<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 2,
                        period = c(30,15)),
    data = data[1:30,])
  
  # Store results and AIC
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # Save Plot for Model 2
  # tiff(paste("Figures/individual_plots/Seattle",i,"_StartTime_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
  
  # --- Model 3: Single component, 15-day period (Semi-lunar) ---
  model3<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 15),
    data = data[1:30,])
  
  # Store results and AIC
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # Save Plot for Model 3
  # tiff(paste("Figures/individual_plots/Seattle",i,"_StartTime_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
}

# Formatting results
results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3)  

# Save all Start Time results to CSV
write.csv(na.omit(results_all1), "results/seattle/individuals/modelResults_individual_30days_StartTime.csv")
write.csv(na.omit(results_all2),"results/seattle/individuals/modelResults_individual_30+15days_StartTime.csv")
write.csv(na.omit(results_all3), "results/seattle/individuals/modelResults_individual_15days_StartTime.csv")
write.csv(na.omit(modelComparison_all), "results/seattle/individuals/ModelComparison_individual_StartTime.csv")

# Generate a manual plot for one specific ID (96866...) using the custom lunarPlot function
results<-results_all2[results_all2$id=="96866_0000000000.cwa_SS2024",]
results<-data.frame("M"=results[1,2], "A"=results[2,2], "A2"=results[3,2], "P"= transformPhase(results[4,2],30), "P2"=transformPhase(results[5,2],15))
data <- reserve[reserve$id=="96866_0000000000.cwa_SS2024",]

tiff(paste("Figures/individual_plots/Seattle_ModelIndividual_30+15.tif", sep=""), width = 1200, height = 900, res = 300,
     compression = "lzw") # file & size
lunarPlot(data, "Sleep onset (hours from midnight)", results, n=2, 30)
dev.off()  



# Smoothed, double-plotted 30+15 fit for every Seattle individual (Fig-1 style)
makeIndivFitPlots(reserve, results_all2, "Sleep onset (hours from midnight)", "Seattle", "Figures/individuals_30_15")

# Fig. 1 Seattle example individual -> big-font VECTOR (axis text/titles ~1.8x)
saveIndivFitSVG(reserve, results_all2, "Sleep onset (hours from midnight)", "Seattle",
                "100737_0000000000.cwa_SS2024",
                "Figures/individuals_30_15/Fig1_Seattle_100737_SS2024_30+15_bigfont.svg")


#### Variable: End Time (Sleep Offset) #####
# Clean data again for End_Time
reserve<-clean_data(sleep, "End_Time", period = 30, n = 7, k = 3, zscore = FALSE)
reserve<-na.omit(reserve)

# Reset result containers
results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")

# Loop through individuals for End Time (same 3 models: 30, 30+15, 15)
for (i in unique(sleep$ID_year)){  
  data<-reserve[reserve$id==i,]
  model1<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 30),
    data = data[1:30,])
  
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_EndTime_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
  model2<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 2,
                        period = c(30,15)),
    data = data[1:30,])
  
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_EndTime_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
  
  model3<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 15),
    data = data[1:30,])
  
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_EndTime_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
}

results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3)  

# Save End Time results
write.csv(na.omit(results_all1), "results/seattle/individuals/modelResults_individual_30days_EndTime.csv")
write.csv(na.omit(results_all2),"results/seattle/individuals/modelResults_individual_30+15days_EndTime.csv")
write.csv(na.omit(results_all3), "results/seattle/individuals/modelResults_individual_15days_EndTime.csv")
write.csv(na.omit(modelComparison_all), "results/seattle/individuals/ModelComparison_individual_EndTime.csv")




#### Variable: Sleep Duration #####
# Clean data for Duration
reserve<-clean_data(sleep, "Duration", period = 30, n = 7, k = 3, zscore = FALSE)
reserve<-na.omit(reserve)

results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")

# Loop through individuals for Duration (same 3 models)
for (i in unique(sleep$ID_year)){  
  data<-reserve[reserve$id==i,]
  model1<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 30),
    data = data[1:30,])
  
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_Duration_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
  model2<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 2,
                        period = c(30,15)),
    data = data[1:30,])
  
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_Duration_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
  
  model3<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 15),
    data = data[1:30,])
  
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Seattle",i,"_Duration_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
}

results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3)

# Save Duration results
write.csv(na.omit(results_all1), "results/seattle/individuals/modelResults_individual_30days_duration.csv")
write.csv(na.omit(results_all2),"results/seattle/individuals/modelResults_individual_30+15days_duration.csv")
write.csv(na.omit(results_all3), "results/seattle/individuals/modelResults_individual_15days_duration.csv")
write.csv(na.omit(modelComparison_all), "results/seattle/individuals/ModelComparison_individual_duration.csv")





##### TOBAS ANALYSIS #####

sleep <- read.csv ("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)


##### Individual analysis #######

#### Variable: Start Time (Sleep Onset) #####

# Clean Data: Note that zscore = TRUE is used here for Tobas dataset
reserve<-clean_data(sleep, "Start_Time", period = 30, n = 7, k = 3, zscore = TRUE)
reserve<-na.omit(reserve)


results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")

# Loop through individuals for Start Time (Tobas)
for (i in unique(sleep$ID_year)){  
  data<-reserve[reserve$id==i,]
  model1<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 30),
    data = data[1:30,])
  
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_StartTime_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
  model2<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 2,
                        period = c(30,15)),
    data = data[1:30,])
  
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_StartTime_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
  
  model3<-cglmm(
    variable ~ amp_acro(DayinCycle,
                        n_components = 1,
                        period = 15),
    data = data[1:30,])
  
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_StartTime_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep onset (hours from midnight)"))
  # dev.off()
  
}

results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3)  


write.csv(na.omit(results_all1), "results/Tobas/individuals/modelResults_individual_30days_StartTime.csv")
write.csv(na.omit(results_all2),"results/Tobas/individuals/modelResults_individual_30+15days_StartTime.csv")
write.csv(na.omit(results_all3), "results/Tobas/individuals/modelResults_individual_15days_StartTime.csv")
write.csv(na.omit(modelComparison_all), "results/Tobas/individuals/ModelComparison_individual_startTime.csv")

# Smoothed, double-plotted 30+15 fit for every Tobas individual (onset is z-scored here)
makeIndivFitPlots(reserve, results_all2, "Sleep onset (z-score)", "Tobas", "Figures/individuals_30_15")

# Fig. 1 Toba example individual -> big-font VECTOR, y-axis in CLOCK TIME
# (re-fit of the SAME 30+15 model on the raw, non z-scored onset -> identical curve
#  shape as the z-scored fit, but readable clock-time units like the Seattle panel)
saveIndivFitSVG_clock(sleep, "Start_Time", "80_Fall-Winter2023",
                      "Sleep onset (hours from midnight)",
                      "Figures/individuals_30_15/Fig1_Tobas_80_FallWinter2023_30+15_bigfont.svg")




#### End Time (Sleep offset) #####
#reserve<-clean_data(sleep, "End_Time", period = 30, n = 7, k = 3, zscore = TRUE)
reserve<-sleep
reserve<-na.omit(reserve)


results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")


for (i in unique(reserve$ID_year)){  
  data<-reserve[reserve$ID_year==i,]
  model1<-cglmm(
    End_Time ~  amp_acro(DayinCycle,
                         n_components = 1,
                         period = 30),
    data = data)
  
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_EndTime_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
  model2<-cglmm(
    End_Time ~  amp_acro(DayinCycle,
                         n_components = 2,
                         period = c(30,15)),
    data = data)
  
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_EndTime_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
  
  model3<-cglmm(
    End_Time ~  amp_acro(DayinCycle,
                         n_components = 1,
                         period = 15),
    data = data)
  
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_EndTime_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep offset (hours from midnight)"))
  # dev.off()
  
}

results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3) 


write.csv(na.omit(results_all1), "results/Tobas/individuals/modelResults_individual_30days_EndTime.csv")
write.csv(na.omit(results_all2),"results/Tobas/individuals/modelResults_individual_30+15days_EndTime.csv")
write.csv(na.omit(results_all3), "results/Tobas/individuals/modelResults_individual_15days_EndTime.csv")
write.csv(na.omit(modelComparison_all), "results/Tobas/individuals/ModelComparison_individual_EndTime.csv")




#### Duration #####
#reserve<-clean_data(sleep, "Duration", period = 30, n = 7, k = 3, zscore = TRUE)
reserve<-sleep
reserve<-na.omit(reserve)


results_all1<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2<-as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3<-as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all<-data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all)<-c("AIC", "sigma","period","id")


for (i in unique(reserve$ID_year)){  
  data<-reserve[reserve$ID_year==i,]
  model1<-cglmm(
    Duration ~  amp_acro(DayinCycle,
                         n_components = 1,
                         period = 30),
    data = data)
  
  results<-summary(model1)$transformed.table
  results$id<-i
  results$period<-"30"
  results_all1<-rbind(results_all1, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model1$fit), "sigma"=sigma(model1$fit),"period"="30", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_Duration_30.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model1, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
  model2<-cglmm(
    Duration ~  amp_acro(DayinCycle,
                         n_components = 2,
                         period = c(30,15)),
    data = data)
  
  results<-summary(model2)$transformed.table
  results$id<-i
  results$period<-"30+15"
  results_all2<-rbind(results_all2, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model2$fit), "sigma"=sigma(model2$fit),"period"="30+15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_Duration_30+15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model2, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
  
  model3<-cglmm(
    Duration ~  amp_acro(DayinCycle,
                         n_components = 1,
                         period = 15),
    data = data)
  
  results<-summary(model3)$transformed.table
  results$id<-i
  results$period<-"15"
  results_all3<-rbind(results_all3, results)
  modelComparison_all<-rbind(modelComparison_all, data.frame("AIC"=AIC(model3$fit), "sigma"=sigma(model3$fit),"period"="15", "id"=i))
  
  # tiff(paste("Figures/individual_plots/Tobas",i,"_Duration_15.tif", sep=""), width = 2000, height = 1200, res = 300,
  #      compression = "lzw") # file & size
  # print(autoplot(model3, superimpose.data = TRUE)+labs(x = "Night in cycle", y = "Sleep duration (minutes)"))
  # dev.off()
  
}

results_all1<-rownames_to_column(results_all1)  
results_all2<-rownames_to_column(results_all2)  
results_all3<-rownames_to_column(results_all3)



write.csv(na.omit(results_all1), "results/Tobas/individuals/modelResults_individual_30days_duration.csv")
write.csv(na.omit(results_all2),"results/Tobas/individuals/modelResults_individual_30+15days_duration.csv")
write.csv(na.omit(results_all3), "results/Tobas/individuals/modelResults_individual_15days_duration.csv")
write.csv(na.omit(modelComparison_all), "results/Tobas/individuals/ModelComparison_individual_duration.csv")




##### TITI MONKEYS ANALYSIS #####

#### Variable: Activity Onset #####
# Per-individual cosinor fits (30, 30+15 and 15 day) for the activity-onset of each
# titi, used by the gravity-referenced circular plots. Titis are sparser than the
# human cohorts, so each fit is wrapped in tryCatch and individuals with fewer than
# 30 cleaned points are skipped. Same cleaning/smoothing as the other populations.

titi <- read.csv("data/Data_S3_all_onsetData_titi.csv", header = TRUE)
titi$ID_year <- titi$id2                 # individual id expected by clean_data
titi$night   <- as.Date(titi$night)

reserve <- clean_data(titi, "time", period = 30, n = 7, k = 3, zscore = FALSE)
reserve <- na.omit(reserve)

results_all1 <- as.data.frame(matrix(nrow = 0, ncol = 5))
results_all2 <- as.data.frame(matrix(nrow = 0, ncol = 5))
results_all3 <- as.data.frame(matrix(nrow = 0, ncol = 5))
modelComparison_all <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(modelComparison_all) <- c("AIC", "sigma", "period", "id")

for (i in unique(titi$ID_year)) {
  data <- reserve[reserve$id == i, ]
  if (nrow(data) < 30) next

  model1 <- tryCatch(cglmm(variable ~ amp_acro(DayinCycle, n_components = 1, period = 30),
                           data = data[1:30, ]), error = function(e) NULL)
  model2 <- tryCatch(cglmm(variable ~ amp_acro(DayinCycle, n_components = 2, period = c(30, 15)),
                           data = data[1:30, ]), error = function(e) NULL)
  model3 <- tryCatch(cglmm(variable ~ amp_acro(DayinCycle, n_components = 1, period = 15),
                           data = data[1:30, ]), error = function(e) NULL)

  if (!is.null(model1)) {
    results <- summary(model1)$transformed.table; results$id <- i; results$period <- "30"
    results_all1 <- rbind(results_all1, results)
    modelComparison_all <- rbind(modelComparison_all, data.frame("AIC" = AIC(model1$fit), "sigma" = sigma(model1$fit), "period" = "30", "id" = i))
  }
  if (!is.null(model2)) {
    results <- summary(model2)$transformed.table; results$id <- i; results$period <- "30+15"
    results_all2 <- rbind(results_all2, results)
    modelComparison_all <- rbind(modelComparison_all, data.frame("AIC" = AIC(model2$fit), "sigma" = sigma(model2$fit), "period" = "30+15", "id" = i))
  }
  if (!is.null(model3)) {
    results <- summary(model3)$transformed.table; results$id <- i; results$period <- "15"
    results_all3 <- rbind(results_all3, results)
    modelComparison_all <- rbind(modelComparison_all, data.frame("AIC" = AIC(model3$fit), "sigma" = sigma(model3$fit), "period" = "15", "id" = i))
  }
}

results_all1 <- rownames_to_column(results_all1)
results_all2 <- rownames_to_column(results_all2)
results_all3 <- rownames_to_column(results_all3)

dir.create("results/Titis/individuals", showWarnings = FALSE, recursive = TRUE)
write.csv(na.omit(results_all1), "results/Titis/individuals/modelResults_individual_30days_ActivityOnset.csv")
write.csv(na.omit(results_all2), "results/Titis/individuals/modelResults_individual_30+15days_ActivityOnset.csv")
write.csv(na.omit(results_all3), "results/Titis/individuals/modelResults_individual_15days_ActivityOnset.csv")
write.csv(na.omit(modelComparison_all), "results/Titis/individuals/ModelComparison_individual_ActivityOnset.csv")

# Smoothed, double-plotted 30+15 fit for every titi individual (Fig-1 style)
makeIndivFitPlots(reserve, results_all2, "Activity onset (hours from midnight)", "Titi", "Figures/individuals_30_15")
