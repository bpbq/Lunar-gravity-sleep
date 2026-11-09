library(dplyr)
library(TTR) #MA function
library (lubridate) #date/Time format
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
library(lme4)
library(lmerTest)
library(lmtest)

#remove environment
rm(list=ls())

setwd(here())


#Open file with gravity data data
#Open file with gravity data data
# Load combined dataset containing gravitational tide predictions (dg) and illuminance
df <- read.csv ("data/Data_S4_meanDaily_Gravitational_variation.csv", header = TRUE)
df<-na.omit(df)

#Choose the data group to work with
# Separate gravitational/environmental data by city
Tobas<-df[df$city=="Ingeniero Juarez",]
seattle<-df[df$city=="Seattle",]


####### Functions #####

# Function: circularShift
# Circularly shifts a vector x by 'shift' positions.
# Equivalent to rotating the time series, preserving its autocorrelation structure
# but destroying the cross-correlation with the other variable.
# Used to generate null distributions in the phase shuffling.
circularShift <- function(x, shift) {
  n <- length(x)
  shift <- shift %% n
  if (shift == 0) return(x)
  c(x[(shift + 1):n], x[1:shift])
}


# Function: runShuffledCCF
# Runs the cross-correlation analysis on data with the sleep series circularly
# shifted. Returns a data.frame in the same format as crossCorrResults but with
# the shuffled series.
# Parameters:
#   sleep_data   : data.frame with columns ID_year, bout_dataCollection,
#                  Start_Time_ma, End_Time_ma, Duration_ma, mean_Total_C, Illu
#   max_lag      : maximum lag for the CCF (default 12)
#   has_bout_col : TRUE for Tobas (has bout_dataCollection column), FALSE for Seattle
runShuffledCCF <- function(sleep_data, max_lag = 12, has_bout_col = TRUE) {
  res <- data.frame()

  if (has_bout_col) {
    id_col <- "ID_year"
  } else {
    id_col <- "ID_year"
  }

  for (id in unique(sleep_data[[id_col]])) {
    sleep_id <- sleep_data[sleep_data[[id_col]] == id, ]

    if (has_bout_col) {
      bouts <- unique(sleep_id$bout_dataCollection)
    } else {
      bouts <- 1  # Seattle has no bout_dataCollection in this context
    }

    for (bout in bouts) {
      if (has_bout_col) {
        sp <- sleep_id[sleep_id$bout_dataCollection == bout, ]
      } else {
        sp <- sleep_id
      }

      N <- nrow(sp)
      if (N < (max_lag + 2)) next

      # Random shift: avoid trivial shifts (1 or N-1) so as not to bias
      shift <- sample(3:(N - 3), 1)

      st_s  <- circularShift(sp$Start_Time_ma,  shift)
      et_s  <- circularShift(sp$End_Time_ma,    shift)
      dur_s <- circularShift(sp$Duration_ma,    shift)

      # CCF: Sleep onset (shuffled) vs Gravity (original)
      cc <- ccf(sp$mean_Total_C, st_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Sleep onset", var2 = "Daily mean total gravity"))
      # CCF: Sleep onset (shuffled) vs Illuminance (original)
      cc <- ccf(sp$Illu, st_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Sleep onset", var2 = "Illuminance"))
      # CCF: Sleep offset (shuffled) vs Gravity
      cc <- ccf(sp$mean_Total_C, et_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Sleep offset", var2 = "Daily mean total gravity"))
      # CCF: Sleep offset (shuffled) vs Illuminance
      cc <- ccf(sp$Illu, et_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Sleep offset", var2 = "Illuminance"))
      # CCF: Duration (shuffled) vs Gravity
      cc <- ccf(sp$mean_Total_C, dur_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Duration", var2 = "Daily mean total gravity"))
      # CCF: Duration (shuffled) vs Illuminance
      cc <- ccf(sp$Illu, dur_s, max_lag, plot = FALSE)
      res <- rbind(res, data.frame(ID = id, ACF = cc$acf, Lag = cc$lag,
                                   var1 = "Duration", var2 = "Illuminance"))
    }
  }
  return(res)
}


# Function: CrosscorPlot_double
# Plots the CCF for two variables (Gravity vs Illuminance) with:
# - Points: mean ACF per lag
# - Asterisks: significant lags (adjusted parametric p-values)
# - null_band (optional): data.frame with columns Lag, lower, upper corresponding
#   to the 95% CI of the phase-shuffling null distribution. Shown as a gray ribbon.
#   If NULL (default), no band is added.
CrosscorPlot_double <- function(data, asteriscks=c(), offset=0.04,
                                data2, asteriscks2=c(), offset2=0.05,
                                null_band = NULL) {
  ypos=c()
  ypos2=c()
  if(length(asteriscks)>0){
    ypos<-max((data %>% group_by(Lag) %>% summarize(ACF=mean(ACF)))$ACF)+offset}
  if(length(asteriscks2)>0){
    ypos2<-max((data2 %>% group_by(Lag) %>% summarize(ACF=mean(ACF)))$ACF)+offset}

  if(length(ypos)>0 & length(ypos2)>0){
    if(ypos>ypos2){
      ypos2=ypos
      ypos=ypos+0.04
    } else{
      ypos=ypos2
      ypos2=ypos2+0.04
    }
  }

  p <- ggplot()

  # Null distribution band (phase shuffling) if provided
  if (!is.null(null_band)) {
    p <- p + geom_ribbon(data = null_band,
                         aes(x = Lag, ymin = lower, ymax = upper),
                         fill = "gray70", alpha = 0.35, inherit.aes = FALSE)
  }

  p +
    stat_summary(data = data, aes(x = Lag, y = ACF), colour = "black", fill = "black",
                 pch = 16, size = 0.65, position = position_dodge(width = .1)) +
    stat_summary(data = data2, aes(x = Lag, y = ACF), colour = "black", fill = "white",
                 pch = 21, stroke = 0.5, size = 0.65, position = position_dodge(width = .1)) +
    labs(x = "Lag", y = "ACF") +
    coord_cartesian(ylim = c(-0.2, 0.28))+
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
          axis.title.x = element_text(),
          axis.text = element_text(size=11))+
    annotate("text", x=asteriscks, y=ypos, label="*", color="black", size=6)+
    annotate("text", x=asteriscks2, y=ypos2, label="+", color="black", size=5)
}




###### Tobas Analysis #####

sleep <- read.csv ("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)
# Exclude Spring 2017 cohort (likely due to sample size/quality issues)
sleep<-sleep[sleep$year!="Spring2017",]
Tobas$night<-Tobas$NightDate

#original<-sleep


nights<-Tobas[,"night"]

IDs<-unique(sleep$ID_year)

# Convert sleep$night to POSIXct before computing per-individual max/min,
# so that all_combinations$max and $min are POSIXct and the filter doesn't fail.
sleep$night <- as.POSIXct(sleep$night, format= "%m/%d/%Y")

# Create a master grid of all ID x Date combinations within the min/max dates for each subject
# This ensures that missing nights are explicit rows in the dataframe (crucial for imputation)
all_combinations <- expand.grid(ID_year = IDs, night = nights)
all_combinations$max<-max(sleep$night)
all_combinations$min<-min(sleep$night)

all_combinations<- all_combinations %>%
  group_by(ID_year) %>%
  mutate(max= max(sleep$night[sleep$ID_year == ID_year[1]]), min =min(sleep$night[sleep$ID_year == ID_year[1]])) %>%
  ungroup()

all_combinations$night<-as.POSIXct(all_combinations$night, format= "%m/%d/%Y")

all_combinations<- all_combinations %>%
  group_by(ID_year) %>%
  filter(night <= max & night >= min)
# Merge actual sleep data into the master grid
sleep <- merge(x=sleep, y=all_combinations, by = c("ID_year","night"), all.y = TRUE)

sleep<-sleep[,c("ID_year", "year","night", "Start_Time", "Duration", "End_Time")]


sleep<-sleep[order(sleep$ID_year, sleep$night),]


### Data Imputation:
# Fill missing sleep values using linear interpolation, but only for gaps <= 3 days.
# Larger gaps are left as NA to avoid creating fake data over long periods.
sleep = sleep %>%
  group_by (ID_year) %>%
  mutate_at(colnames(subset(sleep, select = -c(ID_year,year,night))), na.approx, maxgap=3, na.rm=FALSE)%>% arrange(ID_year, night)%>%
  ungroup()

sleep <- sleep %>%
  distinct(night, ID_year, .keep_all = TRUE)

sleepTobas<-sleep

# Detect continuous data segments ("bouts")
# Uses Run Length Encoding (rle) to identify consecutive non-NA days
rle_result <- sleep %>%
  group_by(ID_year) %>% summarise(rle_len=rle(is.na(Start_Time))$lengths, rle_values=rle(is.na(Start_Time))$values) %>% arrange(ID_year)


# Assign a unique "bout ID" to continuous segments of data that are at least 30 days long.
# This prepares the data for cross-correlation (which needs continuous time series).
sleep$bout_dataCollection<-NA
index=1
bout=1
id=rle_result$ID_year[1]
for (i in c(1:length(rle_result$ID_year))){
  len=rle_result$rle_len[i]
  val=rle_result$rle_values[i]
  if (id==rle_result$ID_year[i] & len>=30 & val==FALSE){
    sleep$bout_dataCollection[index:(index+len-1)]<-bout
    bout=bout+1
  } else if(id!=rle_result$ID_year[i]& len>=30 & val==FALSE){
    bout=1
    id=rle_result$ID_year[i]
    sleep$bout_dataCollection[index:(index+len-1)]<-bout
    bout=bout+1
  }
  index=index+len
}

# Filter out any data not belonging to a valid 30+ day bout
sleep<-sleep[!is.na(sleep$bout_dataCollection),]

###### Moving average ######
# Apply 7-day Simple Moving Average (SMA) to smooth the sleep metrics
sleep=sleep[order(sleep$ID_year, sleep$night),]

sleep = sleep %>%
  group_by (ID_year, bout_dataCollection) %>%
  mutate(Start_Time_ma = SMA(Start_Time, 7), End_Time_ma = SMA(End_Time, 7), Duration_ma = SMA(Duration, 7)) %>% #picking 7 days as window
  ungroup() #picking 7 days as window

# Shift data back by 3 days to center the moving average
sleep$night<-sleep$night-days(3)
#sleep<-na.omit(sleep)
sleep<-sleep[!is.na(sleep$Start_Time_ma),]

# Merge smoothed sleep data with Environmental Data (Gravity/Illuminance)
Tobas$night<-as.POSIXct(Tobas$night, format= "%m/%d/%Y")
sleep <- merge(x=sleep, y=Tobas, by = "night", all.x = TRUE)
sleepTobas <- merge(x=sleepTobas, y=Tobas, by = "night", all.x = TRUE)


linearModRes<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(linearModRes)<-c("year", "sleepVar", "Lag", "estimate", "pval")
crossCorrResults<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(crossCorrResults)=c("ID", "ACF", "Lag", "var1", "var2", "group")


##### Cross-Correlation Analysis (Whole sample) #####

# Loop through each individual and calculates Cross-Correlation Function (CCF)
# Computes CCF for:
# 1. Sleep Onset vs Gravity
# 2. Sleep Onset vs Illuminance
# 3. Sleep Offset vs Gravity
# 4. Sleep Offset vs Illuminance
# 5. Duration vs Gravity
# 6. Duration vs Illuminance
# Stores results in 'crossCorrResults' dataframe.
crossCorrResults<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(crossCorrResults)=c("ID", "ACF", "Lag", "var1", "var2")
for (id in unique(sleep$ID_year)){
  sleep_id<-sleep[sleep$ID_year==id,]
  for (bout in unique(sleep_id$bout_dataCollection)) {
    sleep_partial<-sleep_id[sleep_id$bout_dataCollection==bout,]
    
    tiff(paste("Figures/individual_plots/crossCorr/Tobas_",id,"_StartTime_12_gravity.tif", sep=""), width = 1500, height = 1200, res = 300,
         compression = "lzw") # file & size
    crossCorr<-(ccf(sleep_partial$mean_Total_C, sleep_partial$Start_Time_ma,12))
    dev.off()
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep onset", n), 
                                                         "var2" = rep("Daily mean total gravity", n), "bout"=rep(bout, n)))
    
    tiff(paste("Figures/individual_plots/crossCorr/Tobas_",id,"_StartTime_12_Illu.tif", sep=""), width = 1500, height = 1200, res = 300,
         compression = "lzw") # file & size
    crossCorr<-(ccf(  sleep_partial$Illu, sleep_partial$Start_Time_ma,12))
    dev.off()
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep onset", n), "var2" = rep("Illuminance", n), "bout"=rep(bout, n)))
    
    
    crossCorr<-(ccf(  sleep_partial$mean_Total_C,   sleep_partial$End_Time_ma,12))
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep offset", n), "var2" = rep("Daily mean total gravity", n), "bout"=rep(bout, n)))
    
    
    crossCorr<-(ccf(  sleep_partial$Illu, sleep_partial$End_Time_ma,12))
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep offset", n), "var2" = rep("Illuminance", n), "bout"=rep(bout, n)))
    
    
    crossCorr<-(ccf(  sleep_partial$mean_Total_C,   sleep_partial$Duration_ma,12))
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Duration", n), "var2" = rep("Daily mean total gravity", n), "bout"=rep(bout, n)))
    
    
    crossCorr<-(ccf(  sleep_partial$Illu, sleep_partial$Duration_ma,12))
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Duration", n), "var2" = rep("Illuminance", n), "bout"=rep(bout, n)))
  }
}

write.csv(crossCorrResults, "results/Tobas_WholeSampleCrosscor.csv")

####### Statistical Testing of CCFs ######

# Aggregate individual CCF results to test if the correlation at each lag is significantly different from 0
# Uses t-tests followed by Benjamini-Hochberg (BH) correction for multiple comparisons
tests_Tobas<- crossCorrResults %>% group_by(var1, var2, Lag) %>%
  summarise(pval= t.test(ACF, mu = 0)$p.value, t=t.test(ACF)$statistic, df=t.test(ACF)$parameter, estimate=t.test(ACF)$estimate, ci_Low=t.test(ACF)$conf.int[1], ci_High=t.test(ACF)$conf.int[1])


tests_Tobas<-tests_Tobas %>% group_by(var1, var2) %>%  mutate(pval_adj=p.adjust(pval, method = "BH"))

write.csv(tests_Tobas, "results/test_CrossCor_Tobas.csv")


####### Phase Shuffling: null distribution for the CCF (Tobas) #######
#
# To test whether the observed cross-correlations exceed what is expected by
# chance, we circularly shift each individual's sleep series by a random number
# of days and recompute the CCF. This preserves the internal autocorrelation
# structure of the sleep series but destroys the cross-correlation with gravity.
#
# The permutational p-value is the proportion of permutations where the |mean CCF|
# is >= the |observed CCF|. It complements the parametric t-tests above.
#
# NOTE: n_perm=1000 is adequate for publication. It can be reduced to 100-200 for
# fast development. The set.seed guarantees reproducibility.

n_perm <- 1000
set.seed(42)

cat("Running phase shuffling for Tobas (", n_perm, "permutations)...\n")

# Identificar todos los pares variable1 x variable2
var_pairs_tobas <- unique(crossCorrResults[, c("var1", "var2")])
null_lags       <- sort(unique(crossCorrResults$Lag))
n_lags          <- length(null_lags)

# Structure to store the mean CCF of each permutation by lag and variable pair
# List named by "var1|var2"; each element is a matrix [n_perm x n_lags]
null_dist_tobas <- list()
for (vp in seq_len(nrow(var_pairs_tobas))) {
  key <- paste(var_pairs_tobas$var1[vp], var_pairs_tobas$var2[vp], sep = "|")
  null_dist_tobas[[key]] <- matrix(NA_real_, nrow = n_perm, ncol = n_lags)
}

for (perm in seq_len(n_perm)) {
  # Generate CCF with shuffled series for this permutation
  ccf_perm <- runShuffledCCF(sleep, max_lag = 12, has_bout_col = TRUE)

  # Average the CCF over individuals for each lag and variable pair
  means_perm <- ccf_perm %>%
    group_by(var1, var2, Lag) %>%
    summarise(mean_ACF = mean(ACF), .groups = "drop")

  # Almacena en la estructura null
  for (vp in seq_len(nrow(var_pairs_tobas))) {
    key      <- paste(var_pairs_tobas$var1[vp], var_pairs_tobas$var2[vp], sep = "|")
    idx_rows <- means_perm$var1 == var_pairs_tobas$var1[vp] &
                means_perm$var2 == var_pairs_tobas$var2[vp]
    idx_lag  <- match(means_perm$Lag[idx_rows], null_lags)
    valid    <- !is.na(idx_lag)
    null_dist_tobas[[key]][perm, idx_lag[valid]] <- means_perm$mean_ACF[idx_rows][valid]
  }

  if (perm %% 100 == 0) cat("  Tobas: permutacion", perm, "/", n_perm, "\n")
}

saveRDS(null_dist_tobas, "results/null_dist_CCF_Tobas.rds")
cat("Tobas null distribution saved.\n")

# --- Compute permutational p-values for Tobas ---
# For each lag and variable pair: proportion of permutations where
# |mean_null_CCF| >= |mean_observed_CCF|
obs_means_tobas <- crossCorrResults %>%
  group_by(var1, var2, Lag) %>%
  summarise(obs_mean_ACF = mean(ACF), .groups = "drop")

perm_pvals_tobas <- obs_means_tobas %>% mutate(pval_perm = NA_real_)

for (r in seq_len(nrow(perm_pvals_tobas))) {
  key     <- paste(perm_pvals_tobas$var1[r], perm_pvals_tobas$var2[r], sep = "|")
  lag_idx <- which(null_lags == perm_pvals_tobas$Lag[r])
  if (length(lag_idx) == 0 || is.null(null_dist_tobas[[key]])) next
  null_vals <- null_dist_tobas[[key]][, lag_idx]
  obs_val   <- perm_pvals_tobas$obs_mean_ACF[r]
  perm_pvals_tobas$pval_perm[r] <- mean(abs(null_vals) >= abs(obs_val), na.rm = TRUE)
}

# BH adjustment on permutational p-values within each variable pair
perm_pvals_tobas <- perm_pvals_tobas %>%
  group_by(var1, var2) %>%
  mutate(pval_perm_adj = p.adjust(pval_perm, method = "BH")) %>%
  ungroup()

write.csv(perm_pvals_tobas, "results/perm_pvals_CCF_Tobas.csv")
cat("Tobas permutational p-values saved.\n")


##### Plotting: Tobas CCF with null distribution band #####
#
# This section can be run on its own without repeating the permutations.
# If the objects are not in the environment, they are loaded automatically from disk.

if (!exists("crossCorrResults") || !("bout" %in% colnames(crossCorrResults))) {
  # Note: the Tobas CSV has a 'bout' column; the Seattle one does not.
  # We use that to detect which one is loaded.
  crossCorrResults <- read.csv("results/Tobas_WholeSampleCrosscor.csv",
                                stringsAsFactors = FALSE)
}
if (!exists("tests_Tobas")) {
  tests_Tobas <- read.csv("results/test_CrossCor_Tobas.csv", stringsAsFactors = FALSE)
}
if (!exists("null_dist_tobas")) {
  null_dist_tobas <- readRDS("results/null_dist_CCF_Tobas.rds")
}
if (!exists("null_lags")) {
  null_lags <- sort(unique(crossCorrResults$Lag))
}

# The gray band represents the 95% CI of the phase-shuffling null distribution.
# It is built by taking the 97.5th percentile of |null_ACF| at each lag (symmetric band).
# If the observed mean falls outside the band, the lag is significant by permutation.

# Helper function: extracts the null band for a (var1, var2) pair from the RDS
getNullBand <- function(null_dist, var1_str, var2_str, lag_grid) {
  key <- paste(var1_str, var2_str, sep = "|")
  if (is.null(null_dist[[key]])) return(NULL)
  mat <- null_dist[[key]]
  # 97.5th percentile of |ACF| at each lag -> symmetric band around 0
  upper <- apply(abs(mat), 2, quantile, 0.975, na.rm = TRUE)
  data.frame(Lag = lag_grid, lower = -upper, upper = upper)
}

null_lags_tobas <- sort(unique(crossCorrResults$Lag))

##### Plotting: Sleep Onset #####
graficar  <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Sleep onset",]
asterisks <- tests_Tobas$Lag[tests_Tobas$var2=="Daily mean total gravity" & tests_Tobas$var1=="Sleep onset" & tests_Tobas$pval_adj<=0.05]
graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Sleep onset",]
asterisks2 <- tests_Tobas$Lag[tests_Tobas$var2=="Illuminance" & tests_Tobas$var1=="Sleep onset" & tests_Tobas$pval_adj<=0.05]
nb_onset   <- getNullBand(null_dist_tobas, "Sleep onset", "Daily mean total gravity", null_lags_tobas)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks,0.06, graficar2, asterisks2, 0.06, null_band=nb_onset)

tiff(paste("Figures/CrossCor_WholeSample_Tobas_SleepOnset_BothIlluAndMeand.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()

##### Plotting: Sleep Offset ####
graficar  <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Sleep offset",]
asterisks <- tests_Tobas$Lag[tests_Tobas$var2=="Daily mean total gravity" & tests_Tobas$var1=="Sleep offset" & tests_Tobas$pval_adj<=0.05]
graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Sleep offset",]
asterisks2 <- tests_Tobas$Lag[tests_Tobas$var2=="Illuminance" & tests_Tobas$var1=="Sleep offset" & tests_Tobas$pval_adj<=0.05]
nb_offset  <- getNullBand(null_dist_tobas, "Sleep offset", "Daily mean total gravity", null_lags_tobas)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks, 0.06, graficar2, asterisks2, 0.06, null_band=nb_offset)

tiff(paste("Figures/CrossCor_WholeSample_Tobas_SleepOnffset_BothIlluAndMeand.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()

###### Plotting: Sleep Duration #####
graficar  <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Duration",]
asterisks <- tests_Tobas$Lag[tests_Tobas$var2=="Daily mean total gravity" & tests_Tobas$var1=="Duration" & tests_Tobas$pval_adj<=0.05]

graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Duration",]
asterisks2 <- tests_Tobas$Lag[tests_Tobas$var2=="Illuminance" & tests_Tobas$var1=="Duration" & tests_Tobas$pval_adj<=0.05]
nb_dur     <- getNullBand(null_dist_tobas, "Duration", "Daily mean total gravity", null_lags_tobas)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks,0.06, graficar2, asterisks2, 0.06, null_band=nb_dur)

tiff(paste("Figures/CrossCor_WholeSample_Tobas_Duration_BothIlluAndMeand.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()





###### Seattle Analysis #####
# Repeats identical workflow for the Seattle population
# Key difference: Excludes two specific cohorts (AW2024, SS2023)
setwd(here())

sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)

sleep<-sleep[sleep$year!="AW2024" & sleep$year!="SS2023",]



nightsPerPerson<-sleep %>% group_by(ID_year) %>%
  summarise(nights=length(ID_year))
mean(nightsPerPerson$nights)

seattle$night<-seattle$NightDate

nights<-seattle[,"night"]

IDs<-unique(sleep$ID_year)

# Convertir sleep$night a POSIXct antes de calcular max/min (Seattle: formato YYYY-MM-DD)
sleep$night <- as.POSIXct(sleep$night, format= "%Y-%m-%d")

# Re-generate master date grid for Seattle participants
all_combinations <- expand.grid(ID_year = IDs, night = nights)
all_combinations$max<-max(sleep$night)
all_combinations$min<-min(sleep$night)

all_combinations<- all_combinations %>%
  group_by(ID_year) %>%
  mutate(max= max(sleep$night[sleep$ID_year == ID_year[1]]), min =min(sleep$night[sleep$ID_year == ID_year[1]])) %>%
  ungroup()

all_combinations$night<-as.POSIXct(all_combinations$night, format= "%m/%d/%Y")

all_combinations<- all_combinations %>%
  group_by(ID_year) %>%
  filter(night <= max & night >= min)
sleep <- merge(x=sleep, y=all_combinations, by = c("ID_year","night"), all.y = TRUE)


sleep<-sleep[,c("ID_year", "year","night", "Start_Time", "Duration", "End_Time")]


sleep<-sleep[order(sleep$ID_year, sleep$night),]


### Seattle Data Cleaning & Imputation (Same parameters: maxgap=3)
sleep = sleep %>%
  group_by (ID_year) %>%
  mutate_at(colnames(subset(sleep, select = -c(ID_year, year, night))), na.approx, maxgap=3, na.rm=FALSE)%>%
  ungroup()

sleep <- sleep %>%
  distinct(night, ID_year, .keep_all = TRUE)

sleepSeattle<-sleep


# Identification of 30+ day bouts for Seattle
rle_result <- sleep %>%
  group_by(ID_year) %>% summarise(rle_len=rle(is.na(Start_Time))$lengths, rle_values=rle(is.na(Start_Time))$values) %>% arrange(ID_year)



sleep$bout_dataCollection<-NA
index=1
bout=1
id=rle_result$ID_year[1]
for (i in c(1:length(rle_result$ID_year))){
  len=rle_result$rle_len[i]
  val=rle_result$rle_values[i]
  if (id==rle_result$ID_year[i] & len>=30 & val==FALSE){
    sleep$bout_dataCollection[index:(index+len-1)]<-bout
    bout=bout+1
  } else if(id!=rle_result$ID_year[i]& len>=30 & val==FALSE){
    bout=1
    id=rle_result$ID_year[i]
    sleep$bout_dataCollection[index:(index+len-1)]<-bout
    bout=bout+1
  }
  index=index+len
}

sleep<-sleep[!is.na(sleep$bout_dataCollection),]


###### Moving average (Seattle) ######

sleep=sleep[which(!is.na(sleep$Start_Time)),]
sleep=sleep[order(sleep$ID_year, sleep$night),]


sleep = sleep %>%
  group_by (ID_year) %>%
  mutate(Start_Time_ma = SMA(Start_Time, 7), End_Time_ma = SMA(End_Time, 7), Duration_ma = SMA(Duration, 7)) %>% #picking 7 days as window
  ungroup() #picking 7 days as window


sleep$night<-sleep$night-days(3)
sleep<-sleep[!is.na(sleep$Start_Time_ma),]

seattle$night<-as.POSIXct(seattle$night, format= "%m/%d/%Y")
sleep <- merge(x=sleep, y=seattle, by = "night", all.x = TRUE)
sleepSeattle <- merge(x=sleepSeattle, y=seattle, by = "night", all.x = TRUE)


linearModRes<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(linearModRes)<-c("year", "sleepVar", "Lag", "estimate", "pval")
crossCorrResults<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(crossCorrResults)=c("ID_year", "ACF", "Lag", "var1", "var2", "group")


# Cross-Correlation Loop for Seattle
crossCorrResults<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(crossCorrResults)=c("ID", "ACF", "Lag", "var1", "var2")
for (id in unique(sleep$ID_year)){
  sleep_partial<-sleep[sleep$ID_year==id,]
  tiff(paste("Figures/individual_plots/crossCorr/Seattle",id,"_StartTime_12_gravity.tif", sep=""), width = 1500, height = 1200, res = 300,
       compression = "lzw") # file & size
  crossCorr<-(ccf(sleep_partial$mean_Total_C, sleep_partial$Start_Time_ma,12))
  dev.off()
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep onset", n), "var2" = rep("Daily mean total gravity", n)))
  
  tiff(paste("Figures/individual_plots/crossCorr/Seattle",id,"_StartTime_12_Illu.tif", sep=""), width = 1500, height = 1200, res = 300,
       compression = "lzw") # file & size
  crossCorr<-(ccf(sleep_partial$Illu, sleep_partial$Start_Time_ma,12))
  dev.off()
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep onset", n), "var2" = rep("Illuminance", n)))
  
  
  crossCorr<-(ccf(  sleep_partial$mean_Total_C,   sleep_partial$End_Time_ma,12))
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep offset", n), "var2" = rep("Daily mean total gravity", n)))
  
  
  crossCorr<-(ccf(  sleep_partial$Illu, sleep_partial$End_Time_ma,12))
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Sleep offset", n), "var2" = rep("Illuminance", n)))
  
  
  crossCorr<-(ccf(  sleep_partial$mean_Total_C,   sleep_partial$Duration_ma,12))
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Duration", n), "var2" = rep("Daily mean total gravity", n)))
  
  crossCorr<-(ccf(  sleep_partial$Illu, sleep_partial$Duration_ma,12))
  n = length(crossCorr$acf)
  crossCorrResults<-rbind(crossCorrResults, data.frame("ID"=rep(id, n),  "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var1" = rep("Duration", n), "var2" = rep("Illuminance", n)))
  
  
}

write.csv(crossCorrResults, "results/Seattle_WholeSampleCrosscor.csv")

####### Statistical Testing (Seattle) ######

tests_seattle<- crossCorrResults %>% group_by(var1, var2, Lag) %>%
  summarise(pval= t.test(ACF, mu = 0)$p.value, t=t.test(ACF)$statistic, df=t.test(ACF)$parameter, estimate=t.test(ACF)$estimate, ci_Low=t.test(ACF)$conf.int[1], ci_High=t.test(ACF)$conf.int[1])

tests_seattle <- tests_seattle %>% group_by(var1, var2) %>%  mutate(pval_adj=p.adjust(pval, method = "BH"))

write.csv(tests_seattle, "results/test_CrossCor_Seattle.csv")


####### Phase Shuffling: null distribution for the CCF (Seattle) #######
# Same procedure as Tobas. Seattle has no bout_dataCollection column in the
# context of this loop, since each individual has a single continuous bout.

cat("Running phase shuffling for Seattle (", n_perm, "permutations)...\n")

var_pairs_seattle <- unique(crossCorrResults[, c("var1", "var2")])

null_dist_seattle <- list()
for (vp in seq_len(nrow(var_pairs_seattle))) {
  key <- paste(var_pairs_seattle$var1[vp], var_pairs_seattle$var2[vp], sep = "|")
  null_dist_seattle[[key]] <- matrix(NA_real_, nrow = n_perm, ncol = n_lags)
}

for (perm in seq_len(n_perm)) {
  # For Seattle the individuals have no bouts, so we reuse runShuffledCCF
  # with has_bout_col = FALSE (the inner loop iterates over ID without filtering by bout)
  ccf_perm <- runShuffledCCF(sleep, max_lag = 12, has_bout_col = FALSE)

  means_perm <- ccf_perm %>%
    group_by(var1, var2, Lag) %>%
    summarise(mean_ACF = mean(ACF), .groups = "drop")

  for (vp in seq_len(nrow(var_pairs_seattle))) {
    key      <- paste(var_pairs_seattle$var1[vp], var_pairs_seattle$var2[vp], sep = "|")
    idx_rows <- means_perm$var1 == var_pairs_seattle$var1[vp] &
                means_perm$var2 == var_pairs_seattle$var2[vp]
    idx_lag  <- match(means_perm$Lag[idx_rows], null_lags)
    valid    <- !is.na(idx_lag)
    null_dist_seattle[[key]][perm, idx_lag[valid]] <- means_perm$mean_ACF[idx_rows][valid]
  }

  if (perm %% 100 == 0) cat("  Seattle: permutation", perm, "/", n_perm, "\n")
}

saveRDS(null_dist_seattle, "results/null_dist_CCF_Seattle.rds")

# --- Seattle permutational p-values ---
obs_means_seattle <- crossCorrResults %>%
  group_by(var1, var2, Lag) %>%
  summarise(obs_mean_ACF = mean(ACF), .groups = "drop")

perm_pvals_seattle <- obs_means_seattle %>% mutate(pval_perm = NA_real_)

for (r in seq_len(nrow(perm_pvals_seattle))) {
  key     <- paste(perm_pvals_seattle$var1[r], perm_pvals_seattle$var2[r], sep = "|")
  lag_idx <- which(null_lags == perm_pvals_seattle$Lag[r])
  if (length(lag_idx) == 0 || is.null(null_dist_seattle[[key]])) next
  null_vals <- null_dist_seattle[[key]][, lag_idx]
  obs_val   <- perm_pvals_seattle$obs_mean_ACF[r]
  perm_pvals_seattle$pval_perm[r] <- mean(abs(null_vals) >= abs(obs_val), na.rm = TRUE)
}

perm_pvals_seattle <- perm_pvals_seattle %>%
  group_by(var1, var2) %>%
  mutate(pval_perm_adj = p.adjust(pval_perm, method = "BH")) %>%
  ungroup()

write.csv(perm_pvals_seattle, "results/perm_pvals_CCF_Seattle.csv")
cat("Seattle permutational p-values saved.\n")


##### Plotting: Seattle CCF with null distribution band #####
#
# Automatic load from disk if the objects are not in the environment.

crossCorrResults  <- read.csv("results/Seattle_WholeSampleCrosscor.csv",
                               stringsAsFactors = FALSE)
if (!exists("tests_seattle")) {
  tests_seattle <- read.csv("results/test_CrossCor_Seattle.csv", stringsAsFactors = FALSE)
}
if (!exists("null_dist_seattle")) {
  null_dist_seattle <- readRDS("results/null_dist_CCF_Seattle.rds")
}

null_lags_seattle <- sort(unique(crossCorrResults$Lag))

##### Plotting: Seattle Onset #####
graficar   <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Sleep onset",]
asterisks  <- tests_seattle$Lag[tests_seattle$var2=="Daily mean total gravity" & tests_seattle$var1=="Sleep onset" & tests_seattle$pval_adj<=0.05]
graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Sleep onset",]
asterisks2 <- tests_seattle$Lag[tests_seattle$var2=="Illuminance" & tests_seattle$var1=="Sleep onset" & tests_seattle$pval_adj<=0.05]
nb_onset   <- getNullBand(null_dist_seattle, "Sleep onset", "Daily mean total gravity", null_lags_seattle)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks, 0.06, graficar2, asterisks2, 0.06, null_band=nb_onset)

tiff(paste("Figures/CrossCor_WholeSample_Seattle_SleepOnset_BothIlluAndMean.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()


##### Plotting: Seattle Offset ####
graficar   <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Sleep offset",]
asterisks  <- tests_seattle$Lag[tests_seattle$var2=="Daily mean total gravity" & tests_seattle$var1=="Sleep offset" & tests_seattle$pval_adj<=0.05]
graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Sleep offset",]
asterisks2 <- tests_seattle$Lag[tests_seattle$var2=="Illuminance" & tests_seattle$var1=="Sleep offset" & tests_seattle$pval_adj<=0.05]
nb_offset  <- getNullBand(null_dist_seattle, "Sleep offset", "Daily mean total gravity", null_lags_seattle)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks, 0.06, graficar2, asterisks2, 0.06, null_band=nb_offset)

tiff(paste("Figures/CrossCor_WholeSample_Seattle_SleepOffset_BothIlluAndMean.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()

###### Plotting: Seattle Duration #####
graficar   <- crossCorrResults[crossCorrResults$var2=="Daily mean total gravity" & crossCorrResults$var1=="Duration",]
asterisks  <- tests_seattle$Lag[tests_seattle$var2=="Daily mean total gravity" & tests_seattle$var1=="Duration" & tests_seattle$pval_adj<=0.05]
graficar2  <- crossCorrResults[crossCorrResults$var2=="Illuminance" & crossCorrResults$var1=="Duration",]
asterisks2 <- tests_seattle$Lag[tests_seattle$var2=="Illuminance" & tests_seattle$var1=="Duration" & tests_seattle$pval_adj<=0.05]
nb_dur     <- getNullBand(null_dist_seattle, "Duration", "Daily mean total gravity", null_lags_seattle)

crossCorrPlot_sleepOnset_Total_C = CrosscorPlot_double(graficar,asterisks, 0.06, graficar2, asterisks2, 0.06, null_band=nb_dur)

tiff(paste("Figures/CrossCor_WholeSample_Seattle_Duration_BothIlluAndMean.tif", sep=""), width = 1300, height = 800, res = 300,
     compression = "lzw")
print(crossCorrPlot_sleepOnset_Total_C)
dev.off()
