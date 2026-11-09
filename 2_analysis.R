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
library(GLMMcosinor)
library(tibble)
library(ggpubr)
library(here)


#remove environment
rm(list=ls())

# Set working directory to the project root (using the 'here' package is good practice for reproducibility)
setwd(here())


#### SEATTLE ANALYSIS ####
# Load the cleaned Seattle sleep dataset
sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)


# Exploratory Data Analysis loop:
# Iterates through each unique 'year' (cohort) in the Seattle data.
# Prints:
# 1. The cohort year name.
# 2. Number of unique participants in that cohort.
# 3. First recorded night.
# 4. Last recorded night.
# 5. Average number of recorded nights per person in that cohort.
for(i in unique(sleep$year)){
  print(i)
  print(length(unique(sleep$ID[sleep$year==i])))
  print(min(sleep$night[sleep$year==i]))
  print(max(sleep$night[sleep$year==i]))
  nightsPerPerson<-sleep[sleep$year==i,] %>% group_by(ID_year) %>%
    summarise(nights=length(ID_year))
  print(mean(nightsPerPerson$nights))
  
}

# General descriptive statistics for the whole Seattle dataset
length(unique(sleep$ID)) # Total unique participants
descriptive<-sleep %>% group_by(ID_year) %>% summarize(age=mean(Age), sexF=Sex[1], ID=ID[1])
mean_sd(descriptive$age) # Mean and SD of age
range(descriptive$age) # Age range
descriptive<-sleep %>% group_by(ID) %>% summarize(sex=Sex[1])
count(descriptive$sex=="F") # Count females
count(descriptive$sex=="M") # Count males

###### GLMMcosinor Analysis: Whole Sample (Seattle) #####
# This section fits a Cosinor model to the entire Seattle dataset without distinguishing by cohort.
# It tests for a 30-day (lunar) rhythm while controlling for Type of Day (Weekend/Weekday), Sex, and Age.
# It includes random effects for the intercept and the cosinor parameters (amplitude/acrophase) by ID.

#### Start_Time (Sleep Onset) #####

mixed_mod <- cglmm(
  Start_Time ~ typeOfDay+ Sex+ Age + amp_acro(DayinCycle,
                                              n_components = 1,
                                              period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
# Save model results
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_StartTime_wholeSample.csv")



#### End_Time (Sleep Offset) #####

mixed_mod <- cglmm(
  End_Time ~ typeOfDay+ Sex+ Age + amp_acro(DayinCycle,
                                            n_components = 1,
                                            period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
# Save model results
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_EndTime_wholeSample.csv")



#### Duration #####

mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex + Age + amp_acro(DayinCycle,
                                              n_components = 1,
                                              period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
# Save model results
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_Duration_wholeSample.csv")




###### GLMMcosinor Analysis: Whole Sample, 30+15 model (Seattle) #####
# Combined lunar (30-day) + semilunar (15-day) version of the pooled (whole-sample)
# models above. Random effects on the first (30-day) component only, as in the
# cohort 30+15 models, for reliable convergence.
# Outputs: results/seattle/sample/GLMMCosinor_*_wholeSample_30+15.csv

#### Start_Time (Sleep Onset): 30+15 ####
mixed_mod <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/seattle/sample/GLMMCosinor_StartTime_wholeSample_30+15.csv")

#### End_Time (Sleep Offset): 30+15 ####
mixed_mod <- cglmm(
  End_Time ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/seattle/sample/GLMMCosinor_EndTime_wholeSample_30+15.csv")

#### Duration: 30+15 ####
mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/seattle/sample/GLMMCosinor_Duration_wholeSample_30+15.csv")




###### COHORT ANALYSIS (Seattle) #####
# This section fits separate Cosinor models for specific cohorts (years).
# It filters out years with fewer than 10 unique participants.

sleep<-sleep %>% group_by(year) %>%
  filter(length(unique(ID_year))>=10)

#### Start_Time #####

# Define factor levels to order the cohorts in the output
sleep$year <- factor(
  sleep$year,
  levels = c( "SS2020", "AW2020", "SS2024")
)

# Fit mixed Cosinor model with 'year' (cohort) as a grouping variable for the rhythm parameters.
# This allows the amplitude and acrophase (phase) to vary by cohort.
mixed_mod <- cglmm(
  Start_Time ~ typeOfDay  + Sex+ Age + year + amp_acro(DayinCycle,
                                                       n_components = 1,
                                                       period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_StartTime.csv")

# Test for differences in acrophase (phase) between cohorts
test_cosinor_levels(mixed_mod, x_str = "year", param = "acr")

#### Duration #####

mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex+ Age+ year + amp_acro(DayinCycle,
                                                   n_components = 1,
                                                   period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_Duration.csv")



#### End_Time #####

mixed_mod <- cglmm(
  End_Time ~ typeOfDay + Sex+ Age+ year + amp_acro(DayinCycle,
                                                   n_components = 1,
                                                   period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/seattle/sample/GLMMCosinor_EndTime.csv")



###### COHORT ANALYSIS: 30+15 day model (Seattle) #####
#
# The combined 30+15 model explicitly captures the lunar (~30 days) and
# semi-lunar (~15 days) components. The 30-day model forces a single peak per
# cycle, which can produce wrong phase estimates when the real pattern is bimodal
# (gravity has peaks at both full moon and new moon).
#
# Random effects: only component 1 (amp_acro1) to avoid the convergence problems
# common in 2-component models with full random structure.
# Outputs: results/seattle/sample/GLMMCosinor_*_30+15.csv

sleep$year <- factor(sleep$year, levels = c("SS2020", "AW2020", "SS2024"))

#### Start_Time: 30+15 model (Seattle cohorts) ####
mixed_mod_30_15 <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/seattle/sample/GLMMCosinor_StartTime_30+15.csv")
# Test acrophase differences of the 30-day component across cohorts
test_cosinor_levels(mixed_mod_30_15, x_str = "year", param = "acr")

#### Duration: 30+15 model (Seattle cohorts) ####
mixed_mod_30_15 <- cglmm(
  Duration ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/seattle/sample/GLMMCosinor_Duration_30+15.csv")

#### End_Time: 30+15 model (Seattle cohorts) ####
mixed_mod_30_15 <- cglmm(
  End_Time ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/seattle/sample/GLMMCosinor_EndTime_30+15.csv")



##### Individual Fits Analysis (Seattle) #####
# This section analyzes the results of individual-level cosinor fits (calculated previously and saved to CSVs).
# It calculates the proportion of individuals who show significant lunar (30-day), semi-lunar (15-day), or mixed (30+15) rhythms.

### Start_Time ###

data_30<-read.csv ("results/seattle/individuals/modelResults_individual_30days_StartTime.csv", header = TRUE)
data_15<-read.csv ("results/seattle/individuals/modelResults_individual_15days_StartTime.csv", header = TRUE)
data_30_15<-read.csv ("results/seattle/individuals/modelResults_individual_30+15days_StartTime.csv", header = TRUE)
comparisons<-read.csv ("results/seattle/individuals/ModelComparison_individual_StartTime.csv", header = TRUE)

# Calculate proportion of significant fits (p <= 0.05) for 30-day and 15-day models
prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))

# For the 30+15 model, check if *any* p-value (for either component) is significant
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))

# Calculate proportion of individuals showing ANY significant rhythm across the three models
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

# Compare models using AIC to see which one fits best for each individual
comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC))
table(comparisons$period) # Print counts of best-fit models


### End_Time ###
# (Repeats same logic as Start_Time)
data_30<-read.csv ("results/seattle/individuals/modelResults_individual_30days_EndTime.csv", header = TRUE)
data_15<-read.csv ("results/seattle/individuals/modelResults_individual_15days_EndTime.csv", header = TRUE)
data_30_15<-read.csv ("results/seattle/individuals/modelResults_individual_30+15days_EndTime.csv", header = TRUE)
comparisons<-read.csv ("results/seattle/individuals/ModelComparison_individual_EndTime.csv", header = TRUE)

prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC))
table(comparisons$period)


### Duration ###
# (Repeats same logic as Start_Time)
data_30<-read.csv ("results/seattle/individuals/modelResults_individual_30days_Duration.csv", header = TRUE)
data_15<-read.csv ("results/seattle/individuals/modelResults_individual_15days_Duration.csv", header = TRUE)
data_30_15<-read.csv ("results/seattle/individuals/modelResults_individual_30+15days_Duration.csv", header = TRUE)
comparisons<-read.csv ("results/seattle/individuals/ModelComparison_individual_Duration.csv", header = TRUE)

prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC))
table(comparisons$period)




####### TOBAS ANALYSIS ########
#remove environment
rm(list=ls())

sleep <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)


# ============================================================
# TOBAS DEMOGRAPHIC MERGE — LOCALLY AVAILABLE ONLY
# ============================================================
# Sex and Age were removed from Data_S2 for privacy (public repository).
# The Toba_complete_dataset.csv file contains the full data and must NOT be
# uploaded to the repository.
#
# To run WITHOUT this file: comment out the block between BEGIN and END.
# The models below use Sex + Age; if this block is commented out,
# you must also replace "+ Sex + Age" with "" in all Tobas formulas.
#
# Mapping: uid (complete, with leading zeros) == ID (public, numeric)
#        sexo: "Masculino"/"Femenino" -> Sex: "M"/"F"
#        age: mean per individual (varies slightly across nights due to exact computation)

# --- BEGIN demographic merge (comment out if Toba_complete_dataset.csv unavailable) ---
demo_tobas <- read.csv("data/Toba_complete_dataset.csv", header = TRUE)

demo_tobas <- demo_tobas %>%
  mutate(ID = as.integer(uid),                              # strip leading zeros
         Sex = ifelse(sexo == "Masculino", "M", "F")) %>%  # homogenize with Seattle
  group_by(ID) %>%
  summarise(Sex = Sex[1],
            Age = mean(age, na.rm = TRUE),
            .groups = "drop")

sleep <- sleep %>%
  mutate(ID = as.integer(ID)) %>%
  left_join(demo_tobas, by = "ID")
# --- END demographic merge ---


###### GLMMcosinor Analysis: Whole Group (Tobas) #####
# Same logic as Seattle whole group analysis: fit 30-day model to all data.

#### Start_Time #####

mixed_mod <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age + amp_acro(DayinCycle,
                                                n_components = 1,
                                                period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_StartTime_wholeSample.csv")



#### End_Time #####

mixed_mod <- cglmm(
  End_Time ~ typeOfDay + Sex + Age + amp_acro(DayinCycle,
                                              n_components = 1,
                                              period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_EndTime_wholeSample.csv")



#### Duration #####

mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex + Age + amp_acro(DayinCycle,
                                              n_components = 1,
                                              period = 30) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_Duration_wholeSample.csv")


###### GLMMcosinor Analysis: Whole Sample, 30+15 model (Tobas) #####
# Combined lunar (30-day) + semilunar (15-day) version of the pooled (whole-sample)
# models above. Random effects on the first (30-day) component only, as in the
# cohort 30+15 models, for reliable convergence.
# Outputs: results/Tobas/sample/GLMMCosinor_*_wholeSample_30+15.csv

#### Start_Time: 30+15 ####
mixed_mod <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/Tobas/sample/GLMMCosinor_StartTime_wholeSample_30+15.csv")

#### End_Time: 30+15 ####
mixed_mod <- cglmm(
  End_Time ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/Tobas/sample/GLMMCosinor_EndTime_wholeSample_30+15.csv")

#### Duration: 30+15 ####
mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex + Age +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15)) + (1 + amp_acro1 | ID),
  data = sleep)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table),
          "results/Tobas/sample/GLMMCosinor_Duration_wholeSample_30+15.csv")





###### COHORT ANALYSIS (Tobas) #####
# Fits Cosinor models by cohort (year) for Tobas data.

sleep<-sleep %>% group_by(year) %>%
  filter(length(unique(ID_year))>=10)

unique(sleep$year)

#### Start_Time #####

# Set reference level for the year factor
sleep$year <- relevel(as.factor(sleep$year), ref = "Fall-Winter2016")

# The commented out code below shows an alternative way to set factor levels manually
# sleep$year <- factor(
#   sleep$year,
#   levels = c("Fall-Winter2016", "Spring2018", "Fall-Winter2023", "Spring2023", "Summer 2024", "Fall-Winter2024")
# )

mixed_mod <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age + year + amp_acro(DayinCycle,
                                                       n_components = 1,
                                                       period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

# Test for differences in phases between cohorts
test_cosinor_levels(mixed_mod, x_str = "year", param = "acr")

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_StartTime.csv")



#### Duration #####
mixed_mod <- cglmm(
  Duration ~ typeOfDay + Sex + Age + year + amp_acro(DayinCycle,
                                                     n_components = 1,
                                                     period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_Duration.csv")



#### End_Time #####

mixed_mod <- cglmm(
  End_Time ~ typeOfDay + Sex + Age + year + amp_acro(DayinCycle,
                                                     n_components = 1,
                                                     period = 30, group = year) + (1 + amp_acro1 | ID),
  data = sleep)

summ<-summary(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Tobas/sample/GLMMCosinor_EndTime.csv")



###### COHORT ANALYSIS: 30+15 day model (Tobas) #####
# Same reasoning as for Seattle: the 30+15 model captures the gravitational
# bimodality (peaks at full moon AND new moon) that the 30-day model omits.
# Particularly relevant for Tobas, where the 15-day (semilunar) component is strong.

sleep$year <- relevel(as.factor(sleep$year), ref = "Fall-Winter2016")

#### Start_Time: 30+15 model (Tobas cohorts) ####
mixed_mod_30_15 <- cglmm(
  Start_Time ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/Tobas/sample/GLMMCosinor_StartTime_30+15.csv")
test_cosinor_levels(mixed_mod_30_15, x_str = "year", param = "acr")

#### Duration: 30+15 model (Tobas cohorts) ####
mixed_mod_30_15 <- cglmm(
  Duration ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/Tobas/sample/GLMMCosinor_Duration_30+15.csv")

#### End_Time: 30+15 model (Tobas cohorts) ####
mixed_mod_30_15 <- cglmm(
  End_Time ~ typeOfDay + Sex + Age + year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) +
    (1 + amp_acro1 | ID),
  data = sleep)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/Tobas/sample/GLMMCosinor_EndTime_30+15.csv")



##### Titi Monkeys Analysis ####
setwd(here())
data <- read.csv ("data/Data_S3_all_onsetData_titi.csv", header = TRUE)
data<-data[,c(2:(length(colnames(data))))] # Remove index column if present

# Fit mixed Cosinor model for Titi Monkey activity onset
# Groups the rhythm by 'Year' (Cohort)
# Includes random effects for individual IDs
mixed_mod <- cglmm(
  time ~ Year + amp_acro(DayinCycle,
                         n_components = 1,
                         period = 30, group = Year) + (1 + amp_acro1 | id2),
  data = data)

summ<-summary(mixed_mod)
autoplot(mixed_mod)
write.csv(as.data.frame(summary(mixed_mod)$transformed.table), "results/Titis/GLMMCosinor_ActivityOnset.csv")


#### Titi Monkeys: 30+15 model ####
# Same combined model for the titi monkeys. The argument is identical:
# if gravity modulates behavior, the model must reflect both periodicities.
mixed_mod_30_15 <- cglmm(
  time ~ Year +
    amp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = Year) +
    (1 + amp_acro1 | id2),
  data = data)

summ_30_15 <- summary(mixed_mod_30_15)
write.csv(as.data.frame(summ_30_15$transformed.table),
          "results/Titis/GLMMCosinor_ActivityOnset_30+15.csv")



##### Individual Fits Analysis (Tobas) #####
# Similar to Seattle individual analysis, but calculating proportions of significant fits.

### Start_Time ###

data_30<-read.csv ("results/Tobas/individuals/modelResults_individual_30days_StartTime.csv", header = TRUE)
data_15<-read.csv ("results/Tobas/individuals/modelResults_individual_15days_StartTime.csv", header = TRUE)
data_30_15<-read.csv ("results/Tobas/individuals/modelResults_individual_30+15days_StartTime.csv", header = TRUE)
comparisons<-read.csv ("results/Tobas/individuals/ModelComparison_individual_StartTime.csv", header = TRUE)

prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC))
table(comparisons$period)


### End_Time ###

data_30<-read.csv ("results/Tobas/individuals/modelResults_individual_30days_EndTime.csv", header = TRUE)
data_15<-read.csv ("results/Tobas/individuals/modelResults_individual_15days_EndTime.csv", header = TRUE)
data_30_15<-read.csv ("results/Tobas/individuals/modelResults_individual_30+15days_EndTime.csv", header = TRUE)
comparisons<-read.csv ("results/Tobas/individuals/ModelComparison_individual_EndTime.csv", header = TRUE)

prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

# Identify IDs that had at least one significant fit
ids<-unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05]))

# Filter model comparisons to include only those significant individuals
comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC) & id %in% ids)
table(comparisons$period)


### Duration ###

data_30<-read.csv ("results/Tobas/individuals/modelResults_individual_30days_duration.csv", header = TRUE)
data_15<-read.csv ("results/Tobas/individuals/modelResults_individual_15days_duration.csv", header = TRUE)
data_30_15<-read.csv ("results/Tobas/individuals/modelResults_individual_30+15days_duration.csv", header = TRUE)
comparisons<-read.csv ("results/Tobas/individuals/ModelComparison_individual_duration.csv", header = TRUE)

prop_30<-count(grepl("amp", data_30$rowname) & data_30$p.value<=0.05)/length(unique(data_30$id))
prop_15<-count(grepl("amp", data_15$rowname) & data_15$p.value<=0.05)/length(unique(data_15$id))
data_30_15<-data_30_15[grepl("amp",data_30_15$rowname),]
data_30_15 <- data_30_15 %>% group_by(id) %>% summarise(min_p=min(p.value))
prop_30_15 <- count(data_30_15$min_p<=0.05)/length(unique(data_30_15$id))
prop_any<-length(unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05])))/length(unique(data_30$id))

ids<-unique(c(data_30$id[grepl("amp", data_30$rowname) & data_30$p.value<=0.05], data_15$id[grepl("amp", data_15$rowname) & data_15$p.value<=0.05],data_30_15$id[data_30_15$min_p<=0.05]))


comparisons<-comparisons %>% group_by(id) %>% filter(AIC==min(AIC) & id %in% ids)
table(comparisons$period)
