##############################################################################
#  Author: Jamie Robertson (The Nature Conservancy, Seattle, Washington, USA)
#  Date modified: February 16, 2021 [last modified]
#  Purpose: Calculate median values and confidence intervals for unique
#           combinations of input tables, here with columns named Scenario, 
#           ScenName, County, Activity, Subactivity, and Year.
#  Comments: 
#  Disclaimer: Use at your own risk. Some portions of this code have not been
#              been completed and may not function properly. The authors take
#              no responsiblity and will not be liable for use of this script.
##############################################################################


#load packages that migh be used throughout the entire script
library(ggplot2)
library(ggforce)
library(grid)
library(truncnorm)
library(tidyverse)
library(generics)
library(gridExtra)
library(gtable)
library(ggpubr)

#library( plyr )
library(dplyr)
#install.packages("plyr","phonTools")
library(plyr)
library(phonTools) #for zeros() function

rm(list=ls())

setwd("~/Workspace/WA_NCS/ANALYSIS_OUTPUTS")
out_dir <- ("~/Workspace/WA_NCS/SUMMED_OUTPUTS")

#in_list <- unlist(ldply( .data = list.files(pattern="*.csv"))) #this stopped working
in_list <- list.files(pattern="*.csv")

# Scenario labels in order of scenario indices, where (1,2,3)==("Lim","Amb","Mod")
sce_names <- c("Lim","Amb","Mod")

# Empty lists to add to in "Create lists..." loop
cty_list <- c()
act_list <- c()
sce_list <- c()
sub_list <- c()

# Create lists of counties, activities, scenarios, and sub-activities
for (f in in_list) {
  # Set variables of each filename part (substring)
  cty <- unlist(strsplit(f, "_"))[1]
  act <- unlist(strsplit(f, "_"))[2]
  sce <- unlist(strsplit(f, "_"))[3]
  sub <- str_sub(unlist(strsplit(f, "_"))[4], 1, -5) # strips ".csv" from the last substring
  # Create lists of files by their substring
  if (!(cty %in% cty_list)) {
    cty_list <- c(cty_list, cty)
  }
  if (!(act %in% act_list)) {
    act_list <- c(act_list, act)
  }
  if (!(sce %in% sce_list)) {
    sce_list <- c(sce_list, sce)
  }
  if (!(sub %in% sub_list)) {
    sub_list <- c(sub_list, sub)
  }
}

# Sum statewide by scenario
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "Year", "Median", "CI_High", "CI_Low")
for (scx in sce_list) {
  sum_list <- c()
  for (f in in_list) {
    if (scx == unlist(strsplit(f, "_"))[3]) {
      sum_list <- c(sum_list, f)
    }
  }
  rx <- data.frame(zeros(1000,32))
  for (x in sum_list) {
    rc <- read.csv(x)#[-c(1)]
    if (is.na(rc[1,2]) | is.null(rc[1,2])) {rc <- zeros(rc)} # changes NA or null input values to zeros
    rx <- rx + rc
  }
  summed <- rx
  # Set variables of each filename part (substring)
  Scenario <- scx
  ScenName <- sce_names[strtoi(scx)]
  Year <- (1:31)
  # Summarize RESULTS to median and 95% CI for each year, each county
  CI_High <- apply(summed[-c(1)],c(2),quantile,0.95, na.rm=T) # [-c(1)] drops iteration number column
  Median <- apply(summed[-c(1)],c(2),median, na.rm=T)         # [-c(1)] drops iteration number column
  CI_Low <- apply(summed[-c(1)],c(2),quantile,0.05, na.rm=T)  # [-c(1)] drops iteration number column
  # combine to one dataframe
  interm_df <- t(rbind(Scenario, ScenName, Year, Median, CI_High, CI_Low))
  out_df <-  rbind(out_df,interm_df)      
#  scen_out_name <- paste0("statewide_all_", ScenName)        # EDIT THIS FOR OUTPUT
#  write.csv(out_df, paste0(out_dir, "/", scen_out_name, ".csv"))
}
out_name <- paste0("statewide_all_scen")        # EDIT THIS FOR OUTPUT
write.csv(out_df, paste0(out_dir, "/", out_name, ".csv"))

# Sum statewide by scenario and activity
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "Activity", "Year", "Median", "CI_High", "CI_Low")
for (scx in sce_list) {
  for (acx in act_list) {  #test ag with act_list[2]
    sum_list <- c()
    for (f in in_list) {
      if ((scx == unlist(strsplit(f, "_"))[3]) & (acx == unlist(strsplit(f, "_"))[2])) { #scenario and activity
        sum_list <- c(sum_list, f)
      }
    }
    rx <- data.frame(zeros(1000,32))
    for (x in sum_list) {
      rc <- read.csv(x)#[-c(1)]
      if (is.na(rc[1,2]) | is.null(rc[1,2])) {rc <- zeros(rc)} # changes NA or null input values to zeros
      rx <- rx + rc
    }
    summed <- rx
    # Set variables of each filename part (substring)
    Activity <- acx
    Scenario <- scx
    ScenName <- sce_names[strtoi(scx)]
    Year <- (1:31)
    # Summarize RESULTS to median and 95% CI for each year, each county
    CI_High <- apply(summed[-c(1)],c(2),quantile,0.95, na.rm=T) # [-c(1)] drops iteration number column
    Median <- apply(summed[-c(1)],c(2),median, na.rm=T)         # [-c(1)] drops iteration number column
    CI_Low <- apply(summed[-c(1)],c(2),quantile,0.05, na.rm=T)  # [-c(1)] drops iteration number column
    # combine to one dataframe
    interm_df <- t(rbind(Scenario, ScenName, Activity, Year, Median, CI_High, CI_Low))
    out_df <-  rbind(out_df,interm_df) 
  }
}
out_name <- paste0("statewide_all_scen_activity_results")        # EDIT THIS FOR OUTPUT
write.csv(out_df, paste0(out_dir, "/", out_name,".csv"))

# Sum statewide by scenario and sub-activity
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "Subactivity", "Year", "Median", "CI_High", "CI_Low")
for (scx in sce_list) {
  for (sux in sub_list) {
    sum_list <- c()
    for (f in in_list) {
      if ((scx == unlist(strsplit(f, "_"))[3]) & (sux == str_sub(unlist(strsplit(f, "_"))[4], 1, -5))) { #scenario and sub-activity
        sum_list <- c(sum_list, f)
      }
    }
    rx <- data.frame(zeros(1000,32))
    for (x in sum_list) {
      rc <- read.csv(x)#[-c(1)]
      if (is.na(rc[1,2]) | is.null(rc[1,2])) {rc <- zeros(rc)} # changes NA or null input values to zeros
      rx <- rx + rc
    }
    summed <- rx
    # Set variables of each filename part (substring)
    Scenario <- scx
    ScenName <- sce_names[strtoi(scx)]
    Subactivity <- sux
    Year <- (1:31)
    # Summarize RESULTS to median and 90% CI for each year, each county
    CI_High <- apply(summed[-c(1)],c(2),quantile,0.95, na.rm=T) # [-c(1)] drops iteration number column
    Median <- apply(summed[-c(1)],c(2),median, na.rm=T)         # [-c(1)] drops iteration number column
    CI_Low <- apply(summed[-c(1)],c(2),quantile,0.05, na.rm=T)  # [-c(1)] drops iteration number column
    # combine to one dataframe
    interm_df <- t(rbind(Scenario, ScenName, Subactivity, Year, Median, CI_High, CI_Low))
    out_df <-  rbind(out_df,interm_df) 
  }
}
out_name <- paste0("statewide_all_scen_subactivity_results")        # EDIT THIS FOR OUTPUT
write.csv(out_df, paste0(out_dir, "/", out_name,".csv"))

# Sum by scenario and county
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "County", "Year", "Median", "CI_High", "CI_Low")
for (scx in sce_list) {
  for (ctx in cty_list) {
    sum_list <- c()
    for (f in in_list) {
      if ((scx == unlist(strsplit(f, "_"))[3]) & (ctx == unlist(strsplit(f, "_"))[1])) { #scenario, county
        sum_list <- c(sum_list, f)
      }
    }
    rx <- data.frame(zeros(1000,32))
    for (x in sum_list) {
      rc <- read.csv(x)#[-c(1)]
      if (is.na(rc[1,2]) | is.null(rc[1,2])) {rc <- zeros(rc)} # changes NA or null input values to zeros
      rx <- rx + rc
    }
    summed <- rx
    # Summarize RESULTS to median and 95% CI for each year, each county
    CI_High <- format(round(apply(summed[-c(1)],c(2),quantile,0.95, na.rm=T), 7), nsmall = 7) # [-c(1)] drops iteration number column. Round to 7th dec.
    Median <- format(round(apply(summed[-c(1)],c(2),median, na.rm=T), 7), nsmall = 7)         # [-c(1)] drops iteration number column. Round to 7th dec.
    CI_Low <- format(round(apply(summed[-c(1)],c(2),quantile,0.05, na.rm=T), 7), nsmall = 7)  # [-c(1)] drops iteration number column. Round to 7th dec.
    Scenario <- scx
    ScenName <- sce_names[strtoi(scx)]
    County <- ctx
    Year <- (1:31)
    # combine to one dataframe
    interm_df <- t(rbind(Scenario, ScenName, County, Year, Median, CI_High, CI_Low))
    out_df <- rbind(out_df, interm_df)
  }
  #out_df <- t(out_df)
  #out_name <- paste0("counties_",sce_names[strtoi(scx)])        # EDIT THIS FOR OUTPUT
}
out_name <- "all_scen_county_results"
write.csv(out_df, paste0(out_dir, "/", out_name,".csv"))

# Sum by scenario, county, and activity
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "County", "Activity", "Year", "Median", "CI_High", "CI_Low")
for (scx in sce_list) {
  for (ctx in cty_list) {
    for (acx in act_list) {
      sum_list <- c()
      for (f in in_list) {
        if ((scx == unlist(strsplit(f, "_"))[3]) & (ctx == unlist(strsplit(f, "_"))[1]) & (acx == unlist(strsplit(f, "_"))[2])) { #scenario, county, and activity
          sum_list <- c(sum_list, f)
        }
      }
      rx <- data.frame(zeros(1000,32))
      for (x in sum_list) {
        rc <- read.csv(x)#[-c(1)]
        if (is.na(rc[1,2]) | is.null(rc[1,2])) {rc <- zeros(rc)} # changes NA or null input values to zeros
        rx <- rx + rc
      }
      summed <- rx
      # Summarize RESULTS to median and 95% CI for each year, each county
      CI_High <- format(round(apply(summed[-c(1)],c(2),quantile,0.95, na.rm=T), 7), nsmall = 7) # [-c(1)] drops iteration number column. Round to 7th dec.
      Median <- format(round(apply(summed[-c(1)],c(2),median, na.rm=T), 7), nsmall = 7)         # [-c(1)] drops iteration number column. Round to 7th dec.
      CI_Low <- format(round(apply(summed[-c(1)],c(2),quantile,0.05, na.rm=T), 7), nsmall = 7)  # [-c(1)] drops iteration number column. Round to 7th dec.
      Scenario <- scx
      ScenName <- sce_names[strtoi(scx)]
      County <- ctx
      Activity <- acx
      Year <- (1:31)
      # combine to one dataframe
      interm_df <- t(rbind(Scenario, ScenName, County, Activity, Year, Median, CI_High, CI_Low))
      out_df <- rbind(out_df, interm_df)
    }
  }
  #out_df <- t(out_df)
  #out_name <- paste0("counties_",sce_names[strtoi(scx)])        # EDIT THIS FOR OUTPUT
}
out_name <- "all_scen_county_activity_results"
write.csv(out_df, paste0(out_dir, "/", out_name,".csv"))

# Bind all outputs into single csv (sum by scenario, county, activity, and sub-activity)
# create empty dataframe that will be added to
out_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(out_df) <- c("Scenario", "ScenName", "County", "Activity", "Subactivity", "Year", "Median", "CI_High", "CI_Low")
# create dataframe of same dimensions as simulations files but is zero values (to replace input files in loop that are NA and NAN)
rx <- data.frame(zeros(1000,32))
for (f in in_list) {
  rc <- read.csv(f)
  # Set variables of each filename part (substring)
  County <- unlist(strsplit(f, "_"))[1]
  Activity <- unlist(strsplit(f, "_"))[2]
  Scenario <- unlist(strsplit(f, "_"))[3]
  ScenName <- sce_names[strtoi(Scenario)]
  Subactivity <- str_sub(unlist(strsplit(f, "_"))[4], 1, -5) # strips ".csv" from the last substring
  Year <- (1:31)
  #[-c(1)]
  # Creates a df of zeros when input is NA or null, OR runs the CI and Median functions.
  if (!is.na(rc[1,2]) & !is.null(rc[1,2])) {
    # Calculate median and 95% CI for each subactivity
    Median <- format(round(apply(rc[-c(1)],c(2),median, na.rm=T), 7), nsmall = 7)         # [-c(1)] drops iteration number column. Round to 7th dec.
    CI_High <- format(round(apply(rc[-c(1)],c(2),quantile,0.95, na.rm=T), 7), nsmall = 7) # [-c(1)] drops iteration number column. Round to 7th dec.
    CI_Low <- format(round(apply(rc[-c(1)],c(2),quantile,0.05, na.rm=T), 7), nsmall = 7)  # [-c(1)] drops iteration number column. Round to 7th dec.
    interm_df <- t(rbind(Scenario, ScenName, County, Activity, Subactivity, Year, Median, CI_High, CI_Low))
  } else {
    Median <- 0
    CI_High <- 0
    CI_Low <- 0
    interm_df <- t(rbind(Scenario, ScenName, County, Activity, Subactivity, Year, Median, CI_High, CI_Low))
  }
  # Combine to one dataframe
  out_df <- rbind(out_df, interm_df)
  #out_df <- t(out_df)
  #out_name <- paste0("counties_",sce_names[strtoi(scx)])        # EDIT THIS FOR OUTPUT
}
out_name <- "all_scen_county_activity_subactivity_results"
write.csv(out_df, paste0(out_dir, "/", out_name,".csv"))
