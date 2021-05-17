##############################################################################
#  Authors: Jamie Robertson (The Nature Conservancy, Seattle, Washington, USA)
#              and Kristina Randrup (formerly of TNC and UW).
#  Attribution: This script draws largely from scripts provided from Graves et
#              al.(2020) and Cameron et al. (2017). [see References]
#  Date modified: February 16, 2021 [last modified]
#  Purpose: Calculate CO2e sequestration of Natural Climate Solutions (NCS)
#              pathways by county in Washington, USA.
#  Comments: Outputs of this script are intended for input to a subsequent
#              script, NCS_confidence_intervals.R, also by this author.
#  Disclaimer: Use at your own risk. Some portions of this code have not been
#              been completed and may not function properly. The authors take
#              no responsiblity and will not be liable for use of this script.
#  References: Graves RA, Haugo RD, Holz A, Nielsen-Pincus M, Jones A, Kellogg
#                B, Macdonald C, Popper K, Schindel M. 2020. Potential 
#                greenhouse gas reductions from Natural Climate Solutions in 
#                Oregon, USA. PLOS ONE15:e0230424. 
#                DOI: 10.1371/journal.pone.0230424.
#              Cameron DR, Marvin DC, Remucal JM, Passero MC. 2017. Ecosystem 
#                management and land conservation can substantially contribute
#                to California's climate mitigation goals. Proceedings of the
#                National Academy of Sciences114:12833 LP -12838.
#                DOI: 10.1073/pnas.1707811114.
##############################################################################

######  REQUIRED INPUT FILES  ######
# NCS_model_county_inputs.csv
# county_maximums_all.csv
# IFM_county_inputs.csv
# IFM_inputs_null.csv
# model_carbon_reforestation_low.csv
# model_carbon_reforestation_mod.csv
# model_carbon_reforestation_high.csv
# reforestation_inputs_county.csv
####################################

#load packages that will be throughout the entire script
library(ggplot2)
library(ggforce)
library(grid)
library(truncnorm)
library(tidyverse)
library(generics)
library(gridExtra)
library(gtable)
library(ggpubr)
library(phonTools)

########################
# simulation_functions.R
########################
# calculate the total area of implementation by year, including the mortality rate	
calc_area = function(imp_rate,discount,years){
  acres_in = imp_rate*discount
  acreage_out = acres_in
  for (i in 2:years){
    acreage_out[i] = (acreage_out[i-1] * discount) + acres_in
  }
  return(acreage_out)
}
# calculate the total area of implementation by year, including mortality rate. 
# modified to allow for different implementation each year
calc_area_mod = function(imp_rates,discount,years){
  acres_in = imp_rates[act,1:years,irate]*discount
  acreage_out = acres_in
  for (i in 2:years){
    acreage_out[i]=(acreage_out[i-1]*d_rate) + acres_in[i]
  }
  return(acreage_out)
}
# FOR RESTORATION PATHWAYS WITH A MAXIMUM: calculate the total area of implementation by year, 
# including mortality rate. Modified to allow for different implementation each year and max area
calc_area_mod2 = function(imp_rates,maximum, discount,years){
  acres_in = imp_rates[act,1:years,irate]*discount
  acreage_out = acres_in
  for (i in 2:years){
    acreage_out[i]=ifelse(((acreage_out[i-1]*d_rate) + acres_in[i]) > maximum[act], maximum[act], 
                          (acreage_out[i-1]*d_rate) + acres_in[i])
  }
  indx = which.max(acreage_out) # find which year maximum restoration is reached
  acreage_out[indx:length(acreage_out)] = rep(acreage_out[indx]) # set all following years to the max 
  return(acreage_out)
}
# calculate the total area of implementation by year, including the mortality rate, and store cohort
calc_cohort_area = function(imp_rate,discount,years){
  acreage_out = matrix(NA,years,years)
  cohort = 1:years
  for (c in cohort){
    yrs = seq(c,length(cohort),1) 
    acreage_out[yrs,c] = imp_rate * discount^(yrs-(c-1))
  }
  return(acreage_out)
}
# calculate the total area of implementation by year, including mortality rate, and store cohort; 
# modified to allow for changing imp-rates over years
calc_cohort_area_iter = function(imp_rate,discount,years){
  acreage_out = matrix(NA,years,years)
  cohort = 1:years
  for (c in cohort){
    yrs = seq(c,length(cohort),1) 
    acreage_out[yrs,c] = imp_rate[c] * discount^(yrs-(c-1))
  }
  return(acreage_out)
}
# calculate wildfire area
calc_wildfire_area = function(fire_mean,fire_sd,years){
  acreage_out =matrix(NA,nrow=length(years),ncol=1)
  for(i in 1:length(years)){
    acreage_out[i]=(rtruncnorm(1, mean = fire_mean, sd = fire_sd, a=0))
  }
  return(acreage_out)
}
# grand mean and grand sd
grand.mean <- function(M, N) {weighted.mean(M, N)}
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}
########################
# NCS_ramped_implementation.R
########################
# Overall simulation model parameters #
########################
# set simulation years
years = seq(1,31,1) # represents 2020 - 2051
yr=seq(2020,2050,1)
yrs=length(years)
# set number of monte carlo iterations
iterations = 1000
# set number of scenarios to run
imp_scenarios = 3  # for WA analysis, we ran 3 scenarios (low, moderate, and ambitious)
# set discount rates 
discount_rate_forest = 1 - 0.05 
# (2%, from Cameron et al. 2017, 6% in Fargione 2018, 5% in Latta, 5% in Diaz, 8% in Golub)
discount_rate_reforest = 1 - 0.075 
# high mortality in riparian forest reforestation, 
# 1.5x the discount rate for forests ~7.5% annual discount rate
discount_rate_other = 1 - 0.01 
# 1% background mortality and reversal of other NCS pathways 
########################
# load inputs and maximums and split up by county
########################
setwd("C:/Users/jrobertson/Documents/Workspace/WA_NCS_Rework1/WA County NCS -- user friendly_202005")
ifm_inputs <- read.csv("IFM_county_inputs.csv") 
ifm_inputs$activity <- as.factor(ifm_inputs$activity)
ifm_inputs$type <- as.factor(ifm_inputs$type)
#Adams, Benton, Franklin, and Grant are null -- no timber harvest data
#Chelan, Ferry, and Okanogan not included because >50% high fire risk
#separate out each county
ifm_asotin <-ifm_inputs[1:5,]
ifm_clallam <-ifm_inputs[6:10,]
ifm_clark <-ifm_inputs[11:15,]
ifm_columbia <-ifm_inputs[16:20,]
ifm_cowlitz <-ifm_inputs[21:25,]
ifm_douglas <-ifm_inputs[26:30,]
ifm_garfield <-ifm_inputs[31:35,]
ifm_graysharbor <-ifm_inputs[36:40,]
ifm_island <-ifm_inputs[41:45,]
ifm_jefferson <-ifm_inputs[46:50,]
ifm_king <-ifm_inputs[51:55,]
ifm_kitsap <-ifm_inputs[56:60,]
ifm_kittitas <-ifm_inputs[61:65,]
ifm_klickitat <-ifm_inputs[66:70,]
ifm_lewis <-ifm_inputs[71:75,]
ifm_lincoln <-ifm_inputs[76:80,]
ifm_mason <-ifm_inputs[81:85,]
ifm_pacific <-ifm_inputs[86:90,]
ifm_pendoreille <-ifm_inputs[91:95,]
ifm_pierce <-ifm_inputs[96:100,]
ifm_sanjuan <-ifm_inputs[101:105,]
ifm_skagit <-ifm_inputs[106:110,]
ifm_skamania <-ifm_inputs[111:115,]
ifm_snohomish <-ifm_inputs[116:120,]
ifm_spokane <-ifm_inputs[121:125,]
ifm_stevens <-ifm_inputs[126:130,]
ifm_thurston <-ifm_inputs[131:135,]
ifm_wahkiakum <-ifm_inputs[136:140,]
ifm_wallawalla <-ifm_inputs[141:145,]
ifm_whatcom <-ifm_inputs[146:150,]
ifm_whitman <-ifm_inputs[151:155,]
ifm_yakima <-ifm_inputs[156:160,]
ifm_null <- read.csv("IFM_inputs_null.csv")
ifm_null$activity <- as.factor(ifm_null$activity)
ifm_null$type <- as.factor(ifm_null$type)
# reforestation
reforest_inputs <- read.csv("reforestation_inputs_county.csv")
reforest_inputs$activity <- as.factor(reforest_inputs$activity)
reforest_inputs$activity <- droplevels(reforest_inputs$activity)
reforest_inputs$type <- as.factor(reforest_inputs$type)
#ferry(2), chelan(2), columbia(2), garfield(2), kittitas(1), okanogan(2), yakima(1)
reforest_imprates <- reforest_inputs[,c(2,4,9:14)]
reforest_imprates <- reforest_imprates %>%
  mutate(cv_rates = mean_reveg*2, 
         targetMP2030 = mean_reveg*2, targetMP2040=targetMP2030*2,
         targetMP2050=ifelse(targetMP2040*2>1.0,1.0,targetMP2040*2),
         targetCI2030=mean_reveg*2, targetCI2040=mean_reveg*3)
reforest_ferry <- reforest_imprates[c(1:3),]
reforest_chelan <- reforest_imprates[c(4:6),]
reforest_columbia <- reforest_imprates[c(7:9),]
reforest_garfield <- reforest_imprates[c(10:12),]
reforest_kittitas <- reforest_imprates[c(13:15),]
reforest_okanogan <- reforest_imprates[c(16:18),]
reforest_yakima <- reforest_imprates[c(19:21),]
reforest_null <- reforest_imprates[c(22:24),]
# all other NCS
NCS_inputs = read.csv("NCS_model_county_inputs.csv", stringsAsFactors = F) 
#split up by county
NCS_adams <- NCS_inputs[1:10,]
NCS_asotin <- NCS_inputs[11:20,]
NCS_benton <- NCS_inputs[21:30,]
NCS_chelan <- NCS_inputs[31:40,]
NCS_clallam <- NCS_inputs[41:50,]
NCS_clark <- NCS_inputs[51:60,]
NCS_columbia <- NCS_inputs[61:70,]
NCS_cowlitz <- NCS_inputs[71:80,]
NCS_douglas <- NCS_inputs[81:90,]
NCS_ferry <- NCS_inputs[91:100,]
NCS_franklin <- NCS_inputs[101:110,]
NCS_garfield <- NCS_inputs[111:120,]
NCS_grant <- NCS_inputs[121:130,]
NCS_graysharbor <- NCS_inputs[131:140,]
NCS_island <- NCS_inputs[141:150,]
NCS_jefferson <- NCS_inputs[151:160,]
NCS_king <- NCS_inputs[161:170,]
NCS_kitsap <- NCS_inputs[171:180,]
NCS_kittitas <- NCS_inputs[181:190,]
NCS_klickitat <- NCS_inputs[191:200,]
NCS_lewis <- NCS_inputs[201:210,]
NCS_lincoln <- NCS_inputs[211:220,]
NCS_mason <- NCS_inputs[221:230,]
NCS_okanogan <- NCS_inputs[231:240,]
NCS_pacific <- NCS_inputs[241:250,]
NCS_pendoreille <- NCS_inputs[251:260,]
NCS_pierce <- NCS_inputs[261:270,]
NCS_sanjuan <- NCS_inputs[271:280,]
NCS_skagit <- NCS_inputs[281:290,]
NCS_skamania <- NCS_inputs[291:300,]
NCS_snohomish <- NCS_inputs[301:310,]
NCS_spokane <- NCS_inputs[311:320,]
NCS_stevens <- NCS_inputs[321:330,]
NCS_thurston <- NCS_inputs[331:340,]
NCS_wahkiakum <- NCS_inputs[341:350,]
NCS_wallawalla <- NCS_inputs[351:360,]
NCS_whatcom <- NCS_inputs[361:370,]
NCS_whitman <- NCS_inputs[371:380,]
NCS_yakima <- NCS_inputs[381:390,]
NCS_null <- NCS_inputs[391:400,]
# maximums 
county_maximums <- read.csv("county_maximums_all.csv", stringsAsFactors = F)
########################
# set parameters for county of choice
########################
in_cnty_filename <- ("~/Workspace/WA_NCS_Rework1/WA County NCS -- user friendly_202005/NCS_model_county_inputs.csv")
read_counties <- read.csv(in_cnty_filename, stringsAsFactor=FALSE)
cnty_list <- unlist(lapply(list(read_counties$county), unique))
cnty_i <- 1:length(cnty_list) #set this to the index of the county in cnty_list to restart mid-loop if needed
for (cnty in cnty_list[cnty_i]) {
  print(paste0("working on: ", cnty))
#} #This bracket is for testing the above loop. The full loop ends after writing output csv files later.

county_active <- cnty

#ifm parameters (using null for high fire severity and no-data counties listed)
#ifm_inputs_county <- ifm_null #not using if below conditional statement works
if (county_active %in% list("Adams","Benton","Franklin","Grant","Ferry","Chelan","Okanogan")) {
  ifm_inputs_county <- ifm_null
} else {ifm_inputs_county <- get(paste0("ifm_",casefold(gsub(" ", "", county_active))))
}
#reforest parameters (using null for counties not listed due to no replanting or no-data)
#reforest_imprates_county <- reforest_null #not using if below conditional if-else statement works
if (county_active %in% list("Ferry", "Chelan", "Columbia", "Garfield", "Kittitas", "Okanogan", "Yakima")) {
  reforest_imprates_county <- get(paste0("reforest_",casefold(county_active)))
} else {reforest_imprates_county <- reforest_null
}
reforest_imprates_county$activity <- as.factor(reforest_imprates_county$activity)
reforest_imprates_county$activity <- droplevels(reforest_imprates_county$activity)
reforest_act <- 3
#overall parameters
NCS_inputs_county <- get(paste0("NCS_",casefold(gsub(" ", "", county_active))))
NCS_inputs_county$activity <- as.factor(NCS_inputs_county$activity)
NCS_inputs_county$type <- as.factor(NCS_inputs_county$type)
county_max_active <- county_maximums %>% filter(County == county_active)
########################
# 1. Deferred Harvest (IFM)
########################
ifm_carbonstocks <- ifm_inputs_county[,c(2:4,7:8)] # extract carbon stock data 
ifm_carbonrates <- ifm_inputs_county[,c(2:4,5:6)]
IFM_ImpRates <- ifm_inputs_county[,c(2:3,9:11)] # set up IFM implementation rates
IFM_ImpRates <- IFM_ImpRates %>%
  mutate(startMP = cv/2,
         targetCI2030=baseline*c(.25,.25,.70,.85,.70),
         slopeCI2030=(baseline-targetCI2030)/10)
# set scenario implementation rates
scenario.names <- c("Limited","Moderate","Ambitious")
ifm_imp <- array(NA,dim=c(length(IFM_ImpRates$activity),length(years),imp_scenarios), 
                 dimnames=list(ifm_carbonstocks$activity,years,scenario.names))
# Activities indices: 1 Fed, 2 Other, 3 Private, 4 State, 5 Prvt_seq

# LIMITED
ifm_imp[,10:31,1] <- IFM_ImpRates$cv_scenario # historical variation
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation reduction
  ifm_imp[,cnum,1] <- (IFM_ImpRates$cv_scenario/10)*cnum
}
ifm_imp[,1,2] <- IFM_ImpRates$baseline*IFM_ImpRates$startMP 

# AMBITIOUS
if(county_max_active$EastWest=="West"){
  for (cnum in 2:10){ # ramp up quickly over first 10 years ##WEST AMBITIOUS###
    ifm_imp[c(3,5),cnum,2] <- if_else(IFM_ImpRates$startMP[c(3,5)]*1.2^(cnum-1) > 0.4, 
                                      IFM_ImpRates$baseline[c(3,5)]*0.4, 
                                      IFM_ImpRates$baseline[c(3,5)]*(IFM_ImpRates$startMP[c(3,5)]*1.2^(cnum-1)))
    ifm_imp[4,cnum,2]<-if_else(IFM_ImpRates$startMP[4]*1.2^(cnum-1) > 0.32, 
                               IFM_ImpRates$baseline[4]*0.32, 
                               IFM_ImpRates$baseline[4]*(IFM_ImpRates$startMP[4]*1.2^(cnum-1)))
    ifm_imp[c(1,2),cnum,2] <-if_else(IFM_ImpRates$startMP[c(1,2)]*1.2^(cnum-1) > 0.99, 
                                     IFM_ImpRates$baseline[c(1,2)]*1.0, 
                                     IFM_ImpRates$baseline[c(1,2)]*(IFM_ImpRates$startMP[c(1,2)]*1.2^(cnum-1)))                             
  }
  for (cnum in 11:31){  # maintain constant implementation rate for years 11 - 31
    ifm_imp[c(3,5),cnum,2] <- IFM_ImpRates$baseline[c(3,5)]*0.4
    ifm_imp[4,cnum,2]<- IFM_ImpRates$baseline[4]*0.32
    ifm_imp[c(1,2),cnum,2] <- IFM_ImpRates$baseline[c(1,2)]*1.0
  }
}
if(county_max_active$EastWest=="East"){
  for (cnum in 1:10){ #all east counties to 100% ###EAST AMBITIOUS###
    ifm_imp[c(1:5),cnum,2] <-if_else(IFM_ImpRates$startMP[c(1:5)]*1.2^(cnum-1) > 0.99, 
                                     IFM_ImpRates$baseline[c(1:5)]*1.0, 
                                     IFM_ImpRates$baseline[c(1:5)]*(IFM_ImpRates$startMP[c(1:5)]*1.2^(cnum-1)))
  }
  for (cnum in 11:31){ #maintain implementation 
    ifm_imp[c(1:5),cnum,2] <- IFM_ImpRates$baseline[c(1:5)]*1.0
  }
}

# MODERATE
for (cnum in 1:10){ # constrained implementation 
  ifm_imp[,cnum,3] <- IFM_ImpRates$slopeCI2030*cnum
}
for (cnum in 11:31){
  ifm_imp[,cnum,3] <- IFM_ImpRates$slopeCI2030*10
}
# set up empty data frames to hold simulation results
ifm_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),dim(ifm_carbonstocks)[1]))
ifm_carbon_sum_seq = array(NA, dim=c(imp_scenarios,iterations,length(years),dim(ifm_carbonstocks)[1])) # ADDED FOR TEST
ifm_carbon_sum_sto = array(NA, dim=c(imp_scenarios,iterations,length(years),dim(ifm_carbonstocks)[1])) # ADDED FOR TEST
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:dim(ifm_carbonstocks)[1]) {
    #act = 1
    act_name = ifm_carbonstocks$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = ifm_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = ifm_imp[act,y,2]
      }else {
        imp_rate = ifm_imp[act,y,3]
      }
      d_rate = discount_rate_forest
      ########################
      #	AVOIDED CONVERSION 	  
      ########################
      act_stock_mean = ifm_carbonstocks[which(ifm_carbonstocks$activity %in% act_name), 4]
      act_stock_uc	 = ifm_carbonstocks[which(ifm_carbonstocks$activity %in% act_name), 5]
      act_seq_mean = ifm_carbonrates[which(ifm_carbonrates$activity %in% act_name), 4]
      act_seq_sd	 = ifm_carbonrates[which(ifm_carbonrates$activity %in% act_name), 5]
      imp_per_year = calc_area_mod(ifm_imp, d_rate, yrs)
      # y=1
      
      ### start TEST ###
      ifm_carbon_sum_seq[irate,,y, act] = sapply(imp_per_year[y], function(x) {
        # ongoing sequestration                                      
        (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1)
      })
      ifm_carbon_sum_sto[irate,,y, act] = sapply(imp_per_year[y], function(x) {
        # single event AC portion
        (rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1)
      })
      ### end TEST ###
      
      ifm_carbon_sum[irate,,y, act] = sapply(imp_per_year[y], function(x) {
        # ongoing sequestration                                      
        (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
          # single event AC portion
          (
            rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
          )
      })
    }
  }
  results_out[[irate]] = ifm_carbon_sum
  print(paste('finished IFM implementation scenario ', irate))
} 
# Summarize IFM RESULTS to median and 95% CI for each year, each activity
IFM_CI_High <- apply(ifm_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(IFM_CI_High) <- scenario.names
IFM_median <- apply(ifm_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(IFM_median) <- scenario.names
dimnames(IFM_median)=list(scenario.names,years,ifm_carbonstocks$activity)
IFM_CI_Low <- apply(ifm_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(IFM_CI_Low) <- scenario.names
# reshape data
scen <- rep(c("LowImp","Ambitious", "Moderate"), each = 31)
LowImp_IFM <- cbind(rowSums(cbind(IFM_median[1,1:31,1:5])),
                    rowSums(cbind(IFM_CI_High[1,1:31,1:5])),
                    rowSums(cbind(IFM_CI_Low[1,1:31,1:5])))
colnames(LowImp_IFM)=c("med","upper","lower")
Amb_IFM <- cbind(rowSums(cbind(IFM_median[2,1:31,1:5])),
                 rowSums(cbind(IFM_CI_High[2,1:31,1:5])),
                 rowSums(cbind(IFM_CI_Low[2,1:31,1:5])))
colnames(Amb_IFM)=c("med","upper","lower")
Mod_IFM <- cbind(rowSums(cbind(IFM_median[3,1:31,1:5])),
                 rowSums(cbind(IFM_CI_High[3,1:31,1:5])),
                 rowSums(cbind(IFM_CI_Low[3,1:31,1:5])))
colnames(Mod_IFM)=c("med","upper","lower")
IFM_by_year_all_ownerships <- cbind.data.frame(rbind(LowImp_IFM,Amb_IFM,Mod_IFM),scen)
IFM_by_year_all_ownerships$year <- rep(yr,3)
IFM_by_year_all_ownerships$actname <- rep("Timber Harvest",93)

########################
# 2. Reforestation after wildfires (reforest)
########################
# carbon data
model_carbon_reforest_low <- read.csv("model_carbon_reforestation_low.csv")
model_carbon_reforest_mod <- read.csv("model_carbon_reforestation_mod.csv")
model_carbon_reforest_hi <- read.csv("model_carbon_reforestation_high.csv")
# build array to keep iteration implementation rates for each scenario
iter=1:1000
# first simulate 1000 iterations of fire*30 years
fire_area <- array(NA,dim=c(length(reforest_imprates_county$activity),1000,length(years)), 
                   dimnames=list(reforest_imprates_county$activity,iter,years))
for (i in 1:1000){
  fire_area[1,iter,1:31] <- calc_wildfire_area(reforest_imprates_county$fire_mean[1],
                                               reforest_imprates_county$fire_sd[1],years)
  fire_area[2,iter,1:31] <- calc_wildfire_area(reforest_imprates_county$fire_mean[2],
                                               reforest_imprates_county$fire_sd[2],years)
  fire_area[3,iter,1:31] <- calc_wildfire_area(reforest_imprates_county$fire_mean[3],
                                               reforest_imprates_county$fire_sd[3],years)
}
# set scenario implementation rates
scenario.names <- c("Low","Ambitious","Moderate")
reforest_imp <- array(NA,dim=c(length(reforest_imprates_county$activity),length(years),
                               imp_scenarios), 
                      dimnames=list(reforest_imprates_county$activity,years,scenario.names))
reforest_imp[,10:31,1] <- reforest_imprates_county$cv_rates-reforest_imprates_county$mean_reveg # historical variation
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation (essentially a doubling by 2030)
  reforest_imp[,cnum,1] <- ((reforest_imprates_county$cv_rates-reforest_imprates_county$mean_reveg)/10)*cnum
}
# max potential start year implementation
for (cnum in 1:10){ # ramp up quickly over first 10 years
  reforest_imp[,cnum,2] <- ((reforest_imprates_county$targetMP2030-reforest_imprates_county$mean_reveg)/10)*cnum
}
for (cnum in 11:20){  # maintain constant implementation rate for years 11 - 30
  reforest_imp[,cnum,2] <- (reforest_imprates_county$targetMP2030-reforest_imprates_county$mean_reveg)+
    ((reforest_imprates_county$targetMP2040-reforest_imprates_county$targetMP2030)/10)*(cnum-10)
}
for (cnum in 21:30){  # maintain constant implementation rate for years 11 - 30
  reforest_imp[,cnum,2] <- (reforest_imprates_county$targetMP2040-reforest_imprates_county$mean_reveg)+
    ((reforest_imprates_county$targetMP2050-reforest_imprates_county$targetMP2040)/10)*(cnum-20)
  reforest_imp[,31,2] <- (reforest_imprates_county$targetMP2040-reforest_imprates_county$mean_reveg)+
    ((reforest_imprates_county$targetMP2050-reforest_imprates_county$targetMP2040)/10)*10
}
# constrained implementation scenarios
for (cnum in 1:10){ # ramp up to 2030 target
  reforest_imp[,cnum,3] <- ((reforest_imprates_county$targetCI2030-reforest_imprates_county$mean_reveg)/10)*cnum
}
for (cnum in 11:20){  # ramp up to 2040 
  reforest_imp[,cnum,3] <- (reforest_imprates_county$targetCI2030-reforest_imprates_county$mean_reveg)+
    ((reforest_imprates_county$targetCI2040-reforest_imprates_county$targetCI2030)/10)*(cnum-10)
}
for (cnum in 21:31){  # ramp up to 2050 target
  reforest_imp[,cnum,3] <- (reforest_imprates_county$targetCI2030-reforest_imprates_county$mean_reveg)+
    ((reforest_imprates_county$targetCI2040-reforest_imprates_county$targetCI2030))
}
# set up empty data frames or lists
reforest_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),3))
results_out = list()
# SCENARIO LOOP
for (irate in 1:3) {
  # irate=3
  # ACTIVITY LOOP
  for (act in 1:reforest_act) {
    # act = 1 -- for reforestation rate = 0, don't run that loop
    act_name = reforest_inputs$activity[act]
    imp_rate = reforest_imp[act,1:31,irate]*fire_area[act,,]
    d_rate = discount_rate_forest
    ########################
    # SEQUESTRATION 		 
    ########################
    # get the correct C seq rates
    if (act_name == 'reforestation_low') {
      carbon_rates_reforestation = model_carbon_reforest_low
    }else if (act_name == 'reforestation_moderate') {
      carbon_rates_reforestation = model_carbon_reforest_mod
    }else {
      carbon_rates_reforestation = model_carbon_reforest_hi
    }
    # apply different seq rate and sd depending on age of cohort
    for (i in 1:iterations){
      # i=1
      if (is.na(imp_rate)|is.null(imp_rate)|is.nan(imp_rate)) {
        imp_rate = zeros(imp_rate)
      } 
      i_imp = imp_rate[i,]
      # calculate cohort area
      cohort_area = calc_cohort_area_iter(i_imp,d_rate,length(years))
      cohort_carbon =	apply(cohort_area, 1, function(x) {
        y = x[!is.na(x)]
        act_mean = rev(carbon_rates_reforestation[1:length(y), 2])
        act_sd	 = rev(carbon_rates_reforestation[1:length(y), 3])
        out = sapply(1:length(y), function(p) {
          y[p] * rnorm(1, mean = act_mean[p], sd = act_sd[p])
        })
        return(out)
      })
      reforest_carbon_sum[irate,i, , act] = t(as.numeric(lapply(cohort_carbon,sum))*-1)
    }
  }
  
  results_out[[irate]] = reforest_carbon_sum
  print(paste('finished Reforest implementation scenario ', irate))
}
#####RECIEVING ERROR HERE#####
# "Error in FUN(X[[i]], ...) : invalid 'type' (list) of argument
# In addition: There were 31 warnings (use warnings() to see them)"
##############################

# Summarize post-wildfire reforest RESULTS to median and 95% CI for each year, each activity
reforest_CI_High <- apply(reforest_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(reforest_CI_High) <- scenario.names
reforest_median <- apply(reforest_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(reforest_median) <- scenario.names
reforest_CI_Low <- apply(reforest_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(reforest_CI_Low) <- scenario.names
# reshape data
scen <- rep(c("LowImp","Ambitious", "Moderate"), each = 31)
LowImp_reforest<- cbind(rowSums(cbind(reforest_median[1,1:31,1:3]), na.rm=T),
                        rowSums(cbind(reforest_CI_High[1,1:31,1:3]), na.rm=T),
                        rowSums(cbind(reforest_CI_Low[1,1:31,1:3]), na.rm=T))
colnames(LowImp_reforest)=c("med","upper","lower")
Amb_reforest <- cbind(rowSums(cbind(reforest_median[2,1:31,1:3]), na.rm=T),
                      rowSums(cbind(reforest_CI_High[2,1:31,1:3]), na.rm=T),
                      rowSums(cbind(reforest_CI_Low[2,1:31,1:3]), na.rm=T))
colnames(Amb_reforest)=c("med","upper","lower")
Mod_reforest <- cbind(rowSums(cbind(reforest_median[3,1:31,1:3]), na.rm=T),
                      rowSums(cbind(reforest_CI_High[3,1:31,1:3]), na.rm=T),
                      rowSums(cbind(reforest_CI_Low[3,1:31,1:3]), na.rm=T))
colnames(Mod_reforest)=c("med","upper","lower")
reforest_by_year_all_prodclass <- cbind.data.frame(rbind(LowImp_reforest,
                                                         Amb_reforest,Mod_reforest),scen)
reforest_by_year_all_prodclass$year <- rep(yr,3)
reforest_by_year_all_prodclass$actname <- rep("Replant Wildfires",93)

### OUTPUTS
write.csv(reforest_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_1intensL.csv"))
write.csv(reforest_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_2intensM.csv"))
write.csv(reforest_carbon_sum[1,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_3intensH.csv"))
write.csv(reforest_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_1intensL.csv"))
write.csv(reforest_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_2intensM.csv"))
write.csv(reforest_carbon_sum[2,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_3intensH.csv"))
write.csv(reforest_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_1intensL.csv"))
write.csv(reforest_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_2intensM.csv"))
write.csv(reforest_carbon_sum[3,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_3intensH.csv"))
}


########################
# 3. Avoided Conversion of Forests (ACF)
########################
AC_forests_imprates = NCS_inputs_county[c(1:2),c(1,2,4,5)]
AC_forests_imprates$activity <-droplevels(AC_forests_imprates$activity)
ACF_carbon_stocks = NCS_inputs_county[c(1:2),c(1:2,8:9)]
ACF_carbon_stocks$activity <- droplevels(ACF_carbon_stocks$activity)
ACF_carbon_rates = NCS_inputs_county[c(1:2),c(1:2,6:7)]
ACF_carbon_rates$activity <- droplevels(ACF_carbon_rates$activity)
# set scenario implementation rates
ACF_imp <- array(NA,dim=c(length(AC_forests_imprates$activity),length(years),imp_scenarios), 
                 dimnames=list(AC_forests_imprates$activity,years,scenario.names))
# assume historical variation of 10%, to be reached by 2030
ACF_imp[,10:31,1] <- AC_forests_imprates$baseline*0.10 # historical variation
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  ACF_imp[,cnum,1] <- ((AC_forests_imprates$baseline*0.10)/10)*cnum
}
#max potential
for (cnum in 1:10){ # ramp up quickly over first 10 years to zero conversion
  ACF_imp[,cnum,2] <- (AC_forests_imprates$baseline/10)*cnum
}
ACF_imp[,10:31,2] <- AC_forests_imprates$baseline
# constrained implementation
for (cnum in 1:9){
  ACF_imp[,cnum,3] <- ((AC_forests_imprates$baseline*0.5)/10)*cnum
}
ACF_imp[,10:31,3] <- AC_forests_imprates$baseline*0.5
# set up empty data frames to hold simulation results
ACF_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),
                                 length(AC_forests_imprates$activity)))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:dim(ACF_carbon_stocks)[1]) {
    # act = 1
    act_name = ACF_carbon_stocks$activity[act]
    for (y in 1:length(years)){  ## get implementation rate for each year and activity
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = ACF_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = ACF_imp[act,y,2]
      }else {
        imp_rate = ACF_imp[act,y,3]
      }
      d_rate = discount_rate_forest
      ########################
      #	AVOIDED CONVERSION 	   
      ########################
      act_stock_mean = ACF_carbon_stocks[which(ACF_carbon_stocks$activity %in% act_name), 3]
      act_stock_uc	 = ACF_carbon_stocks[which(ACF_carbon_stocks$activity %in% act_name), 4]
      act_seq_mean = ACF_carbon_rates[act, 3]
      act_seq_sd	 = ACF_carbon_rates[act, 4]
      area_per_year = calc_area_mod(ACF_imp, d_rate, yrs)
      # y=1
      ACF_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
        # ongoing sequestration
        (rtruncnorm(iterations, a=0, b=Inf, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
          # single event AC portion
          (
            rtruncnorm(iterations, a=0,b=Inf,mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
          )
      })
    }
  }
  results_out[[irate]] = ACF_carbon_sum
  print(paste('finished implementation scenario ', irate))
}
# Summarize ACF RESULTS to median and 95% CI for each year, each activity
ACF_CI_High <- apply(ACF_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(ACF_CI_High) <- scenario.names
ACF_median <- apply(ACF_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(ACF_median) <- scenario.names
ACF_CI_Low <- apply(ACF_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(ACF_CI_Low) <- scenario.names
# reshape data
scen <- rep(c("LowImp","Ambitious", "Moderate"), each = 31)
LowImp_ACF <- cbind(rowSums(cbind(ACF_median[1,1:31,1:2])),
                    rowSums(cbind(ACF_CI_High[1,1:31,1:2])),
                    rowSums(cbind(ACF_CI_Low[1,1:31,1:2])))
colnames(LowImp_ACF)=c("med","upper","lower")
Amb_ACF <- cbind(rowSums(cbind(ACF_median[2,1:31,1:2])),
                 rowSums(cbind(ACF_CI_High[2,1:31,1:2])),
                 rowSums(cbind(ACF_CI_Low[2,1:31,1:2])))
colnames(Amb_ACF)=c("med","upper","lower")
Mod_ACF <- cbind(rowSums(cbind(ACF_median[3,1:31,1:2])),
                 rowSums(cbind(ACF_CI_High[3,1:31,1:2])),
                 rowSums(cbind(ACF_CI_Low[3,1:31,1:2])))
colnames(Mod_ACF)=c("med","upper","lower")
ACF_by_year_all <- cbind.data.frame(rbind(LowImp_ACF,Amb_ACF,Mod_ACF),scen)
ACF_by_year_all$year <- rep(yr,3)
ACF_by_year_all$actname <- rep("Forest- Avoided Conversion",93)
########################
# 4. Reforestation of riparian areas (rref)
########################
rref_imprates = NCS_inputs_county[6,c(1:5)]
rref_imprates$activity <-droplevels(rref_imprates$activity)
rref_carbon_stocks = NCS_inputs_county[6,c(1:2,8:9)]
rref_carbon_stocks$activity <- droplevels(rref_carbon_stocks$activity)
rref_carbon_rates = NCS_inputs_county[6,c(1:2,6:7)]
rref_carbon_rates$activity <- droplevels(rref_carbon_rates$activity)
# set scenario implementation rates
scenario.names <- c("HV","MP","CI")
rref_imprates <- rref_imprates %>%
  mutate(target2030MP=baseline*2, target2040MP=target2030MP*2, target2050MP=target2040MP*2)
rref_maxtarget <- county_max_active$rref_max
rref_imp <- array(NA,dim=c(length(rref_imprates$activity),length(years),imp_scenarios), 
                  dimnames=list(rref_imprates$activity,years,scenario.names))
# assume historical variation to be reached by 2030
rref_imp[,10:31,1] <- rref_imprates$baseline*rref_imprates$cv
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  rref_imp[,cnum,1] <- ((rref_imprates$baseline*rref_imprates$cv)/10)*cnum
}
#max potential  ### doubling every 10 years, expressed as + delta from baseline
for (cnum in 1:10){ # ramp up quickly
  rref_imp[,cnum,2] <- ((rref_imprates$target2030MP-rref_imprates$baseline)/10)*cnum
}
for (cnum in 11:20){ 
  rref_imp[,cnum,2] <- ((rref_imprates$target2040MP-rref_imprates$target2030MP)/10)*(cnum-10) + 
    (rref_imprates$target2030MP-rref_imprates$baseline)
}
for (cnum in 21:31){ 
  rref_imp[,cnum,2] <- ((rref_imprates$target2050MP-rref_imprates$target2040MP)/10)*(cnum-20) + 
    (rref_imprates$target2040MP-rref_imprates$baseline)
}
# constrained implementation
for (cnum in 1:10){
  rref_imp[,cnum,3] <- ((rref_imprates$baseline*0.50)/10)*cnum
}
for (cnum in 11:20){
  rref_imp[,cnum,3] <- ((rref_imprates$baseline - rref_imprates$baseline*0.50)/10)*(cnum-10)+
    rref_imprates$baseline*0.50
}
for (cnum in 21:31){
  rref_imp[,cnum,3] <- ((rref_imprates$baseline*1.5 - rref_imprates$baseline)/11)*(cnum-20)+
    rref_imprates$baseline
}
# set up empty data frames or lists
rref_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),1))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:1) {
    # act = 1
    act_name = rref_carbon_rates$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = rref_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = rref_imp[act,y,2]
      }else {
        imp_rate = rref_imp[act,y,3]
      }
      d_rate = discount_rate_reforest
      ########################
      #	SEQUESTRATION 		
      ########################
      seq_mean = rref_carbon_rates[act, 3]
      seq_sd	 = rref_carbon_rates[act, 4]
      area_per_year = calc_area_mod2(rref_imp, rref_maxtarget, d_rate, length(years))
      rref_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
        rnorm(iterations, mean = seq_mean, sd = seq_sd) * x * -1
      })
    }
  }
  results_out[[irate]] = reforest_carbon_sum
  print(paste('finished Reforest implementation scenario ', irate))
} 
# Summarize rref RESULTS to median and 95% CI for each year, each activity
rref_CI_High <- apply(rref_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(rref_CI_High) <- scenario.names
rref_median <- apply(rref_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(rref_median) <- scenario.names
rref_CI_Low <- apply(rref_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(rref_CI_Low) <- scenario.names
# reshape data
scen <- rep(c("LowImp","Ambitious", "Moderate"), each = 31)
LowImp_rref <- cbind(rowSums(cbind(rref_median[1,1:31,])),
                     rowSums(cbind(rref_CI_High[1,1:31,])),
                     rowSums(cbind(rref_CI_Low[1,1:31,])))
colnames(LowImp_rref)=c("med","upper","lower")
Amb_rref <- cbind(rowSums(cbind(rref_median[2,1:31,])),
                  rowSums(cbind(rref_CI_High[2,1:31,])),
                  rowSums(cbind(rref_CI_Low[2,1:31,])))
colnames(Amb_rref)=c("med","upper","lower")
Mod_rref <- cbind(rowSums(cbind(rref_median[3,1:31,])),
                  rowSums(cbind(rref_CI_High[3,1:31,])),
                  rowSums(cbind(rref_CI_Low[3,1:31,])))
colnames(Mod_rref)=c("med","upper","lower")
rref_by_year_all <- cbind.data.frame(rbind(LowImp_rref,Amb_rref,Mod_rref),scen)
rref_by_year_all$year <- rep(yr,3)
rref_by_year_all$actname <- rep("Riparian Reforestation",93)

########################
# 5. Sage-steppe avoided conversion and restoration (sage)
########################
sage_inputs <- NCS_inputs_county[c(3:4),]
sage_inputs$activity <-droplevels(sage_inputs$activity)
sage_carbon_stocks = sage_inputs[,c(1:3,8:9)]
sage_carbon_rates = sage_inputs[,c(1:3,6:7)]
# set scenario implementation rates
sage_max <- c(county_max_active$sage_ac_max,county_max_active$sage_restor_max)
scenario.names <- c("HV","MP","CI")
sage_inputs <- sage_inputs %>%
  mutate(max=sage_max)
sage_imp <- array(NA,dim=c(length(sage_inputs$activity),length(years),imp_scenarios), 
                  dimnames=list(sage_inputs$activity,years,scenario.names))
# assume historical variation of 10%, to be reached by 2030
sage_imp[,10:31,1] <- sage_inputs$baseline*0.10 # historical variation
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  sage_imp[,cnum,1] <- ((sage_inputs$baseline*0.10)/10)*cnum
}
# maximum potential 
for (cnum in 1:9) {
  sage_imp[1,cnum,2] <- sage_inputs$baseline[1]/10*cnum
  sage_imp[2,cnum,2] <- (sage_inputs$baseline[2]*3-sage_inputs$baseline[2])/10*cnum
}
sage_imp[1,10:31,2] <- sage_inputs$baseline[1]
sage_imp[2,10:31,2] <- sage_inputs$baseline[2]*3-sage_inputs$baseline[2]
# constrained implementation
for (cnum in 1:9) {
  sage_imp[1,cnum,3] <- sage_inputs$baseline[1]*0.1
  sage_imp[2,cnum,3] <- (sage_inputs$baseline[2])/10*cnum
}
for (cnum in 10:31) {
  sage_imp[1,cnum,3] <- sage_inputs$baseline[1]*0.2
  sage_imp[2,cnum,3] <- (sage_inputs$baseline[2])/20*(cnum-10)+sage_inputs$baseline[2]
}
# set up empty data frames or lists
sage_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),2))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:1) {
    # act = 1
    act_name =sage_carbon_rates$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = sage_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = sage_imp[act,y,2]
      }else {
        imp_rate = sage_imp[act,y,3]
      }
      d_rate = discount_rate_other
      ########################
      # SEQUESTRATION 		   
      ########################
      if (sage_carbon_rates$type[act] == 'seq') {
        act_mean = sage_carbon_rates[act, 4]
        act_sd	 = sage_carbon_rates[act, 5]
        #     
        area_per_year = calc_area_mod(sage_imp, d_rate, length(years))
        #     
        sage_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          rnorm(iterations, mean = act_mean, sd = act_sd) * x * -1
        })
      }
      ########################
      # AVOIDED CONVERSION 	
      ########################
      if (sage_carbon_rates$type[act] == 'AC') {
        act_stock_mean = sage_carbon_stocks[which(sage_carbon_stocks$activity %in% act_name), 4]
        act_stock_uc	 = sage_carbon_stocks[which(sage_carbon_stocks$activity %in% act_name), 5]
        act_seq_mean = sage_carbon_rates[act, 4]
        act_seq_sd	 = sage_carbon_rates[act, 5]
        area_per_year = calc_area_mod(sage_imp, d_rate, yrs)
        # y=1
        sage_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # single event AC portion
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
            )
        })
      }
      ##### End Sage Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = sage_carbon_sum
  print(paste('finished sage implementation scenario ', irate))
} 
# Summarize sage RESULTS to median and 95% CI for each year, each activity
sage_CI_High <- apply(sage_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(sage_CI_High) <- scenario.names
sage_CI_High[is.na(sage_CI_High)] <- 0

#########ADD BELOW
#sage_carbon_sum[is.na(sage_carbon_sum)] <- 0
#high_sage <- sage_carbon_sum[1,,,1]+sage_carbon_sum[1,,,2]
#filename <- paste("NEWOUTPUTS/highsage",county_active,"csv",sep =".")
#write.csv(high_sage,filename,row.names=FALSE)
#########ADD ABOVE

sage_median <- apply(sage_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(sage_median) <- scenario.names
sage_median[is.na(sage_median)] <- 0
sage_CI_Low <- apply(sage_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(sage_CI_Low) <- scenario.names
sage_CI_Low[is.na(sage_CI_Low)] <- 0
# reshape data
scen <- rep(c("LowImp","Ambitious","Moderate"), each = 31)
LowImp_sage <- cbind(rowSums(cbind(sage_median[1,1:31,1:2])),
                     rowSums(cbind(sage_CI_High[1,1:31,1:2])),
                     rowSums(cbind(sage_CI_Low[1,1:31,1:2])))
colnames(LowImp_sage)=c("med","upper","lower")
Amb_sage <- cbind(rowSums(cbind(sage_median[2,1:31,1:2])),
                  rowSums(cbind(sage_CI_High[2,1:31,1:2])),
                  rowSums(cbind(sage_CI_Low[2,1:31,1:2])))
colnames(Amb_sage)=c("med","upper","lower")
Mod_sage <- cbind(rowSums(cbind(sage_median[3,1:31,1:2])),
                  rowSums(cbind(sage_CI_High[3,1:31,1:2])),
                  rowSums(cbind(sage_CI_Low[3,1:31,1:2])))
colnames(Mod_sage)=c("med","upper","lower")
sage_by_year_all <- cbind.data.frame(rbind(LowImp_sage,Amb_sage,Mod_sage),scen)
sage_by_year_all$year <- rep(yr,3)
sage_by_year_all$actname <- rep("Sagebrush-steppe pathways",93)
########################
# 6. Restoration of tidal wetlands (tide)
########################
tide_inputs <- NCS_inputs_county[5,]
tide_inputs$activity <-droplevels(tide_inputs$activity)
tide_max <- county_max_active$tide_max
# set scenario implementation  rates
scenario.names <- c("HV","MP","CI")
tide_inputs <- tide_inputs %>%
  mutate(max=tide_max)
tide_imp <- array(NA,dim=c(length(tide_inputs$activity),length(years),imp_scenarios), 
                  dimnames=list(tide_inputs$activity,years,scenario.names))
# assume historical variation 
tide_imp[,10:31,1] <- tide_inputs$baseline # increase by 50 ha in 2030 (>100% h.v.)
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  tide_imp[,cnum,1] <- (tide_inputs$baseline*tide_inputs$cv*cnum)
}
# ambitious implementation -- max is not reached -- rate doubled every few years 
for (cnum in 1:5) {
  tide_imp[,cnum,2] <- (tide_inputs$baseline)/5*cnum 
}
for (cnum in 6:15) {
  tide_imp[,cnum,2] <- (tide_inputs$baseline*4-tide_inputs$baseline*2)/10*(cnum-5)+
    tide_inputs$baseline
}
for (cnum in 16:26) {
  tide_imp[,cnum,2] <- (tide_inputs$baseline*8-tide_inputs$baseline*4)/15*(cnum-15)+
    tide_inputs$baseline*3
}
for (cnum in 27:31) {
  tide_imp[,cnum,2] <- (tide_inputs$baseline*16-tide_inputs$baseline*8)/20*(cnum-26)+
    tide_inputs$baseline*7
}
#constrained implementation
for (cnum in 1:10) {
  tide_imp[,cnum,3] <- (2*tide_inputs$baseline/10*cnum)
}
tide_imp[,11:31,3] <- 2*tide_inputs$baseline
# check maximum
sum(tide_imp[,1:31,2]) < tide_max
# set up empty data frames or lists
tide_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),1))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:length(tide_inputs$activity)) {
    # act = 1
    act_name =tide_inputs$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = tide_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = tide_imp[act,y,2]
      }else {
        imp_rate = tide_imp[act,y,3]
      }
      d_rate = discount_rate_other
      ########################
      # AVOIDED CONVERSION 	  
      ########################
      if (tide_inputs$type[act] == 'AC') {
        act_stock_mean = tide_inputs[which(tide_inputs$activity %in% act_name), 8]
        act_stock_uc	 = tide_inputs[which(tide_inputs$activity %in% act_name), 9]
        act_seq_mean = tide_inputs[act, 6]
        act_seq_sd	 = tide_inputs[act, 7]
        area_per_year = calc_area_mod2(tide_imp, tide_max, d_rate, yrs)
        # y=1
        tide_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # for tidal wetlands AC portion is a per ha per year rate
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * x * -1
            )
        })
      }
      ##### End tide Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = tide_carbon_sum
  print(paste('finished tide implementation scenario ', irate))
}
# Summarize tide RESULTS to median and 95% CI for each year, each activity
tide_CI_High <- apply(tide_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(tide_CI_High) <- scenario.names
tide_median <- apply(tide_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(tide_median) <- scenario.names
tide_CI_Low <- apply(tide_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(tide_CI_Low) <- scenario.names
# reshape data
scen <- rep(c("LowImp","Ambitious", "Moderate"), each = 31)
LowImp_tide <- cbind(rowSums(cbind(tide_median[1,1:31,])),
                     rowSums(cbind(tide_CI_High[1,1:31,])),
                     rowSums(cbind(tide_CI_Low[1,1:31,])))
colnames(LowImp_tide)=c("med","upper","lower")
Amb_tide <- cbind(rowSums(cbind(tide_median[2,1:31,])),
                  rowSums(cbind(tide_CI_High[2,1:31,])),
                  rowSums(cbind(tide_CI_Low[2,1:31,])))
colnames(Amb_tide)=c("med","upper","lower")
Mod_tide <- cbind(rowSums(cbind(tide_median[3,1:31,])),
                  rowSums(cbind(tide_CI_High[3,1:31,])),
                  rowSums(cbind(tide_CI_Low[3,1:31,])))
colnames(Mod_tide)=c("med","upper","lower")
tide_by_year_all <- cbind.data.frame(rbind(LowImp_tide,Amb_tide,Mod_tide),scen)
tide_by_year_all$year <- rep(yr,3)
tide_by_year_all$actname <- rep("Tidal Wetlands",93)
########################
# 7. Grassland avoided converstion (grass)
########################
grass_inputs <- NCS_inputs_county[7,]
grass_inputs$activity <-droplevels(grass_inputs$activity)
# set implementation rates
scenario.names <- c("HV","MP","CI")
grass_imp <- array(NA,dim=c(length(grass_inputs$activity),length(years),imp_scenarios), 
                   dimnames=list(grass_inputs$activity,years,scenario.names))
# assume historical variation 
grass_imp[,10:31,1] <- grass_inputs$baseline*0.10 # 10% of baseline
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  grass_imp[,cnum,1] <- (grass_inputs$baseline*0.10)/10*cnum
}
# maximum potential 
for (cnum in 1:10) {
  grass_imp[,cnum,2] <- (grass_inputs$baseline)/10*cnum # decrease by equivalent of baseline by 2030
}
grass_imp[,11:31,2] <- grass_inputs$baseline # 100% of baseline (i.e. no conversion after 2030)
# constrained implementation - decreases conversion by 50% by 2030 and then keeps that
for (cnum in 1:10) {
  grass_imp[,cnum,3] <- ((grass_inputs$baseline)*0.5)/10*cnum 
}
for (cnum in 11:31) {
  grass_imp[,cnum,3] <- ((grass_inputs$baseline)*0.5)/21*(cnum-10) + grass_inputs$baseline*0.5
}
# set up empty data frames or lists
grass_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),1))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:length(grass_inputs$activity)) {
    # act = 1
    act_name =grass_inputs$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = grass_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = grass_imp[act,y,2]
      }else {
        imp_rate = grass_imp[act,y,3]
      }
      d_rate = discount_rate_other
      ########################
      # AVOIDED CONVERSION 	
      ########################
      if (grass_inputs$type[act] == 'AC') {
        act_stock_mean = grass_inputs[which(grass_inputs$activity %in% act_name), 8]
        act_stock_uc	 = grass_inputs[which(grass_inputs$activity %in% act_name), 9]
        act_seq_mean = grass_inputs[act, 6]
        act_seq_sd	 = grass_inputs[act, 7]
        area_per_year = calc_area_mod(grass_imp, d_rate, yrs)
        # y=1
        grass_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # single event AC portion
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
            )
        })
      }
      ##### End grass Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = grass_carbon_sum
  print(paste('finished grass implementation scenario ', irate))
} 
# Summarize grass RESULTS to median and 95% CI for each year, each activity
grass_CI_High <- apply(grass_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(grass_CI_High) <- scenario.names
grass_median <- apply(grass_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(grass_median) <- scenario.names
grass_CI_Low <- apply(grass_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(grass_CI_Low) <- scenario.names
# reshape data
LowImp_grass <- cbind(rowSums(cbind(grass_median[1,1:31,])),
                      rowSums(cbind(grass_CI_High[1,1:31,])),
                      rowSums(cbind(grass_CI_Low[1,1:31,])))
colnames(LowImp_grass)=c("med","upper","lower")
Amb_grass <- cbind(rowSums(cbind(grass_median[2,1:31,])),
                   rowSums(cbind(grass_CI_High[2,1:31,])),
                   rowSums(cbind(grass_CI_Low[2,1:31,])))
colnames(Amb_grass)=c("med","upper","lower")
Mod_grass <- cbind(rowSums(cbind(grass_median[3,1:31,])),
                   rowSums(cbind(grass_CI_High[3,1:31,])),
                   rowSums(cbind(grass_CI_Low[3,1:31,])))
colnames(Mod_grass)=c("med","upper","lower")
grass_by_year_all <- cbind.data.frame(rbind(LowImp_grass,Amb_grass,Mod_grass),scen)
grass_by_year_all$year <- rep(yr,3)
grass_by_year_all$actname <- rep("Grassland - AC",93)
########################
# 8. Agricultural pathways (ag)
########################
ag_inputs <- NCS_inputs_county[c(8:10),]
ag_inputs$activity <-droplevels(ag_inputs$activity)
ag_inputs$max <- c(county_max_active$ag_cc_max, county_max_active$ag_notill_max,NA)
# set scenario implementation rates
scenario.names <- c("HV","MP","CI")
ag_imp <- array(NA,dim=c(length(ag_inputs$activity),length(years),imp_scenarios), 
                dimnames=list(ag_inputs$activity,years,scenario.names))
# assume historical variation 
for (cnum in 1:9){  # ramp up period of 10 years to reach the historical variation
  ag_imp[1,cnum,1] <- (ag_inputs$baseline[1]*0.40)/10*cnum
  ag_imp[2,cnum,1] <- (ag_inputs$baseline[2]*0.30)/10*cnum
  ag_imp[3,cnum,1] <- (ag_inputs$baseline[3]*(ag_inputs$cv[3]*0.5))/10*cnum 
  # for N reduction, we assume HV reduction on 50% of acres
}
for (cnum in 10:31){ # historical variation
  ag_imp[1,cnum,1] <- (ag_inputs$baseline[1]*0.40)
  ag_imp[2,cnum,1] <- (ag_inputs$baseline[2]*0.30)
  ag_imp[3,cnum,1] <- (ag_inputs$baseline[3])*(ag_inputs$cv[3]*0.5)
}
# maximum potential
for (cnum in 1:31 ){ # linear increases in cover crop and no-till to reach max by 2050 
  # (max is 50% of cropland and all tilled cropland)
  ag_imp[1,cnum,2] <- ((as.numeric(ag_inputs$max[1])*0.5)-ag_inputs$baseline[1])/30*cnum
  ag_imp[2,cnum,2] <- (as.numeric(ag_inputs$max[2])-ag_inputs$baseline[2])/30*cnum
}
for (cnum in 1:10) { # linear decrease to reach target of 40% reduction by 2030
  ag_imp[3,cnum,2] <- ((ag_inputs$baseline[3]*0.40)/10)*cnum
}
ag_imp[3,11:31,2] <- (ag_inputs$baseline[3]*0.40)
# constrained implementation
for (cnum in 1:10){ # cover crops increase by 150% and no-till doubled by 2030
  ag_imp[1,cnum,3] <- (ag_inputs$baseline[1]*1.5)/10*cnum
  ag_imp[2,cnum,3] <- ag_inputs$baseline[2]/10*cnum
}
for (cnum in 11:31) {# cover crops quadruple by 2050, no-till increase by 150%
  ag_imp[1,cnum,3] <- ((ag_inputs$baseline[1]*4-ag_inputs$baseline[1]*1.5)/20*(cnum-10))+
    (ag_inputs$baseline[1]*1.5)
  ag_imp[2,cnum,3] <- ag_inputs$baseline[2]*1.5
}
for (cnum in 1:10) {# linear decrease to reach target of 25% on 40% of acres (0.25*0.40) by 2030
  # decrease by 0.25 on 100% of acres by 2050
  ag_imp[3,cnum,3] <- ((ag_inputs$baseline[3])*(0.25*0.4)/10)*cnum
}
for (cnum in 11:31) {
  ag_imp[3,cnum,3] <- (((ag_inputs$baseline[3]*0.25-ag_inputs$baseline[3]*0.25*0.4)/21)*
                         (cnum-10))+(ag_inputs$baseline[3]*0.25*0.40)
}
# set up empty data frames or lists
ag_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),length(ag_inputs$activity)))
results_out = list()
# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:length(ag_inputs$activity)) {
    # act = 2
    act_name =ag_inputs$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = ag_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = ag_imp[act,y,2]
      }else {
        imp_rate = ag_imp[act,y,3]
      }
      d_rate = discount_rate_other
      ########################
      # SEQUESTRATION 	
      ########################
      # Cover Crops 
      if (act_name == 'cover_crop_use') {
        act_mean = ag_inputs[act, 6]
        act_sd	 = ag_inputs[act, 7]
        ag_carbon_sum[irate, ,y, act] = rnorm(iterations, mean = act_mean, sd = act_sd) *imp_rate*d_rate* -1
      }
      # No-till Agriculture 
      if (act_name == 'no_till_ag') {
        act_min = ag_inputs[act,6]
        act_max = ag_inputs[act,7]
        ag_carbon_sum[irate, ,y,act] = runif(iterations,min=act_min, max=act_max)*imp_rate*d_rate* -1 # ongoing sequestration
      }
      ########################
      # AVOIDED CONVERSION 	 
      ########################
      if (act_name == 'N_management') {
        act_min = ag_inputs[act,8]
        act_max = ag_inputs[act,9]
        area_per_year = calc_area_mod(ag_imp, d_rate, length(years))
        #
        ag_carbon_sum[irate, ,y,act] = runif(iterations,min=act_min, max=act_max)*imp_rate * d_rate * -1 # avoided emissions
      }
      ##### End ag Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = ag_carbon_sum
  print(paste('finished ag implementation scenario ', irate))
} 
# Summarize ag RESULTS to median and 95% CI for each year, each activity
ag_CI_High <- apply(ag_carbon_sum,c(1,3,4),quantile,0.95, na.rm=T)/1e6
rownames(ag_CI_High) <- scenario.names
ag_median <- apply(ag_carbon_sum,c(1,3,4),median, na.rm=T)/1e6
rownames(ag_median) <- scenario.names
ag_CI_Low <- apply(ag_carbon_sum,c(1,3,4),quantile,0.05, na.rm=T)/1e6
rownames(ag_CI_Low) <- scenario.names
# Cover Crops
CC_lowimp <- cbind(ag_median[1,1:31,1],ag_CI_High[1,1:31,1],ag_CI_Low[1,1:31,1])
CC_amb <- cbind(ag_median[2,1:31,1],ag_CI_High[2,1:31,1],ag_CI_Low[2,1:31,1])
CC_mod <- cbind(ag_median[3,1:31,1],ag_CI_High[3,1:31,1],ag_CI_Low[3,1:31,1])
CC_by_year <- cbind.data.frame(rbind(CC_lowimp,CC_amb,CC_mod),scen)
CC_by_year$year <- rep(yr,3)
CC_by_year$actname <- rep("Cover Crops",93)
# No-till
NT_lowimp <- cbind(ag_median[1,1:31,2],ag_CI_High[1,1:31,2],ag_CI_Low[1,1:31,2])
NT_amb <- cbind(ag_median[2,1:31,2],ag_CI_High[2,1:31,2],ag_CI_Low[2,1:31,2])
NT_mod <- cbind(ag_median[3,1:31,2],ag_CI_High[3,1:31,2],ag_CI_Low[3,1:31,2])
NT_by_year <- cbind.data.frame(rbind(NT_lowimp,NT_amb,NT_mod),scen)
NT_by_year$year <- rep(yr,3)
NT_by_year$actname <- rep("No-Till",93)
# Nitrogen
Nmgmt_lowimp <- cbind(ag_median[1,1:31,3],ag_CI_High[1,1:31,3],ag_CI_Low[1,1:31,3])
Nmgmt_amb <- cbind(ag_median[2,1:31,3],ag_CI_High[2,1:31,3],ag_CI_Low[2,1:31,3])
Nmgmt_mod <- cbind(ag_median[3,1:31,3],ag_CI_High[3,1:31,3],ag_CI_Low[3,1:31,3])
Nmgmt_by_year <- cbind.data.frame(rbind(Nmgmt_lowimp,Nmgmt_amb,Nmgmt_mod),scen)
Nmgmt_by_year$year <- rep(yr,3)
Nmgmt_by_year$actname <- rep("N Mgmt",93)
ag_by_year_all_act <- rbind(CC_by_year, NT_by_year,Nmgmt_by_year)
colnames(ag_by_year_all_act)[1:3]=c("med","upper","lower")
# reshape data
LowImp_ag <- cbind(rowSums(cbind(ag_median[1,1:31,])),
                   rowSums(cbind(ag_CI_High[1,1:31,])),
                   rowSums(cbind(ag_CI_Low[1,1:31,])))
colnames(LowImp_ag)=c("med","upper","lower")
Amb_ag <- cbind(rowSums(cbind(ag_median[2,1:31,])),
                rowSums(cbind(ag_CI_High[2,1:31,])),
                rowSums(cbind(ag_CI_Low[2,1:31,])))
colnames(Amb_ag)=c("med","upper","lower")
Mod_ag <- cbind(rowSums(cbind(ag_median[3,1:31,])),
                rowSums(cbind(ag_CI_High[3,1:31,])),
                rowSums(cbind(ag_CI_Low[3,1:31,])))
colnames(Mod_ag)=c("med","upper","lower")
ag_by_year_all <- cbind.data.frame(rbind(LowImp_ag,Amb_ag,Mod_ag),scen)
ag_by_year_all$year <- rep(yr,3)
ag_by_year_all$actname <- rep("Agricultural Pathways",93)

########################
# Write out iterations to aggregate for statewide CIs
########################
write.csv(ifm_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_1_1feder.csv"))
write.csv(ifm_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_1_2other.csv"))
write.csv(ifm_carbon_sum[1,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_1_3priva.csv"))
write.csv(ifm_carbon_sum[1,,,4]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_1_4state.csv"))
write.csv(ifm_carbon_sum[1,,,5]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_1_5prseq.csv"))
write.csv(ifm_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_2_1feder.csv"))
write.csv(ifm_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_2_2other.csv"))
write.csv(ifm_carbon_sum[2,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_2_3priva.csv"))
write.csv(ifm_carbon_sum[2,,,4]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_2_4state.csv"))
write.csv(ifm_carbon_sum[2,,,5]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_2_5prseq.csv"))
write.csv(ifm_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_3_1feder.csv"))
write.csv(ifm_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_3_2other.csv"))
write.csv(ifm_carbon_sum[3,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_3_3priva.csv"))
write.csv(ifm_carbon_sum[3,,,4]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_3_4state.csv"))
write.csv(ifm_carbon_sum[3,,,5]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ifm_3_5prseq.csv"))
write.csv(reforest_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_1intensL.csv"))
write.csv(reforest_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_2intensM.csv"))
write.csv(reforest_carbon_sum[1,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_1_3intensH.csv"))
write.csv(reforest_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_1intensL.csv"))
write.csv(reforest_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_2intensM.csv"))
write.csv(reforest_carbon_sum[2,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_2_3intensH.csv"))
write.csv(reforest_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_1intensL.csv"))
write.csv(reforest_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_2intensM.csv"))
write.csv(reforest_carbon_sum[3,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_reforest_3_3intensH.csv"))
write.csv(ACF_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_1_1f2rur.csv"))
write.csv(ACF_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_1_2f2urb.csv"))
write.csv(ACF_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_2_1f2rur.csv"))
write.csv(ACF_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_2_2f2urb.csv"))
write.csv(ACF_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_3_1f2rur.csv"))
write.csv(ACF_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ACF_3_2f2urb.csv"))
write.csv(rref_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_rref_1_1rref.csv"))
write.csv(rref_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_rref_2_1rref.csv"))
write.csv(rref_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_rref_3_1rref.csv"))
write.csv(sage_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_1_1s2ag.csv"))
write.csv(sage_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_1_2ag2s.csv"))
write.csv(sage_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_2_1s2ag.csv"))
write.csv(sage_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_2_2ag2s.csv"))
write.csv(sage_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_3_1s2ag.csv"))
write.csv(sage_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_sage_3_2ag2s.csv"))
write.csv(tide_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_tide_1_1tide.csv"))
write.csv(tide_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_tide_2_1tide.csv"))
write.csv(tide_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_tide_3_1tide.csv"))
write.csv(grass_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_grass_1_1grass.csv"))
write.csv(grass_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_grass_2_1grass.csv"))
write.csv(grass_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_grass_3_1grass.csv"))
write.csv(ag_carbon_sum[1,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_1_1cc.csv"))
write.csv(ag_carbon_sum[1,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_1_2nt.csv"))
write.csv(ag_carbon_sum[1,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_1_3nm.csv"))
write.csv(ag_carbon_sum[2,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_2_1cc.csv"))
write.csv(ag_carbon_sum[2,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_2_2nt.csv"))
write.csv(ag_carbon_sum[2,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_2_3nm.csv"))
write.csv(ag_carbon_sum[3,,,1]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_3_1cc.csv"))
write.csv(ag_carbon_sum[3,,,2]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_3_2nt.csv"))
write.csv(ag_carbon_sum[3,,,3]/1e6, paste0("ANALYSIS_OUTPUTS/", county_active, "_ag_3_3nm.csv"))

#} # closing bracket for loop running through each county at ~line 410

rm(list=ls()) #remove all variables before running new county
