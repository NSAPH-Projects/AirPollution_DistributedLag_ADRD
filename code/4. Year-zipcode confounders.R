# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: gather census and BRFSS data for every year/zipcode     
#' Inputs: denominator data (2000-2016), BRFSS data (exclude b/c no zip after 2012)                             
#' Outputs: year/zip confounders   
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: <32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(dplyr)
library(NSAPHutils)
library(zipcode)
data(zipcode)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"


##### 1. Create zip-year data from National Causal dataset #####
load(file="/nfs/home/D/dam9096/shared_space/ci3_analysis/National_Causal/National_Causal/2016_temp/aggregate_data.RData")
setDT(aggregate_data)
dt <- aggregate_data[, .(mean_bmi = first(mean_bmi), 
                         smoke_rate = first(smoke_rate), 
                         hispanic = first(hispanic),
                         pct_blk = first(pct_blk), 
                         medhouseholdincome = first(medhouseholdincome), 
                         medianhousevalue = first(medianhousevalue),
                         poverty = first(poverty), 
                         education = first(education), 
                         popdensity = first(popdensity),
                         pct_owner_occ = first(pct_owner_occ), 
                         summer_tmmx = first(summer_tmmx), 
                         winter_tmmx = first(winter_tmmx),
                         summer_rmax = first(summer_rmax), 
                         winter_rmax = first(winter_rmax)), by = .(year, zip)]
dt <- merge(dt, zipcode, by = "zip", all.x = TRUE)
setnames(dt, "state", "statecode")
dt$zip <- as.numeric(dt$zip)

write_fst(dt, paste0(dir_data, "denom/year_zip_confounders.fst"))
