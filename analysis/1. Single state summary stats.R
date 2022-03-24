# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: single state distributed lag analysis
#' Inputs: merged denominator and exposure data                                
#' Outputs:       
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
state <- "MA"
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "primary" # primary or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
qid_dat[entry == 2000, sum(dead), by = year]
qid_dat[entry == 2000, uniqueN(qid), by = year]
qid_dat[entry == 2000, mean(first_hosp), by = year]
qid_dat[entry == 2000, sum(first_hosp)/uniqueN(qid)]
qid_dat[entry == 2000, sum(first_hosp), by = year]
qid_dat[entry == 2000 & first_hosp == TRUE, mean(age)]
qid_dat[entry == 2000 & dead == TRUE, mean(age)]
qid_dat[entry <= 2007 & year == 2016, .N]
qid_dat[entry <= 2007 & year == 2016, mean(first_hosp)]
