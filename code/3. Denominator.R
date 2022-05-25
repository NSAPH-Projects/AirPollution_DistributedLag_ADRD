# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: create individual denominator data for years 2009-2016     
#' Inputs: denominator data (2000-2016)                                 
#' Outputs: cleaned yearly qid enrollment and individual data         
#' Author: Daniel Mork                                                         
#' Last updated: Mar 21, 2022                                                  
#' Memory to run: 32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(dplyr)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"
dir_denominator <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/"

##### 1. Extract individual data 2009-2016 #####
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)[11:18]
# example <- read_fst(f[1], from = 1, to = 10000, as.data.table = TRUE)
# names(example)
# example
# > names(example)
# [1] "zip"                          "year"                         "qid"                          "dodflag"                     
# [5] "bene_dod"                     "sex"                          "race"                         "age"                         
# [9] "hmo_mo"                       "hmoind"                       "statecode"                    "latitude"                    
# [13] "longitude"                    "dual"                         "death"                        "dead"                        
# [17] "entry_age"                    "entry_year"                   "entry_age_break"              "followup_year"               
# [21] "followup_year_plus_one"       "pm25_ensemble"                "pm25_no_interp"               "pm25_nn"                     
# [25] "ozone"                        "ozone_no_interp"              "zcta"                         "poverty"                     
# [29] "popdensity"                   "medianhousevalue"             "pct_blk"                      "medhouseholdincome"          
# [33] "pct_owner_occ"                "hispanic"                     "education"                    "population"                  
# [37] "zcta_no_interp"               "poverty_no_interp"            "popdensity_no_interp"         "medianhousevalue_no_interp"  
# [41] "pct_blk_no_interp"            "medhouseholdincome_no_interp" "pct_owner_occ_no_interp"      "hispanic_no_interp"          
# [45] "education_no_interp"          "population_no_interp"         "smoke_rate"                   "mean_bmi"                    
# [49] "smoke_rate_no_interp"         "mean_bmi_no_interp"           "amb_visit_pct"                "a1c_exm_pct"                 
# [53] "amb_visit_pct_no_interp"      "a1c_exm_pct_no_interp"        "tmmx"                         "rmax"                        
# [57] "pr"                           "cluster_cat"                  "fips_no_interp"               "fips"                        
# [61] "summer_tmmx"                  "summer_rmax"                  "winter_tmmx"                  "winter_rmax" 

myvars <- c("qid", "year", "zip", "sex", "race", "age", "dual", "dead", 
            "hmo_mo", "fips")
f_idx <- 1
for (yr in 2009:2016) {
  cat("\nReading", f[f_idx])
  dt <- read_fst(f[f_idx], columns = myvars, as.data.table = TRUE)
  setkey(dt, year, zip, qid)
  dt <- dt[year >= 2009 & # 2009-2016
             hmo_mo == 0 & # FFS only
             zip >= 501 & # remove invalid zips
             race != 0 & # remove unknown race
             age >= 64 & # remove invalid age
             sex != 0 # remove unknown sex
           ][, race2 := ifelse(race == 1, "wht", # recode race
                               ifelse(race == 2, "blk",
                                      ifelse(race == 5, "his", "oth")))
             ][, sexM := ifelse(sex == 1, 1, 0)][] # recode sex
  dt[, race := NULL]
  setnames(dt, "race2", "race")
  dt <- unique(dt[complete.cases(dt)], by = c("year", "zip", "qid"))
  cat("\n", dt[,.N], "records")
  write_fst(dt, paste0(dir_data, "denom/qid_data_", yr, ".fst"))
  rm(dt); gc()
  f_idx <- f_idx + 1
}
