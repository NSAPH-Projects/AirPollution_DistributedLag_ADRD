# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: match exposures to FFS beneficiary zipcode of residence               
#' Inputs: denominator files, QD exposure data, meteorological data              
#' Outputs: exposures in grid by qid (rows) and year (cols)                    
#' Author: Daniel Mork                                                         
#' Last updated: Mar 08, 2022                                                  
#' Memory to run: 96 GB
# ############################################################################ #
rm(list = ls())
gc()

##### Setup #####
library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 24)
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"
dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
dir_pm25_comp <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25_components/"
dir_no2 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/"
dir_ozone <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/seasonal/zipcode/"
dir_temp <- "/nfs/home/D/dam9096/shared_space/ci3_confounders/data_for_analysis/prepped_temperature/annual/"
dir_denominator <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/"

##### Read denom files ##### 
cat("Reading denominator files...")
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip", "hmo_mo", "age", "dead")
dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
setkey(dt, year, zip, qid) # 651,691,916 total person years of data
# ! Check how many people removed
dt[year < 2000 | zip < 501 | hmo_mo != 0, .N] # 162,177,208 person-years

dt <- dt[year >= 2000 &
           zip >= 501 & 
           hmo_mo == 0] # valid zips are >= 00501, retain FFS only (no HMO)
dt <- unique(dt, by = c("year", "zip", "qid"))
cat(dt[,.N], "records\n")
# 489,514,708


##### Enrollment period and continuous enrollment #####
cat("Creating entry/exit data...\n")
setkey(dt, qid, year)
# function to get last year of continuous enrollment in Medicare FFS
last_yr_cont <- function(year) { first((first(year):2017)[!(first(year):2017 %in% year)]) - 1 }
qid_entry_exit <- dt[order(qid, year), # ensure first/last relates to entry/exit year
                     .(entry = first(year), 
                       exit = last_yr_cont(year), # left FFS or died or other censored
                       entry_age = first(age), # get age at entry
                       exit_age = last(age), # get age at exit
                       n_years = .N, 
                       n_unique_zips = uniqueN(zip)), 
                     by = qid]
qid_entry_exit <- merge(qid_entry_exit, dt[, .(qid, year, dead)], 
                        by.x = c("qid", "exit"), by.y = c("qid", "year"),
                        all.x = TRUE)
write_fst(qid_entry_exit, paste0(dir_data, "denom/qid_entry_exit.fst"))

##### Potential sample #####
qid_entry_exit[entry == 2000, .N] # 27,133,737
qid_entry_exit[entry == 2000 & exit < 2010 & !dead, .N] # 3,666,985
qid_entry_exit[entry == 2000 & exit < 2010 & !dead, 
               mean(exit - entry + entry_age)] # 77.7 avg age censor
qid_entry_exit[entry == 2000 & exit < 2010 & dead, .N] # 12,325,236
qid_entry_exit[entry == 2000 & exit < 2010 & dead, 
               mean(exit - entry + entry_age)] # 82.8 avg age death

##### Read exposure data #####
cat("Reading exposure data...\n")
pm25_data <- fread(paste0(dir_pm25, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
no2_data <- fread(paste0(dir_no2, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
ozone_data <- fread(paste0(dir_ozone, "ozone_seasonalavg_zipcode.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
temp_rh_data <- fread(paste0(dir_temp, "temperature_annual_zipcode.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
all_exp <- merge(merge(merge(pm25_data, 
                             no2_data, by = c("zip", "year")), 
                       ozone_data, by = c("zip", "year")),
                 temp_rh_data, by = c("zip", "year"))
rm(pm25_data, no2_data, ozone_data, temp_rh_data)
setkey(all_exp, year, zip)

##### Merge denominator and exposure data based on zipcode and year #####
cat("Merging exposures with denominator...\n")
dt <- merge(dt, all_exp, by = c("year", "zip"))
rm(all_exp)
gc()

##### Exposures for each individual (rows) organized by year (cols) #####
# zipcode
cat("Creating yearly zip of residence file...\n")
qid_yr_zip <- dcast(dt, qid ~ year, value.var = "zip", fill = NA)
write_fst(qid_yr_zip, paste0(dir_data, "qid_yr_exposures/qid_yr_zip.fst"))
rm(qid_yr_zip); gc()
# PM2.5
cat("Creating yearly PM2.5 exposure file...\n")
qid_yr_pm25 <- dcast(dt, qid ~ year, value.var = "pm25", fill = NA)
write_fst(qid_yr_pm25, paste0(dir_data, "qid_yr_exposures/qid_yr_pm25.fst"))
rm(qid_yr_pm25); gc()
# NO2
cat("Creating yearly NO2 exposure file...\n")
qid_yr_no2 <- dcast(dt, qid ~ year, value.var = "no2", fill = NA)
write_fst(qid_yr_no2, paste0(dir_data, "qid_yr_exposures/qid_yr_no2.fst"))
rm(qid_yr_no2); gc()
# Ozone
cat("Creating yearly Ozone exposure file...\n")
qid_yr_ozone <- dcast(dt, qid ~ year, value.var = "ozone_summer", fill = NA)
write_fst(qid_yr_ozone, paste0(dir_data, "qid_yr_exposures/qid_yr_ozone.fst"))
rm(qid_yr_ozone); gc()
# Max temp
cat("Creating yearly max temp exposure file...\n")
qid_yr_tmmx <- dcast(dt, qid ~ year, value.var = "tmmx", fill = NA)
write_fst(qid_yr_tmmx, paste0(dir_data, "qid_yr_exposures/qid_yr_tmmx.fst"))
rm(qid_yr_tmmx); gc()
# Max humiditiy
cat("Creating yearly max humidity exposure file...\n")
qid_yr_rmax <- dcast(dt, qid ~ year, value.var = "rmax", fill = NA)
write_fst(qid_yr_rmax, paste0(dir_data, "qid_yr_exposures/qid_yr_rmax.fst"))
rm(qid_yr_rmax); gc()
# Precip amount
cat("Creating yearly precip exposure file...\n")
qid_yr_pr <- dcast(dt, qid ~ year, value.var = "pr", fill = NA)
write_fst(qid_yr_pr, paste0(dir_data, "qid_yr_exposures/qid_yr_pr.fst"))
rm(qid_yr_pr); gc()



##### PM Components #####
pm_comp <- rbindlist(lapply(2000:2016, function(yr) {
  dt <- fread(paste0(dir_pm25_comp, yr, ".csv"))[, year := yr][]

  dt
}))
setkey(pm_comp, year, ZIP)
comp_names <- names(pm_comp)[3:17]
dt <- merge(dt, pm_comp[, -1], 
            by.x = c("year", "zip"), by.y = c("year", "ZIP"),
            all.x = TRUE)

for (n in comp_names) {
  cat("Creating yearly", n, "exposure file...\n")
  qid_yr_comp <- dcast(dt, qid ~ year, value.var = n, fill = NA)
  write_fst(qid_yr_comp, 
            paste0(dir_data, "qid_yr_exposures/qid_yr_pm25comp_", n, ".fst"))
  rm(qid_yr_comp); gc()
}
