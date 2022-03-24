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

##### 0. Setup #####
library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 24)
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"
dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
dir_no2 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/"
dir_ozone <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/annual/zipcode/requaia_predictions/ywei_aggregation/"
dir_temp <- "/nfs/home/D/dam9096/shared_space/ci3_confounders/data_for_analysis/prepped_temperature/annual/"
dir_denominator <- "/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/"

##### 1. Read denom files ##### 
cat("Reading denominator files...")
f <- list.files(dir_denominator, pattern = "\\.fst", full.names = TRUE)
myvars <- c("qid", "year", "zip", "hmo_mo", "age")
dt <- rbindlist(lapply(f[2:18], read_fst, columns = myvars, as.data.table = TRUE))
setkey(dt, year, zip, qid)
# ! Check how many people removed
dt <- dt[year >= 2000 &
           zip >= 501 & 
           hmo_mo == 0] # valid zips are >= 00501, retain FFS only (no HMO)
dt <- unique(dt, by = c("year", "zip", "qid"))
cat(dt[,.N], "records\n")


##### 2. Read exposure data #####
cat("Reading exposure data...\n")
pm25_data <- fread(paste0(dir_pm25, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
no2_data <- fread(paste0(dir_no2, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
ozone_data <- fread(paste0(dir_ozone, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
temp_rh_data <- fread(paste0(dir_temp, "temperature_annual_zipcode.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
all_exp <- merge(merge(merge(pm25_data, 
                             no2_data, by = c("zip", "year")), 
                       ozone_data, by = c("zip", "year")),
                 temp_rh_data, by = c("zip", "year"))
rm(pm25_data, no2_data, ozone_data, temp_rh_data)
setkey(all_exp, year, zip)

##### 3. Merge denominator and exposure data based on zipcode and year #####
cat("Merging exposures with denominator...\n")
dt <- merge(dt, all_exp, by = c("year", "zip"))
rm(all_exp)
gc()

##### 4. Exposures for each individual (rows) organized by year (cols) #####
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
qid_yr_ozone <- dcast(dt, qid ~ year, value.var = "ozone", fill = NA)
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

##### 5. Enrollment period and continuous enrollment #####
cat("Creating entry/exit data...\n")
setkey(dt, qid, year)
qid_entry_exit <- dt[order(qid, year), # ensure first/last relates to entry/exit year
                     .(entry = first(year), 
                       exit = last(year), 
                       entry_age = first(age), # get age at entry
                       exit_age = last(age), # get age at exit
                       n_years = .N, 
                       cont_enroll = (last(year) - first(year) + 1 == .N),
                       correct_ages = (last(age) - first(age) + 1 == .N),
                       n_unique_zips = uniqueN(zip)), 
                     by = qid]
write_fst(qid_entry_exit, paste0(dir_data, "denom/qid_entry_exit.fst"))

qid_entry_exit[, .N, by = .(cont_enroll, correct_ages)]
# cont_enroll correct_ages        N
# 1:        TRUE         TRUE 55998926
# 2:       FALSE        FALSE  2995822
# 3:        TRUE        FALSE  5271784
# 4:       FALSE         TRUE   134856


qid_entry_exit[cont_enroll & correct_ages, .N, by = n_unique_zips][order(n_unique_zips)]
# n_unique_zips        N
# 1:             1 45533904
# 2:             2  8292390
# 3:             3  1679976
# 4:             4   367085
# 5:             5    90344
# 6:             6    24589
# 7:             7     7261
# 8:             8     2286
# 9:             9      735
# 10:            10+    356

qid_entry_exit[cont_enroll & correct_ages, .N, by = n_years][order(n_years)]
# n_years       N
# 1:       1 6338465
# 2:       2 5812470
# 3:       3 4963635
# 4:       4 4712190
# 5:       5 4271579
# 6:       6 3978715
# 7:       7 3472749
# 8:       8 3184883
# 9:       9 2853617
# 10:      10 2043007
# 11:      11 1802460
# 12:      12 1705681
# 13:      13 1611292
# 14:      14 1534877
# 15:      15 1424418
# 16:      16 1210542
# 17:      17 5078346
