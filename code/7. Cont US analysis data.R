# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: continental US analysis data creation  
#' Inputs: denominator data (2000-2016)                                 
#' Outputs: cleaned yearly qid enrollment and individual data         
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 64 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
# 4 regions based on census?
NE <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
S <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", 
       "AR", "MD", "DE", "OK", "TX", "LA")
MW <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
W <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary, secondary, or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load all data #####

# year/zip confounders
yr_zip_confounders <- read_fst(paste0(dir_data, "denom/year_zip_confounders.fst"), 
                               as.data.table = TRUE)
yr_zip_confounders[, region := ifelse(statecode %in% NE, "NE",
                                      ifelse(statecode %in% S, "S", 
                                             ifelse(statecode %in% MW, "MW",
                                                    ifelse(statecode %in% W, "W", "Other"))))]
setkey(yr_zip_confounders, year, zip)
valid_zips <- sort(unique(yr_zip_confounders[statecode != "Other", zip]))

# hospitalization data
hosp_dat <- read_fst(paste0(dir_data, "hospitalization/First_hosp_", 
                            AD_ADRD, "_", code_type, ".fst"), 
                     as.data.table = TRUE)
hosp_dat[, first_hosp := TRUE]
setkey(hosp_dat, QID)

# FFS beneficiary enrollment history
qid_entry_exit <- read_fst(paste0(dir_data, "denom/qid_entry_exit.fst"), 
                           as.data.table = TRUE)
qid_entry_exit[entry == 2000, .N] # 27,133,737
qid_entry_exit[entry == 2000, mean(entry_age)] # 27,133,737
# Merge with hospitalization records
setkey(qid_entry_exit, qid)
temp <- merge(hosp_dat[year < 2010, .(QID, year, first_hosp)], 
              qid_entry_exit[entry == 2000, .(qid, entry, exit, entry_age, dead)], 
              by.x = 'QID', by.y = "qid", all.y = TRUE)
temp[is.na(first_hosp), first_hosp := FALSE]
temp[is.na(year), year := exit]
temp[entry == 2000 & year < 2010 & first_hosp, .N] # 4,369,407
temp[entry == 2000 & year < 2010 & first_hosp, 
     mean(exit - entry + entry_age)] # 85.8
temp[entry == 2000 & exit < 2010 & dead & !first_hosp, .N] # 9,170,192
temp[entry == 2000 & exit < 2010 & dead & !first_hosp, 
     mean(exit - entry + entry_age)] # 81.8 avg age death
temp[entry == 2000 & exit < 2010 & !dead & !first_hosp, .N] # 3,385,436
temp[entry == 2000 & exit < 2010 & !dead & !first_hosp, 
               mean(exit - entry + entry_age)] # 77.4 avg age censor
rm(temp)


# Restrict to 2010 and later
pre2009_qid <- sort(unique(hosp_dat[year <= 2009, QID]))
hosp_dat <- hosp_dat[year > 2009]
qid_entry_exit <- qid_entry_exit[exit > 2009 & # find all QIDs with exit after 2009
                                   entry == 2000 & # entry 2000
                                   !(qid %in% pre2009_qid)] # eliminate prior ADRD hosp
setkey(qid_entry_exit, qid)



# FFS beneficiary yearly data
qid_dat <- data.table()
for (yr in 2010:2016) {
  cat("\nReading year", yr)
  dt <- read_fst(paste0(dir_data, "denom/qid_data_", yr, ".fst"), 
                 as.data.table = TRUE)
  dt <- dt[zip %in% valid_zips & # restrict to continental US zipcodes
             qid %in% qid_entry_exit$qid & # restrict to exit after 2009, continuous enrollment, no prior ADRD hosp
             !(qid %in% hosp_dat[year < yr, QID])] # eliminate ADRD hosp prior to current year
  cat(" -", dt[,.N], "records")
  qid_dat <- rbindlist(list(qid_dat, dt))
  rm(dt)
}
setkey(qid_dat, year, zip, qid)


##### 2. Merge data sources #####
qid_dat <- merge(qid_dat, hosp_dat[, .(QID, year, first_hosp)],
                 by.x = c("year", "qid"), by.y = c("year", "QID"), 
                 all.x = TRUE)
qid_dat[is.na(first_hosp), first_hosp := FALSE]
qid_dat <- merge(qid_dat, yr_zip_confounders, by = c("year", "zip"), all.x = TRUE)
qid_dat <- merge(qid_dat, qid_entry_exit, by = "qid", all.x = TRUE)

setkey(qid_dat, year, zip, qid)
write_fst(qid_dat, paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                          code_type, "_qid.fst"))
rm(hosp_dat, yr_zip_confounders, qid_entry_exit, pre2009_qid); gc()


##### 3. Create corresponding exposure data #####
exposures <- c("pm25", "no2", "ozone") # from pm25, no2, ozone, tmmx, rmax, pr
for (e in exposures) {
  cat("\nCreating exposure data:", e, "...")
  qid_yr_exp <- read_fst(paste0(dir_data, "qid_yr_exposures/qid_yr_", e, ".fst"), 
                         as.data.table = TRUE)
  qid_yr_exp <- qid_yr_exp[qid %in% qid_dat$qid]
  setkey(qid_yr_exp, qid)
  exp_dat <- rbindlist(lapply(2016:2009, function(y) {
    dt <- merge(qid_dat[year == y, .(qid, year, zip)], 
                qid_yr_exp[, c("qid", as.character(2000:y)), with = FALSE], 
                by = "qid", all.x = TRUE)
    setnames(dt, c("qid", "year", "zip", paste0("lag", (y - 2000):0)))
    dt
  }), use.names = TRUE, fill = TRUE)
  setkey(exp_dat, year, zip, qid)
  write_fst(exp_dat, paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                            code_type, "_", e, ".fst"))
  rm(exp_dat)
  cat(" complete.")
}
