# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: region analysis data creation  
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


region <- NE
region_name <- "NE"
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary, secondary, or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load all data for single region #####

# year/zip confounders
yr_zip_confounders <- read_fst(paste0(dir_data, "denom/year_zip_confounders.fst"), 
                               as.data.table = TRUE)
yr_zip_confounders <- yr_zip_confounders[statecode %in% region]
setkey(yr_zip_confounders, year, zip)
region_zips <- sort(unique(yr_zip_confounders$zip))

# hospitalization data
hosp_dat <- read_fst(paste0(dir_data, 
                            "hospitalization/First_hosp_", AD_ADRD, "_", code_type, ".fst"), 
                     as.data.table = TRUE)
hosp_dat[, first_hosp := TRUE]
pre2009_qid <- hosp_dat[year < 2009, QID]
hosp_dat <- hosp_dat[year >= 2009]

# FFS beneficiary enrollment history
qid_entry_exit <- read_fst(paste0(dir_data, "denom/qid_entry_exit.fst"), 
                           as.data.table = TRUE)
qid_entry_exit <- qid_entry_exit[exit >= 2009 & # find all QIDs with exit after 2009
                                   entry == 2000 &
                                   !(qid %in% pre2009_qid) & # eliminate prior ADRD hosp
                                   cont_enroll == TRUE] # restrict to continuously enrolled
setkey(qid_entry_exit, qid)

# FFS beneficiary yearly data
qid_dat <- data.table()
for (yr in 2009:2016) {
  cat("\nReading year", yr)
  dt <- read_fst(paste0(dir_data, "denom/qid_data_", yr, ".fst"), 
                 as.data.table = TRUE)
  dt <- dt[zip %in% region_zips & # restrict to region zipcodes
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
write_fst(qid_dat, paste0(dir_data, "analysis/region", region_name, "_", AD_ADRD, "_", 
                          code_type, "_qid.fst"))
rm(hosp_dat, yr_zip_confounders, qid_entry_exit, pre2009_qid); gc()


##### 3. Create corresponding exposure data #####
# exposures <- c("pm25", "no2", "ozone", "tmmx", "rmax", "pr") # from pm25, no2, ozone, tmmx, rmax, pr
exposures <- paste0("pm25comp_", c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
                                   "no3", "oc", "pb", "si", "so4", "v", "z"))
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
  write_fst(exp_dat, paste0(dir_data, "analysis/region", region_name, "_", AD_ADRD, "_", 
                            code_type, "_", e, ".fst"))
  rm(exp_dat)
  cat(" complete.")
}
