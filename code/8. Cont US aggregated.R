# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: continental US analysis aggregated by entry/exit, sex, age, dual, etc. 
#' Inputs: denominator data (2000-2016)                                 
#' Outputs: cleaned yearly qid enrollment and individual data         
#' Author: Daniel Mork                                                         
#' Last updated: Jul. 21, 2022                                               
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
setDTthreads(threads = 4)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

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
qid_entry_exit[exit - entry >= 10, .N] # 17841844
qid_entry_exit[exit - entry >= 10, .N, by = entry][order(entry)]
# entry        N
# 1:  2000 11141516
# 2:  2001  1241729
# 3:  2002  1264695
# 4:  2003  1122231
# 5:  2004  1051211
# 6:  2005  1010477
# 7:  2006  1009985
qid_entry_exit[exit - entry >= 10, mean(entry_age)] # 70.0974
# Merge with hospitalization records
# setkey(qid_entry_exit, qid)
# temp <- merge(hosp_dat[, .(QID, year, first_hosp)], 
#               qid_entry_exit[, .(qid, entry, exit, entry_age, dead)], 
#               by.x = 'QID', by.y = "qid", all.y = TRUE)
# temp[is.na(first_hosp), first_hosp := FALSE]
# temp[is.na(year), year := exit]
# temp[year - entry < 10 & first_hosp, .N] # 1159263
# temp[year - entry < 10 & first_hosp, 
#      mean(exit - entry + entry_age)] # 85.8
# temp[exit - entry < 10 & dead & !first_hosp, .N] # 9,170,192
# temp[exit - entry < 10 & dead & !first_hosp, 
#      mean(exit - entry + entry_age)] # 81.8 avg age death
# temp[exit - entry < 10 & !dead & !first_hosp, .N] # 3,385,436
# temp[exit - entry < 10 & !dead & !first_hosp, 
#                mean(exit - entry + entry_age)] # 77.4 avg age censor
# rm(temp)


# Restrict to 2010 and later
qid_entry_exit <- merge(qid_entry_exit, hosp_dat,
                        by.x = "qid", by.y = "QID", all.x = TRUE)
qid_entry_exit[is.na(first_hosp), c("first_hosp", "year") := list(FALSE, exit)]
qid_entry_exit <- qid_entry_exit[year - entry >= 10] # eliminate prior ADRD hosp
eligible_qid <- sort(unique(qid_entry_exit[, qid]))
setkey(qid_entry_exit, qid)



# FFS beneficiary yearly data
qid_dat <- data.table()
for (yr in 2010:2016) {
  cat("\nReading year", yr)
  dt <- read_fst(paste0(dir_data, "denom/qid_data_", yr, ".fst"), 
                 as.data.table = TRUE)
  dt <- dt[zip %in% valid_zips & # restrict to continental US zipcodes
             qid %in% eligible_qid & # restrict to exit after 2009, continuous enrollment, no prior ADRD hosp
             !(qid %in% hosp_dat[year < yr, QID])] # eliminate ADRD hosp prior to current year
  cat(" -", dt[,.N], "records")
  qid_dat <- rbindlist(list(qid_dat, dt))
  rm(dt)
}
setkey(qid_dat, year, zip, qid)


##### 2. Merge data sources #####
qid_dat <- merge(qid_dat, yr_zip_confounders, by = c("year", "zip"), all.x = TRUE)
qid_dat <- merge(qid_dat, qid_entry_exit, by = "qid", all.x = TRUE)
qid_dat[, year.y := NULL]
qid_dat[, dead.y := NULL]
setnames(qid_dat, c("year.x", "dead.x", "zip.x"), c("year", "dead", "zip"))
qid_dat <- qid_dat[n_unique_zips == 1 & (year - entry >= 10)]

setkey(qid_dat, year, zip, qid)
write_fst(qid_dat, paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                          code_type, "_qid.fst"))
rm(hosp_dat, yr_zip_confounders, qid_entry_exit, eligible_qid); gc()


##### 3. Create corresponding exposure data #####
exposures <- c("pm25") # from pm25, no2, ozone, tmmx, rmax, pr
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
  write_fst(exp_dat, paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                            code_type, "_", e, ".fst"))
  rm(exp_dat)
  cat(" complete.")
}





##### 4. Aggregate for mixtures data #####
qid_dat[, PIR := medianhousevalue / (medhouseholdincome + 0.001)]
qid_dat[, age_corrected := entry_age + year - entry]
qid_dat <- qid_dat[!is.infinite(qid_dat$PIR) &
                     complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                               region, statecode,
                                               education, poverty, pct_blk, hispanic,
                                               popdensity, pct_owner_occ, PIR,
                                               winter_tmmx, winter_rmax,
                                               summer_tmmx, summer_rmax)])]
qid_dat[, .N]
qid_dat[, dual := dual == "1"]
qid_dat[, dual_any := ifelse(any(dual), "Yes", "No"), by = qid]
qid_dat[, sex := ifelse(sexM, "M", "F"), by = qid]
qid_dat[, age_cat := cut(age_corrected, breaks = c(74, 84, 114), labels = c("75-84", "85+"))]
qid_dat_agg_zip <- qid_dat[,.(at_risk = .N, n_hosp = sum(first_hosp), 
                              rate = sum(first_hosp / (year - entry)) / .N,
                              avg_age = mean(age_corrected),
                              his = mean(hispanic), blk = mean(pct_blk),
                              PIR = mean(PIR), pov = mean(poverty), ed = mean(education),
                              oo = mean(pct_owner_occ), dens = mean(popdensity),
                              s_tmmx = mean(summer_tmmx), w_tmmx = mean(winter_tmmx),
                              s_rmax = mean(summer_rmax), w_rmax = mean(winter_rmax),
                              lat = mean(latitude), long = mean(longitude), 
                              state = first(statecode), region = first(region)),
                           by = .(year, entry, zip, sex, age_cat, dual_any, race)]
qid_dat_agg_zip[rate == 0, rate := min(qid_dat_agg_zip$rate[qid_dat_agg_zip$rate > 0])]
qid_dat_agg_zip[, lograte := log(rate)]
write_fst(qid_dat_agg_zip, paste0(dir_data, "analysis/DLMM_agg/contUS_agg_", 
                                  AD_ADRD, "_", code_type, "_agg_zip.fst"))

##### 5. Create corresponding exposure data #####
exposures <- paste0("pm25comp_", c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
                                   "no3", "oc", "pb", "si", "so4", "v", "z"))
for (e in exposures) {
  cat("\nCreating exposure data:", e, "...")
  qid_yr_exp <- read_fst(paste0(dir_data, "qid_yr_exposures/qid_yr_", e, ".fst"), 
                         as.data.table = TRUE)
  qid_yr_exp <- qid_yr_exp[qid %in% qid_dat_agg_zip$qid]
  setkey(qid_yr_exp, qid)
  exp_dat <- rbindlist(lapply(2016:2009, function(y) {
    dt <- merge(qid_dat_agg_zip[year == y, .(qid, year, zip)], 
                qid_yr_exp[, c("qid", as.character(2000:y)), with = FALSE], 
                by = "qid", all.x = TRUE)
    setnames(dt, c("qid", "year", "zip", paste0("lag", (y - 2000):0)))
    dt
  }), use.names = TRUE, fill = TRUE)
  setkey(exp_dat, year, zip, qid)
  write_fst(exp_dat, paste0(dir_data, "analysis/DLMM_agg/contUS_agg_", AD_ADRD, "_", 
                            code_type, "_", e, ".fst"))
  rm(exp_dat)
  cat(" complete.")
}
