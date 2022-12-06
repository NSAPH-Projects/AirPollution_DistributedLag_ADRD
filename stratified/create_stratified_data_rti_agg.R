# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis - Stratified
#' Code: create data
#' Inputs: merged denominator and exposure data                                
#' Outputs: subset data
#' Author: Daniel Mork            
# ############################################################################ #
rm(list = ls())
gc()

#
# ----- Setup -----
#
n_threads = 48
Sys.setenv(OMP_NUM_THREADS = n_threads)

library(data.table)
library(fst)
options(stringsAsFactors = FALSE)
setDTthreads(threads = n_threads)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
dir_strat <- paste0(dir_data, "analysis/DLM_Strat/stratified_data/")

# RTI race
rti <- fread(paste0("~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/",
                    "data/auxiliary_medicare_cols/rti_race_2010.csv"))
rti[, rti_race_cd := as.integer(rti_race_cd)]
rti[rti_race_cd != 0, 
    race_rti := c("wht", "blk", "oth", "api", "his", "nat")[rti_race_cd]]
setkey(rti, qid)


# Load and merge data
qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_qid_clean.fst"),
                    as.data.table = TRUE)
pm25 <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_pm25_clean.fst"),
                 as.data.table = TRUE)
setnames(pm25, paste0("pm25_lag", 1:10))
no2 <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_no2_clean.fst"),
                 as.data.table = TRUE)
setnames(no2, paste0("no2_lag", 1:10))
ozone <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_ozone_clean.fst"),
                 as.data.table = TRUE)
setnames(ozone, paste0("ozone_lag", 1:10))
qid_dat <- cbind(qid_dat, pm25, no2, ozone)
setkey(qid_dat, qid, year)
qid_dat <- merge(qid_dat, rti[, .(qid, race_rti)], by = c("qid"), all.x = TRUE)
qid_dat[!is.na(race_rti), race := race_rti]
qid_dat[, n_hosp := first_hosp + 0L][, n_pers := 1][, age := age_corrected]



# Aggregate data for people who haven't moved
qid_dat[, age_grp := cut(age_corrected, c(74, 79, 84, 89, 94, 99))]
agg_dat <- qid_dat[n_unique_zips == 1, 
                   .(n_hosp = sum(first_hosp), 
                     n_pers = .N, 
                     age = mean(age_corrected)),
                   by = .(year, zip, age_grp, sexM, race, dual)]
cols <- c("statecode", "education", "poverty", "PIR", "pct_blk", "hispanic", 
          "popdensity", "region", "pct_owner_occ", "smoke_rate", "mean_bmi",
          "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax",
          paste0("pm25_lag", 1:10), paste0("no2_lag", 1:10), paste0("ozone_lag", 1:10))
agg_zip <- qid_dat[n_unique_zips == 1, 
                   lapply(.SD, first), .SDcols = cols,
                   by = .(year, zip, age_grp, sexM, race, dual)]
agg_dat <- agg_dat[agg_zip, on = .(year, zip, age_grp, sexM, race, dual)]
dat <- rbind(agg_dat, qid_dat[n_unique_zips > 1, .SD, .SDcols = names(agg_dat)])
rm(agg_dat, agg_zip, qid_dat, pm25, no2, ozone, rti); gc()


#
# ----- Create stratified data -----
#
for (r in c("nat", "api", "oth", "his", "blk", "wht")) {
for (d in 0:1) {
for (s in 0:1) {
  cat("\n", "race", r, "dual", d, "sexM", s, "young", "n =", 
      dat[race == r & dual == d & sexM == s & age < 81, .N])
  write_fst(dat[race == r & dual == d & sexM == s & age < 81], 
            paste0(dir_strat, "racerti-", r, "_dual-", d, 
                   "_sexM-", s, "_age-young.fst"))
  cat("\n", "race", r, "dual", d, "sexM", s, "old", "n =", 
      dat[race == r & dual == d & sexM == s & age >= 81, .N])
  write_fst(dat[race == r & dual == d & sexM == s & age >= 81], 
            paste0(dir_strat, "racerti-", r, "_dual-", d, 
                   "_sexM-", s, "_age-old.fst"))
} # end for s
} # end for d
} # end for r

rm(list=ls())
gc()
