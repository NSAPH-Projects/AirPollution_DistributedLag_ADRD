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
AD_ADRD <- "ADRD" # AD or ADRD outcome
code_type <- "any" # primary or any outcome
n_threads = 48
Sys.setenv(OMP_NUM_THREADS = n_threads)

library(data.table)
library(fst)
options(stringsAsFactors = FALSE)
setDTthreads(threads = n_threads)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
dir_strat <- paste0(dir_data, "analysis/DLM_Strat/stratified_data/")

# RTI race
rti <- rbindlist(lapply(c(2010), function(y) {
  d <- fread(paste0("~/nsaph_projects/pm_no2_o3-adrd_hosp-medicare-causalgps/",
                    "data/auxiliary_medicare_cols/rti_race_", y, ".csv"))[,year := y][]
}), fill = TRUE)
rti[, rti_race_cd := as.integer(rti_race_cd)]
rti[rti_race_cd != 0, race_rti := c("wht", "blk", "oth", "api", "his", "nat")[rti_race_cd]]

#
# ----- Create stratified data -----
#
for (exp in c("pm25", "no2", "ozone")) {
  qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_qid_clean.fst"),
                      as.data.table = TRUE)
  exp_dat <- read_fst(paste0(dir_data, "analysis/contUS_ADRD_any_", exp, "_clean.fst"),
                       as.data.table = TRUE)
  qid_dat <- cbind(qid_dat, exp_dat)
  rm(exp_dat)
  qid_dat[, dual_any := any(dual == "1"), by = qid]
  qid_dat <- merge(qid_dat, rti[, .(qid, year, race_rti)], 
                   by = c("qid", "year"), all.x = TRUE)
  qid_dat[!is.na(race_rti), race := race_rti]
  
  for (r in c("nat", "api", "oth", "his", "blk", "wht")) {
  for (d in 0:1) {
  for (s in 0:1) {
    qid_y <- qid_dat[race == r & dual_any == d & sexM == s & 
                       age_corrected < 81 & year == 2010, qid]
    qid_o <- qid_dat[race == r & dual_any == d & sexM == s & 
                       age_corrected > 80 & year == 2010, qid]
    idx_y <- which(qid_dat$year == 2010 & qid_dat$qid %in% qid_y)
    idx_o <- which(qid_dat$year == 2010 & qid_dat$qid %in% qid_o)
    for (yr in 2011:2016) {
      idx_y <- c(idx_y, which(qid_dat$year == yr & qid_dat$qid %in% qid_y))
      idx_o <- c(idx_o, which(qid_dat$year == yr & qid_dat$qid %in% qid_o))
      qid_y <- qid_dat[idx_y][year == yr, qid]
      qid_o <- qid_dat[idx_o][year == yr, qid]
    }
    cat("\n", exp, "race", r, "dual", d, "sexM", s, "young", "n =", length(idx_y))
    cat("\n", exp, "race", r, "dual", d, "sexM", s, "old", "n =", length(idx_o))
    write_fst(qid_dat[idx_y], 
              paste0(dir_strat, "exp-", exp, "_racerti-", r, "_dual-", d, 
                     "_sexM-", s, "_age-young.fst"))
    write_fst(qid_dat[idx_o], 
              paste0(dir_strat, "exp-", exp, "_racerti-", r, "_dual-", d, 
                     "_sexM-", s, "_age-old.fst"))
  } # end for s
  } # end for d
  } # end for r
} # end for exp

rm(list=ls())
gc()
