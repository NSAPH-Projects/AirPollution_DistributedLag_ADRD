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
n_threads = 24
Sys.setenv(OMP_NUM_THREADS = n_threads)

library(data.table)
library(fst)
library(mgcv)
library(dlmtree)
options(stringsAsFactors = FALSE)
setDTthreads(threads = n_threads)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
dir_strat <- paste0(dir_data, "analysis/DLM_Strat/stratified_data/")
dir_results <- paste0(dir_data, "analysis/DLM_Strat/stratified_analysis/")
file_list <- list.files(dir_strat)
file_list <- file_list[!grepl("exp-", file_list)]
#
# ----- Send jobs to workers -----
#
if (FALSE) {
library(parallel)
library(doParallel)
library(foreach)

cl <- makeCluster(min(length(file_list), 4), "PSOCK")
registerDoParallel(cl)
getDoParWorkers()
ret <- foreach(file = file_list,
               .packages = c("data.table", "fst", "mgcv", "dlmtree"),
               .export = c("dir_strat", "dir_results", "dir_data"),
               .errorhandling = "remove",
               .verbose = TRUE) %dopar% 
  # Begin worker job -----
  { 
    sink(paste0(dir_data, "analysis/DLM_Strat/log_var.txt"), append = TRUE)
    file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
    desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                   function(i) strsplit(i, "-", TRUE)[[1]])
    names(desc) <- sapply(desc, function(i) i[1])
    desc <- lapply(desc, function(i) i[2])
    
    cat("\nReading file", file)
    dat <- read_fst(paste0(dir_strat, file), as.data.table = TRUE)
    desc$n <- dat[, .N]
    
    # estimate initial parameters 
    form <- as.formula(cbind(n_hosp, n_pers - n_hosp) ~
                         factor(year) + factor(statecode) + 
                         age + education + poverty + PIR +
                         pct_blk + hispanic + popdensity + pct_owner_occ +
                         smoke_rate + mean_bmi +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax)
    coef <- coef(bam(form, data = dat, family = binomial))
    
    form <- as.formula(n_hosp ~
                         factor(year) + factor(statecode) + 
                         age + education + poverty + PIR +
                         pct_blk + hispanic + popdensity + pct_owner_occ +
                         smoke_rate + mean_bmi +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax)
    for (exp in c("pm25", "no2", "ozone")) {
      desc$exp <- exp
      # skip file if results complete
      if (!(paste0("exp-", exp, "_", file_start, ".rda") %in% 
            list.files(dir_results))) {
        cat("\nRunning model: exp =", exp, "file = ", file_start, "n =", desc$n, "\n")
        
        m <- tdlnm(form, data = dat,
                   exposure.data = as.matrix(dat[, .SD, .SDcols = paste0(exp, "_lag", 1:10)]),
                   exposure.splits = 0,
                   n.trees = 10, n.burn = 1000, n.iter = 2000, n.thin = 5,
                   family = "logit", binomial.size = dat[, n_pers],
                   initial.params = coef)
        s <- summary(m)
        save(desc, s, file = paste0(dir_results, "exp-", exp, "_", file_start, ".rda"))
      }
    }
    TRUE
  } # End worker job -----
stopCluster(cl)
}


#
# ----- Run remaining results (if parallel failed) -----
#

for (file in file_list) {
  file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
  cat("\n", file_start)
  desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                 function(i) strsplit(i, "-", TRUE)[[1]])
  names(desc) <- sapply(desc, function(i) i[1])
  desc <- lapply(desc, function(i) i[2])
  
  dat <- read_fst(paste0(dir_strat, file), as.data.table = TRUE)
  desc$n <- dat[, .N]
  cat(" (n = ", dat[,.N], ")")
  form <- as.formula(cbind(n_hosp, n_pers - n_hosp) ~
                       factor(year) + factor(statecode) + 
                       age + education + poverty + PIR +
                       pct_blk + hispanic + popdensity + pct_owner_occ +
                       smoke_rate + mean_bmi +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax)
  coef <- coef(bam(form, data = dat, family = binomial))
  
  form <- as.formula(n_hosp ~
                       factor(year) + factor(statecode) + 
                       age + education + poverty + PIR +
                       pct_blk + hispanic + popdensity + pct_owner_occ +
                       smoke_rate + mean_bmi +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax)
  for (exp in c("pm25", "no2", "ozone")) {
    # skip file if results complete
    if (!(paste0("exp-", exp, "_", file_start, ".rda") %in% 
          list.files(dir_results))) {
      cat("\nRunning model: exp =", exp, "file = ", file_start, "n =", desc$n, "\n")
      desc$exp <- exp
      cat("", exp)
      m <- tdlnm(form, data = dat,
                 exposure.data = as.matrix(dat[, .SD, .SDcols = paste0(exp, "_lag", 1:10)]),
                 exposure.splits = 0,
                 n.trees = 10, n.burn = 1000, n.iter = 2000, n.thin = 5,
                 family = "logit", binomial.size = dat[, n_pers],
                 initial.params = coef)
      s <- summary(m)
      
      save(desc, s, file = paste0(dir_results, "exp-", exp, "_", file_start, ".rda"))
    }
  }
}
