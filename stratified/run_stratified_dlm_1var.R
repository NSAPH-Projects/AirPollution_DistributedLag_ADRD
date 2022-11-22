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
dir_results <- paste0(dir_data, "analysis/DLM_Strat/stratified_analysis_1var/")
file_list <- list.files(dir_strat)

#
# ----- Create list of stratifications -----
#
desc <- list()
for (file in (file_list)) {
  file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
  desc[[file]] <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                         function(i) strsplit(i, "-", TRUE)[[1]])
  names(desc[[file]]) <- sapply(desc[[file]], function(i) i[1])
  desc[[file]] <- lapply(desc[[file]], function(i) i[2])
  desc[[file]]$file <- file
}
desc <- rbindlist(desc)


#
# ----- Send jobs to workers -----
#
if (FALSE) {
library(parallel)
library(doParallel)
library(foreach)

cl <- makeCluster(2, "FORK")
registerDoParallel(cl)
getDoParWorkers()
ret <- foreach(exp = c("no2", "ozone", "pm25"),
               .packages = c("data.table", "fst", "mgcv", "dlmtree"),
               .export = c("dir_strat", "dir_results", "desc"),
               .errorhandling = "remove",
               .verbose = TRUE) %:%
  foreach(var = c("racerti", "dual", "sexM", "age")) %dopar% 
  # Begin worker job -----
  { 
    sink(paste0(dir_data, "analysis/DLM_Strat/log_1var.txt"), append = TRUE)
    n_threads = 2
    Sys.setenv(OMP_NUM_THREADS = n_threads)
    options(stringsAsFactors = FALSE)
    
    for (level in unique(desc[[var]])) {
      if (!(paste0("exp-", exp, "_var-", var, "_level-", level, ".rda") %in% 
            list.files(dir_results))) {
        qid_dat <- 
          rbindlist(lapply(which(desc[[var]] == level & desc$exp == exp), function(x) {
          read_fst(paste0(dir_strat, desc$file[x]), as.data.table = TRUE)
        }))
        n <- qid_dat[, .N]
        cat("Running model: exp =", exp, "var = ", var, "level =", level,
            "n =", n, "\n")
        
        form <- as.formula(first_hosp ~
                             factor(year) + factor(statecode) + 
                             age_corrected + education + poverty + PIR +
                             pct_blk + hispanic + popdensity + pct_owner_occ)
        coef <- coef(bam(form, data = qid_dat, family = binomial))
        m <- tdlnm(form, data = qid_dat,
                   exposure.data = as.matrix(qid_dat[,paste0("lag", 1:10)]),
                   exposure.splits = 0,
                   n.trees = 5, n.burn = 1000, n.iter = 2000, n.thin = 5,
                   tree.params = c(.5, 2),
                   family = "logit", initial.params = coef, verbose = FALSE)
        s <- summary(m)
        
        save(var, level, exp, n, s,
             file = paste0(dir_results, 
                           "exp-", exp, "_var-", var, "_level-", level, ".rda"))
        TRUE
      } else {
        FALSE
      }
    } # end for over levels of var
    
  } # End worker job -----
stopCluster(cl)
}


#
# ----- Run remaining results (if parallel failed) -----
#
# remain <- c()
# for (file in file_list) {
#   file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
#   if (!(paste0(file_start, ".rda") %in% list.files(dir_results))) {
#     remain <- c(remain, file)
#   }
# }
# 
# for (file in (remain)) {
#   file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
#   cat("\n", file_start)
#   if (!(paste0(file_start, ".rda") %in% list.files(dir_results))) {
#     desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
#                    function(i) strsplit(i, "-", TRUE)[[1]])
#     names(desc) <- sapply(desc, function(i) i[1])
#     desc <- lapply(desc, function(i) i[2])
#     
#     qid_dat <- read_fst(paste0(dir_strat, file), as.data.table = TRUE)
#     desc$n <- qid_dat[, .N]
#     cat(" (n = ", qid_dat[,.N], ")")
#     form <- as.formula(first_hosp ~
#                          factor(year) + factor(statecode) + 
#                          age_corrected + education + poverty + PIR +
#                          pct_blk + hispanic + popdensity + pct_owner_occ)
#     coef <- coef(bam(form, data = qid_dat, family = binomial))
#     m <- tdlnm(form, data = qid_dat,
#                exposure.data = as.matrix(qid_dat[,paste0("lag", 1:10)]),
#                exposure.splits = 0,
#                n.trees = 10, n.burn = 1000, n.iter = 2000, n.thin = 5,
#                tree.params = c(.5, 2),
#                family = "logit", initial.params = coef)
#     s <- summary(m)
#     
#     save(desc, s, file = paste0(dir_results, file_start, ".rda"))
#     cat(" complete")
#   } else {
#     cat(" complete")
#   }
# }
