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
file_list <- file_list[grepl("ozone", file_list)]
#
# ----- Send jobs to workers -----
#
if (FALSE) {
library(parallel)
library(doParallel)
library(foreach)

cl <- makeCluster(min(length(file_list), 6), "FORK")
registerDoParallel(cl)
getDoParWorkers()
ret <- foreach(file = file_list,
               .packages = c("data.table", "fst", "mgcv", "dlmtree"),
               .export = c("dir_strat", "dir_results"),
               .errorhandling = "remove",
               .verbose = TRUE) %dopar% 
  # Begin worker job -----
  { 
    n_threads = 2
    Sys.setenv(OMP_NUM_THREADS = n_threads)
    options(stringsAsFactors = FALSE)
    
    file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
    # skip file if results complete
    if (!(paste0(file_start, ".rda") %in% list.files(dir_results))) {
      desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                     function(i) strsplit(i, "-", TRUE)[[1]])
      names(desc) <- sapply(desc, function(i) i[1])
      desc <- lapply(desc, function(i) i[2])

      qid_dat <- read_fst(paste0(dir_strat, file), as.data.table = TRUE)
      desc$n <- qid_dat[, .N]
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
                 family = "logit", initial.params = coef)
      s <- summary(m)

      save(desc, s, file = paste0(dir_results, file_start, ".rda"))
      TRUE
    } else {
      FALSE
    }
  } # End worker job -----
stopCluster(cl)
}


#
# ----- Run remaining results (if parallel failed) -----
#
remain <- c()
for (file in file_list) {
  file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
  if (!(paste0(file_start, ".rda") %in% list.files(dir_results))) {
    remain <- c(remain, file)
  }
}

for (file in (remain)) {
  file_start <- strsplit(file, ".", fixed = TRUE)[[1]][1]
  cat("\n", file_start)
  if (!(paste0(file_start, ".rda") %in% list.files(dir_results))) {
    desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                   function(i) strsplit(i, "-", TRUE)[[1]])
    names(desc) <- sapply(desc, function(i) i[1])
    desc <- lapply(desc, function(i) i[2])
    
    qid_dat <- read_fst(paste0(dir_strat, file), as.data.table = TRUE)
    desc$n <- qid_dat[, .N]
    cat(" (n = ", qid_dat[,.N], ")")
    form <- as.formula(first_hosp ~
                         factor(year) + factor(statecode) + 
                         age_corrected + education + poverty + PIR +
                         pct_blk + hispanic + popdensity + pct_owner_occ)
    coef <- coef(bam(form, data = qid_dat, family = binomial))
    m <- tdlnm(form, data = qid_dat,
               exposure.data = as.matrix(qid_dat[,paste0("lag", 1:10)]),
               exposure.splits = 0,
               n.trees = 10, n.burn = 1000, n.iter = 2000, n.thin = 5,
               tree.params = c(.5, 2),
               family = "logit", initial.params = coef)
    s <- summary(m)
    
    save(desc, s, file = paste0(dir_results, file_start, ".rda"))
    cat(" complete")
  } else {
    cat(" complete")
  }
}
