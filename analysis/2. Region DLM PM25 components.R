# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: region distributed lag analysis
#' Inputs: merged denominator and exposure data                                
#' Outputs:       
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 64 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
region <- "NE"
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
exp_names <- c(paste0("pm25comp_", c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
                                   "no3", "oc", "pb", "si", "so4", "v", "z")),
                    "no2", "ozone")
lag_cols <- paste0("lag", 1:10)
exposures <- list()
for (e in exp_names) {
  exposures[[e]] <- as.matrix(
    read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, 
                    "_", code_type, "_", e, ".fst"), columns = lag_cols,
             as.data.table = TRUE))
}


##### 2. Restrict to complete follow-up 2010 through 2016, no deaths #####
idx <- which(qid_dat$entry == 2000 & # follow up from 2000
               qid_dat$year > 2009) # year 2010 and later for 10 lag years exposures
qid_dat <- qid_dat[idx]
exposures <- lapply(names(exposures), function(e) exposures[[e]][idx,])
names(exposures) <- exp_names

##### Remove rows with missing exposure data #####
idx <- 1:nrow(qid_dat)
for (e in exp_names) {
  idx <- intersect(idx, which(rowSums(is.na(exposures[[e]])) == 0))
}
qid_dat <- qid_dat[idx]
exposures <- lapply(names(exposures), function(e) exposures[[e]][idx,])
names(exposures) <- exp_names

##### Data fixes #####
qid_dat[, age_corrected := entry_age + year - entry]
qid_dat[, PIR := medianhousevalue/medhouseholdincome]
# clear missing data, infinite PIR
idx2 <- which(!is.infinite(qid_dat$PIR) &
                complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                        education, poverty, pct_blk, hispanic,
                                        popdensity, pct_owner_occ, PIR)]))
qid_dat <- qid_dat[idx2]
exposures <- lapply(names(exposures), function(e) exposures[[e]][idx2,])
names(exposures) <- exp_names

##### Remove QID after missing year of data #####
idx3 <- which(qid_dat$year == 2010)
qid_rm <- qid_dat[idx3, qid]
for (yr in 2011:2016) {
  idx3 <- c(idx3, which(qid_dat$year == yr & qid_dat$qid %in% qid_rm))
  qid_rm <- qid_dat[idx3][year == yr, qid]
}
qid_dat <- qid_dat[idx3]
exposures <- lapply(names(exposures), function(e) exposures[[e]][idx3,])
names(exposures) <- exp_names
rm(idx, idx2, idx3, qid_rm)

##### Summary stats #####
qid_dat[, .N]
qid_dat[, uniqueN(qid)]
qid_dat[, uniqueN(qid), by = year]
qid_dat[, sum(first_hosp)]
qid_dat[year == 2000, mean(age_corrected)]
qid_dat[, .(sum(first_hosp), sum(dead)), by = year]
qid_dat[, sum(first_hosp)/uniqueN(qid)]
qid_dat[, sum(dead)/uniqueN(qid)]
qid_dat[, .(sum(first_hosp)/uniqueN(qid), sum(dead)/uniqueN(qid)), by = year]


#### Lag 1 correlations ####
l1exp <- sapply(exposures, function(e) e[,1])
l1cor <- cor(l1exp)
corrplot::corrplot(l1cor, type = "lower", method = "square")


##### 4. DLMM analysis - TDLMM #####
library(dlmtree)
library(mgcv)

set.seed(8372)
qid_samp <- qid_dat[year == 2010, qid][sample(qid_dat[year == 2010, .N], 500000)]
qid_dat[qid %in% qid_samp, .N] #2484214

form <- as.formula(first_hosp ~ factor(year) +
                     age_corrected + I(age_corrected^2) + 
                     factor(dual) + factor(race) + factor(sexM) +
                     education + poverty + PIR +
                     pct_blk + hispanic + popdensity + pct_owner_occ)
init.params <- bam(form,
                   data = qid_dat[qid %in% qid_samp], 
                   family = binomial,
                   nthreads = 1, samfrac = 0.05, chunk.size = 5000,
                   control = gam.control(trace = TRUE))$coefficients

mm <- tdlmm(form,
            data = qid_dat,
            exposure.data = exposures, # 15 pm components, no2, ozone
            mixture.interactions = "noself",
            family = "logit", shrinkage = "exposures",
            subset = which(qid_dat$qid %in% qid_samp),
            mix.prior = 1.404, # prior inc prob = 0.9
            n.trees = 50, 
            n.burn = 1000, n.iter = 2000, n.thin = 4,
            initial.params = init.params)
save(mm, file = paste0(dir_data, "analysis/tdlmm500k_50_0.9_region", 
                       region, "_", AD_ADRD, "_", 
                       code_type, "_pm25comp.rda"))
(ss <- summary(mm))
plot.summary.tdlmm(ss)

