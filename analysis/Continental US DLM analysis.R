# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: continental US distributed lag analysis
#' Inputs: merged denominator and exposure data                                
#' Outputs: distributed lag results
#' Author: Daniel Mork            
# ############################################################################ #
rm(list = ls())
gc()

##### Setup #####
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
n_threads = 10
Sys.setenv(OMP_NUM_THREADS = n_threads)


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = n_threads)
threads_fst(n_threads)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

##### Load and clean relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
pm25_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                            code_type, "_pm25.fst"), as.data.table = TRUE)
no2_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_",
                           code_type, "_no2.fst"), as.data.table = TRUE)
ozone_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_",
                             code_type, "_ozone.fst"), as.data.table = TRUE)
pm25_dat <- as.matrix(pm25_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
no2_dat <- as.matrix(no2_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
ozone_dat <- as.matrix(ozone_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0


##### + Restrict to complete follow-up 2010 through 2016, no deaths #####
idx <- which(rowSums(is.na(pm25_dat[, 1:10])) == 0 & # no missing exposures lags 1-10
               rowSums(is.na(no2_dat[, 1:10])) == 0 &
               rowSums(is.na(ozone_dat[, 1:10])) == 0 &
               qid_dat$entry == 2000 & # follow up from 2000
               qid_dat$year > 2009) # year 2010 and later for 10 lag years exposures
qid_dat <- qid_dat[idx]
pm25_dat <- pm25_dat[idx,]
no2_dat <- no2_dat[idx,]
ozone_dat <- ozone_dat[idx,]

##### + Data fixes #####
qid_dat[, age_corrected := entry_age + year - entry]
qid_dat[, PIR := medianhousevalue/medhouseholdincome]
# clear missing data, zero exposure, infinite PIR
idx2 <- which(rowSums(pm25_dat[, 1:10] == 0) == 0 &
                rowSums(no2_dat[, 1:10] == 0) == 0 &
                rowSums(ozone_dat[, 1:10] == 0) == 0 &
                !is.infinite(qid_dat$PIR) &
                complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                          region, statecode,
                                          education, poverty, pct_blk, hispanic,
                                          popdensity, pct_owner_occ, PIR)]))
qid_dat <- qid_dat[idx2]
pm25_dat <- pm25_dat[idx2,]
no2_dat <- no2_dat[idx2,]
ozone_dat <- ozone_dat[idx2,]

##### + Remove QID if missing earlier year of data #####
idx3 <- which(qid_dat$year == 2010)
qid_rm <- qid_dat[idx3, qid]
for (yr in 2011:2016) {
  idx3 <- c(idx3, which(qid_dat$year == yr & qid_dat$qid %in% qid_rm))
  qid_rm <- qid_dat[idx3][year == yr, qid]
}
qid_dat <- qid_dat[idx3]
pm25_dat <- pm25_dat[idx3,]
no2_dat <- no2_dat[idx3,]
ozone_dat <- ozone_dat[idx3,]
rm(idx, idx2, idx3, qid_rm)
qid_dat[,dead.y:=NULL]
setnames(qid_dat,"dead.x", "dead")

write_fst(qid_dat, paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                          code_type, "_qid_clean.fst"))
write_fst(as.data.frame(pm25_dat), paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                          code_type, "_pm25_clean.fst"))
write_fst(as.data.frame(no2_dat), paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                         code_type, "_no2_clean.fst"))
write_fst(as.data.frame(ozone_dat), paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                           code_type, "_ozone_clean.fst"))


##### Summary stats #####
(n <- qid_dat[year == 2010, .N]) # unique individuals
qid_dat[year == 2010, mean(age_corrected)]
qid_dat[, .N] # person years = 42288065
qid_dat[year == 2010, .(.N, .N/n), by = race] # breakdown by race
qid_dat[year == 2010, .(.N, .N/n), by = sexM] # breakdown by sex
qid_dat[year == 2010, .(.N, .N/n), by = dual] # breakdown by dual eligible
qid_dat[year == 2010, .(age = mean(age_corrected), sd = sd(age_corrected),
                        min = min(age_corrected), max = max(age_corrected))] # age
qid_dat[year == 2010, .(age7584 = sum(age_corrected >=75 & age_corrected < 80),
                        age8594 = sum(age_corrected >=80 & age_corrected < 85),
                        age8594 = sum(age_corrected >=85 & age_corrected < 90),
                        age8594 = sum(age_corrected >=90 & age_corrected < 95),
                        age95 = sum(age_corrected >= 95),
                        age7584p = sum(age_corrected >=75 & age_corrected < 80)/.N,
                        age8594p = sum(age_corrected >=80 & age_corrected < 85)/.N,
                        age8594p = sum(age_corrected >=85 & age_corrected < 90)/.N,
                        age8594p = sum(age_corrected >=90 & age_corrected < 95)/.N,
                        age95p = sum(age_corrected >= 95)/.N)] # age

qid_dat[, .(n_events = sum(first_hosp),
            prop_events = sum(first_hosp)/uniqueN(qid), 
            n_deaths = sum(dead.x),
            prop_deaths = sum(dead.x)/uniqueN(qid))] # prop ADRD events or death
qid_dat[, .(N = uniqueN(qid), # proportions by year
            n_events = sum(first_hosp), prop_events = sum(first_hosp)/uniqueN(qid), 
            n_deaths = sum(dead), prop_deaths = sum(dead)/uniqueN(qid)), by = year]
qid_dat[first_hosp == 1, .(age = mean(age_corrected), sd = sd(age_corrected))]
qid_dat[dead.x == 1, .(age = mean(age_corrected), sd = sd(age_corrected))]

##### Exposure distribution #####
pm25_c <- c(pm25_dat[which(qid_dat$year == 2010), 1:10])
no2_c <- c(no2_dat[which(qid_dat$year == 2010), 1:10])
oz_c <- c(ozone_dat[which(qid_dat$year == 2010), 1:10])
for (yr in 2011:2016) {
  id_yr <- which(qid_dat$year == yr)
  pm25_c <- c(pm25_c, pm25_dat[id_yr, 1])
  no2_c <- c(no2_c, no2_dat[id_yr, 1])
  oz_c <- c(oz_c, ozone_dat[id_yr, 1])
}
IQR(pm25_c)
# [1] 4.103418
IQR(no2_c)
# [1] 14.98296
IQR(oz_c)
# [1] 10.08973
quantile(pm25_c, c(0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995))
# 0.5%       10%       25%       50%       75%       90%     99.5% 
# 2.927854  6.699503  8.634108 10.595629 12.737525 14.571118 20.149561 
quantile(no2_c, c(0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995))
# 0.5%       10%       25%       50%       75%       90%     99.5% 
# 3.928043  8.468661 12.064160 18.257298 27.047125 35.314847 51.711734 
quantile(oz_c, c(0.005, 0.1, 0.25, 0.5, 0.75, 0.9, 0.995))
# 0.5%      10%      25%      50%      75%      90%    99.5% 
# 23.92878 35.54242 41.53470 46.68333 51.62444 55.60097 68.38790 

##### Lagged exposure averages #####
exp_avg <- rbindlist(
  lapply(2010:2016, function(yr) {
    w <- which(qid_dat$year == yr)
    data.table(year = yr, lag = 1:10, 
               t(sapply(1:10, function(l) {
                 c(mean(pm25_dat[w,l]), quantile(pm25_dat[w, l], c(.25, .75)))
               })))
  })
)
setnames(exp_avg, c("year", "lag", "mean", "q1", "q3"))
library(ggplot2)
ggplot(exp_avg, aes(x = year - lag, y = mean, ymin = q1, ymax = q3, fill = factor(year),
                    color = factor(year), group = factor(year))) +
  geom_ribbon(alpha = 0.2) +
  geom_line(size = 2) +
  facet_wrap("year", nrow = 2) +
  theme_bw(base_size = 24) + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(2000, 2004, 2008, 2012)) +
  labs(x = "Year", y = "Exposure Concentration")




##### DLM #####
library(mgcv)
library(dlnm)
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                           code_type, "_qid_clean.fst"),
                    as.data.table = TRUE)
id10 <- which(qid_dat$year == 2010)

for (exp in c("pm25", "no2", "ozone")) {
  cat("\n", exp)
  exp_dat <- as.matrix(read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                       code_type, "_", exp, "_clean.fst")))[,1:10]
  exp_iqr <- IQR(exp_dat[id10, 1:10])
  cb_dlm <- crossbasis(exp_dat[, 1:10], c(1, 10), 
                       arglag = list(fun = "ns", knots = exp(log(10/2)), 
                                     intercept = TRUE), 
                       argvar = list(fun = "lin"))
  
  prD <- bam(dead ~ cb_dlm + factor(year),
             data = qid_dat, family = binomial,
             chunk.size = 5000,
             control = gam.control(trace = TRUE))
  prD_L <- bam(dead ~ cb_dlm +
                 factor(year) + factor(statecode) +
                 age_corrected + I(age_corrected^2) +
                 factor(dual) + factor(race) + factor(sexM) +
                 education + poverty + PIR +
                 pct_blk + hispanic + popdensity + pct_owner_occ,
               data = qid_dat, family = binomial,
               chunk.size = 5000,
               control = gam.control(trace = TRUE))
  qid_dat[, IPWcr := (1 - prD$fitted.values) / (1 - prD_L$fitted.values)]
  qid_dat[IPWcr > 10, IPWcr := 10]
  rm(prD, prD_L)
  
  m3 <- bam(first_hosp ~ cb_dlm +
              factor(year) + factor(statecode) +
              factor(dual) + factor(race) + factor(sexM) +
              age_corrected + I(age_corrected^2) + 
              education + poverty + PIR +
              pct_blk + hispanic + popdensity + pct_owner_occ,
            data = qid_dat, weights = qid_dat$IPWcr, 
            family = binomial,
            chunk.size = 5000,
            control = gam.control(trace = TRUE))
  cp_dlm <- crosspred(cb_dlm, m3, cen = 0, at = 0:1, bylag = 0.2)
  plot(cp_dlm, "slices", var = 1, 
       xlab = "Years prior", ylab = "Hazard Odds")
  res <- list(cp = cp_dlm)
  rm(cb_dlm, m3); gc()
  
  save(res, exp_iqr, knots,
       file = paste0(dir_data, "analysis/gamns3log_dlm_comprisk_contUS_", 
                     AD_ADRD, "_", code_type, "_", exp, ".rda"))
}


##### DLNM #####
library(mgcv)
library(dlnm)
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                           code_type, "_qid_clean.fst"), as.data.table = TRUE)
for (exp in c("pm25", "no2", "ozone")) {
  cat(exp)
  exp_dat <- as.matrix(read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                       code_type, "_", exp, "_clean.fst")))[,1:10]
  # id10 <- which(qid_dat$year == 2010)
  exp_c <- c(exp_dat[which(qid_dat$year == 2010), 1:10])
  for (yr in 2011:2016) {
    id_yr <- which(qid_dat$year == yr)
    exp_c <- c(exp_c, exp_dat[id_yr, 1])
  }
  
  k <- quantile(exp_c, c(.25, .5, .75))
  cb_temp <- crossbasis(exp_dat[1:10000, 1:10], c(1, 10), 
                        arglag = list(fun = "ns", knots = exp(log(10)/2), intercept = TRUE), 
                        argvar = list(fun = "ns", knots = k))
  # create crossbasis in parts, data too large for dlnm package
  cat(", cb")
  s <- floor(seq(0, nrow(exp_dat), length.out = 11))
  cb_dlnm <- do.call(rbind, lapply(1:10, function(i) {
    crossbasis(exp_dat[(s[i] + 1):s[i + 1], 1:10], c(1, 10),
               arglag = list(fun = "ns", knots = exp(log(10)/2), intercept = TRUE),
               argvar = list(fun = "ns", knots = k))
  }))
  class(cb_dlnm) <- c("crossbasis", "matrix")
  attributes(cb_dlnm)$range <- range(exp_dat[, 1:10])
  attributes(cb_dlnm)$df <- attributes(cb_temp)$df
  attributes(cb_dlnm)$lag <- attributes(cb_temp)$lag
  attributes(cb_dlnm)$argvar <- attributes(cb_temp)$argvar
  attributes(cb_dlnm)$arglag <- attributes(cb_temp)$arglag
  attributes(cb_dlnm)$dimnames <- list()
  attributes(cb_dlnm)$dimnames[[1]] <- NULL
  attributes(cb_dlnm)$dimnames[[2]] <- attributes(cb_temp)$dimnames[[2]]
  
  cat(", ipw")
  prD <- bam(dead ~ cb_dlnm + factor(year),
             data = qid_dat, family = binomial,
             nthreads = n_threads, chunk.size = 5000, 
             control = gam.control(trace = TRUE))
  prD_L <- bam(dead ~ cb_dlnm +
                 factor(year) + factor(statecode) +
                 age_corrected + I(age_corrected^2) +
                 factor(dual) + factor(race) + factor(sexM) +
                 education + poverty + PIR +
                 pct_blk + hispanic + popdensity + pct_owner_occ,
               data = qid_dat, family = binomial,
               nthreads = n_threads, chunk.size = 5000, 
               control = gam.control(trace = TRUE))
  qid_dat[, IPWcr := (1 - prD$fitted.values) / (1 - prD_L$fitted.values)]
  qid_dat[IPWcr > 10, IPWcr := 10]
  rm(prD, prD_L)
  
  m3 <- bam(first_hosp ~ cb_dlnm +
              factor(year) + factor(statecode) +
              factor(dual) + factor(race) + factor(sexM) +
              age_corrected + I(age_corrected^2) + 
              education + poverty + PIR +
              pct_blk + hispanic + popdensity + pct_owner_occ,
            data = qid_dat, weights = qid_dat[,IPWcr],
            family = binomial,
            nthreads = n_threads, chunk.size = 5000, 
            control = gam.control(trace = TRUE))
  summary(m3)
  
  pred_k <- sort(c(pm25 = 5, no2 = 10, ozone = 60, epa_pm25 = 12, epa_no2 = 53, epa_ozone = 70,
                   quantile(exp_c, c(0, .005, .1, .25, .5, .75, .9, .995, 1))))
  cp_dlnm_50 <- crosspred(cb_dlnm, m3, cen = pred_k["50%"], at = pred_k, bylag = 0.2)
  cp_dlnm_10 <- crosspred(cb_dlnm, m3, cen = pred_k["10%"], at = pred_k, bylag = 0.2)
  cp_dlnm_005 <- crosspred(cb_dlnm, m3, cen = pred_k["0.5%"], at = pred_k, bylag = 0.2)
  cp_dlnm_who <- crosspred(cb_dlnm, m3, cen = pred_k[exp], at = pred_k, bylag = 0.2)
  cp_dlnm_epa <- crosspred(cb_dlnm, m3, cen = pred_k[paste0("epa_", exp)], at = pred_k, bylag = 0.2)
  plot(cp_dlnm_005, "contour", main = exp, xlab = "Exposure concentration", ylab = "Years prior")
  
  save(cp_dlnm_epa, cp_dlnm_who, cp_dlnm_005, cp_dlnm_10, cp_dlnm_50, pred_k,
       file = paste0(dir_data, "analysis/gamns_dlnmlog35_expdist15_comprisk_who-epa_contUSsamp_", AD_ADRD, "_", 
                     code_type, "_", exp, ".rda"))
  
}








##### Sensitivity: 1 yr exposure #####
library(mgcv)
library(dlnm)
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                           code_type, "_qid_clean.fst"),
                    as.data.table = TRUE)
id10 <- which(qid_dat$year == 2010)
for (exp in c("pm25", "no2", "ozone", "oxidant")) {
  cat("\n", exp)
  exp_dat <- as.matrix(read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                                       code_type, "_", exp, "_clean.fst")))[,1:10]
  exp_iqr <- IQR(exp_dat[id10, 1:10])
  
  prD <- bam(dead ~ exp_dat[,1] + factor(year),
             data = qid_dat, family = binomial,
             nthreads = 1, samfrac = 0.05, chunk.size = 5000,
             control = gam.control(trace = TRUE))
  prD_L <- bam(dead ~ exp_dat[,1] +
                 factor(year) + factor(statecode) +
                 age_corrected + I(age_corrected^2) +
                 factor(dual) + factor(race) + factor(sexM) +
                 education + poverty + PIR +
                 pct_blk + hispanic + popdensity + pct_owner_occ,
               data = qid_dat, family = binomial,
               nthreads = 1, samfrac = 0.05, chunk.size = 5000,
               control = gam.control(trace = TRUE))
  qid_dat[, IPWcr := (1 - prD$fitted.values) / (1 - prD_L$fitted.values)]
  qid_dat[IPWcr > 10, IPWcr := 10]
  rm(prD, prD_L)
  
  m3 <- bam(first_hosp ~ exp_dat[,1] +
              factor(year) + factor(statecode) +
              factor(dual) + factor(race) + factor(sexM) +
              age_corrected + I(age_corrected^2) + 
              education + poverty + PIR +
              pct_blk + hispanic + popdensity + pct_owner_occ,
            data = qid_dat, weights = qid_dat$IPWcr, 
            family = binomial,
            nthreads = 1, samfrac = 0.05, chunk.size = 5000,
            control = gam.control(trace = TRUE))
  s3 <- summary(m3)
  res <- list(est = s3$p.coeff, se = s3$se, iqr = exp_iqr)
  cat(", HO = ", exp(m3$coefficients[2] * exp_iqr))
  save(res, exp_iqr,
       file = paste0(dir_data, "analysis/gam_1yr_comprisk_contUS_", 
                     AD_ADRD, "_", code_type, "_", exp, ".rda"))
  rm(cb_dlm, m3); gc()
}


##### Sensitivity: TDLM #####
library(mgcv)
# library(devtools)
# install_github("danielmork/dlmtree")
library(dlmtree)
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

qid_dat <- read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                           code_type, "_qid_clean.fst"),
                    as.data.table = TRUE)

set.seed(5165)
n <- qid_dat[year == 2010, .N]
samp_size <- floor(.25 * n)
qid_samp <- qid_dat[year == 2010, qid
                    ][sample(n, samp_size)]
idx <- which(qid_dat$qid %in% qid_samp)
length(idx) # 10569701
#qid_dat <- qid_dat[idx]
gc()

ip <- bam(first_hosp ~ 
            factor(year) + factor(statecode) +
            factor(dual) + factor(race) + factor(sexM) +
            age_corrected +
            education + poverty + PIR +
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
          family = binomial,
          control = gam.control(trace = TRUE))$coef

for (exp in c("pm25", "no2", "ozone")) {
  exp_dat <- as.matrix(
    read_fst(paste0(dir_data, "analysis/contUS_", AD_ADRD, "_", 
                    code_type, "_", exp, "_clean.fst")))[idx, 1:10]
  m4 <- tdlnm(first_hosp ~ 
                factor(year) + factor(statecode) +
                factor(dual) + factor(race) + factor(sexM) +
                age_corrected + 
                education + poverty + PIR +
                pct_blk + hispanic + popdensity + pct_owner_occ,
              data = qid_dat[idx], 
              exposure.data = exp_dat, 
              exposure.splits = 0,
              n.trees = 10, n.burn = 500, n.iter = 1000, n.thin = 4,
              family = "logit", initial.params = ip)
  save(ip, m4,
       file = paste0(dir_data, "analysis/tdlm_contUS_samp25per_", 
                     AD_ADRD, "_", code_type, "_", exp, ".rda"))
  print(s4 <- summary(m4))
  print(plot(s4))
}








