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
setDTthreads(threads = 16)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
pm25_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_", 
                            code_type, "_pm25.fst"), as.data.table = TRUE)
no2_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
                           code_type, "_no2.fst"), as.data.table = TRUE)
ozone_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
                             code_type, "_ozone.fst"), as.data.table = TRUE)
tmmx_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
                             code_type, "_tmmx.fst"), as.data.table = TRUE)
rmax_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
                             code_type, "_rmax.fst"), as.data.table = TRUE)
pr_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
                             code_type, "_pr.fst"), as.data.table = TRUE)
pm25_dat <- as.matrix(pm25_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
no2_dat <- as.matrix(no2_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
ozone_dat <- as.matrix(ozone_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
tmmx_dat <- as.matrix(tmmx_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
rmax_dat <- as.matrix(rmax_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
pr_dat <- as.matrix(pr_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0


##### 2. Restrict to complete follow-up 2010 through 2016, no deaths #####
idx <- which(rowSums(is.na(pm25_dat[, 1:10])) == 0 & # no missing exposures lags 1-10
               # rowSums(is.na(no2_dat[, 1:10])) == 0 &
               # rowSums(is.na(ozone_dat[, 1:10])) == 0 &
               # rowSums(is.na(tmmx_dat[, 1:10])) == 0 &
               # rowSums(is.na(rmax_dat[, 1:10])) == 0 &
               # rowSums(is.na(pr_dat[, 1:10])) == 0 &
               qid_dat$entry == 2000 & # follow up from 2000
               qid_dat$year > 2009) # year 2010 and later for 10 lag years exposures
qid_dat <- qid_dat[idx]
pm25_dat <- pm25_dat[idx,]
no2_dat <- no2_dat[idx,]
ozone_dat <- ozone_dat[idx,]
tmmx_dat <- tmmx_dat[idx,]
rmax_dat <- rmax_dat[idx,]
pr_dat <- pr_dat[idx,]

##### Data fixes #####
qid_dat[, age_corrected := entry_age + year - entry]
qid_dat[, PIR := medianhousevalue/medhouseholdincome]
# clear missing data, infinite PIR
idx2 <- which(rowSums(pm25_dat[, 1:10] == 0) == 0 &
                !is.infinite(qid_dat$PIR) &
                complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                        education, poverty, pct_blk, hispanic,
                                        popdensity, pct_owner_occ, PIR)]))
qid_dat <- qid_dat[idx2]
pm25_dat <- pm25_dat[idx2,]
no2_dat <- no2_dat[idx2,]
ozone_dat <- ozone_dat[idx2,]
tmmx_dat <- tmmx_dat[idx2,]
rmax_dat <- rmax_dat[idx2,]
pr_dat <- pr_dat[idx2,]

##### Remove QID after missing year of data #####
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
tmmx_dat <- tmmx_dat[idx3,]
rmax_dat <- rmax_dat[idx3,]
pr_dat <- pr_dat[idx3,]
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

##### 3. DLM analysis - TDLM #####
library(dlmtree)
library(mgcv)
form <- as.formula(first_hosp ~ factor(year) +
                     age_corrected + I(age_corrected^2) + 
                     latitude * longitude +
                     factor(dual) + factor(race) + factor(sexM) +
                     education + poverty + PIR +
                     pct_blk + hispanic + popdensity + pct_owner_occ)
init.params <- bam(form,
                   data = qid_dat, 
                   family = binomial,
                   nthreads = 8, samfrac = 0.05, chunk.size = 5000,
                   control = gam.control(trace = TRUE))$coef
pm_range <- quantile(pm25_dat[,1:10], c(0.005, 0.995))
(splits <- seq(from = pm_range[1], to = pm_range[2], length.out = 20))
set.seed(38173)
m <- tdlnm(form,
           data = qid_dat,
           exposure.data = pm25_dat[, 1:10], 
           exposure.splits = splits,
           exposure.se = sd(pm25_dat[, 1:10]) / 2,
           family = "logit", 
           n.trees = 10, n.burn = 1000, n.iter = 2500, n.thin = 2,
           initial.params = init.params)
save(m, file = paste0(dir_data, "analysis/tdlnm_region", region, "_", AD_ADRD, "_", 
                      code_type, "_pm25.rda"))
(s <- summary(m))#, cenval = 5))
plot(s)


##### 4. DLMM analysis - TDLMM #####
library(dlmtree)
library(mgcv)
form <- as.formula(first_hosp ~ factor(year) +
                     age_corrected + I(age_corrected^2) + 
                     factor(dual) + factor(race) + factor(sexM) +
                     education + poverty + PIR +
                     pct_blk + hispanic + popdensity + pct_owner_occ)
init.params <- bam(form,
                   data = qid_dat, 
                   family = binomial,
                   nthreads = 8, samfrac = 0.1, chunk.size = 5000,
                   control = gam.control(trace = TRUE))
set.seed(8372)
qid_samp <- qid_dat[year == 2010, qid][sample(qid_dat[year == 2010, .N], 200000)]
qid_dat[qid %in% qid_samp, .N] #994428
mm <- tdlmm(form,
            data = qid_dat,
            exposure.data = list(pm25 = pm25_dat[,1:10], 
                                 no2 = no2_dat[,1:10], 
                                 ozone = ozone_dat[,1:10],
                                 tmmx = tmmx_dat[,1:10],
                                 rmax = rmax_dat[,1:10],
                                 pr = pr_dat[,1:10]), 
            mixture.interactions = "all",
            family = "logit",
            subset = which(qid_dat$qid %in% qid_samp),
            mix.prior = 1.263, # prior inc prob = 0.9
            n.trees = 20, n.burn = 2000, n.iter = 5000, n.thin = 5,
            initial.params = init.params$coef)
save(mm, file = paste0(dir_data, "analysis/tdlmm_region", region, "_", AD_ADRD, "_", 
                       code_type, "samp200kqid_pm25_no2_ozone_tmmx_rmax_pr_allint.rda"))
(ss <- summary(mm))
plot(ss)

##### Exposure correlations #####
lag1_exp <- data.table(pm25 = pm25_dat[,1], no2 = no2_dat[,1], ozone = ozone_dat[,1])
cor(lag1_exp)


##### 5. DLM analysis - GAM #####
library(dlnm)
library(mgcv)
cb_pm <- crossbasis(pm25_dat[, 1:10], c(1, 10),
                    argvar = list(fun = "lin"), 
                    arglag = list(fun = "ps"))
cb_no <- crossbasis(no2_dat[, 1:10], c(1, 10),
                    argvar = list(fun = "lin"), 
                    arglag = list(fun = "ps"))
cb_oz <- crossbasis(ozone_dat[, 1:10], c(1, 10),
                    argvar = list(fun = "lin"), 
                    arglag = list(fun = "ps"))
cb_pm_pen <- cbPen(cb_pm)
cb_no_pen <- cbPen(cb_no)
cb_oz_pen <- cbPen(cb_oz)
m2 <- bam(first_hosp ~ factor(year) +
            cb_pm + #cb_no + cb_oz +
            age_corrected + I(age_corrected^2) + 
            factor(dual) + factor(race) + factor(sexM) +
            education + poverty + PIR +
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
         paraPen = list(cb_pm = cb_pm_pen),
                         # cb_no = cb_no_pen, 
                         # cb_oz = cb_oz_pen),
          family = binomial,
          nthreads = 16, samfrac = 0.1, chunk.size = 5000,
          control = gam.control(trace = TRUE))
summary(m2)
cp_pm <- crosspred(cb_pm, m2, cen = 5, at = 5:15, bylag = 0.2)
cp_no <- crosspred(cb_no, m2, cen = 6, at = 6:47, bylag = 0.2)
cp_oz <- crosspred(cb_oz, m2, cen = 28, at = 28:41, bylag = 0.2)
save(cp_pm, #cp_no, cp_oz, 
     file = paste0(dir_data, "analysis/gam-ps_region", region, "_", AD_ADRD, "_", 
                       code_type, "_pm25.rda"))
plot(cp_pm, "slices", var = 10, xlab = "Years prior", ylab = "Hazard Odds (PM2.5)")
plot(cp_no, "slices", var = 20, xlab = "Years prior", ylab = "Hazard Odds (NO2)")
plot(cp_oz, "slices", var = 36, xlab = "Years prior", ylab = "Hazard Odds (Ozone)")
plot(cp_pm, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")
plot(cp_no, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")
plot(cp_oz, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")


##### 5.2 DLNM analysis - GAM #####
library(dlnm)
library(mgcv)
cb_pm <- crossbasis(pm25_dat[, 1:10], c(1, 10),
                    argvar = list(fun = "ps", df = 5), 
                    arglag = list(fun = "ps", df = 5))
# cb_no <- crossbasis(no2_dat[, 1:10], c(1, 10),
#                     argvar = list(fun = "lin"), 
#                     arglag = list(fun = "ps"))
# cb_oz <- crossbasis(ozone_dat[, 1:10], c(1, 10),
#                     argvar = list(fun = "lin"), 
#                     arglag = list(fun = "ps"))
cb_pm_pen <- cbPen(cb_pm)
# cb_no_pen <- cbPen(cb_no)
# cb_oz_pen <- cbPen(cb_oz)
m2 <- bam(first_hosp ~ factor(year) +
            cb_pm + #cb_no + cb_oz +
            age_corrected + I(age_corrected^2) + 
            factor(dual) + factor(race) + factor(sexM) +
            education + poverty + PIR +
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
          #paraPen = list(cb_pm = cb_pm_pen), 
                         # cb_no = cb_no_pen, 
                         # cb_oz = cb_oz_pen),
          family = binomial,
          nthreads = 16, samfrac = 0.1, chunk.size = 5000,
          control = gam.control(trace = TRUE))
summary(m2)
cp_pm <- crosspred(cb_pm, m2, cen = 5, at = 5:15, bylag = 0.2)
# cp_no <- crosspred(cb_no, m2, cen = 6, at = 6:47, bylag = 0.2)
# cp_oz <- crosspred(cb_oz, m2, cen = 28, at = 28:41, bylag = 0.2)
save(cp_pm, #cp_no, cp_oz, 
     file = paste0(dir_data, "analysis/gam-ps_dlnm_region", region, "_", AD_ADRD, "_", 
                   code_type, "_pm25.rda"))
plot(cp_pm, "slices", var = 10, xlab = "Years prior", ylab = "Hazard Odds (PM2.5)")
# plot(cp_no, "slices", var = 20, xlab = "Years prior", ylab = "Hazard Odds (NO2)")
# plot(cp_oz, "slices", var = 36, xlab = "Years prior", ylab = "Hazard Odds (Ozone)")
plot(cp_pm, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")
# plot(cp_no, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")
# plot(cp_oz, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")


##### 6. Comp risk DLM analysis #####
library(dlnm)
library(mgcv)
cb <- crossbasis(pm25_dat[, 1:10], c(1, 10),
                 argvar = list(fun = "lin"), 
                 arglag = list(fun = "cr", df = 5))
cb_pen <- cbPen(cb)
prD <- bam(dead ~ factor(year) - 1 + cb,
           data = qid_dat, 
           paraPen = list(cb = cb_pen),
           family = binomial,
           nthreads = 16,
           control = gam.control(trace = TRUE))
prD_L <- bam(dead ~ factor(year) - 1 + cb +
               age_corrected + I(age_corrected^2) + 
               factor(dual) + factor(race) + factor(sexM) +
               education + poverty +
               I(medianhousevalue/medhouseholdincome) +
               pct_blk + hispanic + popdensity + pct_owner_occ,
             data = qid_dat, 
             paraPen = list(cb = cb_pen),
             family = binomial,
             nthreads = 16,
             control = gam.control(trace = TRUE))
qid_dat[, pD := 1 / (1-prD_L$fitted.values)]
m3 <- bam(first_hosp ~ factor(year) - 1 + 
            factor(region) +
            cb + 
            age_corrected + I(age_corrected^2) + 
            factor(dual) + factor(race) + factor(sexM) +
            education + poverty + PIR +
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, weights = qid_dat$pD,
          paraPen = list(cb = cb_pen),
          family = binomial,
          nthreads = 16,
          control = gam.control(trace = TRUE))
summary(m3)
cp <- crosspred(cb, m3, cen = 5, at = 2:25, bylag = 0.2)
save(cp, file = paste0(dir_data, "analysis/gamcr_comp-risk_contUS_", AD_ADRD, "_", 
                       code_type, ".rda"))
plot(cp, "slices", var = 10, xlab = "Years prior", ylab = "Hazard Odds")
plot(cp, "overall", xlab = "Exposure Concentration", ylab = "Cumulative Hazard")