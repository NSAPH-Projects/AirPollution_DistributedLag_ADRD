# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: single state distributed lag analysis
#' Inputs: merged denominator and exposure data                                
#' Outputs:       
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
state <- "MA"
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary, secondary, or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 4)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
pm25_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_pm25.fst"), as.data.table = TRUE)
no2_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_no2.fst"), as.data.table = TRUE)
ozone_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_ozone.fst"), as.data.table = TRUE)
tmmx_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_tmmx.fst"), as.data.table = TRUE)
rmax_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_rmax.fst"), as.data.table = TRUE)
pr_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_pr.fst"), as.data.table = TRUE)
pm25_dat <- as.matrix(pm25_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys
no2_dat <- as.matrix(no2_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys
ozone_dat <- as.matrix(ozone_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys
tmmx_dat <- as.matrix(tmmx_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys
rmax_dat <- as.matrix(rmax_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys
pr_dat <- as.matrix(pr_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys



##### 2. Restrict to complete follow-up 2009 through 2016, no deaths #####
idx <- which(rowSums(is.na(pm25_dat[, 1:10])) == 0 & # no missing exposures
               rowSums(is.na(no2_dat[, 1:10])) == 0 &
               rowSums(is.na(ozone_dat[, 1:10])) == 0 &
               rowSums(is.na(tmmx_dat[, 1:10])) == 0 &
               rowSums(is.na(rmax_dat[, 1:10])) == 0 &
               rowSums(is.na(pr_dat[, 1:10])) == 0 &
               qid_dat$entry == 2000) # follow up from 2000
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
idx2 <- which(complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                        education, poverty, medianhousevalue,
                                        medhouseholdincome, pct_blk, hispanic,
                                        popdensity, pct_owner_occ, PIR)]) &
                !is.infinite(qid_dat$PIR))
qid_dat <- qid_dat[idx2]
pm25_dat <- pm25_dat[idx2,]
no2_dat <- no2_dat[idx2,]
ozone_dat <- ozone_dat[idx2,]
tmmx_dat <- tmmx_dat[idx2,]
rmax_dat <- rmax_dat[idx2,]
pr_dat <- pr_dat[idx2,]

##### Remove QID after missing year of data #####
idx3 <- which(qid_dat$year == 2010)
qid_rm <- qid_dat[year == 2010, qid]
for (yr in 2011:2016) {
  idx3 <- c(idx3, which(qid_dat$year == yr & 
                          qid_dat$qid %in% qid_rm))
  qid_rm <- qid_dat[idx3][year == yr, qid]
}
qid_dat <- qid_dat[idx3]
pm25_dat <- pm25_dat[idx3,]
no2_dat <- no2_dat[idx3,]
ozone_dat <- ozone_dat[idx3,]
tmmx_dat <- tmmx_dat[idx3,]
rmax_dat <- rmax_dat[idx3,]
pr_dat <- pr_dat[idx3,]
setkey(qid_dat, year, zip, qid)
rm(idx, idx2, idx3, qid_rm)

##### Summary stats #####
qid_dat[, .N]
qid_dat[, uniqueN(qid)]
qid_dat[, uniqueN(qid), by = year]
qid_dat[, sum(first_hosp)]
qid_dat[year==2010,mean(age)]
qid_dat[(first_hosp),mean(age)]
qid_dat[, .(sum(first_hosp), sum(dead)), by = year]
qid_dat[, sum(first_hosp)/uniqueN(qid)]
qid_dat[, sum(dead)/uniqueN(qid)]
qid_dat[, .(sum(first_hosp)/uniqueN(qid), sum(dead)/uniqueN(qid)), by = year]

##### 3. DLM analysis - TDLM #####
library(dlmtree)
pm_range <- range(pm25_dat[,1:10])
(splits <- seq(from = 4, to = 25, length.out = 32))
m <- tdlnm(first_hosp ~ factor(year) - 1 +
             age_corrected + I(age_corrected^2) + 
             factor(dual) + factor(race) + factor(sexM) +
             education + poverty + PIR +
             pct_blk + hispanic + popdensity + pct_owner_occ,
           data = qid_dat,
           exposure.data = pm25_dat[,1:10], 
           exposure.splits = splits,
           exposure.se = sd(pm25_dat[,1:10]) / 2,
           family = "logit",
           n.trees = 10, n.burn = 2000, n.iter = 5000, n.thin = 5, 
           max.threads = 2)
save(m, file = paste0(dir_data, "analysis/tdlnm_state", state, "_", AD_ADRD, "_", 
                      code_type, "_pm25.rda"))
(s <- summary(m, cenval = 5))#, cenval = splits[3]))
plot(s)
plot(s, "slice", time = 9)
plot(s, "slice", val = 14)
plot(s, "cumulative")


##### 4. DLMM analysis - TDLMM #####
library(dlmtree)
library(mgcv)
prD_L <- bam(dead ~ factor(year) - 1 +
               age_corrected + I(age_corrected^2) + 
               factor(dual) + factor(race) + factor(sexM) +
               education + poverty + PIR +
               pct_blk + hispanic + popdensity + pct_owner_occ,
             data = qid_dat, 
             family = binomial)
prD <- bam(dead ~ factor(year) - 1,
           data = qid_dat, 
           family = binomial)
qid_dat[, pD := predict(prD, newdata = .SD, type = "response") /
          predict(prD_L, newdata = .SD, type = "response")]
mm <- tdlmm(first_hosp ~ factor(year) - 1 + pD +
              age_corrected + I(age_corrected^2) + 
              factor(dual) + factor(race) + factor(sexM) +
              education + poverty +
              I(medianhousevalue/medhouseholdincome) +
              pct_blk + hispanic + popdensity + pct_owner_occ,
            data = qid_dat,
            exposure.data = list(pm25 = pm25_dat[,1:10], 
                                 no2 = no2_dat[,1:10], 
                                 ozone = ozone_dat[,1:10],
                                 tmmx = tmmx_dat[,1:10],
                                 rmax = rmax_dat[,1:10],
                                 pr = pr_dat[,1:10]), 
            mixture.interactions = "noself",
            family = "logit",
            mix.prior = 0.167, # Prior inc prob = 0.5
            n.trees = 20, n.burn = 5000, n.iter = 10000, n.thin = 5)
save(mm, file = paste0(dir_data, "analysis/tdlmm_state", state, "_", AD_ADRD, "_", 
                       code_type, "_pm25_no2_ozone_tmmx_rmax_pr.rda"))
(ss <- summary(mm, conf.level = .95))
plot(ss)
mix_cor <- cor(data.frame(pm25 = pm25_dat[,1], no2 = no2_dat[,1], ozone = ozone_dat[,1],
                      tmax = tmmx_dat[,1], rmax = rmax_dat[,1], pr = pr_dat[,1]))
corrplot::corrplot(mix_cor, type = "upper")

##### 5. DLM analysis - GAM #####
library(dlnm)
library(mgcv)
cb <- crossbasis(pm25_dat[,1:10], c(0, 9),
                 argvar = list(fun = "lin"), 
                 arglag = list(fun = "ps"))
cb_pen <- cbPen(cb)
m2 <- bam(first_hosp ~ factor(year) - 1 + 
            s(fips, bs = "re") + # random effect of county
            age + I(age^2) + 
            factor(dual) + factor(race) + factor(sexM) +
            cb + 
            pD +
            education + poverty + PIR + # price to income ratio
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
          paraPen = list(cb = cb_pen),
          family = binomial,
          nthreads = 8,
          control = gam.control(trace = TRUE))
summary(m2)
cp <- crosspred(cb, m2, cen = 5, at = seq(2, 25, by = 0.5), bylag = 0.2)
plot(cp, "slices", var = 10)
plot(cp, "overall")




## 