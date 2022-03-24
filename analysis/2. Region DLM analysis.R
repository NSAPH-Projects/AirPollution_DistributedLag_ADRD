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
# no2_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
#                            code_type, "_no2.fst"), as.data.table = TRUE)
# ozone_dat <- read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, "_",
#                              code_type, "_ozone.fst"), as.data.table = TRUE)
# all.equal(qid_dat$qid, pm25_dat$qid)

##### 2. Restrict to complete follow-up 2010 through 2016, no deaths #####
pm25_dat <- as.matrix(pm25_dat[, -c(1:3, 20)])[, 16:1] # reorder exposures, remove keys and lag0
idx <- which(rowSums(is.na(pm25_dat[, 1:10])) == 0 & # no missing exposures lags 1-10
               rowSums(pm25_dat[, 1:10] == 0) == 0 & # no zero exposures
               qid_dat$entry == 2000 & # follow up from 2000
               qid_dat$year > 2009) # year 2010 and later for 10 lag years exposures
qid_dat <- qid_dat[idx]
pm25_dat <- pm25_dat[idx,]

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

##### Remove QID after missing year of data #####
idx3 <- which(qid_dat$year == 2010)
qid_rm <- qid_dat[idx3, qid]
for (yr in 2011:2016) {
  idx3 <- c(idx3, which(qid_dat$year == yr & qid_dat$qid %in% qid_rm))
  qid_rm <- qid_dat[idx3][year == yr, qid]
}
qid_dat <- qid_dat[idx3]
pm25_dat <- pm25_dat[idx3,]
setkey(qid_dat, year, zip, qid)
rm(idx, idx2, idx3, qid_rm)

##### Summary stats #####
qid_dat[, .N]
qid_dat[, uniqueN(qid)]
qid_dat[, uniqueN(qid), by = year]
qid_dat[, sum(first_hosp)]
qid_dat[, .(sum(first_hosp), sum(dead)), by = year]
qid_dat[, sum(first_hosp)/uniqueN(qid)]
qid_dat[, sum(dead)/uniqueN(qid)]
qid_dat[, .(sum(first_hosp)/uniqueN(qid), sum(dead)/uniqueN(qid)), by = year]

##### 3. DLM analysis - TDLM #####
library(dlmtree)
(splits <- quantile(pm25_dat[, 1:10], 1:19/20))
m <- tdlnm(first_hosp ~ factor(year) - 1 +
             age_corrected + I(age_corrected^2) + 
             factor(dual) + factor(race) + factor(sexM) +
             education + poverty + PIR +
             pct_blk + hispanic + popdensity + pct_owner_occ,
           data = qid_dat,
           exposure.data = pm25_dat[, 1:10], 
           exposure.splits = splits,
           exposure.se = sd(pm25_dat[, 1:10]) / 2,
           family = "logit",
           n.trees = 10, n.burn = 2000, n.iter = 5000, n.thin = 5, 
           max.threads = 8)
save(m, file = paste0(dir_data, "analysis/tdlnm_region", region, "_", AD_ADRD, "_", 
                      code_type, "_pm25.rda"))
(s <- summary(m, cenval = m$Xsplits[1]))
plot(s)


##### 4. DLMM analysis - TDLMM #####
mm <- tdlmm(first_hosp ~ age + dual + race + sexM, 
            data = qid_dat,
            exposure.data = list(pm25 = pm25_dat, no2 = no2_dat, ozone = ozone_dat), 
            mixture.interactions = "noself",
            family = "logit",
            subset = which(qid_dat$year == 2016), mix.prior = 0.132,
            n.trees = 10, n.burn = 5000, n.iter = 10000, n.thin = 5)
save(mm, file = paste0(dir_data, "analysis/tdlmm_region", region, "_", AD_ADRD, "_", 
                       code_type, "_pm25_no2_ozone.rda"))
(ss <- summary(mm))
plot(ss)


##### 5. DLM analysis - GAM #####
library(dlnm)
library(mgcv)
cb <- crossbasis(pm25_dat[, 1:10], c(1, 10),
                 argvar = list(fun = "lin"), 
                 arglag = list(fun = "cr"))
cb_pen <- cbPen(cb)
m2 <- bam(first_hosp ~ factor(year) - 1 +
            cb + 
            age_corrected + I(age_corrected^2) + 
            factor(dual) + factor(race) + factor(sexM) +
            education + poverty + PIR +
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
          paraPen = list(cb = cb_pen),
          family = binomial,
          nthreads = 16,
          control = gam.control(trace = TRUE))
summary(m2)
cp <- crosspred(cb, m2, cen = 5, at = 5:15, bylag = 0.2)
plot(cp, "slices", var = 6)
plot(cp, "overall")
