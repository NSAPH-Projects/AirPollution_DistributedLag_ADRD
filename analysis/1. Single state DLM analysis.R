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
setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
pm25_dat <- read_fst(paste0(dir_data, "analysis/state", state, "_", AD_ADRD, "_", 
                            code_type, "_pm25.fst"), as.data.table = TRUE)

##### 2. Restrict to complete follow-up 2009 through 2016, no deaths #####
pm25_dat <- as.matrix(pm25_dat[, -c(1:3)])[, 17:1] # reorder exposures, remove keys
idx <- which(rowSums(is.na(pm25_dat[, 1:10])) == 0 & # no missing exposures
               rowSums(pm25_dat[, 1:10] == 0) == 0 & # no zero exposures
               qid_dat$entry == 2000) # follow up from 2000
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
idx3 <- which(qid_dat$year == 2009)
qid_rm <- qid_dat[year == 2009, qid]
for (yr in 2010:2016) {
  idx3 <- c(idx3, which(qid_dat$year == yr & qid_dat$qid %in% qid_rm))
  qid_rm <- qid_dat[idx3][year == yr, qid]
}
qid_dat <- qid_dat[idx3]
pm25_dat <- pm25_dat[idx3, ]
setkey(qid_dat, year, zip, qid)
rm(idx, idx2, idx3, qid_rm)

##### Summary stats #####
qid_dat[, .N]
qid_dat[, .N, by = year]
qid_dat[, uniqueN(qid)]
qid_dat[, uniqueN(qid), by = year]
qid_dat[, sum(first_hosp)]
qid_dat[year==2009,mean(age)]
qid_dat[(first_hosp),mean(age)]
qid_dat[, .(sum(first_hosp), sum(dead)), by = year]
qid_dat[, sum(first_hosp)/uniqueN(qid)]
qid_dat[, sum(dead)/uniqueN(qid)]
qid_dat[, .(sum(first_hosp)/uniqueN(qid), sum(dead)/uniqueN(qid)), by = year]

##### 3. DLM analysis - TDLM #####
library(dlmtree)
(splits <- quantile(pm25_dat[,1:10], c(0.01, 0.05, 1:9/10, 0.95, 0.99)))
m <- tdlnm(first_hosp ~ factor(year) - 1 + factor(fips) +
             age_corrected + I(age_corrected^2) + dual + race + sexM +
             education + poverty +
             I(medianhousevalue/medhouseholdincome) +
             pct_blk + hispanic + popdensity + pct_owner_occ,
           data = qid_dat,
           exposure.data = pm25_dat[,1:10], 
           exposure.splits = splits,
           exposure.se = sd(pm25_dat[,1:10]) / 2,
           family = "logit",
           n.trees = 10, n.burn = 2000, n.iter = 5000, n.thin = 5, 
           max.threads = 6)
save(m, file = paste0(dir_data, "analysis/tdlnm_state", state, "_", AD_ADRD, "_", 
                      code_type, "_pm25.rda"))
(s <- summary(m, cenval = 5))#, cenval = splits[3]))
plot(s)
plot(s, "slice", time = 1)
plot(s, "slice", val = 15)
plot(s, "cumulative")


##### 4. DLMM analysis - TDLMM #####
mm <- tdlmm(first_hosp ~ age + I(age^2) + dual + race + sexM +
              education + poverty +
              I(medianhousevalue/medhouseholdincome) +
              pct_blk + hispanic + popdensity + pct_owner_occ,
            data = qid_dat,
            exposure.data = list(pm25 = pm25_dat, no2 = no2_dat, ozone = ozone_dat), 
            mixture.interactions = "noself",
            family = "logit",
            subset = which(complete.cases(qid_dat)),
            mix.prior = 0.132,
            n.trees = 10, n.burn = 2000, n.iter = 5000, n.thin = 2)
save(mm, file = paste0(dir_data, "analysis/tdlmm_state", state, "_", AD_ADRD, "_", 
                       code_type, "_pm25_no2_ozone.rda"))
(ss <- summary(mm))
plot(ss)


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
            education + poverty + PIR + # price to income ratio
            pct_blk + hispanic + popdensity + pct_owner_occ,
          data = qid_dat, 
          paraPen = list(cb = cb_pen),
          family = binomial,
          nthreads = 8,
          control = gam.control(trace = TRUE))
summary(m2)
cp <- crosspred(cb, m2, cen = 5, at = seq(5, 15, by = 0.1), bylag = 0.2)
plot(cp, "slices", var = 10)
plot(cp, "overall")




## 