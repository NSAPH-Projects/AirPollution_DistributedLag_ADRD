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
                    "no2", "ozone", "tmmx", "rmax")
lag_cols <- paste0("lag", 1:10)
exposures <- list()
for (e in exp_names) {
  exposures[[e]] <- as.matrix(
    read_fst(paste0(dir_data, "analysis/region", region, "_", AD_ADRD, 
                    "_", code_type, "_", e, ".fst"), columns = lag_cols,
             as.data.table = TRUE))
}
exp_names <- c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
               "no3", "oc", "pb", "si", "so4", "v", "z",
               "no2", "ozone", "tmmx", "rmax")
names(exposures) <- exp_names


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
qid_samp <- qid_dat[year == 2010, qid][sample(qid_dat[year == 2010, .N], 1000000)]
qid_dat[qid %in% qid_samp, .N] # 4972500

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
            exposure.data = exposures, # 15 pm components, no2, ozone, tmmx, rmax
            mixture.interactions = "noself",
            family = "logit", shrinkage = "exposures",
            subset = which(qid_dat$qid %in% qid_samp),
            mix.prior = 0.347, # prior inc prob = 0.5
            n.trees = 20, 
            n.burn = 1000, n.iter = 2000, n.thin = 4,
            initial.params = init.params)
save(mm, file = paste0(dir_data, "analysis/tdlmm1M_20_0.5_region", 
                       region, "_", AD_ADRD, "_", 
                       code_type, "_pm25comp.rda"))
(ss <- summary(mm))
plot.summary.tdlmm(ss)

##### Adjusting for changes in co-exposures ######
idx <- which(qid_dat$qid %in% qid_samp & qid_dat$year == 2010)
exposureDat <- do.call(cbind.data.frame, 
                       lapply(exposures, function(e) c(e[idx,])))
expMean <- colMeans(exposureDat) # Empirical means
expIQRLevels <-                   # 25/75 percentiles of each exposure
  lapply(exposureDat, function(i) quantile(i, probs = c(0.25, 0.75)))
do.call(cbind.data.frame, expIQRLevels)

# predLevels[[predictor]][[response]] =
#  use level of 'predictor' to estimate 'response'.
predLevels <- lapply(exposureDat, function(i) list())
# Loop over pairs of exposures
for (predExp in names(predLevels)) {
  for (respExp in names(predLevels)) {
    cat("\nPredicting", respExp, "at 25/75 percentiles of", predExp)
    if (predExp == respExp) { 
      # If predictor/response exposure same, set to 25/75 percentiles.
      predLevels[[predExp]][[respExp]] <- expIQRLevels[[predExp]]
      
    } else { 
      # Predict `response` (respExp) at 25/75 percentiles of `predictor` (predExp).
      # Other methods could be substituted here to achieve prediction
      formula <- as.formula(paste(respExp, "~s(", predExp, ", k = 5, bs = 'cr')"))
      model <- bam(formula, dat = exposureDat)
      newdat <- data.frame(expIQRLevels[[predExp]])
      names(newdat) <- predExp
      predLevels[[predExp]][[respExp]] <- predict(model, newdata = newdat)
    }
    cat(":", predLevels[[predExp]][[respExp]])
  }
}


modelSummary <- summary(mm, keep.mcmc = TRUE, verbose = FALSE)
expDLM <- list() # For output
max_t <- 10
for (exposure in names(modelSummary$DLM)) {
  # MCMC results for main effect of interest
  mcmc <- modelSummary$DLM[[exposure]]$mcmc * diff(expIQRLevels[[exposure]])
  # Add MCMC results for expected changes in co-exposures
  for (coexposure in names(modelSummary$DLM)) {
    if (exposure != coexposure)
      mcmc <- mcmc + 
        modelSummary$DLM[[coexposure]]$mcmc * diff(predLevels[[exposure]][[coexposure]])
  }
  # Add MCMC results due to interactions: 
  # 1) same time interactions, use expected changes in main and co-exposure
  # 2) diff time interactions, use expected change in main, mean in co-exposure
  # 3) same time interactions (b/t co-exp), use expected change in co-exposures
  # 4) diff time interactions (b/t co-exp), cancels (set at mean)
  for (mix in 1:length(modelSummary$MIX)) {
    # Interaction between main and co-exposure
    if (modelSummary$MIX[[mix]]$rows == exposure) {
      mcmc <- mcmc +
        # add same time interactions using expected changes
        t(sapply(1:max_t, function(k) {
          modelSummary$MIX[[mix]]$mcmc[k,k,]
        })) * 
        (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1]) +
        # add different time interactions using IQR change in primary exposure
        # (other exposures set to mean)
        t(sapply(1:max_t, function(k) {
          colSums(modelSummary$MIX[[mix]]$mcmc[k,-k,]) # interaction sum removing time k
        })) * 
        diff(predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]]) *
        expMean[modelSummary$MIX[[mix]]$cols]
      # Interaction between main and co-exposure
    } else if (modelSummary$MIX[[mix]]$cols == exposure) {
      mcmc <- mcmc +
        # add same time interactions using expected changes
        t(sapply(1:max_t, function(k) {
          modelSummary$MIX[[mix]]$mcmc[k,k,]
        })) * 
        (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1]) +
        # add different time interactions using IQR change in primary exposure
        # (other exposure set to mean)
        t(sapply(1:max_t, function(k) {
          colSums(modelSummary$MIX[[mix]]$mcmc[-k,k,]) # interaction sum removing time k
        })) * 
        diff(predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]]) *
        expMean[modelSummary$MIX[[mix]]$rows]
      # Interaction between two co-exposures
    } else {
      mcmc <- mcmc +
        # add same time interactions using expected changes
        t(sapply(1:max_t, function(k) {
          modelSummary$MIX[[mix]]$mcmc[k,k,]
        })) * 
        (predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][2] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][2] -
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$rows]][1] * 
           predLevels[[exposure]][[modelSummary$MIX[[mix]]$cols]][1])
      # no different time interactions (they are set at mean, so equals zero)
    }
  }
  # Summarize main effects, adjusted for changes in co-exposures
  expDLM[[exposure]] <- 
    data.frame("Name" = exposure,
               "Week" = 1:modelSummary$nLags,
               "Effect" = apply(mcmc, 1, mean),
               # Credible interval at levels 0.025 and 0.975
               "Lower" =  apply(mcmc, 1, quantile, probs = .025),
               "Upper" =  apply(mcmc, 1, quantile, probs = .975))
  expDLM[[exposure]]$CW <- # Critical window where CI does not contain zero
    (expDLM[[exposure]]$Lower > 0 | expDLM[[exposure]]$Upper < 0)
}
rm(modelSummary, tdlmmResult)




#### plot ####
library(ggplot2)
expDLM <- do.call(rbind.data.frame, expDLM)
ggplot(expDLM, aes(x = Week, y = exp(Effect), ymin = exp(Lower), ymax = exp(Upper))) +
  geom_hline(yintercept = 1, color = "red", size = 1) +
  geom_ribbon(fill = "grey", col = "grey", size = 2) +
  geom_line(size = 1) +
  facet_wrap(~Name, ncol = 6) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#DDDDDD")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Change in outcome for IQR change in exposure \nadjusting for expected changes in co-exposures") +
  xlab("Years prior") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


setDT(expDLM)
cumDLM <- expDLM[, .(Effect = sum(Effect), Lower = sum(Lower), Upper = sum(Upper)), by = Name]
ggplot(cumDLM, aes(x = Name, y = exp(Effect), ymin = exp(Lower), ymax = exp(Upper))) +
  geom_hline(yintercept = 1, color = "red") +
  geom_point(size = 2) +
  geom_errorbar(size = 1, width = 0.5) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0.2)) +
  labs(x = "Years prior", y = "Cumulative hazard odds")
