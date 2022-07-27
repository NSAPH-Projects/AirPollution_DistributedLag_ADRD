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
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any


library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 2)
Sys.setenv("OMP_NUM_THREADS" = 2)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

##### 1. Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/DLMM_agg/contUS_agg_", 
                           AD_ADRD, "_", code_type, "_agg_zip.fst"), as.data.table = TRUE)
exp_names <- c(paste0("pm25comp_", c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
                                   "no3", "oc", "pb", "si", "so4", "v", "z")))
lag_cols <- paste0("lag", 1:10)
exposures <- list()
for (e in exp_names) {
  exposures[[e]] <- as.matrix(
    read_fst(paste0(dir_data, "analysis/DLMM_agg/contUS_agg_", AD_ADRD, "_", 
                    code_type, "_", e, ".fst"), columns = lag_cols,
             as.data.table = TRUE))
}
exp_names <- c("br", "ca", "cu", "ec", "fe", "k", "nh4", "ni",
               "no3", "oc", "pb", "si", "so4", "v", "z")
names(exposures) <- exp_names

##### Remove rows with missing exposure data #####
idx <- 1:nrow(qid_dat)
for (e in exp_names) {
  idx <- intersect(idx, which(rowSums(is.na(exposures[[e]])) == 0))
}
qid_dat <- qid_dat[idx]
exposures <- lapply(names(exposures), function(e) exposures[[e]][idx,])
names(exposures) <- exp_names

##### Summary stats #####
qid_dat[, .N]
hist(qid_dat$lograte)

##### State rate data #####
qid_dat[, wt := at_risk / sum(at_risk), by = .(year, state, sexM, age_cat, dual_any, race)]
state_dat <- qid_dat[, .(at_risk = sum(at_risk), n_hosp = sum(n_hosp), 
                         rate = sum(n_hosp) / sum(at_risk),
                         avg_age = sum(wt * avg_age),
                         his = sum(wt * his), blk = sum(wt * blk),
                         PIR = sum(wt * PIR), pov = sum(wt * pov), ed = sum(wt * ed),
                         oo = sum(wt * oo), dens = sum(wt * dens),
                         s_tmmx = sum(wt * s_tmmx), w_tmmx = sum(wt * w_tmmx),
                         s_rmax = sum(wt * s_rmax), w_rmax = sum(wt * w_rmax),
                         lat = sum(wt * lat), long = sum(wt * long), 
                         state = first(state), region = first(region)),
                     by = .(year, state, sexM, age_cat, dual_any, race)]
state_dat[rate == 0, rate := min(state_dat$rate[state_dat$rate > 0])]
state_dat[, lograte := log(rate)]
state_exp <- lapply(exposures, function(e) matrix(0, state_dat[,.N], 10))
for (i in 1:nrow(state_dat)) {
  if (i %% 1000 == 0) cat(".")
  idx <- which(qid_dat$state == state_dat$state[i] &
                 qid_dat$year == state_dat$year[i] &
                 qid_dat$sexM == state_dat$sexM[i] &
                 qid_dat$age_cat == state_dat$age_cat[i] &
                 qid_dat$dual_any == state_dat$dual_any[i] &
                 qid_dat$race == state_dat$race[i])
  for (e in names(exposures)) {
    state_exp[[e]][i, ] <- crossprod(qid_dat$wt[idx], exposures[[e]][idx,,drop = F])
  }
}


#### Lag 1 correlations ####
l1exp <- sapply(exposures, function(e) e[,1])
l1cor <- cor(l1exp)
corrplot::corrplot(l1cor, type = "lower", method = "square")


##### 4. DLMM analysis - TDLMM #####
library(dlmtree)
form <- as.formula(lograte ~ factor(age_cat) + sexM + dual_any + factor(race) + 
                     factor(year) + factor(state) + his + blk + PIR + pov + ed + 
                     oo + dens + avg_age)
mm <- tdlmm(form,
            data = qid_dat,
            exposure.data = exposures, # 15 pm components
            mixture.interactions = "noself",
            shrinkage = "exposures",
            mix.prior = 1, # prior inc prob = 0.5
            n.trees = 20, # limited to 20 exposures, 10 interactions
            n.burn = 200, n.iter = 500, n.thin = 5)
save(mm, file = paste0(dir_data, "analysis/DLMM_agg/tdlmm_100atrisk_", AD_ADRD, "_",
                       code_type, "_pm25comp.rda"))
(ss <- summary(mm))
plot(rowMeans(mm$termNodes))
plot(rowMeans(mm$termNodes2))
plot(ss)

##### Adjusting for changes in co-exposures ######
# idx <- which(qid_dat$qid %in% qid_samp & qid_dat$year == 2010)
exposureDat <- do.call(cbind.data.frame, 
                       lapply(exposures, function(e) c(e[,1])))
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
      formula <- as.formula(paste(respExp, "~s(", predExp, ", k = 3, bs = 'cr')"))
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
  if (length(modelSummary$MIX) > 0) {
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
#rm(modelSummary, tdlmmResult)




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
