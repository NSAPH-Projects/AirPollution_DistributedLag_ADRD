# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: DLNM plots
#' Inputs: merged denominator and exposure data                                
#' Outputs:       
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 96 GB
# ############################################################################ #
rm(list = ls())
gc()

##### Setup #####
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any


library(data.table)
library(fst)
library(NSAPHutils)
library(mgcv)
library(dlnm)
library(ggplot2)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 1)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

##### Load DLNM results #####
load(paste0(dir_data, "analysis/gamns_dlnmlog34_comprisk_contUSsamp_", AD_ADRD, "_", 
            code_type, "_pm25.rda"))
pm_dlnm <- cp_dlnm
pm_levs <- k

load(paste0(dir_data, "analysis/gamns_dlnmlog34_comprisk_contUSsamp_", AD_ADRD, "_", 
            code_type, "_no2.rda"))
no_dlnm <- cp_dlnm
no_levs <- k

load(paste0(dir_data, "analysis/gamns_dlnmlog34_comprisk_contUSsamp_", AD_ADRD, "_", 
            code_type, "_ozone.rda"))
oz_dlnm <- cp_dlnm
oz_levs <- k


##### Create data table of results #####
lag_vals <- seq(1, 10, by = 0.2)
exp_dt <- rbindlist(list(
  rbindlist(lapply(2:6, function(q) {
    data.table(exp = "PM2.5",
               lag = lag_vals,
               quantile = names(pm_levs)[q],
               val = pm_levs[q],
               RR = pm_dlnm$matRRfit[paste0(pm_levs[q]),],
               RRlow = pm_dlnm$matRRlow[paste0(pm_levs[q]),],
               RRhigh = pm_dlnm$matRRhigh[paste0(pm_levs[q]),],
               cumRR = pm_dlnm$allRRfit[paste0(pm_levs[q])],
               cumRRlow = pm_dlnm$allRRlow[paste0(pm_levs[q])],
               cumRRhigh = pm_dlnm$allRRhigh[paste0(pm_levs[q])]) 
  })),
  rbindlist(lapply(2:6, function(q) {
    data.table(exp = "NO2",
               lag = lag_vals,
               quantile = names(no_levs)[q],
               val = no_levs[q],
               RR = no_dlnm$matRRfit[paste0(no_levs[q]),],
               RRlow = no_dlnm$matRRlow[paste0(no_levs[q]),],
               RRhigh = no_dlnm$matRRhigh[paste0(no_levs[q]),],
               cumRR = no_dlnm$allRRfit[paste0(no_levs[q])],
               cumRRlow = no_dlnm$allRRlow[paste0(no_levs[q])],
               cumRRhigh = no_dlnm$allRRhigh[paste0(no_levs[q])]) 
  })),
  rbindlist(lapply(2:6, function(q) {
    data.table(exp = "Summer Ozone",
               lag = lag_vals,
               quantile = names(oz_levs)[q],
               val = oz_levs[q],
               RR = oz_dlnm$matRRfit[paste0(oz_levs[q]),],
               RRlow = oz_dlnm$matRRlow[paste0(oz_levs[q]),],
               RRhigh = oz_dlnm$matRRhigh[paste0(oz_levs[q]),],
               cumRR = oz_dlnm$allRRfit[paste0(oz_levs[q])],
               cumRRlow = oz_dlnm$allRRlow[paste0(oz_levs[q])],
               cumRRhigh = oz_dlnm$allRRhigh[paste0(oz_levs[q])]) 
  }))))


##### DLNM Plot #####
exp_dt$exp <- factor(exp_dt$exp,
                     levels = c("PM2.5", "NO2", "Summer Ozone"))
ggplot(exp_dt) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  # geom_ribbon(aes(x = lag, ymin = RRlow, ymax = RRhigh, fill = quantile), alpha = 0.3) +
  geom_line(aes(x = lag, y = RR, color = quantile, linetype = quantile), size = 2) +
  facet_grid( ~ exp) +
  theme_bw(base_size = 32) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm")) +
  scale_x_continuous(breaks = (1:5)*2, minor_breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.19)) +
  labs(x = "Years prior", y = "Hazard odds", 
       color = "", linetype = "", fill = "") +
  theme(legend.position = c(0.93, 0.82), 
        legend.background = element_blank(),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2.5, "cm"))

##### Cumulative effect plot #####
ggplot(exp_dt[lag == 1]) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(aes(x = quantile, ymin = cumRRlow, ymax = cumRRhigh), width = 0.5, size = 2) +
  geom_point(aes(x = quantile, y = cumRR), size = 3) +
  facet_grid(~exp) + 
  theme_bw(base_size = 32) +
  labs(x = "Exposure percentile", y = "Cumulative hazard odds")
