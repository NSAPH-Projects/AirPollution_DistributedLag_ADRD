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
library(mgcv)
library(dlnm)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = FALSE)
setDTthreads(threads = 1)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

##### Load DLNM results #####
perc <- c("0.5%", "10%", "25%", "50%", "75%", "90%", "99.5%")
load(paste0(dir_data, "analysis/gamns_dlnmlog35_expdist15_comprisk_who-epa_contUSsamp_ADRD_any_pm25.rda"))
pm_dlnm <- cp_dlnm_005
pm_levs <- pred_k[c(perc, "epa_pm25", "pm25")]

load(paste0(dir_data, "analysis/gamns_dlnmlog35_expdist15_comprisk_who-epa_contUSsamp_ADRD_any_no2.rda"))
no_dlnm <- cp_dlnm_005
no_levs <- pred_k[c(perc, "epa_no2", "no2")]

load(paste0(dir_data, "analysis/gamns_dlnmlog35_expdist15_comprisk_who-epa_contUSsamp_ADRD_any_ozone.rda"))
oz_dlnm <- cp_dlnm_005
oz_levs <- pred_k[c(perc, "epa_ozone", "ozone")]

names(pm_levs) <- names(no_levs) <- names(oz_levs) <- c(perc, "NAAQS", "AQG")


##### Create data table of results #####
lag_vals <- seq(1, 10, by = 0.2)
exp_dt <- rbindlist(list(
  rbindlist(lapply(1:length(pm_levs), function(q) {
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
  rbindlist(lapply(1:length(pm_levs), function(q) {
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
  rbindlist(lapply(1:length(pm_levs), function(q) {
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
ggplot(exp_dt[quantile != "NAAQS" & quantile != "AQG" & lag %in% c(1) & quantile != "99.5%"]) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_ribbon(aes(x = val, ymin = RRlow, ymax = RRhigh), alpha = 0.3) +
  geom_line(aes(x = val, y = RR, group = lag), size = 2) +
  facet_grid( ~ exp, scales = "free_x") +
  theme_bw(base_size = 32) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm")) +
  # scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0.05, 0)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.13)) +
  labs(x = "Exposure concentration", y = "Odds ratio", 
       color = "Lag", linetype = "Lag", fill = "") +
  theme(legend.position = "none") +
  labs(title = "Exposure-response: 1 year prior")
  # theme(legend.position = c(0.93, 0.82), 
  #       legend.background = element_blank(),
  #       legend.key.height = unit(1, "cm"),
  #       legend.key.width = unit(2.5, "cm"))

##### Cumulative effect plot #####
ggplot(exp_dt[lag == 1 & quantile != "0.5%" & quantile != "99.5%"], 
       aes(reorder_within(quantile, val, exp), cumRR, 
           ymin = cumRRlow, ymax = cumRRhigh)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(width = 0.5, size = 2) +
  geom_point(size = 3) +
  facet_grid(~exp, scales = "free_x") + 
  scale_x_reordered() +
  theme_bw(base_size = 32) +
  labs(x = "Exposure concentration", y = "Cumulative odds ratio\nrelative to 0.5%") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

exp_dt[lag == 1 & quantile != "0.5%" & quantile != "99.5%", .(exp, quantile, cumRR, cumRRlow, cumRRhigh)]
