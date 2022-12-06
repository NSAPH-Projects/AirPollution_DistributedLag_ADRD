library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dlmtree)

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
who_epa_diff <- c(pm25 = 7, no2 = 43, ozone = 10)

plot_dat <- list()
for (exp in c("pm25", "no2", "ozone")) {
  rm(res)
  load(paste0(dir_data, "analysis/gamns3log_dlm_comprisk_contUS_", AD_ADRD, "_", 
              code_type, "_", exp, ".rda"))
  fit <- res$cp$matfit[2,] * who_epa_diff[exp]#exp_iqr
  se <- res$cp$matse[2,] * 1.96 * who_epa_diff[exp]#exp_iqr
  fit_cum <- res$cp$allfit[2] * who_epa_diff[exp]#exp_iqr
  se_cum <- res$cp$allse[2] * 1.96 * who_epa_diff[exp]#exp_iqr
  plot_dat[[exp]] <- data.table(lag = seq(1, 10, by = 0.2),
                                est = exp(fit),
                                low = exp(fit - se),
                                high = exp(fit + se),
                                cum = exp(fit_cum),
                                cum_low = exp(fit_cum - se_cum),
                                cum_high = exp(fit_cum + se_cum),
                                exp = exp)
}
plot_dat <- rbindlist(plot_dat)
plot_dat$model <- "Main analysis"

plot_dat_eco <- list()
for (exp in c("pm25", "no2", "ozone")) {
  rm(res)
  load(paste0(dir_data, "analysis/gamns3log_dlm_comprisk_contUS_", AD_ADRD, "_", 
              code_type, "_", exp, "_v2.rda"))
  fit <- res$cp$matfit[2,] * who_epa_diff[exp]#exp_iqr
  se <- res$cp$matse[2,] * 1.96 * who_epa_diff[exp]#exp_iqr
  fit_cum <- res$cp$allfit[2] * who_epa_diff[exp]#exp_iqr
  se_cum <- res$cp$allse[2] * 1.96 * who_epa_diff[exp]#exp_iqr
  plot_dat_eco[[exp]] <- data.table(lag = seq(1, 10, by = 0.2),
                                est = exp(fit),
                                low = exp(fit - se),
                                high = exp(fit + se),
                                cum = exp(fit_cum),
                                cum_low = exp(fit_cum - se_cum),
                                cum_high = exp(fit_cum + se_cum),
                                exp = exp)
}
plot_dat_eco <- rbindlist(plot_dat_eco)
plot_dat_eco$model <- "Ecological risk factors"


plot_dat2 <- list()
rm(res)
for (exp in c("pm25", "no2", "ozone")) {
  load(paste0(dir_data, "analysis/tdlm_contUS_samp25per_", 
              AD_ADRD, "_", code_type, "_", exp, ".rda"))
  s <- summary(m4)
  df <- 3
  fit <- s$matfit * who_epa_diff[exp]#exp_iqr
  fit_cum <- sum(s$matfit) * who_epa_diff[exp]#exp_iqr
  plot_dat2[[exp]] <- data.table(lag = seq(1, 10, by = 1),
                                 est = exp(fit),
                                 low = exp(s$cilower * who_epa_diff[exp]),
                                 high = exp(s$ciupper * who_epa_diff[exp]),
                                 cum = exp(fit_cum),
                                 cum_low = exp(s$cumulative.effect[1,3] / s$cumulative.effect[1,1] * who_epa_diff[exp]),
                                 cum_high = exp(s$cumulative.effect[1,4] / s$cumulative.effect[1,1] * who_epa_diff[exp]),
                                 exp = exp)
}
plot_dat2 <- rbindlist(plot_dat2)
plot_dat2$model <- "TDLM"

plot_dat3 <- list()
for (exp in c("pm25", "no2", "ozone")) {
  rm(res)
  load(paste0(dir_data, "analysis/gamns3log_dlm5yr_comprisk_contUS_", 
              AD_ADRD, "_", code_type, "_", exp, ".rda"))
  fit <- res$cp$matfit[2,] * who_epa_diff[exp]#exp_iqr
  se <- res$cp$matse[2,] * 1.96 * who_epa_diff[exp]#exp_iqr
  fit_cum <- res$cp$allfit[2] * who_epa_diff[exp]#exp_iqr
  se_cum <- res$cp$allse[2] * 1.96 * who_epa_diff[exp]#exp_iqr
  plot_dat3[[exp]] <- data.table(lag = seq(1, 5, by = 0.2),
                                est = exp(fit),
                                low = exp(fit - se),
                                high = exp(fit + se),
                                cum = exp(fit_cum),
                                cum_low = exp(fit_cum - se_cum),
                                cum_high = exp(fit_cum + se_cum),
                                exp = exp)
}
plot_dat3 <- rbindlist(plot_dat3)
plot_dat3$model <- "5-year"

pd <- rbindlist(list(plot_dat, plot_dat2, plot_dat3, plot_dat_eco))

#### ggplot ####
pd$exp <- factor(pd$exp, levels = c("pm25", "no2", "ozone"),
                 labels = c("PM[2.5]", "NO[2]", "Ozone"))
pd$model <- factor(pd$model, levels = c("Main analysis", "TDLM", "5-year",
                                        "Ecological risk factors"))
ggplot(pd, aes(x = lag, y = est, ymin = low, ymax = high, 
               color = model, linetype = model, fill = model)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_ribbon(alpha = 0.5, size = 0) +
  geom_line(size = 2) +
  facet_grid(~exp, labeller = label_parsed) +
  theme_bw(base_size = 32) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0.05, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.19)) +
  labs(x = "Years prior", y = "Odds ratio (NAAQS vs AQG)", 
       color = "", fill = "", linetype = "") +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") #+
  # guides(color = guide_legend(ncol = 2, byrow = TRUE),
  #        linetype = guide_legend(ncol = 2, byrow = TRUE),
  #        fill = guide_legend(ncol = 2, byrow = TRUE))

ggplot(pd[lag == 1], 
       aes(x = exp, y = cum, ymin = cum_low, ymax = cum_high,
           color = model)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar( width = 0.5, size = 2, position = position_dodge(width = 1)) +
  geom_point(size = 3, position = position_dodge(width = 1)) +
  theme_bw(base_size = 32) +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels = parse(text = levels(pd$exp))) +
  labs(x = "", y = "Cumulative odds ratio", color = "") +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(ncol = 2, byrow = TRUE))


plot_dat[lag == 1]
