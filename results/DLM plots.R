library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
who_epa_diff <- c(pm25 = 7, no2 = 43, ozone = 10)

plot_dat <- list()
for (exp in c("pm25", "no2", "ozone")) {
  load(paste0(dir_data, "analysis/gamns3log_dlm_comprisk_contUS_", AD_ADRD, "_", 
              code_type, "_", exp, ".rda"))
  if (!("cp" %in% names(res))) {
    res <- res[[3]]
  }
  df <- 3
  fit <- res$cp$matfit[2,] * who_epa_diff[exp]#exp_iqr
  se <- res$cp$matse[2,] * 1.96 * who_epa_diff[exp]#exp_iqr
  fit_cum <- res$cp$allfit[2] * who_epa_diff[exp]#exp_iqr
  se_cum <- res$cp$allse[2] * 1.96 * who_epa_diff[exp]#exp_iqr
  plot_dat[[exp]] <- data.table(lag = seq(1, 10, by = 0.2),
                                est = exp(fit),
                                low = exp(fit - se),
                                high = exp(fit + se),
                                cum = exp(fit_cum),
                                cum_low = exp(fit_cum - se_cum ),
                                cum_high = exp(fit_cum + se_cum),
                                df = df,
                                exp = exp)
}
plot_dat <- rbindlist(plot_dat)

plot_dat[lag %in% 1:10, .(lag, est = paste0(round(est,4), 
                                            " (", round(low,4), ", ", 
                                            round(high,4), ")"), exp)] %>%
  pivot_wider(c("lag", "exp"), names_from = "exp", values_from = "est")

#### ggplot: df = 3 ####
plot_dat$exp <- factor(plot_dat$exp, levels = c("pm25", "no2", "ozone"),
                       labels = c("PM[2.5]", "NO[2]", "Summer~Ozone"))
ggplot(plot_dat, aes(x = lag, y = est, 
                     ymin = low, ymax = high, group = df)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_ribbon(fill = "grey", color = "grey", size = 1) +
  geom_line(size = 2) +
  facet_grid(~exp, labeller = label_parsed) +
  theme_bw(base_size = 32) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0.05, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.19)) +
  labs(x = "Years prior", y = "Odds ratio (NAAQS vs AQG)")

ggplot(plot_dat[lag == 1], 
       aes(x = exp, y = cum, ymin = cum_low, ymax = cum_high)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar( width = 0.5, size = 2) +
  geom_point(size = 3) +
  scale_x_discrete(labels = parse(text = levels(plot_dat$exp))) +
  theme_bw(base_size = 32) +
  labs(x = "", y = "Cumulative odds ratio\n(NAAQS vs AQG)")


plot_dat[lag == 1]