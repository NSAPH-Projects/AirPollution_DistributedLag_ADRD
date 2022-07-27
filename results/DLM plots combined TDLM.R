library(data.table)
library(dlmtree)
library(ggplot2)

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
who_epa_diff <- c(pm25 = 7, no2 = 43, ozone = 10)

plot_dat2 <- list()
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
                                df = df,
                                exp = exp)
}
plot_dat2 <- rbindlist(plot_dat2)

#### ggplot: df = 3 ####
plot_dat2$exp <- factor(plot_dat2$exp, levels = c("pm25", "no2", "ozone"),
                       labels = c("PM2.5", "NO2", "Summer Ozone"))
ggplot(plot_dat2, aes(x = lag, y = est, 
                     ymin = low, ymax = high, group = df)) +
  geom_hline(yintercept = 1, color = "red", size = 1) +
  geom_ribbon(fill = "grey", color = "grey", size = 0) +
  geom_line(size = 2) +
  facet_grid(~exp) +
  theme_bw(base_size = 32) +
  scale_x_continuous(breaks = 1:9, minor_breaks = NULL, expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.19)) +
  labs(x = "Years prior", y = "Hazard odds")

ggplot(plot_dat2[lag == 1], 
       aes(x = exp, y = cum, ymin = cum_low, ymax = cum_high)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar( width = 0.5, size = 2) +
  geom_point(size = 3) +
  theme_bw(base_size = 32) +
  labs(x = "Years prior", y = "Cumulative hazard odds")


plot_dat2[lag == 1]


plot_dat2[, model := "TDLM"]
plot_dat[, model := "GAM"]
plot_dat_comb <- rbindlist(list(plot_dat, plot_dat2))


ggplot(plot_dat_comb[lag == 1], 
       aes(x = exp, y = cum, ymin = cum_low, ymax = cum_high,
           color = model)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar( width = 0.5, size = 2, position = position_dodge(width = 0.8)) +
  geom_point(size = 3, position = position_dodge(width = 0.8)) +
  theme_bw(base_size = 32) +
  labs(x = "Years prior", y = "Cumulative hazard odds", color = "")
