library(data.table)
library(ggplot2)

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"

plot_dat <- list()
for (exp in c("pm25", "no2", "ozone")) {
  load(paste0(dir_data, "analysis/gamns3log_dlm_comprisk_contUS_", AD_ADRD, "_", 
              code_type, "_", exp, ".rda"))
  df <- 3
  fit <- res[[df]]$cp$matfit[2,] * exp_iqr
  se <- res[[df]]$cp$matse[2,] * 1.96 * exp_iqr
  fit_cum <- res[[df]]$cp$allfit[2] * exp_iqr
  se_cum <- res[[df]]$cp$allse[2] * exp_iqr * 1.96
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

#### ggplot: df = 3 ####
plot_dat$exp <- factor(plot_dat$exp, levels = c("pm25", "no2", "ozone"),
                       labels = c("PM2.5", "NO2", "Summer Ozone"))
ggplot(plot_dat, aes(x = lag, y = est, 
                              ymin = low, ymax = high, group = df)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_ribbon(fill = "grey", color = "grey", size = 2) +
  geom_line(size = 2) +
  facet_grid(~exp) +
  theme_bw(base_size = 32) +
  scale_x_continuous(breaks = 1:5 * 2, minor_breaks = NULL, expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.96, 1.19)) +
  labs(x = "Years prior", y = "Hazard odds")

ggplot(plot_dat[lag == 1], 
       aes(x = exp, y = cum, ymin = cum_low, ymax = cum_high)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar( width = 0.5, size = 2) +
  geom_point(size = 3) +
  theme_bw(base_size = 32) +
  labs(x = "Years prior", y = "Cumulative hazard odds")


plot_dat[lag == 1]

# 
# load(paste0(dir_data, "analysis/gamns_state_dlmbic_boots_contUS_",
#             AD_ADRD, "_", code_type, "_", exp, ".rda"))
# est_bic <- exp(res[[df_bic]]$cp$matfit[q3,] - res[[df_bic]]$cp$matfit[q1,])
# se_bic <- apply(boot_mat, 2, sd) * exp_iqr * sqrt(boot_size / n)
#### BIC m/n boot plot data ####
plot_dat_boot_bic <-
  data.table(lag = seq(1, 10, by = 0.2),
             est = est_bic,
             low = est_bic - 1.96 * se_bic,
             high = est_bic + 1.96 * se_bic)
