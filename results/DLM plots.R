library(data.table)
library(ggplot2)

AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"
boot_size <- 533165
n <- 8507437

exp <- "no2"
load(paste0(dir_data, "analysis/gamns_state_dlmaic_contUS_", AD_ADRD, "_", 
            code_type, "_", exp, ".rda"))
load(paste0(dir_data, "analysis/gamcrpen-state_contUS_", AD_ADRD, "_",
                   code_type, "_", exp, ".rda"))
q1 <- which(aic_res[[2]]$cp$predvar == k[3])
q3 <- which(aic_res[[2]]$cp$predvar == k[5])
exp_iqr <- k[5] - k[3]
df_bic <- which.min(sapply(1:10, function(df) aic_res[[df]]$bic))
df_aic <- which.min(sapply(1:10, function(df) aic_res[[df]]$aic))


load(paste0(dir_data, "analysis/gamns_state_dlmbic_boots_contUS_",
            AD_ADRD, "_", code_type, "_", exp, ".rda"))
est_bic <- exp(aic_res[[df_bic]]$cp$matfit[q3,] - aic_res[[df_bic]]$cp$matfit[q1,])
se_bic <- apply(boot_mat, 2, sd) * exp_iqr * sqrt(boot_size / n)

#### df plot data ####
plot_dat <- rbindlist(lapply(2:10, function(df) {
  data.table(lag = seq(1, 10, by = 0.2),
             est = exp(aic_res[[df]]$cp$matfit[q3,] - aic_res[[df]]$cp$matfit[q1,]),
             low = exp(log(aic_res[[df]]$cp$matRRlow[q3,]) - log(aic_res[[df]]$cp$matRRlow[q1,])),
             high = exp(log(aic_res[[df]]$cp$matRRhigh[q3,]) - log(aic_res[[df]]$cp$matRRhigh[q1,])),
             df = df)
}))
#### AIC/BIC/Pen plot data ####
plot_dat_aic_bic_pen <- rbind.data.frame(
  data.frame(lag = seq(1, 10, by = 0.2),
             est = exp(aic_res[[df_aic]]$cp$matfit[q3,] - aic_res[[df_aic]]$cp$matfit[q1,]),
             low = exp(log(aic_res[[df_aic]]$cp$matRRlow[q3,]) - log(aic_res[[df_aic]]$cp$matRRlow[q1,])),
             high = exp(log(aic_res[[df_aic]]$cp$matRRhigh[q3,]) - log(aic_res[[df_aic]]$cp$matRRhigh[q1,])),
             group = "AIC"),
  data.frame(lag = seq(1, 10, by = 0.2),
             est = exp(aic_res[[df_bic]]$cp$matfit[q3,] - aic_res[[df_bic]]$cp$matfit[q1,]),
             low = exp(log(aic_res[[df_bic]]$cp$matRRlow[q3,]) - log(aic_res[[df_bic]]$cp$matRRlow[q1,])),
             high = exp(log(aic_res[[df_bic]]$cp$matRRhigh[q3,]) - log(aic_res[[df_bic]]$cp$matRRhigh[q1,])),
             group = "BIC"),
  data.frame(lag = seq(1, 10, by = 0.2),
             est = exp(cp_dlm$matfit[q3,] - cp_dlm$matfit[q1,]),
             low = exp(log(cp_dlm$matRRlow[q3,]) - log(cp_dlm$matRRlow[q1,])),
             high = exp(log(cp_dlm$matRRhigh[q3,]) - log(cp_dlm$matRRhigh[q1,])),
             group = "Pen")
)
#### BIC m/n boot plot data ####
plot_dat_boot_bic <-
  data.table(lag = seq(1, 10, by = 0.2),
             est = est_bic,
             low = est_bic - 1.96 * se_bic,
             high = est_bic + 1.96 * se_bic)

#### ggplot: BIC m/n boot CI ####
ggplot(plot_dat_boot_bic, aes(x = lag, y = est, ymin = low, ymax = high)) +
  geom_hline(yintercept = 1) +
  geom_ribbon(fill = "grey") +
  geom_line() +
  theme_bw(base_size = 24) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0)) +
  labs(x = "Years prior", y = "Hazard odds, IQR increase",
       title = paste0("Distributed lag: ", exp, ", IQR = ", round(exp_iqr, 2)), 
       subtitle = paste0("BIC selected df = ", df_bic, ", m/n bootstrap CI"))

#### ggplot: AIC/BIC/Pen ####
ggplot(plot_dat_aic_bic_pen,
       aes(x = lag, y = est, ymin = low, ymax = high, group = group, color = group, fill = group)) +
  geom_hline(yintercept = 1) +
  geom_ribbon(alpha = 0.2) + 
  geom_line(size = 2) +
  theme_bw(base_size = 24) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0)) +
  labs(x = "Years prior", y = "Hazard odds, IQR increase",
       title = paste0("Distributed lag: ", exp, ", IQR = ", round(exp_iqr, 2)))

#### ggplot: df ####
ggplot(plot_dat, aes(x = lag, y = est, ymin = low, ymax = high, 
                     group = df, color = factor(df))) +
  geom_hline(yintercept = 1) +
  geom_line(size = 2, alpha = 0.5) +
  theme_bw(base_size = 24) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0)) +
  labs(x = "Years prior", y = "Hazard odds, IQR increase", color = "df",
       title = paste0("Distributed lag: ", exp, ", IQR = ", round(exp_iqr, 2)))

#### ggplot: df = 3 ####
ggplot(plot_dat[df == 3], aes(x = lag, y = est, 
                               ymin = low, ymax = high, group = df)) +
  geom_hline(yintercept = 1) +
  geom_ribbon(alpha = 0.5) +
  geom_line(size = 2) +
  theme_bw(base_size = 24) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0)) +
  labs(x = "Years prior", y = "Hazard odds, IQR increase", color = "df",
       title = paste0("Distributed lag: ", exp, ", IQR = ", round(exp_iqr, 2)))



which.min(sapply(1:10, function(df) aic_res[[df]]$aic))
plot(1:10, sapply(1:10, function(df) aic_res[[df]]$aic))
which.min(sapply(1:10, function(df) aic_res[[df]]$bic))
plot(1:10, sapply(1:10, function(df) aic_res[[df]]$bic))
