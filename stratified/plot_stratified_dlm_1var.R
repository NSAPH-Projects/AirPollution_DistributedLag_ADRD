library(ggplot2)
library(data.table)
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
dir_results <- paste0(dir_data, "analysis/DLM_Strat/stratified_analysis_1var/")
exp_change <- c(pm25 = 4.11, no2 = 14.99, ozone = 10)
files <- list.files(dir_results)

plot_dat <- list()
for (f in files) {
  load(paste0(dir_results, f))
  plot_dat[[f]] <- cbind.data.frame(
    n, var, level, exp,
    data.frame(change = exp_change[exp][[1]],
               lag = 1:10,
               est = exp(s$matfit * exp_change[exp]),
               low = exp(s$cilower * exp_change[exp]),
               high = exp(s$ciupper * exp_change[exp])),
    exp(s$cumulative.effect[1,] * exp_change[exp] / 
          s$cumulative.effect$vals[1])[2:4]
  )
}

plot_dat <- rbindlist(plot_dat)
plot_dat$eff <- ifelse(plot_dat$low > 1 | plot_dat$high < 1, TRUE, FALSE)
# NO2
ggplot(plot_dat[exp == "no2"], 
       aes(x = lag, y = est, ymin = low, ymax = high)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  geom_point(aes(x = lag, y = ifelse(eff, est, NA)), size = 1) +
  theme_minimal(base_size = 16) +
  facet_grid(~var + level) +
  scale_x_continuous(expand = c(0.02, 0), breaks = 1:10, minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "NO2",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "no2" & lag == 1],
       aes(x = level, y = mean, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~var, scales = "free_x") +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "NO2")

# PM2.5
ggplot(plot_dat[exp == "pm25"], 
       aes(x = lag, y = est, ymin = low, ymax = high)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  geom_point(aes(x = lag, y = ifelse(eff, est, NA)), size = 1) +
  theme_minimal(base_size = 16) +
  facet_grid(~var + level) +
  scale_x_continuous(expand = c(0.02, 0), breaks = 1:10, minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "PM2.5",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "pm25" & lag == 1],
       aes(x = level, y = mean, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~var, scales = "free_x") +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "PM2.5")

# Ozone
ggplot(plot_dat[exp == "ozone"], 
       aes(x = lag, y = est, ymin = low, ymax = high)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  geom_point(aes(x = lag, y = ifelse(eff, est, NA)), size = 1) +
  theme_minimal(base_size = 16) +
  facet_grid(~var + level) +
  scale_x_continuous(expand = c(0.02, 0), breaks = 1:10, minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "Ozone",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "ozone" & lag == 1],
       aes(x = level, y = mean, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~var, scales = "free_x") +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "Ozone")
