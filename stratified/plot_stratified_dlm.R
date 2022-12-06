library(ggplot2)
library(data.table)
dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"
dir_results <- paste0(dir_data, "analysis/DLM_Strat/stratified_analysis/")
exp_change <- c(pm25 = 4.11, no2 = 14.99, ozone = 10)
files <- list.files(dir_results)

plot_dat <- list()
for (f in files) {
  load(paste0(dir_results, f))
  file_start <- strsplit(f, ".", fixed = TRUE)[[1]][1]
  desc <- lapply(strsplit(file_start, "_", TRUE)[[1]],
                 function(i) strsplit(i, "-", TRUE)[[1]])
  names(desc) <- sapply(desc, function(i) i[1])
  desc <- lapply(desc, function(i) i[2])
  # desc[['n']] <- NULL
  plot_dat[[f]] <- cbind.data.frame(
    desc,
    data.frame(change = exp_change[desc$exp][[1]],
               lag = 1:10,
               est = exp(s$matfit * exp_change[desc$exp]),
               low = exp(s$cilower * exp_change[desc$exp]),
               high = exp(s$ciupper * exp_change[desc$exp])),
    exp(s$cumulative.effect[1,] * exp_change[desc$exp] / 
          s$cumulative.effect$vals[1])[2:4]
  )
}

plot_dat <- rbindlist(plot_dat)
plot_dat$sexM <- factor(plot_dat$sexM, 0:1, labels = c("Female", "Male"))
plot_dat$race <- factor(plot_dat$race,
                        levels = c("nat", "api", "oth", "his", "blk", "wht"),
                        labels = c("Nat. Amer.", "Asian/PI", "Other", "Hisp.",
                                   "Black", "White"))
plot_dat$dual <- factor(plot_dat$dual,
                        levels = 0:1,
                        labels = c("Medicare FFS", "Medicaid eligible"))
plot_dat$age <- factor(plot_dat$age,
                       levels = c("young", "old"),
                       labels = c("75-80", "81+"))
plot_dat$eff <- plot_dat$low > 1 | plot_dat$high < 1





# NO2 -----
ggplot(plot_dat[exp == "no2"], 
       aes(x = lag, y = est, 
           ymin = low, ymax = high, 
           linetype = sexM, color = sexM, fill = sexM, group = sexM)) +
  # credible intervals
  geom_ribbon(alpha = 0.25) +
  # periods of susceptibility
  geom_point(aes(x = lag, y = ifelse(eff == 1, est, NA)), size = 1) +
  # reference line
  geom_hline(yintercept = 1, color = "black", size = 0.5, linetype = 1) +
  # theme and display
  theme_minimal(base_size = 16) +
  facet_grid(race ~ dual + age) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10), minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "NO2",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "no2"],
       aes(x = race, y = mean, ymin = lower, ymax = upper,
           color = sexM)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~dual + age) +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "NO2",
       color = "")





# PM2.5 -----
ggplot(plot_dat[exp == "pm25"], 
       aes(x = lag, y = est, ymin = low, ymax = high, linetype = sexM,
           color = sexM, fill = sexM, group = sexM)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(alpha = 0.1) +
  # geom_line(size = 1) +
  geom_point(aes(x = lag, y = ifelse(eff == 1, est, NA)), size = 1) +
  theme_minimal(base_size = 16) +
  facet_grid(race ~ dual + age, scales = "free_y") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10), minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "PM2.5",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "pm25"],
       aes(x = race, y = mean, ymin = lower, ymax = upper,
           color = sexM)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~dual + age) +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "PM2.5",
       color = "")

# Ozone
ggplot(plot_dat[exp == "ozone"], 
       aes(x = lag, y = est, ymin = low, ymax = high, linetype = sexM,
           color = sexM, fill = sexM, group = sexM)) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(alpha = 0.1) +
  # geom_line(size = 1) +
  geom_point(aes(x = lag, y = ifelse(eff == 1, est, NA)), size = 1) +
  theme_minimal(base_size = 16) +
  facet_grid(race ~ dual + age, scales = "free_y") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10), minor_breaks = NULL) +
  scale_color_brewer(type = "div", palette = "Set1", 
                     aesthetics = c("color", "fill")) +
  labs(x = "Years Prior", y = "Odds ratio (IQR increase)", title = "Ozone",
       color = "", fill = "", linetype = "")

ggplot(plot_dat[exp == "ozone"],
       aes(x = race, y = mean, ymin = lower, ymax = upper,
           color = sexM)) +
  geom_hline(yintercept = 1, color = "black", size = 1) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~dual + age) +
  scale_color_brewer(type = "div", palette = "Set1") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Cumulative odds ratio (IQR increase)", title = "Ozone",
       color = "")
