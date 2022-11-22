rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 2)
dir_pm25 <- "/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/"
pm25_data <- fread(paste0(dir_pm25, "all_years.csv"))[, zip := ZIP][, ZIP := NULL][zip >= 501]
pm25_data <- dcast(pm25_data, zip ~ year, value.var = 'pm25')
pm25_data <- as.matrix(pm25_data[complete.cases(pm25_data)][,-1])

exp_dat <- rbind.data.frame(
  data.frame(year = 2000:2006, exp = pm25_data[100, 1:7], grp = 1),
  data.frame(year = 2000:2010, exp = pm25_data[500, 1:11], grp = 2),
  data.frame(year = 2000:2012, exp = pm25_data[1000, 1:13], grp = 3),
  data.frame(year = 2000:2013, exp = pm25_data[1200, 1:14], grp = 4),
  data.frame(year = 2000:2016, exp = pm25_data[2000, 1:17], grp = 5)
)

ggplot(exp_dat, aes(x = year, y = exp, group = grp, color = factor(grp))) +
  geom_line(size = 2) +
  geom_vline(xintercept = 2009.5, linetype = 2, size = 2) +
  geom_vline(xintercept = 2016.5, linetype = 2, size = 2) +
  theme_bw(base_size = 24) +
  labs(x = "Year", y = "Exposure") +
  scale_x_continuous(breaks = 2000 + 4 * 0:4) +
  theme(legend.position = 'none')
