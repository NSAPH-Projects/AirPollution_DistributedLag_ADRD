rm(list = ls())
gc()

##### Setup #####
AD_ADRD <- "ADRD" # AD or ADRD
code_type <- "any" # primary or any
n_threads = 2
Sys.setenv(OMP_NUM_THREADS = n_threads)


library(data.table)
library(fst)
library(ggplot2)

options(stringsAsFactors = FALSE)
setDTthreads(threads = n_threads)
threads_fst(n_threads)

dir_data <- "/n/home_fasse/dmork/projects/adrd_dlm/data/"

##### Load relevant data #####
qid_dat <- read_fst(paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                           code_type, "_qid.fst"), as.data.table = TRUE)
pm25_dat <- read_fst(paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                            code_type, "_pm25.fst"), as.data.table = TRUE)

##### Merge PM/qid data #####
qid_dat <- merge(qid_dat, pm25_dat, by = c("qid", "year"), all.x = TRUE)
rm(pm25_dat)
setnames(qid_dat, "zip.x", "zip")


##### Create aggregated datasets #####
qid_dat[, PIR := medianhousevalue / (medhouseholdincome + 0.001)]
qid_dat[, age_corrected := entry_age + year - entry]
qid_dat <- qid_dat[!is.infinite(qid_dat$PIR) &
                     complete.cases(qid_dat[,.(year, age_corrected, dual, race, sexM,
                                               region, statecode,
                                               education, poverty, pct_blk, hispanic,
                                               popdensity, pct_owner_occ, PIR,
                                               winter_tmmx, winter_rmax,
                                               summer_tmmx, summer_rmax,
                                               lag1, lag2, lag3, lag4, lag5,
                                               lag6, lag7, lag8, lag9, lag10)])]
qid_dat <- qid_dat[lag1 > 0 & lag2 > 0 & lag3 > 0 & lag4 > 0 & lag5 > 0 &
                     lag6 > 0 & lag7 > 0 & lag8 > 0 & lag9 > 0 & lag10 > 0]
qid_dat[, .N]
qid_dat[, dual := dual == "1"]
qid_dat[, dual_any := ifelse(any(dual), "Yes", "No"), by = qid]
qid_dat[, sex := ifelse(sexM, "M", "F"), by = qid]
qid_dat[, age_cat := cut(age_corrected, breaks = c(74, 84, 114), labels = c("75-84", "85+"))]
qid_dat_agg_zip <- qid_dat[,.(at_risk = .N, n_hosp = sum(first_hosp), 
                              rate = sum(first_hosp / (year - entry)) / .N,
                              avg_age = mean(age_corrected),
                              his = mean(hispanic), blk = mean(pct_blk),
                              PIR = mean(PIR), pov = mean(poverty), ed = mean(education),
                              oo = mean(pct_owner_occ), dens = mean(popdensity),
                              s_tmmx = mean(summer_tmmx), w_tmmx = mean(winter_tmmx),
                              s_rmax = mean(summer_rmax), w_rmax = mean(winter_rmax),
                              lat = mean(latitude), long = mean(longitude), 
                              state = first(statecode), region = first(region),
                              pm1 = mean(lag1), pm2 = mean(lag2), pm3 = mean(lag3),
                              pm4 = mean(lag4), pm5 = mean(lag5), pm6 = mean(lag4),
                              pm7 = mean(lag7), pm8 = mean(lag8), pm9 = mean(lag9),
                              pm10 = mean(lag10)),
                           by = .(year, entry, zip, sex, age_cat, dual_any, race)]
qid_dat_agg_zip[rate == 0, rate := min(qid_dat_agg_zip$rate[qid_dat_agg_zip$rate > 0])]
qid_dat_agg_zip[, lograte := log(rate)]
write_fst(qid_dat_agg_zip, paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                  code_type, "_agg_zip.fst"))
qid_dat_agg_fips <- qid_dat[,.(at_risk = .N, n_hosp = sum(first_hosp), 
                               rate = sum(first_hosp / (.N * (year - entry))),
                               avg_age = mean(age_corrected),
                               his = mean(hispanic), blk = mean(pct_blk),
                               PIR = mean(PIR), pov = mean(poverty), ed = mean(education),
                               oo = mean(pct_owner_occ), dens = mean(popdensity),
                               s_tmmx = mean(summer_tmmx), w_tmmx = mean(winter_tmmx),
                               s_rmax = mean(summer_rmax), w_rmax = mean(winter_rmax),
                               lat = mean(latitude), long = mean(longitude), 
                               state = first(statecode), region = first(region),
                               pm1 = mean(lag1), pm2 = mean(lag2), pm3 = mean(lag3),
                               pm4 = mean(lag4), pm5 = mean(lag5), pm6 = mean(lag4),
                               pm7 = mean(lag7), pm8 = mean(lag8), pm9 = mean(lag9),
                               pm10 = mean(lag10)),
                            by = .(year, entry, fips, sex, age_cat, dual_any, race)]
qid_dat_agg_fips[rate == 0, rate := min(qid_dat_agg_fips$rate[qid_dat_agg_fips$rate > 0])]
qid_dat_agg_fips[, lograte := log(rate)]
write_fst(qid_dat_agg_fips, paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                   code_type, "_agg_fips.fst"))
qid_dat_agg_state <- qid_dat[,.(at_risk = .N, n_hosp = sum(first_hosp), 
                                rate = sum(first_hosp / (.N * (year - entry))),
                                avg_age = mean(age_corrected),
                                his = mean(hispanic), blk = mean(pct_blk),
                                PIR = mean(PIR), pov = mean(poverty), ed = mean(education),
                                oo = mean(pct_owner_occ), dens = mean(popdensity),
                                s_tmmx = mean(summer_tmmx), w_tmmx = mean(winter_tmmx),
                                s_rmax = mean(summer_rmax), w_rmax = mean(winter_rmax),
                                lat = mean(latitude), long = mean(longitude), 
                                state = first(statecode), region = first(region),
                                pm1 = mean(lag1), pm2 = mean(lag2), pm3 = mean(lag3),
                                pm4 = mean(lag4), pm5 = mean(lag5), pm6 = mean(lag4),
                                pm7 = mean(lag7), pm8 = mean(lag8), pm9 = mean(lag9),
                                pm10 = mean(lag10)),
                             by = .(year, entry, statecode, sex, age_cat, dual_any, race)]
qid_dat_agg_state[rate == 0, rate := min(qid_dat_agg_state$rate[qid_dat_agg_state$rate > 0])]
qid_dat_agg_state[, lograte := log(rate)]
write_fst(qid_dat_agg_state, paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                    code_type, "_agg_state.fst"))


##### Analysis #####
qid_dat_agg_zip <- read_fst(paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                     code_type, "_agg_zip.fst"), as.data.table = TRUE)
qid_dat_agg_fips <- read_fst(paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                     code_type, "_agg_fips.fst"), as.data.table = TRUE)
qid_dat_agg_state <- read_fst(paste0(dir_data, "analysis/HDLM/contUS_agg_", AD_ADRD, "_", 
                                     code_type, "_agg_state.fst"), as.data.table = TRUE)

qid_dat_agg <- qid_dat_agg_state
exp_dat <- as.matrix(qid_dat_agg[, paste0("pm", 1:10)])
exp_iqr <- IQR(exp_dat)
exp_mean <- mean(exp_dat)
exp_dat <- (exp_dat - exp_mean) / exp_iqr

library(dlnm)
library(mgcv)
cb_dlm <- crossbasis(as.matrix(exp_dat), c(1, 10), 
                     arglag = list(fun = "ns", knots = exp(seq(0, log(10), length = 3))[2], 
                                   intercept = TRUE), 
                     argvar = list(fun = "lin"))
m3 <- gam(lograte ~ cb_dlm + factor(age_cat) + sex + dual_any + factor(state) +
            factor(race) + factor(year) + factor(entry) + 
            his + blk + PIR + pov + ed + oo + dens + avg_age,
          data = qid_dat_agg)
cp_dlm <- crosspred(cb_dlm, m3, cen = 0, at = 0:1, bylag = 0.2)
plotdat <- data.frame(lag = seq(1, 10, by = 0.2), est = 100*(exp(cp_dlm$matfit[2,])-1), 
                      low = 100*(exp(cp_dlm$matlow[2,])-1), high = 100*(exp(cp_dlm$mathigh[2,])-1))
ggplot(plotdat) + 
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(aes(x = lag, ymin = low, ymax = high), fill = "grey") +
  geom_line(aes(x = lag, y = est)) + theme_bw(base_size = 24) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Years Prior", y = "% Change in Hospitalization Rate")

library(dlmtree)
hdlm <- dlmtree(lograte ~ factor(age_cat) + sex + dual_any + factor(state) +
                  factor(race) + factor(year) + factor(entry) + 
                  his + blk + PIR + pov + ed + oo + dens + avg_age, 
                data = qid_dat_agg, 
                exposure.data = exp_dat,
                dlm.type = "shared", 
                tree.modifiers = c("sex", "dual_any", "race", "age_cat", "region"), 
                modifier.splits = 2,
                n.trees = 20, n.burn = 5000, n.iter = 10000, n.thin = 5, 
                save.data = FALSE)
save(hdlm, file = paste0(dir_data, "analysis/HDLM/contUS_", AD_ADRD, "_", 
                         code_type, "_hdlm_agg_state.rda"))

plot(rowMeans(hdlm$termNodesMod))
plot(rowMeans(hdlm$termNodesDLM))
colMeans(hdlm$modCount > 0)
boxplot(hdlm$modInf)

plotDLM <- function(estDLM_data, groups = 1, trans = "log") {
  if (trans == "log") {
    estDLM_data$plotData$est <- 100*(exp(estDLM_data$plotData$est) - 1)
    estDLM_data$plotData$upper <- 100*(exp(estDLM_data$plotData$upper) - 1)
    estDLM_data$plotData$lower <- 100*(exp(estDLM_data$plotData$lower) - 1)
  }
  if (groups == 2) {
    estDLM_data$plotData$grp1 <- sapply(strsplit(estDLM_data$plotDat$group, " & "), function(g) g[1])
    estDLM_data$plotData$grp2 <- sapply(strsplit(estDLM_data$plotDat$group, " & "), function(g) g[2])
  }
  p <- ggplot(estDLM_data$plotData) +
    geom_hline(yintercept = 0, color = "red") +
    geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = "grey") +
    geom_line(aes(x = time, y = est)) +
    facet_wrap(~group)
  if (groups == 2) {
    p <- p+ facet_grid(grp1 ~ grp2)
  } else {
    p <- p + facet_wrap(~group)
  }
  p <- p + theme_bw(base_size = 24) +
    scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0)) +
    labs(x = "Years Prior", y = "% change in Hospitalization Rate")
  p
}

createGrpIdx <- function(model, data, mod1) {
  m1 <- which(model$modNames == mod1)
  m1sp <- model$modSplitValRef[[m1]]
  m1num <- model$modIsNum[m1]
  if (model$modIsNum[m1]) {
    grp1_idx <- lapply(1:(length(m1sp) + 1), function(i) {
      if (i == 1) {
        which(data[[mod1]] < m1sp[i])
      } else if (i <= length(m1sp)) {
        which(data[[mod1]] >= m1sp[i-1] & data[[mod1]] < m1sp[i])
      } else {
        which(data[[mod1]] >= m1sp[i-1])
      }
    })
    names(grp1_idx) <- sapply(1:(length(m1sp) + 1), function(i) {
      if (i == 1) {
        paste0(mod1, " < ", round(m1sp[i], 2))
      } else if (i <= length(m1sp)) {
        paste0(mod1, " in [", round(m1sp[i-1], 2), ", ", round(m1sp[i], 2),")")
      } else {
        paste0(mod1, " >= ", round(m1sp[i-1], 2))
      }
    })
  } else {
    grp1_idx <- lapply(m1sp, function(s) {
      which(data[[mod1]] == s)
    })
    names(grp1_idx) <- sapply(m1sp, function(s) {
      paste0(mod1, " = ", s)
    })
  }
  return(grp1_idx)
}
create2GrpIdx <- function(model, data, mod1, mod2) {
  grp1_idx <- createGrpIdx(model, data, mod1)
  grp2_idx <- createGrpIdx(model, data, mod2)
  comb_idx <- list()
  for (i in 1:length(grp1_idx)) {
    for (j in 1:length(grp2_idx)) {
      comb_idx[[paste0(names(grp1_idx)[i], " & ", names(grp2_idx)[j])]] <-
        intersect(grp1_idx[[i]], grp2_idx[[j]])
    }
  }
  return(comb_idx)
}



allDLM <- estDLM(hdlm, qid_dat_agg, group.index = list("all" = 1:qid_dat_agg[,.N]))
plotDLM(allDLM) + theme_bw(base_size = 24)

ageDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "age_cat"))
plotDLM(ageDLM)

raceDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "race"))
plotDLM(raceDLM)

dualDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "dual_any"))
plotDLM(dualDLM)

sexDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "sex"))
plotDLM(sexDLM)

#densDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "dens"))
#plotDLM(densDLM)

regDLM <- estDLM(hdlm, qid_dat_agg, createGrpIdx(hdlm, qid_dat_agg, "region"))
plotDLM(regDLM)



raceageDLM <- estDLM(hdlm, qid_dat_agg, create2GrpIdx(hdlm, qid_dat_agg, "age_cat", "race"))
plotDLM(raceageDLM, 2) + theme_bw(base_size = 24)

racedualDLM <- estDLM(hdlm, qid_dat_agg, create2GrpIdx(hdlm, qid_dat_agg, "dual_any", "race"))
plotDLM(racedualDLM, 2) + theme_bw(base_size = 24)

agedualDLM <- estDLM(hdlm, qid_dat_agg, create2GrpIdx(hdlm, qid_dat_agg, "dual_any", "age_cat"))
plotDLM(agedualDLM, 2) + theme_bw(base_size = 24)

raceregionDLM <- estDLM(hdlm, qid_dat_agg, create2GrpIdx(hdlm, qid_dat_agg, "race", "region"))
plotDLM(raceregionDLM, 2) + theme_bw(base_size = 24)

ageregionDLM <- estDLM(hdlm, qid_dat_agg, create2GrpIdx(hdlm, qid_dat_agg, "age_cat", "region"))
plotDLM(ageregionDLM, 2) + theme_bw(base_size = 24)
