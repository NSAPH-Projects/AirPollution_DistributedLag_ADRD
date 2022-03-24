# ############################################################################ #
#' Project: ADRD Distributed Lag Analysis                                      
#' Code: extract first hospitalization for AD/ADRD based on CCW algorithm      
#' Inputs: ADRD hospitalization data (2000-2016)                                 
#' Outputs: first AD/ADRD primary and secondary hospitalization                
#' Author: Daniel Mork                                                         
#' Last updated: Mar 09, 2022                                                  
#' Memory to run: 32 GB
# ############################################################################ #
rm(list = ls())
gc()

##### 0. Setup #####
library(data.table)
library(fst)
library(NSAPHutils)
options(stringsAsFactors = FALSE)

setDTthreads(threads = 8)

dir_data <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/dmork/Data/DLM_ADRD/"
dir_hosp <- "/nfs/home/D/dam9096/shared_space/ci3_analysis/data_ADRDhospitalization/ADRDhospitalization_CCWlist/"

##### 1. Load hospitalization data #####
f <- list.files(dir_hosp, pattern = "\\.fst", full.names = TRUE)
ADRDhosp <- rbindlist(lapply(f, read_fst, as.data.table = TRUE))
setkey(ADRDhosp, year, QID)


##### 2. Generate first hospitalization #####
firstPrimaryAD <- ADRDhosp[AD_primary == TRUE, 
                           .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstPrimaryAD, paste0(dir_data, "hospitalization/First_hosp_AD_primary.fst"))
# year     N
# 1: 2000 31828
# 2: 2001 35022
# 3: 2002 36785
# 4: 2003 37389
# 5: 2004 37905
# 6: 2005 39923
# 7: 2006 39978
# 8: 2007 39690
# 9: 2008 37435
# 10: 2009 34748
# 11: 2010 33606
# 12: 2011 31292
# 13: 2012 27135
# 14: 2013 23689
# 15: 2014 21824
# 16: 2015 22186
# 17: 2016 20701

firstPrimaryADRD <- ADRDhosp[ADRD_primary == TRUE, 
                             .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstPrimaryADRD, paste0(dir_data, "hospitalization/First_hosp_ADRD_primary.fst"))
# year     N
# 1: 2000 77270
# 2: 2001 76896
# 3: 2002 76369
# 4: 2003 76083
# 5: 2004 73663
# 6: 2005 72002
# 7: 2006 66686
# 8: 2007 64397
# 9: 2008 60622
# 10: 2009 57688
# 11: 2010 58709
# 12: 2011 54904
# 13: 2012 53489
# 14: 2013 49532
# 15: 2014 48167
# 16: 2015 49417
# 17: 2016 46964

firstAnyAD <- ADRDhosp[AD_any == TRUE, 
                       .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstAnyAD, paste0(dir_data, "hospitalization/First_hosp_AD_any.fst"))
# year      N
# 1: 2000 260335
# 2: 2001 233988
# 3: 2002 226679
# 4: 2003 221818
# 5: 2004 221262
# 6: 2005 216285
# 7: 2006 200953
# 8: 2007 195254
# 9: 2008 185836
# 10: 2009 173329
# 11: 2010 172537
# 12: 2011 168330
# 13: 2012 151921
# 14: 2013 136947
# 15: 2014 128238
# 16: 2015 132184
# 17: 2016 127486

firstAnyADRD <- ADRDhosp[ADRD_any == TRUE, 
                         .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstAnyADRD, paste0(dir_data, "hospitalization/First_hosp_ADRD_any.fst"))
# year      N
# 1: 2000 687372
# 2: 2001 583704
# 3: 2002 544365
# 4: 2003 523328
# 5: 2004 497407
# 6: 2005 482350
# 7: 2006 458650
# 8: 2007 454507
# 9: 2008 443028
# 10: 2009 422160
# 11: 2010 439483
# 12: 2011 450657
# 13: 2012 428407
# 14: 2013 415478
# 15: 2014 414892
# 16: 2015 453065
# 17: 2016 454502


firstSecondaryAD <- ADRDhosp[AD_any == TRUE & AD_primary == FALSE, 
                       .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstSecondaryAD, paste0(dir_data, "hospitalization/First_hosp_AD_secondary.fst"))
# year      N
# 1: 2000 237205
# 2: 2001 214253
# 3: 2002 208526
# 4: 2003 204702
# 5: 2004 204375
# 6: 2005 198349
# 7: 2006 182726
# 8: 2007 177195
# 9: 2008 169367
# 10: 2009 157853
# 11: 2010 157165
# 12: 2011 154122
# 13: 2012 139400
# 14: 2013 125644
# 15: 2014 117397
# 16: 2015 120933
# 17: 2016 116927


firstSecondaryADRD <- ADRDhosp[ADRD_any == TRUE & ADRD_primary == FALSE, 
                             .(ADATE = min(ADATE), year = min(year)), by = QID]
write_fst(firstSecondaryADRD, paste0(dir_data, "hospitalization/First_hosp_ADRD_secondary.fst"))
# year      N
# 1: 2000 638973
# 2: 2001 552818
# 3: 2002 520616
# 4: 2003 502729
# 5: 2004 478663
# 6: 2005 465266
# 7: 2006 444205
# 8: 2007 440764
# 9: 2008 431309
# 10: 2009 410215
# 11: 2010 426739
# 12: 2011 439544
# 13: 2012 416575
# 14: 2013 404357
# 15: 2014 403531
# 16: 2015 441486
# 17: 2016 443060