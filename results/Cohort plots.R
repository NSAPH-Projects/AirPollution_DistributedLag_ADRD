n_exit <- qid_entry_exit[entry==2000, .N, by =exit][order(exit)]
ggplot(n_exit[exit>=2009], aes(x = exit, y = N/1000000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  labs(x = "Year exited cohort", y = "# enrolled since 2000 without\n ADRD hospitalization (Millions)")

n_exit2 <- n_exit[exit >= 2009, sum(N)] - c(0, cumsum(n_exit[exit>=2009, N]))
ggplot(data = data.frame(exit = 2009:2016, N = n_exit2[1:8]), 
       aes(x = exit, y = N/1000000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  labs(x = "Year", 
       y = "# enrollees since 2000 remaining\nin study (Millions)")

qid_entry_exit[entry==2000 & exit >= 2009, .N, by = exit
               ][, .(N = sum(N * (exit - 2008)))]

n_hosp<-hosp_dat[,.N,by=year]
ggplot(n_hosp[year>=2009], aes(x = year, y = N/1000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  labs(x = "Year", 
       y = "# of first ADRD hospitalizations (thousands)")
