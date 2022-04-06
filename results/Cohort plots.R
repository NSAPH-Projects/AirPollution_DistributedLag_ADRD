library(ggplot2)
n_exit <- qid_dat[, .N, by = year][order(year)]
ggplot(n_exit[year>=2010], aes(x = year, y = N/1000000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  scale_x_continuous(breaks = 2010:2016) +
  labs(x = "Year", 
       y = "# enrollees since 2000 remaining\nin study (Millions)")

n_hosp <- qid_dat[, .(N = sum(first_hosp)), by = year][order(year)]
ggplot(n_hosp[year>=2010], aes(x = year, y = N/1000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  scale_x_continuous(breaks = 2010:2016) +
  labs(x = "Year", 
       y = "# of first ADRD hospitalizations (thousands)")


n_dead <- qid_dat[first_hosp==0, .(N = sum(dead)), by = year][order(year)]
ggplot(n_hosp[year>=2009], aes(x = year, y = N/1000)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw(base_size = 20) +
  labs(x = "Year", 
       y = "# of first ADRD hospitalizations (thousands)")
