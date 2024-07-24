library(tidyverse)
library(ggpubr)

# pCO2 ----
load("out/clp1e4_ms_zjc.rda")
zhaojiachuan = post.clp1
zjc.age = read.csv("data/data.csv") %>% 
  filter(site == "Zhaojiachuan" & age < 2.6)
zjc.age = zjc.age$age
zjc.co2 = data.frame(cbind("yes", zjc.age, 
                          t(apply(zhaojiachuan$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.co2) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")

load("out/clp1e4_ms_zjc_nod18O.rda")
zhaojiachuan2 = post.clp1
zjc.age = read.csv("data/data.csv") %>% 
  filter(site == "Zhaojiachuan" & age < 2.6)
zjc.age = zjc.age$age
zjc.co2.w = data.frame(cbind("no", zjc.age, 
                           t(apply(zhaojiachuan2$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.co2.w) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")

post.co2 = rbind(zjc.co2, zjc.co2.w) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p1 = ggplot(post.co2, aes(x = age, y = median, fill = d18O)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2, color = "ivory3") +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = expression("CO"[2]))

# MAP ----
zjc.map = data.frame(cbind("yes", zjc.age, 
                          t(apply(zhaojiachuan$BUGSoutput$sims.list$MAP, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.map) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
zjc.map.w = data.frame(cbind("no", zjc.age, 
                           t(apply(zhaojiachuan2$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.map.w) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
post.map = rbind(zjc.map, zjc.map.w) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p2 = ggplot(post.map, aes(x = age, y = median, fill = d18O)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2, color = "ivory3") +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = "MAP")

# MAT ----
zjc.mat = data.frame(cbind("yes", zjc.age, 
                           t(apply(zhaojiachuan$BUGSoutput$sims.list$MAT, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.mat) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
zjc.mat.w = data.frame(cbind("no", zjc.age, 
                           t(apply(zhaojiachuan2$BUGSoutput$sims.list$MAT, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.mat.w) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
post.mat = rbind(zjc.mat, zjc.mat.w) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p3 = ggplot(post.mat, aes(x = age, y = median, fill = d18O)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = "MAT")

# Sz ----
zjc.sz = data.frame(cbind("yes", zjc.age, 
                           t(apply(zhaojiachuan$BUGSoutput$sims.list$S_z, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.sz) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
zjc.sz.w = data.frame(cbind("no", zjc.age, 
                          t(apply(zhaojiachuan2$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.sz.w) = c("d18O", "age", "x5", "x25", "median", "x75", "x95")
post.sz = rbind(zjc.sz, zjc.sz.w) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p4 = ggplot(post.sz, aes(x = age, y = median, fill = d18O)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  scale_y_log10() +
  labs(x = "age (Ma)", y = expression("S"[(z)]))

ggarrange(p1, p2, p3, p4, nrow = 4, ncol = 1, common.legend = TRUE, align = "hv")
ggsave("figure/post.ms.zjc_d18O.jpg", height = 9, width = 5)


# iteration ----
source("code/constructors.R")
source("code/helpers.R")
load("out/clp1e4_ms.rda")

post = post.clp1
ages = age$age
ages = unique(ages)
plot.jpi(ages, post$BUGSoutput$sims.list$pCO2, n = 100)
plot.jpi(ages, post$BUGSoutput$sims.list$MAP, n = 100)
plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_pf, n = 100)
plot.jpi(ages, post$BUGSoutput$sims.list$MAT)
plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_to)
plot.jpi(ages, post$BUGSoutput$sims.list$S_z)
plot.jpi(ages, post$BUGSoutput$sims.list$f_R)
plot.jpi(ages, post$BUGSoutput$sims.list$d18.p)

