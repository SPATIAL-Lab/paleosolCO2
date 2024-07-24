library(tidyverse)
library(ggpubr)

# pCO2 ----
load("out/clp1e4_ms_fx.rda")
fuxian = post.clp1
fx.age = read.csv("data/data.csv") %>% 
  filter(site == "Fuxian" & age < 2.6) %>%
  drop_na(d13Co)
fx.age = fx.age$age
fx.co2 = data.frame(cbind("fuxian", fx.age, 
                          t(apply(fuxian$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/clp1e4_ms_zjc.rda")
zhaojiachuan = post.clp1
zjc.age = read.csv("data/data.csv") %>% 
  filter(site == "Zhaojiachuan" & age < 2.6)
zjc.age = zjc.age$age
zjc.co2 = data.frame(cbind("zhaojiachuan", zjc.age, 
                          t(apply(zhaojiachuan$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/clp1e4_ms_lc.rda")
luochuan = post.clp1
lc.age = read.csv("data/data.csv") %>% 
  filter(site == "Luochuan" & age < 2.6)
lc.age = lc.age$age
lc.co2 = data.frame(cbind("luochuan", lc.age, 
                           t(apply(luochuan$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

post.co2 = rbind(fx.co2, zjc.co2, lc.co2) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p1 = ggplot(post.co2, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = expression("CO"[2]))

# MAP ----
fx.map = data.frame(cbind("fuxian", fx.age, 
                          t(apply(fuxian$BUGSoutput$sims.list$MAP, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")
zjc.map = data.frame(cbind("zhaojiachuan", fx.age, 
                          t(apply(zhaojiachuan$BUGSoutput$sims.list$MAP, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")
lc.map = data.frame(cbind("luochuan", lc.age, 
                          t(apply(luochuan$BUGSoutput$sims.list$MAP, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")
post.map = rbind(fx.map, zjc.map, lc.map) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p2 = ggplot(post.map, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = "MAP")

# MAT ----
fx.mat = data.frame(cbind("fuxian", fx.age, 
                          t(apply(fuxian$BUGSoutput$sims.list$MAT, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.mat) = c("site", "age", "x5", "x25", "median", "x75", "x95")
zjc.mat = data.frame(cbind("zhaojiachuan", fx.age, 
                           t(apply(zhaojiachuan$BUGSoutput$sims.list$MAT, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.mat) = c("site", "age", "x5", "x25", "median", "x75", "x95")
lc.mat = data.frame(cbind("luochuan", lc.age, 
                          t(apply(luochuan$BUGSoutput$sims.list$MAT, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.mat) = c("site", "age", "x5", "x25", "median", "x75", "x95")
post.mat = rbind(fx.mat, zjc.mat, lc.mat) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p3 = ggplot(post.mat, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "age (Ma)", y = "MAT")

# Sz ----
fx.sz = data.frame(cbind("fuxian", fx.age, 
                          t(apply(fuxian$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")
zjc.sz = data.frame(cbind("zhaojiachuan", fx.age, 
                           t(apply(zhaojiachuan$BUGSoutput$sims.list$S_z, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")
lc.sz = data.frame(cbind("luochuan", lc.age, 
                          t(apply(luochuan$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")
post.sz = rbind(fx.sz, zjc.sz, lc.sz) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p4 = ggplot(post.sz, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 3)) +
  scale_y_log10() +
  labs(x = "age (Ma)", y = expression("S"[(z)]))

ggarrange(p1, p2, p3, p4, nrow = 4, ncol = 1, common.legend = TRUE, align = "hv")
ggsave("figure/post.ms.jpg", height = 9, width = 5)


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

