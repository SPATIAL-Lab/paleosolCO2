library(tidyverse)
library(ggpubr)

# pCO2 ----
load("out/ms_rc_lt_1e3.rda")
lantian = post.clp
lt.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Lantian")
lt.age = lt.age$age
lt.co2 = data.frame(cbind("Lantian", lt.age, 
                          t(apply(lantian$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lt.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_sl_1e3.rda")
shilou = post.clp
sl.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Shilou")
sl.age = sl.age$age
sl.co2 = data.frame(cbind("Shilou", sl.age, 
                           t(apply(shilou$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(sl.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_jx_1e3.rda")
jiaxian = post.clp
jx.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Jiaxian")
jx.age = jx.age$age
jx.co2 = data.frame(cbind("Jiaxian", jx.age, 
                          t(apply(jiaxian$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(jx.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")

post.co2 = rbind(lt.co2, sl.co2, jx.co2) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p1 = ggplot(post.co2, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(2, 8)) +
  labs(x = "age (Ma)", y = expression("CO"[2]))
p1

# S(z) ----
load("out/ms_rc_lt_1e3.rda")
lantian = post.clp
lt.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Lantian")
lt.age = lt.age$age
lt.sz = data.frame(cbind("Lantian", lt.age, 
                          t(apply(lantian$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lt.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_sl_1e3.rda")
shilou = post.clp
sl.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Shilou")
sl.age = sl.age$age
sl.sz = data.frame(cbind("Shilou", sl.age, 
                          t(apply(shilou$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(sl.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_jx_1e3.rda")
jiaxian = post.clp
jx.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Jiaxian")
jx.age = jx.age$age
jx.sz = data.frame(cbind("Jiaxian", jx.age, 
                          t(apply(jiaxian$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(jx.sz) = c("site", "age", "x5", "x25", "median", "x75", "x95")

post.sz = rbind(lt.sz, sl.sz, jx.sz) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p2 = ggplot(post.sz, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(2, 8)) +
  scale_y_continuous(limits = c(0, 5000)) +
  labs(x = "age (Ma)", y = expression("S"[(z)]))
p2

# MAP ----
load("out/ms_rc_lt_1e3.rda")
lantian = post.clp
lt.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Lantian")
lt.age = lt.age$age
lt.map = data.frame(cbind("Lantian", lt.age, 
                         t(apply(lantian$BUGSoutput$sims.list$MAP, 2, quantile, 
                                 c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lt.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_sl_1e3.rda")
shilou = post.clp
sl.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Shilou")
sl.age = sl.age$age
sl.map = data.frame(cbind("Shilou", sl.age, 
                         t(apply(shilou$BUGSoutput$sims.list$MAP, 2, quantile, 
                                 c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(sl.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")

load("out/ms_rc_jx_1e3.rda")
jiaxian = post.clp
jx.age = read.csv("data/red_clay_data.csv") %>% 
  filter(site == "Jiaxian")
jx.age = jx.age$age
jx.map = data.frame(cbind("Jiaxian", jx.age, 
                         t(apply(jiaxian$BUGSoutput$sims.list$MAP, 2, quantile, 
                                 c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(jx.map) = c("site", "age", "x5", "x25", "median", "x75", "x95")

post.map = rbind(lt.map, sl.map, jx.map) %>%
  mutate(across(c("age", "x5", "x25", "median", "x75", "x95"), as.numeric))

p3 = ggplot(post.map, aes(x = age, y = median, fill = site)) +
  geom_errorbar(aes(ymin = x25, ymax = x75), size = 0.2) +
  geom_point(size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  scale_x_continuous(limits = c(2, 8)) +
  scale_y_continuous(limits = c(100, 800)) +
  labs(x = "age (Ma)", y = "MAP (mm)")
p3
