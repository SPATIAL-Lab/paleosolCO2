library(tidyverse)
library(ggpubr)
source("code/helpers.R")
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

# read data ----
load("out/ms_lc_2e5.rda")
ms = read.csv("data/data.csv") %>%
  filter(site == "Zhaojiachuan" & age < 2.6) %>% drop_na(d13Co)
age = unique(ms$age)
pco2 = post.clp$BUGSoutput$mean$pCO2
MAP = post.clp$BUGSoutput$mean$MAP
PPCQ = post.clp$BUGSoutput$mean$PPCQ
PCQ_pf = post.clp$BUGSoutput$mean$PCQ_pf
MAT = post.clp$BUGSoutput$mean$MAT
Tsoil = post.clp$BUGSoutput$mean$Tsoil
PCQ_to = post.clp$BUGSoutput$mean$PCQ_to
Sz = post.clp$BUGSoutput$mean$S_z
z_m = post.clp$BUGSoutput$mean$z_m
AET = post.clp$BUGSoutput$mean$AET_PCQ
ms = data.frame(age, pco2, MAP, MAT, Sz, PPCQ, PCQ_to, PCQ_pf, Tsoil, z_m, AET)

load("out/ts_zjc_1e5_30ppm_normal.rda")
age = seq(-2.6, 0, by = 0.1)
pco2 = post.clp$BUGSoutput$mean$pCO2
MAP = post.clp$BUGSoutput$mean$MAP
PPCQ = post.clp$BUGSoutput$mean$PPCQ
PCQ_pf = post.clp$BUGSoutput$mean$PCQ_pf
MAT = post.clp$BUGSoutput$mean$MAT
Tsoil = post.clp$BUGSoutput$mean$Tsoil
Sz = post.clp$BUGSoutput$mean$S_z
Tsoil = post.clp$BUGSoutput$mean$Tsoil
PCQ_to = post.clp$BUGSoutput$mean$PCQ_to
z_m = post.clp$BUGSoutput$mean$z_m
AET = post.clp$BUGSoutput$mean$AET_PCQ
ts = data.frame(age, pco2, MAP, MAT, Sz, PPCQ, PCQ_to, PCQ_pf, Tsoil, z_m, AET)
ts$age = -ts$age

# plot ----
p1 = ggplot() +
  geom_point(data = ms, aes(x = age, y = pco2), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = pco2), size = 3, color = pal[6]) +
  theme_bw()
p1
p2 = ggplot() +
  geom_point(data = ms, aes(x = age, y = MAP), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = MAP), size = 3, color = pal[6]) +
  theme_bw()
p2
p3 = ggplot() +
  geom_point(data = ms, aes(x = age, y = PPCQ), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = PPCQ), size = 3, color = pal[6]) +
  theme_bw()
p3
p4 = ggplot() +
  geom_point(data = ms, aes(x = age, y = MAT), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = MAT), size = 3, color = pal[6]) +
  theme_bw()
p4
p5 = ggplot() +
  geom_point(data = ms, aes(x = age, y = Tsoil), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = Tsoil), size = 3, color = pal[6]) +
  theme_bw()
p5
p6 = ggplot() +
  geom_point(data = ms, aes(x = age, y = Sz), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = Sz), size = 3, color = pal[6]) +
  scale_y_continuous(limits = c(0, 5000)) +
  theme_bw()
p6
p7 = ggplot() +
  geom_point(data = ms, aes(x = age, y = z_m), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = z_m), size = 3, color = pal[6]) +
  theme_bw()
p7

p8 = ggplot() +
  geom_point(data = ms, aes(x = age, y = PCQ_pf), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = PCQ_pf), size = 3, color = pal[6]) +
  theme_bw()
p8
p9 = ggplot() +
  geom_point(data = ms, aes(x = age, y = PCQ_to), size = 3, color = pal[2]) +
  geom_point(data = ts, aes(x = age, y = PCQ_to), size = 3, color = pal[6]) +
  theme_bw()
p9

ggarrange(p1, p2, p3, p4, p5, p6, p8, p9, nrow = 4, ncol = 2)
ggsave("figure/ms_ts_comparison_zjc_normal.jpg", width = 8.6, height = 9.4)
``