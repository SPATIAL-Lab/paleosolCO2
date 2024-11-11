library(tidyverse)
library(ggpubr)
library(readxl)
source("code/constructors.R")

theme = theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 10),
              plot.title = element_text(hjust = 0.1, vjust = -10))

# inverted vs bayes w/o age model ----
# glacial data
inv = read_xlsx("data/Dataset S1.xlsx", sheet = 2)
inv = inv[1:141, c(1, 3, 11:13)]
names(inv) = c("section", "age", "co2", "low", "high")
inv$age = inv$age / 1000
zjc.inv = inv %>% filter(section == "Zhaojiachuan")
fx.inv = inv %>% filter(section == "Fuxian")

load("out/ms_zjc_1e4.rda")
zjc = post.clp
zjc.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
zjc.co2 = data.frame(cbind("zhaojiachuan", zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")
zjc.co2[,2:7] = lapply(zjc.co2[,2:7], as.numeric)

p1 = ggplot(zjc.co2, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = zjc.inv, aes(x = age, y = co2,
                                  ymin = co2 - low, ymax = co2 + high), fill = "royalblue1", alpha = 0.3) +
  geom_line(color = "firebrick2") +
  geom_line(data = zjc.inv, aes(x = age, y = co2), color = "royalblue") +
  theme_bw() + theme +
  ggtitle("Zhaojiachuan") +
  labs(x = "Age (Ma)", y = expression(italic(p)*"CO"[2]*" (ppm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p1

load("out/ms_fx_1e4.rda")
fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.co2 = data.frame(cbind("fuxian", fx.age$age, 
                           t(apply(fx$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")
fx.co2[,2:7] = lapply(fx.co2[,2:7], as.numeric)

p2 = ggplot(fx.co2, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = fx.inv, aes(x = age, y = co2,
                                  ymin = co2 - low, ymax = co2 + high), fill = "royalblue1", alpha = 0.3) +
  geom_line(color = "firebrick2") +
  geom_line(data = fx.inv, aes(x = age, y = co2), color = "royalblue") +
  theme_bw() + theme +
  ggtitle("Fuxian") +
  labs(x = "Age (Ma)", y = expression(italic(p)*"CO"[2]*" (ppm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p2
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv")

# interglacial data 
inv = read_xlsx("data/lc_inv.xlsx", sheet = 2)
inv = inv[, c("age", "CO2", "CO2.low", "CO2.high")]
load("out/ms_lc_1e4.rda")
lc = post.clp
lc.age = read.csv("data/loess_interglacial.csv")
lc.co2 = data.frame(cbind("luochuan", lc.age$age, 
                           t(apply(lc$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.co2) = c("site", "age", "x5", "x25", "median", "x75", "x95")
lc.co2[,2:7] = lapply(lc.co2[,2:7], as.numeric)
p3 = ggplot(lc.co2, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = inv, aes(x = age, y = CO2,
                                 ymin = CO2 - CO2.low, ymax = CO2 + CO2.high), fill = "royalblue1", alpha = 0.3) +
  geom_line(color = "firebrick2") +
  geom_line(data = inv, aes(x = age, y = CO2), color = "royalblue") +
  theme_bw() + theme +
  ggtitle("Luochuan") +
  labs(x = "Age (Ma)", y = expression(italic(p)*"CO"[2]*" (ppm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p3
ggarrange(p1, p2, p3, nrow = 1, ncol = 3, align = "hv")
ggsave("figure/ms_inv_comparison.jpg", width = 12, height = 3.5)

# d18c model ----
load("out/ms_zjc_1e4.rda")
zjc = post.clp
zjc.sz = post.clp$BUGSoutput$mean$S_z
zjc.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
zjc.age = zjc.age[order(zjc.age$age),]
zjc.dat = data.frame("Zhaojiachuan", zjc.age$d18c, zjc.sz)
names(zjc.dat) = c("site", "d18c", "Sz")
load("out/ms_fx_1e4.rda")
fx = post.clp
fx.sz = post.clp$BUGSoutput$mean$S_z
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.age = fx.age[order(fx.age$age),]
fx.dat = data.frame("Fuxian", fx.age$d18c, fx.sz)
names(fx.dat) = c("site", "d18c", "Sz")
dat = rbind(zjc.dat, fx.dat)
m1 = nls(data = dat, Sz ~ a*exp(-b*d18c), start = list(a = 30, b = 0.3))
dat$Sz_pred = predict(m1, newdata = dat)
ggplot(dat, aes(x = d18c, y = Sz)) +
  geom_point(aes(fill = site), size = 4, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  geom_line(aes(x = d18c, y = Sz_pred), linetype = "dashed", size = 1) +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O"[c]*" (\u2030)"),
       y = expression("S"[(z)]*" (ppm)"))

# other parameters ----
load("out/ms_zjc_1e4.rda")
zjc = post.clp
zjc.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
load("out/ms_fx_1e4.rda")
fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
load("out/ms_lc_1e4.rda")
lc = post.clp
lc.age = read.csv("data/loess_interglacial.csv")

zjc.MAP = data.frame(cbind("Zhaojiachuan", zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
fx.MAP = data.frame(cbind("Fuxian", fx.age$age, 
                           t(apply(fx$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
lc.MAP = data.frame(cbind("Luochuan", lc.age$age, 
                          t(apply(lc$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
MAP = rbind(zjc.MAP, fx.MAP, lc.MAP)
names(MAP) = c("site", "age", "x5", "x25", "median", "x75", "x95")
MAP[,2:7] = lapply(MAP[,2:7], as.numeric)

ggplot(MAP, aes(x = age, y = median, group = site, fill = site)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), alpha = 0.3) +
  geom_line(aes(color = site)) +
  scale_fill_manual(values = c("firebrick2", "royalblue", "gray")) +
  scale_color_manual(values = c("firebrick2", "royalblue", "black")) +
  theme_bw() + theme +
  ggtitle("Tsoil") +
  labs(x = "Age (Ma)", y = expression("Tsoil (degC)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))


## Fuxian + D47 ----
# CO2
load("out/ms_fx_1e4.rda")
fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.age = fx.age[order(fx.age$age),]
load("out/ms_fx_1e4_D47.rda")
fx.47 = post.clp
fx.age1 = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.age2 = read.csv("data/data.csv") %>% drop_na(D47)
ages = ts(fx.age1$age, fx.age2$age)
ai = ages$ts
fx1 = data.frame(cbind("no", fx.age$age, 
                          t(apply(fx$BUGSoutput$sims.list$DIFC, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx1) = c("D47", "age", "x5", "x25", "median", "x75", "x95")
fx2 = data.frame(cbind("yes", ai, 
                           t(apply(fx.47$BUGSoutput$sims.list$DIFC, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx2) = c("D47", "age", "x5", "x25", "median", "x75", "x95")
params = rbind(fx1, fx2)
params[,2:7] = lapply(params[,2:7], as.numeric)

ggplot(params, aes(x = age, y = median, fill = D47)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), alpha = 0.3) +
  geom_line(aes(color = D47)) +
  geom_point(aes(color = D47), shape = 21, size = 3, fill = "white") +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  # ggtitle("Zhaojiachuan") +
  labs(x = "Age (Ma)",
       # y = expression(italic(p)*"CO"[2]*" (ppm)")
       # y = expression(paste("T"[soil]*" (", degree, "C)"))
       # y = expression("S"[z]*" (ppmv)")
       y = expression("R (molC/cm"^"2"*"/s)")
       ) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5)) +
  scale_y_continuous()

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

