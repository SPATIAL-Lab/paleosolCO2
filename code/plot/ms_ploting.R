library(tidyverse)
library(ggpubr)
library(readxl)
source("code/constructors.R")
theme = theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 10),
              plot.title = element_text(hjust = 0.1, vjust = -10))
# function to create data frame used to draw parameter curves
fm = function(site, var, param) {
  if(site == "Luochuan") {
    age = read.csv("data/loess_interglacial.csv")
  } else {
    age = read.csv("data/loess_glacial.csv") %>% filter(section == site)
  }
  post.var = data.frame(cbind(site, age$age,
                              t(apply(param$BUGSoutput$sims.list[[var]], 2, quantile, 
                                      c(0.05, 0.25, 0.5, 0.75, 0.95)))))
  names(post.var) = c("site", "age", "x5", "x25", "median", "x75", "x95")
  post.var[, 2:7] = lapply(post.var[, 2:7], as.numeric)
  results = post.var
}

# input params constraints ----
load("out/ms_fx_1e4.rda")
base = post.clp
load("out/ms_params_constraint/ms_fx_MAP_1e4.rda")
param = post.clp
dat.base = fm("Fuxian", "Tsoil", base)
dat.param = fm("Fuxian", "Tsoil", param)
ggplot(dat.base, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = dat.param, aes(x = age, y = median, ymin = x25, ymax = x75), fill = "royalblue1", alpha = 0.3) +
  geom_point(color = "firebrick2", shape = 21, size = 3) +
  geom_point(data = dat.param, aes(x = age, y = median), color = "royalblue", shape = 21, size = 3) +
  theme_bw() + theme +
  ggtitle("MAP") + 
  theme(plot.title = element_text(hjust = 0.9)) +
  labs(x = "Age (Ma)", 
       # y = expression(italic(p)*"CO"[2]*" (ppm)")
       # y = "MAP (mm)"
       y = expression(paste("T"[soil]*" (", degree, "C)"))
       ) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
ggsave("figure/prior_constraints_ms/MAP_MAT.jpg", width = 4.9, height = 3.5)

# inverted vs bayes w/o age model ----
# glacial data
inv = read_xlsx("data/Dataset S1.xlsx", sheet = 2)
inv = inv[1:141, c(1, 3, 11:13)]
names(inv) = c("section", "age", "co2", "low", "high")
inv$age = inv$age / 1000
zjc.inv = inv %>% filter(section == "Zhaojiachuan")
fx.inv = inv %>% filter(section == "Fuxian")

load("out/ms_fx_1e4.rda")
fx = post.clp
fx.ms = fm("Fuxian", "pCO2", fx)
p1 = ggplot(fx.ms, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = fx.inv, aes(x = age, y = co2,
                                  ymin = co2 - low, ymax = co2 + high), fill = "royalblue1", alpha = 0.3) +
  geom_line(color = "firebrick2") +
  geom_line(data = fx.inv, aes(x = age, y = co2), color = "royalblue") +
  theme_bw() + theme +
  ggtitle("Fuxian") +
  labs(x = "Age (Ma)", y = expression(italic(p)*"CO"[2]*" (ppm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p1

load("out/ms_zjc_1e4.rda")
zjc = post.clp
zjc.ms = fm("Zhaojiachuan", "pCO2", zjc)
p2 = ggplot(zjc.ms, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "salmon", alpha = 0.3) +
  geom_ribbon(data = zjc.inv, aes(x = age, y = co2,
                                  ymin = co2 - low, ymax = co2 + high), fill = "royalblue1", alpha = 0.3) +
  geom_line(color = "firebrick2") +
  geom_line(data = zjc.inv, aes(x = age, y = co2), color = "royalblue") +
  theme_bw() + theme +
  ggtitle("Zhaojiachuan") +
  labs(x = "Age (Ma)", y = expression(italic(p)*"CO"[2]*" (ppm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p2

inv = read_xlsx("data/lc_inv.xlsx", sheet = 2)
inv = inv[, c("age", "CO2", "CO2.low", "CO2.high")]
load("out/ms_lc_1e4.rda")
lc = post.clp
lc.ms = fm("Luochuan", "pCO2", lc)
p3 = ggplot(lc.ms, aes(x = age, y = median)) +
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
zjc.ms = fm("Zhaojiachuan", "S_z", zjc)
zjc.d18 = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
zjc.sz = data.frame(cbind("Zhaojiachuan", zjc.d18$d18c, zjc.ms$median, zjc.ms$x25, zjc.ms$x75))
load("out/ms_fx_1e4.rda")
fx = post.clp
fx.ms = fm("Fuxian", "S_z", fx)
fx.d18 = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.sz = data.frame(cbind("Fuxian", fx.d18$d18c, fx.ms$median, fx.ms$x25, fx.ms$x75))
dat = rbind(zjc.sz, fx.sz)
names(dat) = c("site", "d18c", "Sz", "Sz.low", "Sz.high")
dat[2:5] = lapply(dat[, 2:5], as.numeric)
m1 = nls(data = dat, Sz ~ a*exp(-b*d18c), start = list(a = 30, b = 0.3))
dat$Sz_pred = predict(m1, newdata = dat)
p1 = ggplot(dat, aes(x = d18c, y = Sz)) +
  geom_errorbar(aes(ymin = Sz.low, ymax = Sz.high), size = 0.2, width = 0, color = "ivory3") +
  geom_point(aes(fill = site), size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  geom_line(aes(x = d18c, y = Sz_pred), linetype = "dashed", linewidth = 1) +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O"[c]*" (\u2030)"),
       y = expression("S"[(z)]*" (ppm)"), fill = "") +
  scale_y_continuous(limits = c(200, 1100))
p1
# inverse model
inv = read_xlsx("data/Dataset S1.xlsx")
inv = inv[1:40, c(1,6,11,12)]
names(inv) = c("site", "d18c", "Sz", "Sz.sd")
m2 = nls(data = inv, Sz ~ a*exp(-b*d18c), start = list(a = 30, b = 0.3))
inv$Sz_pred = predict(m2, newdata = inv)
p2 = ggplot(inv, aes(x = d18c, y = Sz)) +
  geom_errorbar(aes(ymin = Sz - Sz.sd, ymax = Sz + Sz.sd), size = 0.2, width = 0, color = "ivory3") +
  geom_point(aes(fill = site), size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  geom_line(aes(x = d18c, y = Sz_pred), linetype = "dashed", linewidth = 1) +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O"[c]*" (\u2030)"),
       y = expression("S"[(z)]*" (ppm)"), fill = "") +
  scale_y_continuous(limits = c(200, 1100))
ggarrange(p2, p1, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figure/d18c_Sz_model.jpeg", width = 5.1, height = 3.1)

# other parameters ----
load("out/ms_zjc_1e4_v2.rda")
zjc = post.clp
zjc.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
load("out/ms_fx_1e4_v2.rda")
fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
load("out/ms_lc_1e4_v2.rda")
lc = post.clp
lc.age = read.csv("data/loess_interglacial.csv")

zjc.MAP = data.frame(cbind("Zhaojiachuan", zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
fx.MAP = data.frame(cbind("Fuxian", fx.age$age, 
                           t(apply(fx$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
lc.MAP = data.frame(cbind("Luochuan", lc.age$age, 
                          t(apply(lc$BUGSoutput$sims.list$MAP, 2, quantile, 
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
  ggtitle("MAP") +
  labs(x = "Age (Ma)", y = expression("MAP (mm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))

## Fuxian + D47 ----
load("out/ms_fx_1e4_v2.rda")
fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
fx.age = fx.age[order(fx.age$age),]
load("out/ms_fx_1e4_D47.v2.rda")
fx47 = post.clp
fx47.age = read.csv("data/D47.csv")

fx1 = data.frame(cbind("no", fx.age$age, 
                          t(apply(fx$BUGSoutput$sims.list$MAP, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx1) = c("D47", "age", "x5", "x25", "median", "x75", "x95")
fx2 = data.frame(cbind("yes", fx47.age$age, 
                           t(apply(fx47$BUGSoutput$sims.list$MAP, 2, quantile, 
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
       fill = expression(Delta[47]),
       color = expression(Delta[47]),
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
load("out/ms_fx_1e4_v2.rda")

plot.jpi(ai, post.clp$BUGSoutput$sims.list$pCO2, n = 100)
lines(ai, post.clp$BUGSoutput$median$pCO2, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$S_z, n = 100, ylim = c(0, 3000))
lines(ai, post.clp$BUGSoutput$median$S_z, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$Tsoil, n = 100)
lines(ai, post.clp$BUGSoutput$median$Tsoil, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$L, n = 100)
lines(ai, post.clp$BUGSoutput$median$L, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$MAP, n = 100)
lines(ai, post.clp$BUGSoutput$median$MAP, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$d18p, n = 100)
lines(ai, post.clp$BUGSoutput$median$d18p, col="red", lwd = 5)
plot.jpi(ai, post.clp$BUGSoutput$sims.list$pore, n = 100)
lines(ai, post.clp$BUGSoutput$median$pore, col="red", lwd = 5)

