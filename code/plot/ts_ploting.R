library(tidyverse)
library(ggpubr)
library(readxl)
source("code/constructors.R")
source("code/helpers.R")
# function to create data frame used to plot parameter curves
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
fm.ts = function(site, var, param) {
  age = seq(-2.6, 0, by = 0.1)
  post.var = data.frame(cbind(site, age,
                              t(apply(param$BUGSoutput$sims.list[[var]], 2, quantile, 
                                      c(0.05, 0.25, 0.5, 0.75, 0.95)))))
  names(post.var) = c("site", "age", "x5", "x25", "median", "x75", "x95")
  post.var[, 2:7] = lapply(post.var[, 2:7], as.numeric)
  post.var$age = -post.var$age
  results = post.var
}

theme = theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 10),
              plot.title = element_text(hjust = 0.1, vjust = -10))

## plot proxies ----
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
par(mar = c(4, 4, 1, 4))
plot(0.2, 0, xlim = c(-3, 0), ylim = c(0, 4), axes = FALSE, xlab = "", ylab = "")

yext = range(clp$d13C)
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 2)
clp.d13Crs = cbind(-clp$age,
                  3 + (clp$d13C - min(tix)) / diff(range(tix)))
points(clp.d13Crs[, 1], clp.d13Crs[, 2], col = "black", bg = pal[1], pch = 21, cex = 1)
axis(2, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[c]*" (\u2030)"), 2, line = 2.5, at = 3.5)

yext = range(clp$d18O)
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 2)
clp.d18Ors = cbind(-clp$age,
                   2 + (clp$d18O - min(tix)) / diff(range(tix)))
points(clp.d18Ors[, 1], clp.d18Ors[, 2], col = "black", bg = pal[2], pch = 21, cex = 1)
axis(4, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[c]*" (\u2030)"), 4, line = 2.5, at = 2.5)

yext = range(clp$d13Co, na.rm = TRUE)
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 1)
clp.d13Cors = cbind(-clp$age,
                   1 + (clp$d13Co - min(tix)) / diff(range(tix)))
points(clp.d13Cors[, 1], clp.d13Cors[, 2], col = "black", bg = pal[3], pch = 21, cex = 1)
axis(2, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[o]*" (\u2030)"), 2, line = 2.5, at = 1.5)

yext = range(clp$D47, na.rm = TRUE)
tix = seq(ceiling(max(yext)*100), 
          floor(min(yext)*100), by = -1) / 100
clp.D47rs = cbind(-clp$age,
                    1 - (clp$D47 - min(tix)) / diff(range(tix)))
points(clp.D47rs[, 1], clp.D47rs[, 2], col = "black", bg = pal[5], pch = 21, cex = 1)
axis(4, 1 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta*"47 (\u2030)"), 4, line = 2.5, at = 0.5)

axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()
## comparison with published records ----
load("out/ts_fx_1e4_normal.rda")
ts.fx = post.clp
load("out/ts_zjc_1e5_normal.rda")
ts.zjc = post.clp
load("out/ts_lc_1e4_normal.rda")
ts.lc = post.clp
param = "MAP"
fx.ts = fm.ts("Fuxian", param, ts.fx)
zjc.ts = fm.ts("Zhaojiachuan", param, ts.zjc)
lc.ts = fm.ts("Luochuan", param, ts.lc)
dat = rbind(fx.ts, zjc.ts, lc.ts)
pco2 = read.csv("data/age_co2_data.csv") 
pco2 = pco2[, c(3,4,7:9)] 
names(pco2) = c("proxy", "age", "CO2", "high", "low")
MAP = read_xlsx("data/EASM_meng2018.xlsx", sheet = 4)
MAP = MAP %>% drop_na()
names(MAP) = c("age", "unit", "MIS", "index")
p1 = ggplot(dat, aes(x = age, y = median, fill = site)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = site)) +
  scale_fill_manual(values = c("firebrick2", "royalblue", "azure4")) +
  scale_color_manual(values = c("firebrick2", "royalblue", "azure4")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("JPI") + theme(plot.title = element_text(hjust = 0.9)) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p1
# p2 = ggplot(pco2, aes(x = age/1000, y = CO2, fill = proxy)) +
#   # geom_errorbar(aes(ymin = low, ymax = high), linewidth = 0.2, width = 0, color = "azure3") +
#   geom_point(shape = 21, size = 3) +
#   theme_bw() + theme +
#   labs(x = "Age (Ma)", y = param) +
#   ggtitle("inverted") + 
#   scale_x_continuous(breaks = seq(0, 2.5, 0.5)) + 
#   scale_y_continuous(limits = c(0, 500))
p2 = ggplot(MAP, aes(x = age, y = index)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("inverted") + theme(plot.title = element_text(hjust = 0.9)) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p2

ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figure/MAP_comparison.jpg", width = 7.7, height = 3.5)

## plot w/ uncertainties ----
# load data
# w/o age model
load("out/ms_fx_1e4.rda")
fx = post.clp
load("out/ms_zjc_1e4.rda")
zjc = post.clp
load("out/ms_lc_1e4.rda")
lc = post.clp
# w/ age model
load("out/ts_fx_1e5_normal.rda")
ts.fx = post.clp
load("out/ts_zjc_1e5_normal.rda")
ts.zjc = post.clp
load("out/ts_lc_1e5_normal.rda")
ts.lc = post.clp

# plot
param = "pCO2"
fx.ms = fm("Fuxian", param, fx)
fx.ts = fm.ts("Fuxian", param, ts.fx)
fx.ms$ts = "No"
fx.ts$ts = "Yes"
fx.param = rbind(fx.ms, fx.ts)
p1 = ggplot(fx.param, aes(x = age, y = median, fill = ts)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = ts)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("Fuxian") +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p1

zjc.ms = fm("Zhaojiachuan", param, zjc)
zjc.ts = fm.ts("Zhaojiachuan", param, ts.zjc)
zjc.ms$ts = "No"
zjc.ts$ts = "Yes"
zjc.param = rbind(zjc.ms, zjc.ts)
p2 = ggplot(zjc.param, aes(x = age, y = median, fill = ts)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = ts)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("Zhaojiachuan") +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p2

lc.ms = fm("Luochuan", param, lc)
lc.ts = fm.ts("Luochuan", param, ts.lc)
lc.ms$ts = "No"
lc.ts$ts = "Yes"
lc.param = rbind(lc.ms, lc.ts)
p3 = ggplot(lc.param, aes(x = age, y = median, fill = ts)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = ts)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("Luochuan") +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p3

ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE)
ggsave(paste("figure/ms_ts_", param, ".jpg", sep = ""), width = 9.5, height = 3.4)

## plot w/ and w/o D47 data ----
load("out/ts_fx_1e4_normal.rda")
fx = post.clp
load("out/ts_fx_1e4_D47_normal.rda")
fx.47 = post.clp
param = "pCO2"
fx.ts = fm.ts("Fuxian", param, fx)
fx.ts$D47 = "no"
fx.ts.47 = fm.ts("Fuxian", param, fx.47)
fx.ts.47$D47 = "yes"
dat = rbind(fx.ts, fx.ts.47)
p4 = ggplot(dat, aes(x = age, y = median, fill = D47)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = D47)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = param) +
  ggtitle("Fuxian") +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
p4

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figure/D47.jpg", width = 7.3, height = 5.5)
## plot the effect of tau ----
load("out/ts_params_1e4_normal/ts_fx_10ppm.rda")
small = post.clp
load("out/ts_params_1e4_normal/ts_fx_50ppm.rda")
median = post.clp
load("out/ts_params_1e4_normal/ts_fx_100ppm.rda")
large = post.clp
param = "pCO2"
post.s = fm.ts("Fuxian", param, small) %>% mutate(tau = "10 ppm")
post.m = fm.ts("Fuxian", param, median) %>% mutate(tau = "50 ppm")
post.l = fm.ts("Fuxian", param, large) %>% mutate(tau = "100 ppm")
dat = rbind(post.s, post.m, post.l)
dat$tau = factor(dat$tau, levels = c("10 ppm", "50 ppm", "100 ppm"))
ggplot(dat, aes(x = age, y = median, fill = tau)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = tau)) +
  scale_fill_manual(values = c("firebrick2", "royalblue", "azure4")) +
  scale_color_manual(values = c("firebrick2", "royalblue", "azure4")) +
  theme_bw() + theme +
  ggtitle("Fuxian") + theme(plot.title = element_text(hjust = 0.9)) +
  labs(x = "Age (Ma)", 
       y = expression(italic(p)*"CO"[2]*" (ppmv)"), 
       fill = expression(tau), color = expression(tau)) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))

# plot.jpg ----
load("out/ts_fx_1e5.rda")
param = "pCO2"
plot.jpi(ages, post.clp$BUGSoutput$sims.list[[param]], n = 1000, ylab = param)
lines(ages, post.clp$BUGSoutput$median[[param]], col="red", lwd = 5)



