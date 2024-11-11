library(tidyverse)
library(ggpubr)
library(readxl)
source("code/constructors.R")
source("code/helpers.R")

theme = theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 10),
              plot.title = element_text(hjust = 0.1, vjust = -10))
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

# View(post$BUGSoutput$summary)

## plot proxies ----
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
## plot w/ iterations ----
# for output without age model
load("out/ts_lc_1e5_30ppm_MS.rda")
dt = 0.1
ages = seq(-2.6, 0, by = dt)
# ages = -ages
# dat.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan") %>% filter(age < 2600)
# dat.age = - dat.age$age / 1000

# ms = read.csv("data/data.csv") %>%
#   filter(site == "Luochuan" & age < 2.6) %>% drop_na(d13Co)
# ages = sort(-unique(ms$age), decreasing = TRUE)
# ages = ages$ts

plot(ages, post.clp$BUGSoutput$sims.list$pCO2[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(italic(p)*"CO"[2]), xlim = range(ages), ylim = c(100,500), col=rgb(red=0, green=0, blue=0, alpha=0.1), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$pCO2[i,], col=rgb(red=0, green=0, blue=0, alpha=0.3), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$pCO2, col="red", lwd = 5)
# dat.rs = cbind(dat.age, 100)
# points(dat.rs[, 1], dat.rs[, 2], pch = 24, cex = 1)
# points(ages, post.clp$BUGSoutput$median$pCO2, col="red")
axis(2)
axis(1)

plot(ages, post.clp$BUGSoutput$sims.list$Tsoil[1,], type="l", axes = FALSE,
     xlab = "Age (Ma)", ylab = expression("Tsoil"), 
     xlim = range(ages), ylim = c(10, 25), 
     col=rgb(red=0, green=0, blue=0, alpha=0.1), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$Tsoil[i,], col=rgb(red=0, green=0, blue=0, alpha=0.3), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$Tsoil, col="palegreen", lwd=5)
axis(1)
axis(2)
# dat.rs = cbind(dat.age, 10)
# points(dat.rs[, 1], dat.rs[, 2], pch = 24, cex = 1)


plot(ages, post.clp$BUGSoutput$sims.list$MAP[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("MAP")), xlim = range(ages), ylim = c(100, 750), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$MAP[i,], col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$MAP, col="deepskyblue2", lwd=5)
axis(1)
axis(2)

plot(ages, post.clp$BUGSoutput$sims.list$PPCQ[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("PPCQ")), xlim = range(ages), ylim = c(0, 500), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$PPCQ[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$PPCQ, col="deepskyblue2", lwd=2)
axis(1)
axis(2)
# dat.rs = cbind(dat.age, 0)
# points(dat.rs[, 1], dat.rs[, 2], pch = 24, cex = 1)

plot(ages, post.clp$BUGSoutput$sims.list$S_z[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("S(z)")), xlim = range(ages), ylim = c(0, 5000), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$S_z[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$S_z, col="darkgoldenrod2", lwd=2)
axis(1)
axis(2)

plot(ages, post.clp$BUGSoutput$sims.list$MAT[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("MAT")), xlim = range(ages), ylim = c(0, 20), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$MAT[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$MAT, col="darkgoldenrod2", lwd=2)
axis(1)
axis(2)

plot(ages, post.clp$BUGSoutput$sims.list$PCQ_to[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("PCQ_to")), xlim = range(ages), ylim = c(10, 16), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$PCQ_to[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$PCQ_to, col="palegreen", lwd=2)
axis(1)
axis(2)




## plot w/ uncertainties ----
# load data
# w/o age model
load("out/ms_fx_1e4.rda")
ms.fx = post.clp
fx.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
load("out/ms_zjc_1e4.rda")
ms.zjc = post.clp
zjc.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
load("out/ms_lc_1e4.rda")
ms.lc = post.clp
lc.age = read.csv("data/loess_interglacial.csv")

# w/ age model
load("out/ts_fx_1e5_30ppm.rda")
ts.fx = post.clp
load("out/ts_zjc_1e5_30ppm.rda")
ts.zjc = post.clp
load("out/ts_lc_1e5_30ppm.rda")
ts.lc = post.clp
age = seq(-2.6, 0, by = 0.1)

# Any observations
fx.ms.param = data.frame(cbind("no", fx.age$age, 
                           t(apply(ms.fx$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.ms.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
fx.ms.param[,2:7] = lapply(fx.ms.param[,2:7], as.numeric)
fx.ts.param = data.frame(cbind("yes", age, 
                           t(apply(ts.fx$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.ts.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
fx.ts.param[,2:7] = lapply(fx.ts.param[,2:7], as.numeric)
fx.ts.param$age = - fx.ts.param$age
fx.co2 = rbind(fx.ms.param, fx.ts.param)
p1 = ggplot(fx.co2, aes(x = age, y = median, fill = age_model)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = age_model)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = expression("MAP (mm)")) +
  ggtitle("Fuxian") +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
zjc.ms.param = data.frame(cbind("no", zjc.age$age, 
                             t(apply(ms.zjc$BUGSoutput$sims.list$MAP, 2, quantile, 
                                     c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.ms.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
zjc.ms.param[,2:7] = lapply(zjc.ms.param[,2:7], as.numeric)
zjc.ts.param = data.frame(cbind("yes", age, 
                             t(apply(ts.zjc$BUGSoutput$sims.list$MAP, 2, quantile, 
                                     c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(zjc.ts.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
zjc.ts.param[,2:7] = lapply(zjc.ts.param[,2:7], as.numeric)
zjc.ts.param$age = - zjc.ts.param$age
zjc.co2 = rbind(zjc.ms.param, zjc.ts.param)
p2 = ggplot(zjc.co2, aes(x = age, y = median, fill = age_model)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = age_model)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  ggtitle("Zhaojiachuan") +
  labs(x = "Age (Ma)", y = expression("MAP (mm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))

ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figure/ms_ts_comparison_MAP_normal.jpg", width = 7.3, height = 3.5)

# interglacial
lc.ms.param = data.frame(cbind("no", lc.age$age, 
                                t(apply(ms.lc$BUGSoutput$sims.list$MAP, 2, quantile, 
                                        c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.ms.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
lc.ms.param[,2:7] = lapply(lc.ms.param[,2:7], as.numeric)
lc.ts.param = data.frame(cbind("yes", age, 
                                t(apply(ts.lc$BUGSoutput$sims.list$MAP, 2, quantile, 
                                        c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(lc.ts.param) = c("age_model", "age", "x5", "x25", "median", "x75", "x95")
lc.ts.param[,2:7] = lapply(lc.ts.param[,2:7], as.numeric)
lc.ts.param$age = -lc.ts.param$age
lc.co2 = rbind(lc.ms.param, lc.ts.param)
ggplot(lc.co2, aes(x = age, y = median, fill = age_model)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = age_model)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  ggtitle("Luochuan") +
  labs(x = "Age (Ma)", y = expression("MAP (mm)")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5)) +
  guides(color = "none", fill = "none")

# lc.co2 = data.frame(cbind(lc.age$age, 
#                           t(apply(lc$BUGSoutput$sims.list$pCO2, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
# lc.co2$site = "Luochuan"

fx.co2 = data.frame(cbind(fx.age$age, 
                          t(apply(fx$BUGSoutput$sims.list$pCO2, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
fx.co2$site = "Fuxian"
zjc.co2 = data.frame(cbind(zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$pCO2, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
zjc.co2$site = "Zhaojiachuan"
clp.co2 = rbind(lc.co2, fx.co2, zjc.co2)

lc.mat = data.frame(cbind(lc.age$age, 
                          t(apply(lc$BUGSoutput$sims.list$MAT, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
lc.mat$site = "Luochuan"
fx.mat = data.frame(cbind(fx.age$age, 
                          t(apply(fx$BUGSoutput$sims.list$MAT, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
fx.mat$site = "Fuxian"
zjc.mat = data.frame(cbind(zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$MAT, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
zjc.mat$site = "Zhaojiachuan"
clp.mat = rbind(lc.mat, fx.mat, zjc.mat)

lc.map = data.frame(cbind(lc.age$age, 
                          t(apply(lc$BUGSoutput$sims.list$MAP, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
lc.map$site = "Luochuan"
fx.map = data.frame(cbind(fx.age$age, 
                          t(apply(fx$BUGSoutput$sims.list$MAP, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
fx.map$site = "Fuxian"
zjc.map = data.frame(cbind(zjc.age$age, 
                           t(apply(zjc$BUGSoutput$sims.list$MAP, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
zjc.map$site = "Zhaojiachuan"
clp.map = rbind(lc.map, fx.map, zjc.map)

site = pal[factor(clp.co2$site, levels = c("Luochuan", "Fuxian", "Zhaojiachuan"))]
#png("out/Curves.png", 7, 6, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(0, 3), ylim = c(0, 3.5), axes = FALSE, xlab = "", ylab = "")
legend(x = 2, y = 3.5, legend = c("Luochuan", "Fuxian", "Zhaojiachuan"),
       col = pal, pch = 16, cex = 0.8, pt.cex = 1.5)

yext = range(clp.co2[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 1) * 100
clp.co2rs = cbind(clp.co2[, 1], clp.co2$site,
                  1.8 + (clp.co2[, 2:6] - min(tix)) / diff(range(tix)))
arrows(clp.co2rs[, 1], clp.co2rs[, 4], clp.co2rs[, 1], clp.co2rs[, 6], col = "ivory2",
       angle=90, length=0, code = 0)
points(clp.co2rs[, 1], clp.co2rs[, 5], col = "black", bg = site, pch = 21, cex = 1)
axis(2, 1.8 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("pCO"[2]*" (ppmv)"), 2, line = 2.5, at = 2.2)

yext = range(clp.mat[, 2:6])
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 5)
clp.matrs = cbind(clp.mat[, 1], clp.mat$site,
                  1 + (clp.mat[, 2:6] - min(tix)) / diff(range(tix)))
arrows(clp.matrs[, 1], clp.matrs[, 4], clp.matrs[, 1], clp.matrs[, 6], col = "ivory2",
       angle=90, length=0, code = 0)
points(clp.matrs[, 1], clp.matrs[, 5], col = "black", bg = site, pch = 21, cex = 1)
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("MAT"), 4, line = 2.5, at = 1.5)

yext = range(clp.map[, 2:6])
tix = seq(floor(min(yext)-10), 
          ceiling(max(yext)), by = 100)
clp.maprs = cbind(clp.map[, 1], clp.map$site,
                  0 + (clp.map[, 2:6] - min(tix)) / diff(range(tix)))
arrows(clp.maprs[, 1], clp.maprs[, 4], clp.maprs[, 1], clp.maprs[, 6], col = "ivory2",
       angle=90, length=0, code = 0)
points(clp.maprs[, 1], clp.maprs[, 5], col = "black", bg = site, pch = 21, cex = 1)
axis(2, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("MAT"), 2, line = 2.5, at = 0.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

## plot w/ and w/o D47 ----
load("out/ts_fx_1e5_30ppm_normal_D47.rda")
fx.47 = post.clp
load("out/ts_fx_1e5_30ppm_normal.rda")
fx = post.clp
dt = 0.1
age = seq(-2.6, 0, by = dt)
fx.47.param = data.frame(cbind("yes", age, 
                               t(apply(fx.47$BUGSoutput$sims.list$S_z, 2, quantile, 
                                       c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.47.param) = c("D47", "age", "x5", "x25", "median", "x75", "x95")
fx.47.param[,2:7] = lapply(fx.47.param[,2:7], as.numeric)
fx.no47.param = data.frame(cbind("no", age, 
                               t(apply(fx$BUGSoutput$sims.list$S_z, 2, quantile, 
                                       c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(fx.no47.param) = c("D47", "age", "x5", "x25", "median", "x75", "x95")
fx.no47.param[,2:7] = lapply(fx.no47.param[,2:7], as.numeric)
fx.no47.param$age = fx.no47.param$age
fx.param = rbind(fx.47.param, fx.no47.param)
p1 = ggplot(fx.param, aes(x = age, y = median, fill = D47)) +
  geom_ribbon(aes(ymin = x25, ymax = x75), alpha = 0.3) +
  geom_line(aes(color = D47)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = expression("S(z) (ppmv)"),
       fill = expression(Delta[47]),
       color = expression(Delta[47])) +
  ggtitle("Fuxian") +
  scale_x_continuous(breaks = seq(-2.5, 0, 0.5))
p1

## plot.jpg ---- 

post = post.clp
plot.jpi(ages, post$BUGSoutput$sims.list$pCO2, ylim = c(0, 500))

splot.jpi(ages, log10(post$BUGSoutput$sims.list$S_z))

plot.jpi(ages, post$BUGSoutput$sims.list$f_R)

plot.jpi(ages, post$BUGSoutput$sims.list$MAT)

plot.jpi(ages, post$BUGSoutput$sims.list$MAP)
# points(d13Cc$age, d13Cc$d13Cc)

# plot.jpi(ages, post$BUGSoutput$sims.list$GMT)

plot.jpi(ages, post$BUGSoutput$sims.list$d13Ca)

plot.jpi(ages, post$BUGSoutput$sims.list$MAP)

plot.jpi(ages, post$BUGSoutput$sims.list$tsc)

plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_to)

plot.jpi(ages, post$BUGSoutput$sims.list$Tsoil, n = 500)

plot.jpi(ages, post$BUGSoutput$sims.list$MAP)

plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_pf)

plot.jpi(ages, post$BUGSoutput$sims.list$ha)

plot.jpi(ages, post$BUGSoutput$sims.list$z_m)

# Parameters ----
co2.pri = density(runif(1e6, 100, 400))
co2 = post.clp$BUGSoutput$sims.list$pCO2
co2.post = density(post.clp$BUGSoutput$sims.list$pCO2[,20])

# png("out/Parms.png", 9, 5, "in", res = 300)
# layout(matrix(c(1, 2), nrow = 1))
par(mai = c(1, 0.2, 0.2, 0.2))
plot(co2.pri, xlim = range(co2.pri$x, co2.post$x),
     ylim = range(co2.pri$y, co2.post$y), main = "", axes = FALSE,
     xlab = expression("CO"[2]), lty = 2, lwd = 2)
axis(1)
axis(2, labels = FALSE)
box()
lines(co2.post, lwd = 2)

# plot(ages, post.clp$BUGSoutput$mean$MAP, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$PPCQ, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$ha, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$MAT, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$Tsoil, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$S_z, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$z_m, type = "l")
# plot(ages, post.clp$BUGSoutput$mean$d13Ca, type = "l")
