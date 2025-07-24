rm(list = ls())
source("code/helpers.R")
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

## plot w/ iterations ----
# for output without time-series model
load("out/ts_fx_1e5_30ppm_0.05Ma.rda")
dt = 0.05
ages = seq(-2.6, 0, by = dt)
dat.age = read.csv("data/loess_glacial.csv") %>% filter(section == "Fuxian") %>% filter(age < 2.6)
dat.age$age = -dat.age$age

# ms = read.csv("data/loess_interglacial.csv") %>%
#   filter(age < 2.6)
# ages = sort(-unique(ms$age), decreasing = TRUE)
ages = ages$ts

plot(ages, post.clp$BUGSoutput$sims.list$pCO2[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(italic(p)*"CO"[2]), xlim = range(ages), ylim = c(100,500), col=rgb(red=0, green=0, blue=0, alpha=0.1), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$pCO2[i,], col=rgb(red=0, green=0, blue=0, alpha=0.3), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$pCO2, col="red", lwd = 5)
# dat.rs = cbind(dat.age$age, 100)
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
lines(ages, post.clp$BUGSoutput$median$Tsoil, col="palegreen", lwd=2)
axis(1)
axis(2)

plot(ages, post.clp$BUGSoutput$sims.list$PPCQ[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("PPCQ")), xlim = range(ages), ylim = c(0, 500), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$PPCQ[i,], col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$PPCQ, col="deepskyblue2", lwd=2)
axis(1)
axis(2)

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

plot(ages, post.clp$BUGSoutput$sims.list$MAP[1,], type="l", axes = FALSE, xlab = "Age (Ma)", ylab = expression(paste("MAP")), xlim = range(ages), ylim = c(200, 800), col=rgb(red=0, green=0, blue=0, alpha=0.05), lwd=0.3)
for (i in 2:500) {
  lines(ages, post.clp$BUGSoutput$sims.list$MAP[i,], col=rgb(red=0, green=0, blue=0, alpha=0.2), lwd=0.3)
}
lines(ages, post.clp$BUGSoutput$median$MAP, col="deepskyblue2", lwd=2)
axis(1)
axis(2)



## plot w/ multi-section ----
# dt = 0.02
# ages = seq(-3, 0.1, by = dt)

load("out/ms_fx_1e3.rda")
fx = post.clp
load("out/ms_zjc_1e3.rda")
zjc = post.clp
load("out/ms_lc_1e3.rda")
lc = post.clp
lc.age = read.csv("data/data.csv") %>% filter(site == "Luochuan") %>% filter(age < 2.6)
fx.age = read.csv("data/data.csv") %>% filter(site == "Fuxian") %>% filter(age < 2.6)
zjc.age = read.csv("data/data.csv") %>% filter(site == "Zhaojiachuan") %>% filter(age < 2.6)

lc.co2 = data.frame(cbind(lc.age$age, 
                          t(apply(lc$BUGSoutput$sims.list$pCO2, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
lc.co2$site = "Luochuan"
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
