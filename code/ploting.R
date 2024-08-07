source("code/helpers.R")
load("bigout/post.Rda")

View(post$BUGSoutput$summary)

plot.jpi(ages, post$BUGSoutput$sims.list$pCO2)
points(ages, rep(min(post$BUGSoutput$sims.list$pCO2), 
                 length(ages)), col = si)

plot.jpi(ages, log(post$BUGSoutput$sims.list$S_z))
points(ages, rep(0, length(ages)), col = si)

plot.jpi(ages, post$BUGSoutput$sims.list$DIFC)

plot.jpi(ages, post$BUGSoutput$sims.list$d13Cr, ylim = c(-35, 0))
points(d13Cc$age, d13Cc$d13Cc, col = d13Cc$si)

plot.jpi(ages, post$BUGSoutput$sims.list$GMT)

plot.jpi(ages, post$BUGSoutput$sims.list$d13Ca)

plot.jpi(ages, post$BUGSoutput$sims.list$MAT[, 1,])

plot.jpi(ages, post$BUGSoutput$sims.list$MAP[, 1,])

plot.jpi(ages, post$BUGSoutput$sims.list$ha[, 1,])

plot.jpi(ages, post$BUGSoutput$sims.list$z_m)

plot(post$BUGSoutput$sims.list$ha[, 1,], post$BUGSoutput$sims.list$PPCQ[, 1,])
plot(post$BUGSoutput$sims.list$ha[, 1,], post$BUGSoutput$sims.list$MAP[, 1,])



pal = c("#0099FF", "#5Dc863")

png("out/DataFig.png", 8, 6, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(-70, -52), ylim = c(0, 7), axes = FALSE,
     xlab = "", ylab = "")

### d13Cc
tix = seq(floor(min(td$d13C, na.rm = TRUE)), 
          ceiling(max(td$d13C, na.rm = TRUE)), by = 2)
points(td$Age, 5 + (td$d13C - min(tix)) / diff(range(tix)),
       pch = 21, bg = pal[2])
axis(2, c(5 + (tix - min(tix)) / diff(range(tix))), labels = tix)
mtext(expression(delta^{13}*"C soil"), 2, at = 5.5, 
      line = 2.5)

### D47c
tix = seq(floor(min(td$D47, na.rm = TRUE) * 100), 
          ceiling(max(td$D47, na.rm = TRUE) * 100), by = 3) / 100
points(td$Age, 4 - (td$D47 - min(tix)) / diff(range(tix)),
       pch = 22, bg = pal[2])
axis(2, rev(c(3 + (tix - min(tix)) / diff(range(tix)))), labels = tix)
mtext(expression(Delta^{47}*" soil"), 2, at = 3.5, 
      line = 2.5)

### d18Oc
tix = seq(floor(min(td$d18O, na.rm = TRUE)), 
          ceiling(max(td$d18O, na.rm = TRUE)), by = 2)
points(td$Age, 1 + (td$d18O - min(tix)) / diff(range(tix)),
       pch = 23, bg = pal[2])
axis(2, c(1 + (tix - min(tix)) / diff(range(tix))), labels = tix)
mtext(expression(delta^{18}*"O soil"), 2, at = 1.5, 
      line = 2.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

## Parse data into series
d13Cc = na.exclude(td[c("Age", "d13C", "d13C.stdev")])
d18Oc = na.exclude(td[c("Age", "d18O", "d18O.stdev")])
D47c = na.exclude(td[c("Age", "D47", "D47.stderr")])

## Ages ----
dt = 0.25
ages = seq(-72, -52, by = dt)

#load("bigout/tsf4e3.rda")
load("bigout/tsd1e5_250.rda")
load("bigout/tsmd1e5_250.rda")

tsd.co2 = cbind(ages + dt / 2, 
                t(apply(post.tsd$BUGSoutput$sims.list$pCO2, 2, quantile, 
                c(0.05, 0.25, 0.5, 0.75, 0.95))))
tsd.MAT = cbind(ages + dt / 2, 
                t(apply(post.tsd$BUGSoutput$sims.list$MAT, 2, quantile, 
                        c(0.05, 0.25, 0.5, 0.75, 0.95))))

#png("out/Curves.png", 7, 6, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(-70, -52), ylim = c(0, 2), axes = FALSE,
     xlab = "", ylab = "")
yext = range(tsd.co2[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 5) * 100
tsd.co2rs = cbind(tsd.co2[, 1], 
                  1 + (tsd.co2[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsd.co2rs, "red")
axis(2, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("pCO"[2]*" (ppmv)"), 2, line = 2.5, at = 1.3)

yext = range(tsd.MAT[, 2:6])
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 5)
tsd.MATrs = cbind(tsd.MAT[, 1], 
                  (tsd.MAT[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsd.MATrs, "blue")
axis(4, (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("MAT"), 4, line = 2.5, at = 0.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

# TSMD environmental ----
tsmd.MAT = cbind(-(ages) - dt / 2, 
               t(apply(post.tsmd$BUGSoutput$sims.list$MAT, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))
tsmd.MAP = cbind(-(ages) - dt / 2, 
               t(apply(post.tsmd$BUGSoutput$sims.list$MAP, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))
tsmd.d13Ca = cbind(-(ages) - dt / 2, 
               t(apply(post.tsmd$BUGSoutput$sims.list$d13Ca, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))

png("out/Environmental.png", 7, 6, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(65, 55), ylim = c(0, 5), axes = FALSE,
     xlab = "", ylab = "")

yext = range(tsmd.co2[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 3) * 100
tsmd.co2rs = cbind(tsmd.co2[, 1], 
                  4 + (tsmd.co2[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsmd.co2rs, pal[1])
axis(2, 4 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("pCO"[2]), 2, 2.5, at = 4.5)

yext = range(tsmd.tempC[, 2:6])
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 2)
tsmd.tempCrs = cbind(tsmd.tempC[, 1], 
                   3 + (tsmd.tempC[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsmd.tempCrs, pal[1])
axis(4, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("SST", 4, 2.5, at = 3.5)

yext = range(tsmd.d13Ca[, 2:6])
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 3)
tsmd.d13Cars = cbind(tsmd.d13Ca[, 1], 
                     2 + (tsmd.d13Ca[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsmd.d13Cars, "#9dc3e6")
axis(2, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^{13}*"C"[atm]), 2, 2.5, at = 2.5)

yext = range(tsmd.MAT[, 2:6])
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 4)
tsmd.MATrs = cbind(tsmd.MAT[, 1], 
                     1 + (tsmd.MAT[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsmd.MATrs, pal[2])
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("MAT", 4, 2.5, at = 1.5)

yext = range(tsmd.MAP[, 2:6])
tix = seq(floor(min(yext)/100), ceiling(max(yext)/100), by = 2) * 100
tsmd.MAPrs = cbind(tsmd.MAP[, 1], 
                   (tsmd.MAP[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(tsmd.MAPrs, pal[2])
axis(2, (tix - min(tix)) / diff(range(tix)), tix)
mtext("MAP", 2, 2.5, at = 0.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

# Parameters ----
sens.pri = density(rgamma(1e6, 22.5, 5))
sens.post = density(post.tsmd$BUGSoutput$sims.list$sens)
oe.pri = density(rnorm(1e6, 2e4, sqrt(2e6)) * 1e-6)
oe.post = density(post.tsmd$BUGSoutput$sims.list$orgEff * 1e-6)

png("out/Parms.png", 9, 5, "in", res = 300)
layout(matrix(c(1, 2), nrow = 1))
par(mai = c(1, 0.2, 0.2, 0.2))
plot(sens.pri, xlim = range(sens.pri$x, sens.post$x),
     ylim = range(sens.pri$y, sens.post$y), main = "", axes = FALSE,
     xlab = "Site 1209 climate sensitivity", lty = 2, lwd = 2)
axis(1)
axis(2, labels = FALSE)
box()
lines(sens.post, lwd = 2)
plot(oe.pri, xlim = range(oe.pri$x, oe.post$x),
     ylim = range(oe.pri$y, oe.post$y), main = "", axes = FALSE,
     xlab = expression(delta^{13}*"C"[atm]*" sensitivity (ppt/ppm)"), lty = 2, lwd = 2)
axis(1)
axis(2, labels = FALSE)
box()
lines(oe.post, lwd = 2)
dev.off()

plot(ages, post.tsmd$BUGSoutput$mean$MAP, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$PPCQ, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$ha, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$MAT, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$TmPCQ, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$S_z, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$z_m, type = "l")
plot(ages, post.tsmd$BUGSoutput$mean$d13Ca, type = "l")

plot(ages, post.tsmd$BUGSoutput$mean$d18Of, type = "l")
points(d18Of$age, d18Of$d18O)
plot(ages, post.tsmd$BUGSoutput$mean$d11BGrub, type = "l")
points(-d11BGrub$age, d11BGrub$d11B)
plot(ages, post.tsmd$BUGSoutput$mean$d11BTsac, type = "l")
points(-d11BTsac$age, d11BTsac$d11B)
plot(ages, post.tsmd$BUGSoutput$mean$mgcaf, type = "l")
points(-mgcaf$age, mgcaf$MgCa)
plot(ages, post.tsmd$BUGSoutput$mean$d13Cf, type = "l")
points(-d13Cf$age, d13Cf$d13C)
plot(ages, post.tsmd$BUGSoutput$mean$d13Cc, type = "l")
points(-d13Cc$Age, d13Cc$d13C)
plot(ages, post.tsmd$BUGSoutput$mean$d18Oc, type = "l")
points(-d18Oc$Age, d18Oc$d18O)
plot(ages, post.tsmd$BUGSoutput$mean$D47c, type = "l", ylim = range(D47c$D47))
points(-D47c$Age, D47c$D47)


x = ts$ts[ts$ts_ind[[9]]]
points(x, rep(150, length(x)))
x = ts$ts[ts$ts_ind[[10]]]
points(x, rep(200, length(x)), col = "red")

plot.jpi(ts$ts, post.tsd$BUGSoutput$sims.list$pCO2, xlim = c(-65, -53))
plot.jpi(ts$ts, post.tsm$BUGSoutput$sims.list$pCO2, xlim = c(-60, -53))

plot.jpi(ts$ts, post.ts$BUGSoutput$sims.list$pCO2, xlim = c(-60, -53))
par(new = TRUE)
plot(d11BGrub$age, d11BGrub$d11B, xlim = c(-60, -53), ylim = c(14, 17),
       axes = FALSE, col = "blue")
points(d11BTsac$age, d11BTsac$d11B, col = "red")
