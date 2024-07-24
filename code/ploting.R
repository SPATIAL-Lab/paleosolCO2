source("code/helpers.R")
load("out/clp1e4_ms.rda")

# View(post$BUGSoutput$summary)

## plot.ms ----
# for output without age model
post.clp = post.clp1
ages = ages$ts

clp.co2 = data.frame(cbind(ages, 
                           t(apply(post.clp$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.co2) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.co2, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "ivory") +
  geom_ribbon(aes(ymin = x25, ymax = x75), fill = "ivory2") +
  geom_path(linewidth = 1)

clp.mat = data.frame(cbind(ages, 
                           t(apply(post.clp$BUGSoutput$sims.list$MAT, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.mat) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.mat, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.st = data.frame(cbind(ages, 
                          t(apply(post.clp$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.st) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.st, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.map = data.frame(cbind(ages, 
                           t(apply(post.clp$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.map) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.map, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.pcqpf = data.frame(cbind(ages2, 
                            t(apply(post.clp$BUGSoutput$sims.list$PCQ_pf, 2, quantile, 
                                    c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.ppcq) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.ppcq, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.d18p = data.frame(cbind(ages, 
                            t(apply(post.clp$BUGSoutput$sims.list$d18.p, 2, quantile, 
                                    c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.d18p) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.d18p, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.sz = data.frame(cbind(ages, 
                          t(apply(post.clp$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.sz) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.sz, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.fr = data.frame(cbind(ages, 
                          t(apply(post.clp$BUGSoutput$sims.list$f_R, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.fr) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.fr, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.d13Cc = data.frame(cbind(ages,
                             t(apply(post.clp1$BUGSoutput$sims.list$d13Cc, 2, quantile, 
                                     c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.d13Cc) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.d13Cc, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

clp.ha = data.frame(cbind(ages,
                             t(apply(post.clp1$BUGSoutput$sims.list$ha, 2, quantile, 
                                     c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.ha) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.ha, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), fill = "grey") +
  geom_path()

## plot.jpg ---- 

post = post.clp
plot.jpi(ages, post$BUGSoutput$sims.list$pCO2, ylim = c(0, 500))
# points(ages, rep(min(post$BUGSoutput$sims.list$pCO2), 
#                  length(ages)))

plot.jpi(ages, post$BUGSoutput$sims.list$S_z)
plot.jpi(ages, log10(post$BUGSoutput$sims.list$S_z))
# points(ages, rep(0, length(ages)))
plot.jpi(ages, post$BUGSoutput$sims.list$f_R)

# plot.jpi(ages, post$BUGSoutput$sims.list$DIFC)

plot.jpi(ages, post$BUGSoutput$sims.list$d18.p)
# points(d13Cc$age, d13Cc$d13Cc)

# plot.jpi(ages, post$BUGSoutput$sims.list$GMT)

plot.jpi(ages, post$BUGSoutput$sims.list$d13Ca)

plot.jpi(ages, post$BUGSoutput$sims.list$MAT)

plot.jpi(ages, post$BUGSoutput$sims.list$tsc)

plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_to)

plot.jpi(ages, post$BUGSoutput$sims.list$Tsoil, n = 500)

plot.jpi(ages, post$BUGSoutput$sims.list$MAP)

plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_pf)

plot.jpi(ages, post$BUGSoutput$sims.list$ha)

plot.jpi(ages, post$BUGSoutput$sims.list$z_m)

plot(post$BUGSoutput$sims.list$ha, post$BUGSoutput$sims.list$PPCQ)
plot(post$BUGSoutput$sims.list$ha, post$BUGSoutput$sims.list$MAP)



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
dt = 0.02
ages = seq(-3, 0.1, by = dt)

#load("bigout/tsf4e3.rda")
load("out/clp1e4_ms.rda")

clp.co2 = cbind(ages + dt / 2, 
                t(apply(post.clp$BUGSoutput$sims.list$pCO2, 2, quantile, 
                c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.MAT = cbind(ages + dt / 2, 
                t(apply(post.clp$BUGSoutput$sims.list$MAT, 2, quantile, 
                        c(0.05, 0.25, 0.5, 0.75, 0.95))))

#png("out/Curves.png", 7, 6, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(-3, 0), ylim = c(0, 2), axes = FALSE,
     xlab = "", ylab = "")
yext = range(clp.co2[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 1) * 100
clp.co2rs = cbind(clp.co2[, 1], 
                  1 + (clp.co2[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.co2rs, "red")
axis(2, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("pCO"[2]*" (ppmv)"), 2, line = 2.5, at = 1.5)

yext = range(clp.MAT[, 2:6])
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 2)
clp.MATrs = cbind(clp.MAT[, 1], 
                  (clp.MAT[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.MATrs, "blue")
axis(4, (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("MAT"), 4, line = 2.5, at = 0.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

# clp environmental ----
clp.co2 = cbind(-(ages) - dt / 2, 
                t(apply(post.clp$BUGSoutput$sims.list$pCO2, 2, quantile, 
                        c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.MAT = cbind(-(ages) - dt / 2, 
               t(apply(post.clp$BUGSoutput$sims.list$MAT, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.st = cbind(-(ages) - dt / 2, 
                t(apply(post.clp$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                        c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.MAP = cbind(-(ages) - dt / 2, 
               t(apply(post.clp$BUGSoutput$sims.list$MAP, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.sz = cbind(-(ages) - dt / 2, 
               t(apply(post.clp$BUGSoutput$sims.list$S_z, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))
clp.ppcq = cbind(-(ages) - dt / 2, 
               t(apply(post.clp$BUGSoutput$sims.list$PPCQ, 2, quantile, 
                       c(0.05, 0.25, 0.5, 0.75, 0.95))))

png("out/Environmental_ms.png", width = 5, height = 10, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(3, 0), ylim = c(0, 6), axes = FALSE,
     xlab = "", ylab = "")

yext = range(clp.co2[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 1) * 100
clp.co2rs = cbind(clp.co2[, 1], 
                  5 + (clp.co2[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.co2rs, pal[1])
axis(2, 5 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("pCO"[2]), 2, 2.5, at = 5.5)

yext = range(clp.MAT[, 2:6])
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 2)
clp.MATrs = cbind(clp.MAT[, 1], 
                   4 + (clp.MAT[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.MATrs, pal[1])
axis(4, 4 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("MAT", 4, 2.5, at = 4.5)

yext = range(clp.st[, 2:6])
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 2)
clp.strs = cbind(clp.st[, 1], 
                     3 + (clp.st[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.strs, pal[1])
axis(2, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("T"[soil]), 2, 2.5, at = 3.5)

yext = range(clp.MAP[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 2) * 100
clp.MAPrs = cbind(clp.MAP[, 1], 
                     2 + (clp.MAP[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.MAPrs, pal[1])
axis(4, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("MAP", 4, 2.5, at = 2.5)

yext = range(clp.ppcq[, 2:6])
tix = seq(floor(min(yext) / 100), 
          ceiling(max(yext) / 100), by = 2) * 100
clp.ppcqrs = cbind(clp.ppcq[, 1], 
                  1 + (clp.ppcq[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.ppcqrs, pal[1])
axis(2, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("PPCQ", 2, 2.5, at = 1.5)

yext = range(clp.sz[, 2:6])
tix = seq(floor(min(yext) / 1e4), 
          ceiling(max(yext) / 1e4), by = 1) * 1e4
clp.szrs = cbind(clp.sz[, 1], 
                   (clp.sz[, 2:6] - min(tix)) / diff(range(tix)))
tsdens(clp.szrs, pal[1])
axis(4, (tix - min(tix)) / diff(range(tix)), tix)
mtext("S(z)", 4, 2.5, at = 0.5)

### X axis
axis(1)
mtext("Age (Ma)", 1, line = 2)

dev.off()

# Parameters ----
co2.pri = density(runif(1e6, 150, 400))
co2.post = density(post.clp$BUGSoutput$sims.list$pCO2)
oe.pri = density(rnorm(1e6, 2e4, sqrt(2e6)) * 1e-6)
oe.post = density(post.clp$BUGSoutput$sims.list$orgEff * 1e-6)

png("out/Parms.png", 9, 5, "in", res = 300)
layout(matrix(c(1, 2), nrow = 1))
par(mai = c(1, 0.2, 0.2, 0.2))
plot(co2.pri, xlim = range(co2.pri$x, co2.post$x),
     ylim = range(co2.pri$y, co2.post$y), main = "", axes = FALSE,
     xlab = expression("CO"[2]), lty = 2, lwd = 2)
axis(1)
axis(2, labels = FALSE)
box()
lines(co2.post, lwd = 2)
plot(oe.pri, xlim = range(oe.pri$x, oe.post$x),
     ylim = range(oe.pri$y, oe.post$y), main = "", axes = FALSE,
     xlab = expression(delta^{13}*"C"[atm]*" co2itivity (ppt/ppm)"), lty = 2, lwd = 2)
axis(1)
axis(2, labels = FALSE)
box()
lines(oe.post, lwd = 2)
dev.off()

plot(ages, post.clp$BUGSoutput$mean$MAP, type = "l")
plot(ages, post.clp$BUGSoutput$mean$PPCQ, type = "l")
plot(ages, post.clp$BUGSoutput$mean$ha, type = "l")
plot(ages, post.clp$BUGSoutput$mean$MAT, type = "l")
plot(ages, post.clp$BUGSoutput$mean$Tsoil, type = "l")
plot(ages, post.clp$BUGSoutput$mean$S_z, type = "l")
plot(ages, post.clp$BUGSoutput$mean$z_m, type = "l")
plot(ages, post.clp$BUGSoutput$mean$d13Ca, type = "l")

plot(ages, post.clp$BUGSoutput$mean$d13Cc, type = "l")
points(-d13Cc$Age, d13Cc$d13C)
plot(ages, post.clp$BUGSoutput$mean$d18Oc, type = "l")
points(-d18Oc$Age, d18Oc$d18O)
plot(ages, post.clp$BUGSoutput$mean$D47c, type = "l", ylim = range(D47c$D47))
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
