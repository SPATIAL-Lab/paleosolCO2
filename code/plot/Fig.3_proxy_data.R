rm(list = ls())
pacman::p_load(tidyverse, readxl)

# load and groom data ----
fuxian = read.csv("data/CLP_data/loess_glacial.csv") |>
  filter(section == "Fuxian") |>
  select(age, d13c, d18c, d13o, MS)
fuxian_D47 = read.csv("data/CLP_data/D47.csv") |>
  select(age, D47, D47.sd) |>
  mutate(D47_low = D47 - D47.sd,
         D47_high = D47 + D47.sd)

# plot ----
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
png("figure/Fig.3_fuxian_proxy_records.png", width = 4.5, height = 5.5, units = "in", res = 500)
par(mar = c(4, 4, 4, 4))
plot(-1, 0, xlim = c(0, 2.6), ylim = c(0, 5), axes = FALSE,
     xlab = "", ylab = "")

axis(3, mgp = c(2, 1.2, .5))
mtext("Age (Ma)", 3, line = 2.5)

yext = range(fuxian$d13c)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 1)
d13c.rs = cbind(fuxian$age,
                4 + (fuxian$d13c - min(tix)) / diff(range(tix)))
points(d13c.rs[, 1], d13c.rs[, 2], col = "black", bg = pal[1], pch = 21, cex = 1.1)
axis(4, 4 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[c]*" (\u2030)"), 4, line = 2.5, at = 4.5)

yext = range(fuxian$d18c)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 1)
d18c.rs = cbind(fuxian$age,
                3 + (fuxian$d18c - min(tix)) / diff(range(tix)))
points(d18c.rs[, 1], d18c.rs[,2], col = "black", bg = pal[3], pch = 21, cex = 1.1)
axis(2, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"C"[c]*" (\u2030)"), 2, line = 2.5, at = 3.5)

yext = range(fuxian$d13o)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 1)
d13o.rs = cbind(fuxian$age,
                2 + (fuxian$d13o - min(tix)) / diff(range(tix)))
points(d13o.rs[, 1], d13o.rs[, 2], col = "black", bg = pal[2], pch = 21, cex = 1.1)
axis(4, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"13"*"C"[o]*" (\u2030)"), 4, line = 2.5, at = 2.5)

yext = range(fuxian_D47$D47_low, fuxian_D47$D47_high)
tix = seq(floor(min(yext * 100)), ceiling(max(yext * 100)), by = 1) / 100
D47c.rs = cbind(fuxian_D47$age,
                1 + (fuxian_D47$D47 - min(tix)) / diff(range(tix)),
                1 + (fuxian_D47$D47_low - min(tix)) / diff(range(tix)),
                1 + (fuxian_D47$D47_high - min(tix)) / diff(range(tix)))
arrows(D47c.rs[, 1], D47c.rs[, 3], D47c.rs[, 1], D47c.rs[, 4], 
       col = "grey", angle = 90, length = 0, code = 0)
points(D47c.rs[, 1], D47c.rs[, 2], col = "black", bg = pal[4], pch = 21, cex = 1.3)
axis(2, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta[47]*" (\u2030)"), 2, line = 2.5, at = 1.5)

yext = range(fuxian$MS)
tix = seq(floor(min(yext)-1), ceiling(max(yext)+1), by = 10)
MS.rs = cbind(fuxian$age,
              0 + (fuxian$MS - min(tix)) / diff(range(tix)))
points(MS.rs[, 1], MS.rs[, 2], col = "black", bg = pal[1], pch = 21, cex = 1.1)
axis(4, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(italic(chi)[lf]* " (10"^"-8"*"m"^"3"*"kg"^"-1"*")"), 4, line = 2.5, at = 0.5)

axis(1)
mtext("Age (Ma)", 1, line = 2.2)
text(x = 0.1, y = 5, label = "a", cex = 1.3, font = 1)
text(x = 0.1, y = 4, label = "b", cex = 1.3, font = 1)
text(x = 0.1, y = 2.2, label = "c", cex = 1.3, font = 1)
text(x = 0.1, y = 1.8, label = "d", cex = 1.3, font = 1)
text(x = 0.1, y = .9, label = "e", cex = 1.3, font = 1)

dev.off()
