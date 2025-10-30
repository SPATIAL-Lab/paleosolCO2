rm(list = ls())
library(tidyverse)
source('code/helpers.R')
parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf",
          "Tsoil", "S_z", "f_R", "spre", "pore", "temp_diff")
# load data ----
# CS versions 
load("out/ms_fuxian_D47_MS_ECS_2e5.rda")
post_ms_cs = post.ms
load("out/ts_fuxian_CS_1e5.rda")
post_ts_cs = post.ts

# MS versions
load("out/ms_fuxian_D47_MS_2e5.rda")
load("out/ts_fuxian_MS_8e5.rda")

clp = read.csv("data/CLP_data/loess_glacial.csv") |> 
  filter(section == "Fuxian") |>
  select(age)
ms_age = clp$age

ts_age = seq(-3.0, 0, by = .05)

plot_ts_ms = function(n){
  par(mfrow = c(1, 2), mar = margin(3, 3, 1, 1))
  plot.jpi(ms_age, post.ms$BUGSoutput$sims.list[[parms[n]]], n = 100, xlab = "Age", ylab = parms[n],
           mgp = c(2,1,0))
  plot.jpi(-ts_age, post.ts$BUGSoutput$sims.list[[parms[n]]], n = 5e2, 
           xlab = "Age (Ma)", ylab = parms[n], mgp = c(2, .8, 0))
}

# plot
png("figure/Fig.S5_ms_ts_comparison.png", width = 6, height = 9.6, units = "in", res = 500)
par(mfrow = c(4, 2), mar = margin(3, 3, 1, 1))
plot.jpi(ms_age, post.ms$BUGSoutput$sims.list$MAP, n = 100, 
         xlab = "Age (Ma)", ylim = c(100, 900), ylab = "MAP (mm)",
         xlim = c(0, 2.6), mgp = c(2, .8, 0))
text(2.5, 850, "a", cex = 1.3)
plot.jpi(-ts_age, post.ts$BUGSoutput$sims.list$MAP, n = 200, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         xlim = c(0, 2.6), ylim = c(100, 900), mgp = c(2, .8, 0))
text(0.1, 850, "b", cex = 1.3)
plot.jpi(ms_age, post.ms$BUGSoutput$sims.list$PCQ_pf, n = 100, 
         xlab = "Age (Ma)", ylab = expression(italic(f)[PPCQ]),
         xlim = c(0, 2.6), ylim = c(0, 1), mgp = c(2, .8, 0))
text(2.5, 0.1, "c", cex = 1.3)
plot.jpi(-ts_age, post.ts$BUGSoutput$sims.list$PCQ_pf, n = 200, 
         xlab = "Age (Ma)", ylab = expression(italic(f)[PPCQ]),
         xlim = c(0, 2.6), ylim = c(0, 1), mgp = c(2, .8, 0))
text(2.5, 0.1, "d", cex = 1.3)
plot.jpi(ms_age, post.ms$BUGSoutput$sims.list$S_z, n = 100, 
         xlab = "Age (Ma)", ylab = "S(z) (ppmv)",
         xlim = c(0, 2.6), ylim = c(0, 2e3), mgp = c(2, .8, 0))
text(0.1, 100, "e", cex = 1.3)
plot.jpi(-ts_age, post.ts$BUGSoutput$sims.list$S_z, n = 200, 
         xlab = "Age (Ma)", ylab = "S(z) (ppmv)",
         xlim = c(0, 2.6), ylim = c(0, 2e3), mgp = c(2, .8, 0))
text(0.1, 100, "f", cex = 1.3)

plot.jpi(ms_age, post_ms_cs$BUGSoutput$sims.list$ECS / (5.35*0.64*log(2)), n = 300, 
         xlab = "Age (Ma)", ylab = expression("S"["CO2,LI"]*" (K/W/m"^"2"*")"),
         ylim = c(.8, 3.5),
         mgp = c(2, .8, 0))
text(0.2, 1, "g", cex = 1.3)
text(2, 3.2, "JPI-CS", cex = 1.3)

plot.jpi(-seq(-3, -.1, .05), 0.32 + post_ts_cs$BUGSoutput$sims.list$ECS / (5.35*0.64*log(2)), n = 500, 
         xlab = "Age (Ma)", ylab = expression("S"["CO2,LI"]*" (K/W/m"^"2"*")"),
         ylim = c(.8, 3.5),
         mgp = c(2, .8, 0))
text(2.7, 1, "h", cex = 1.2, font = 2)
text(2, 3.2, "JPI-CS_ts", cex = 1.3)

dev.off()

