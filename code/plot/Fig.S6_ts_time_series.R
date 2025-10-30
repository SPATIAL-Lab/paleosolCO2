rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/helpers.R")

# load Bayesian data ----
load("out/ts_fuxian_8e5.rda")
post = post.ts
load("out/ts_fuxian_D47_8e5.rda")
post_47 = post.ts
load("out/ts_fuxian_MS_8e5.rda")
post_ms = post.ts
load("out/ts_fuxian_D47_MS_ECS_1e5.rda")
post_ecs = post.ts
cat("\014")
# mean(post_47_ms_ecs$BUGSoutput$sd$pCO2)

# time series w/ iterations ----
png("figure/Fig.S6_ts_time_series.png", width = 6.5, height = 7.2, units = "in", res = 500)
par(mfrow = c(3, 2), mar = margin(3, 3, 1, 1))
plot.jpi(seq(-3, -.05, .05), post$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         xlim = c(-3, 0), ylim = c(-5, 25),
         mgp = c(2, .8, 0))
text(-2.7, -2, "a", cex = 1.5, font = 2)
text(-2, 22, labels = "JPI-base_ts", cex = 1.2)

plot.jpi(seq(-3, -.05, .05), post_47$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         xlim = c(-3, 0), ylim = c(-5, 25),
         mgp = c(2, .8, 0))
text(-2.7, -2, "b", cex = 1.5, font = 2)
text(-2, 22, labels = expression("JPI-"*Delta[47]*"_ts"), cex = 1.2)

plot.jpi(seq(-3, -.05, .05), post_47$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         xlim = c(-3, 0), ylim = c(0, 3e3),
         mgp = c(2, .8, 0))
text(-2.7, 2700, "c", cex = 1.5, font = 2)
text(-2, 2500, labels = expression("JPI-"*Delta[47]*"_ts"), cex = 1.2)

plot.jpi(seq(-3, 0, .05), post_ms$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         xlim = c(-3, 0), ylim = c(0, 1e3),
         mgp = c(2, .8, 0))
text(-0.2, 100, "d", cex = 1.5, font = 2)
text(-1.5, 900, labels = "JPI-MS_ts", cex = 1.2)

plot.jpi(seq(-3, 0, .05), post_ms$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         xlim = c(-3, 0), ylim = c(-5, 25),
         mgp = c(2, .8, 0))
text(-2.7, -2, "e", cex = 1.5, font = 2)
text(-2, 22, labels = "JPI-MS_ts", cex = 1.2)

plot.jpi(seq(-3, 0, .05), post_ecs$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         xlim = c(-3, 0), ylim = c(-10, 30),
         mgp = c(2, .8, 0))
text(-2.7, -6, "f", cex = 1.5, font = 2)
text(-2, 25, labels = "JPI-CS_ts", cex = 1.2)
dev.off()
