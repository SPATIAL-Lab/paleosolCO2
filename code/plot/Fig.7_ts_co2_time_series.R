rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/helpers.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(hjust = 0.1, vjust = -10))
# parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf", "f_R", "spre")

# load data ----
# JPI posteriors
load("out/ts_fuxian_8e5.rda")
post = post.ts
load("out/ts_fuxian_D47_8e5.rda")
post_47 = post.ts
load("out/ts_fuxian_MS_8e5.rda")
post_ms = post.ts
load("out/ts_fuxian_CS_1e5.rda")
post_ecs = post.ts
# sample resolution
clp = read.csv("data/CLP_data/loess_glacial.csv") |> 
  filter(section == "Fuxian" & age > 0.05) |>
  select(age)
clp$index = 0
cat("\014")

ages = seq(-3, -0.05, .05)
png("figure/Fig.7_ts_co2_time_series.png", width = 6.4, height = 5.4, units = "in", res = 500)
par(mfrow = c(2, 2), mar = margin(3, 3, 1, 1))
plot.jpi(-ages, post$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         xlim = c(0, 3), ylim = c(-30, 1e3), mgp = c(2, .8, 0))
points(clp$age, clp$index)
text(2.9, 50, "a", cex = 1.3)
text(2.5, 900, "JPI-base", cex = 1)
plot.jpi(-ages, post_47$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         xlim = c(0, 3), ylim = c(-20, 1e3), mgp = c(2, .8, 0))
points(clp$age, clp$index)
text(2.9, 50, "b", cex = 1.3)
text(2.5, 900, expression("JPI-"*Delta[47]), cex = 1)
plot.jpi(-seq(-3, 0, .05), post_ms$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         xlim = c(0, 3), ylim = c(-20, 800), mgp = c(2, .8, 0))
points(clp$age, clp$index)
text(2.9, 50, "c", cex = 1.3)
text(2.5, 750, "JPI-MS", cex = 1)
plot.jpi(-seq(-3, -.1, .05), post_ecs$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         xlim = c(0, 3), ylim = c(-20, 550), mgp = c(2, .8, 0))
points(clp$age, clp$index)
text(2.9, 30, "d", cex = 1.3)
text(1, 500, "JPI-CS", cex = 1)
dev.off()



