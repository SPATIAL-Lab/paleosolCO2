rm(list = ls())
source('code/helpers.R')
parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf",
          "Tsoil", "S_z", "f_R", "spre", "pore", "temp_diff")

load("out/ms_fuxian_D47_MS_1e5.rda")
load("out/ts_fuxian_D47_MS_1e4.rda")

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
for (i in 1:length(parms)) {
  plot_ts_ms(i)
}

