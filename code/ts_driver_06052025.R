rm(list = ls())
pacman::p_load(R2jags, readxl, tidyverse)
source("code/constructors.R")
source("code/helpers.R")

### load and groom data ----
clp = read_csv("data/CLP_data/loess_glacial.csv") |> 
  filter(section == "Fuxian")
D47 = read_csv("data/CLP_data/D47.csv")
d13a = read_csv("data/global_data/d13Ca_tipple.csv") # Tipple et al. (2010)
ice_sheet = read_xlsx("data/global_data/ice_forcing_stap_2018.xlsx", sheet = 2) # Stap et al. (2018)
ice_sheet = ice_sheet[2:nrow(ice_sheet),1:2]
names(ice_sheet) = c("age", "ice_forcing")
ice_sheet = ice_sheet |>
  mutate(across(everything(), as.numeric))

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13c.sd = 0.2
clp$d18c.sd = 0.2
clp$d13o.sd = 0.2
clp$MS.sd = 10

## Make ages negative
clp$age = -clp$age
D47$age = -D47$age
d13a$age = -d13a$age
ice_sheet$age = -ice_sheet$age
## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13c", "d13c.sd")])
d18Oc = na.exclude(clp[c("age", "d18c", "d18c.sd")])
d13Co = na.exclude(clp[c("age", "d13o", "d13o.sd")])
D47c = na.exclude(D47[c("age", "D47", "D47.sd")])
MS = na.exclude(clp[c("age", "MS", "MS.sd")])

# ages = clp$age
dt = 0.05
ages = seq(-3.0, -0.05, by = dt)
d18Oc.ai = get.ind(d18Oc$age, ages)
d13Cc.ai = get.ind(d13Cc$age, ages)
d13Co.ai = get.ind(d13Co$age, ages)
D47c.ai = get.ind(D47c$age, ages)
MS.ai = get.ind(MS$age, ages)
d13Ca = approx(x = d13a$age, y = d13a$d13C, xout = ages)$y
R_LI = approx(x = ice_sheet$age, y = ice_sheet$ice_forcing, xout = ages)$y

d = list(ai = ages, dt = dt, d13Ca = d13Ca, R_ice = R_LI,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = d18Oc.ai,
         d13Co.obs = d13Co[, 2:3], d13Co.ai = d13Co.ai,
         D47c.obs = D47c[, 2:3], D47c.ai = D47c.ai,
         MS.obs = MS[, 2:3], MS.ai = MS.ai)

parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf",
          "Tsoil", "S_z", "f_R", "spre", "pore", "temp_diff")

system.time({post.ts = jags.parallel(d, NULL, parms, "code/models/time_series_06052025.R", 
                        n.iter = 1e5, n.chains = 3, n.burnin = 1e4)})

View(post.ts$BUGSoutput$summary)
save(post.ts, file = "out/ts_fuxian_D47_MS_ECS_1e5.rda")
for (i in 1:length(parms)) {
  name = parms[i]
  plot.jpi(ai, post.ts$BUGSoutput$sims.list[[name]], n = 5e2, 
           xlab = "Age (Ma)", ylab = name, mgp = c(2, .8, 0))
}
