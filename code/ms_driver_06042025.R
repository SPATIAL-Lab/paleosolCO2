rm(list = ls())
pacman::p_load(R2jags, readxl)
source("code/constructors.R")
source("code/helpers.R")

## Read and groom data ----
clp = read.csv("data/CLP_data/loess_glacial.csv") |> 
  filter(section == "Fuxian") |>
  select(age, d13c, d18c, d13o, MS)
D47 = read.csv("data/CLP_data/D47.csv") |>
  select(age, D47, D47.sd)
clp$D47 = approx(x = D47$age, y = D47$D47, xout = clp$age)$y
clp$D47_se = approx(x = D47$age, y = D47$D47.sd, xout = clp$age)$y

d13a = read.csv("data/global_data/d13Ca_tipple.csv") |> 
  filter(age <= 2.6)
clp$d13a = approx(x = d13a$age, y = d13a$d13C, xout = clp$age)$y
R_ice = read.csv("data/global_data/ice_forcing.csv")
clp$R_ice = approx(x = R_ice$age, y = R_ice$ice_forcing, xout = clp$age)$y

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13c_stdev = 0.2
clp$d18c_stdev = 0.2
clp$d13o_stdev = 0.2
clp$MS_stdev = 10

## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13c", "d13c_stdev")])
d18Oc = na.exclude(clp[c("age", "d18c", "d18c_stdev")])
d13Co = na.exclude(clp[c("age", "d13o", "d13o_stdev")])
D47c = na.exclude(clp[c("age", "D47", "D47_se")])
MS = na.exclude(clp[c("age", "MS", "MS_stdev")])

ages = ts(d13Cc$age, d18Oc$age, d13Co$age
          , D47c$age
          , MS$age
          )
tsi = ages$ts_ind
ai = ages$ts

## MCMC ----
d = list(ai = ages$ts, d13Ca = clp$d13a, R_ice = clp$R_ice,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = tsi[[1]],
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = tsi[[2]],
         d13Co.obs = d13Co[, 2:3], d13Co.ai = tsi[[3]]
         , D47c.obs = D47c[, 2:3], D47c.ai = tsi[[4]]
         , MS.obs = MS[, 2:3], MS.ai = tsi[[5]]
)

parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf",
          "Tsoil", "S_z", "f_R", "spre", "pore", "temp_diff")

system.time({post.ms = jags.parallel(d, NULL, parms, "code/models/multi_sample_06042025.R",
                                      n.iter = 1e5, n.chains = 3, n.burnin = 1e4)})

View(post.ms$BUGSoutput$summary)
save(post.ms, file = "out/ms_fuxian_D47_MS_ECS_1e5.rda")

# load("out/ms_fuxian_1e5.rda")
post_data = data.frame(age = ai)
for (i in 1:length(parms)) {
  name = parms[i]
  name_sd = paste0(parms[i], "_sd")
  name_eff = paste0(parms[i], "_eff")
  for (p in 1:nrow(post_data)) {
    post_data[[name]][p] = mean(post.ms$BUGSoutput$sims.list[[name]][, p])
    post_data[[name_sd]][p] = sd(post.ms$BUGSoutput$sims.list[[name]][, p])
    index = nrow(post_data) * (i - 1) + p
    rhat = post.ms$BUGSoutput$summary[index, "Rhat"]
    n_eff = post.ms$BUGSoutput$summary[index, "n.eff"]
    if (rhat < 1.01) {
      post_data[[name_eff]][p] = "positive"
    } else {
      post_data[[name_eff]][p] = "negative"
    }
  }
}
write.csv(post_data, file = "out/ms_fuxian_D47_MS_ECS_1e5.csv")

for (i in 1:length(parms)) {
  name = parms[i]
  plot.jpi(ai, post.ms$BUGSoutput$sims.list[[name]], n = 100, ylab = name)
}


