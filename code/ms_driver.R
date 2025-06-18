rm(list = ls())
pacman::p_load(rjags, R2jags, tidyverse, readxl)
source("code/constructors.R")
source("code/helpers.R")

## Read and groom data ----
clp = read.csv("data/loess_glacial.csv") %>% filter(section == "Zhaojiachuan")
clp = clp[,1:8]
# clp = read.csv("data/loess_interglacial.csv") %>% filter(age < 2.6)
# clp = clp[,1:9]
# clp = read.csv("data/D47.csv")
# clp = clp[order(clp$age),]

md = read.csv("data/d13Ca_tipple.csv") %>% filter(age <= max(clp$age) + 0.1)
clp$d13a = approx(x = md$age, y = md$d13C, xout = clp$age)$y

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13c.stdev = 0.2
clp$d18c.stdev = 0.2
clp$d13o.stdev = 0.2
clp$d13a.stdev = 0.2
# clp$MS.stdev = 10

plot(clp$age, clp$d13o)

## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13c", "d13c.stdev")])
d18Oc = na.exclude(clp[c("age", "d18c", "d18c.stdev")])
d13Co = na.exclude(clp[c("age", "d13o", "d13o.stdev")])
d13Ca = na.exclude(clp[c("age", "d13a", "d13a.stdev")])
# D47c = na.exclude(clp[c("age", "D47", "D47.sd")])
# MS = na.exclude(clp[c("age", "MS", "MS.stdev")])

ages = ts(d13Cc$age, d18Oc$age, d13Co$age, d13Ca$age) # , D47c$age , MS$age
tsi = ages$ts_ind
ai = ages$ts

## MCMC ----
d = list(ai = ages$ts, 
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = tsi[[1]],
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = tsi[[2]],
         # d18Oc.obs2 = d18Oc[, 2:3], d18Oc.ai2 = tsi[[2]],
         d13Co.obs = d13Co[, 2:3], d13Co.ai = tsi[[3]],
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = tsi[[4]]
         # , D47c.obs = D47c[, 2:3], D47c.ai = tsi[[5]]
         # , MS.obs = MS[, 2:3], MS.ai = tsi[[5]]
)

parms = c("pCO2", "MAT", "PCQ_to", "Tsoil", "tsc", "MAP", "PCQ_pf", "PPCQ", "d18p", 
          "d18O.s", "pore", "z_m", "f_R", "S_z", "R_PCQ_S_0", "DIFC", "L", "AI", "Ratio")

system.time({post.clp = jags.parallel(d, NULL, parms, "code/models/multi_sample_v2.R",
                                      n.iter = 1e4, n.chains = 3, n.burnin = 3e3)})

sum.clp = post.clp$BUGSoutput$summary
save(post.clp, file = "out/ms_zjc_1e4.rda")
load("out/ms_zjc_1e4_d18c.rda")
# traceplot(post.clp, varname = "MAP")
