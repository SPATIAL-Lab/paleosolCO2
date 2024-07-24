library(R2jags)
source("code/constructors.R")
source("code/helpers.R")

## Read data
# clp = read.csv("data/data.csv") %>% filter(age < 2.6)
clp = read.csv("data/data.csv") %>% filter(site == "Zhaojiachuan") %>% filter(age < 2.6)
md = read.csv("data/d13Ca_tipple.csv") %>% filter(age < 2.6)

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13C.stdev = 0.2
clp$d18O.stdev = 0.2
clp$d13Co.stdev = 0.2
md$d13Ca.stdev = 0.2

## Make ages negative
clp$age = -clp$age
md$age = -md$age

## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13C", "d13C.stdev")])
d18Oc = na.exclude(clp[c("age", "d18O", "d18O.stdev")])
d13Co = na.exclude(clp[c("age", "d13Co", "d13Co.stdev")])
D47c = na.exclude(clp[c("age", "D47", "D47.sd")])
d13Ca = na.exclude(md[c("age", "d13C", "d13Ca.stdev")])

dt = 0.02
ages = seq(-3, 0.1, by = dt)
d18Oc.ai = get.ind(d18Oc$age, ages)
d13Cc.ai = get.ind(d13Cc$age, ages)
d13Co.ai = get.ind(d13Co$age, ages)
D47c.ai = get.ind(D47c$age, ages)
d13Ca.ai = get.ind(d13Ca$age, ages)

d = list(ai = ages, dt = dt,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = d18Oc.ai,
         d13Co.obs = d13Co[, 2:3], d13Co.ai = d13Co.ai,
         D47c.obs = D47c[, 2:3], D47c.ai = D47c.ai,
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = d13Ca.ai)

parms = c("pCO2", "MAT", "PCQ_to", "Tsoil", "tsc", "MAP", "PCQ_pf", "PPCQ", "d18.p", 
          "d18O.s", "AET_PCQ", "z_m", "f_R", "S_z", "d13Cr", 
          "d13Cc", "d18Oc", "D47c", "d13Ca", "ha")

system.time({post.clp = jags.parallel(d, NULL, parms, "code/models/time_series_discrete.R", 
                        n.iter = 1e4, n.chains = 3)})

sum_clp = post.clp$BUGSoutput$summary
save(post.clp, file = "out/clp1e4_ts_zjc.rda")
# traceplot(post.clp, varname = "pCO2")
# traceplot(post.clp, varname = "MAT")
# traceplot(post.clp, varname = "MAP")
