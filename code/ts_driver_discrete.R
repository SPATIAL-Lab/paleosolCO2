library(R2jags)
source("code/constructors.R")
source("code/helpers.R")

## Read data
td = read.csv("data/BBNP_data.csv")
md = read.csv("data/Shatsky_data.csv")

## Propagate measurement, within outcrop, and within age bin uncertainties
td$d18O.stdev = sqrt(td$d18O.stdev ^ 2 + 0.4 ^ 2 + 1 ^ 1)
td$d13C.stdev = sqrt(td$d13C.stdev ^ 2 + 0.4 ^ 2 + 1 ^ 1)
md$d13C.stdev = rep(1)

## Remove PETM
md = md[md$age <= 55.741 | md$age >= 55.942,]

## Make ages negative
td$Age = -td$Age
md$age = -md$age

## Parse data into series
d13Cc = na.exclude(td[c("Age", "d13C", "d13C.stdev")])
d18Oc = na.exclude(td[c("Age", "d18O", "d18O.stdev")])
D47c = na.exclude(td[c("Age", "D47", "D47.stderr")])
d13Ca = na.exclude(md[c("age", "d13C", "d13C.stdev")])
d13Ca$d13C = d13Ca$d13C - 8

dt = 0.25
ages = seq(-72, -52, by = dt)
d18Oc.ai = get.ind(d18Oc$Age, ages)
d13Cc.ai = get.ind(d13Cc$Age, ages)
D47c.ai = get.ind(D47c$Age, ages)
d13Ca.ai = get.ind(d13Ca$age, ages)

d = list(ai = ages, dt = dt,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = d18Oc.ai,
         D47c.obs = D47c[, 2:3], D47c.ai = D47c.ai,
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = d13Ca.ai)

parms = c("pCO2", "MAT", "MAP", "TmPCQ", "PPCQ", "d18.p", 
          "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Cr", 
          "d13Cc", "d18Oc", "D47c", "d13Ca", "ha")

system.time({post.tsd = jags.parallel(d, NULL, parms, "code/models/time_series_discrete.R", 
                        n.iter = 1e4, n.chains = 3)})

