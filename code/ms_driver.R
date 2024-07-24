library(R2jags)
library(tidyverse)
source("code/constructors.R")
source("code/helpers.R")

## Read data
clp = read.csv("data/data.csv") %>% filter(site == "Zhaojiachuan") %>% filter(age < 2.6)
# clp = read.csv("data/data.csv") %>% filter(age < 2.6)
md = read.csv("data/d13Ca_tipple.csv") %>% filter(age < 2.6)
clp = na.exclude(clp[c("age", "d13C", "d18O", "d13Co")])
clp$d13Ca = approx(x = md$age, y = md$d13C, xout = clp$age)$y

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13C.stdev = 0.2
clp$d18O.stdev = 0.2
clp$d13Co.stdev = 0.2
clp$d13Ca.stdev = 0.2

## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13C", "d13C.stdev")])
d18Oc = na.exclude(clp[c("age", "d18O", "d18O.stdev")])
d13Co = na.exclude(clp[c("age", "d13Co", "d13Co.stdev")])
# D47c = na.exclude(clp[c("age", "D47", "D47.sd")])
d13Ca = na.exclude(clp[c("age", "d13C", "d13Ca.stdev")])

ages = ts(d13Cc$age, d18Oc$age, d13Co$age, d13Ca$age)
tsi = ages$ts_ind

d = list(ai = ages$ts, 
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = tsi[[1]],
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = tsi[[2]],
         # d13Co.obs = d13Co[, 2:3], d13Co.ai = tsi[[3]],
         # D47c.obs = D47c[, 2:3], D47c.ai = tsi[[4]],
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = tsi[[4]]
         )

parms = c("pCO2", "MAT", "PCQ_to", "Tsoil", "tsc", "MAP", "PCQ_pf", "PPCQ", "d18.p", 
          "d18O.s", "AET_PCQ", "z_m", "f_R", "S_z", "d13Cr", 
          "d13Cc", "d18Oc", "D47c", "d13Ca", "ha")

post.clp1 = jags.parallel(d, NULL, parms, "code/models/multi_sample.R",
                         n.iter = 1e4, n.chains = 3)

# sum1 = post.clp1$BUGSoutput$summary
save(post.clp1, file = "out/clp1e4_ms_zjc_nod13Co.rda")
# traceplot(post.clp1, varname = "pCO2")



