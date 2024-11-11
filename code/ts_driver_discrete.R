library(tidyverse)
library(R2jags)
library(readxl)
source("code/constructors.R")
source("code/helpers.R")

### loess ----
# Read data
# Glacial data
clp = read.csv("data/loess_glacial.csv") %>% filter(age < 2.6) %>% filter(section == "Fuxian")
clp = clp[,1:8]
D47 = read.csv("data/data.csv") %>% filter(site == "Fuxian") %>% drop_na(D47)
# # Interglacial data
# clp = read.csv("data/loess_interglacial.csv") %>% filter(age < 2.6)
# clp = clp[,1:9]

# randomly sampling of subset
# t.step = 0.05
# n = 2.6 / t.step + 1
# clp.sub = data.frame()
# for (i in 1:n) {
#   begin = t.step * (i-1)
#   end = t.step * i - 0.01
#   subset = clp %>% filter(age > begin & age < end)
#   sample = subset %>% slice_sample(n = 1)
#   clp.sub = rbind(clp.sub, sample)
# }
# clp = clp.sub

md = read.csv("data/d13Ca_tipple.csv") %>% filter(age <= max(clp$age) + 0.1) # Tipple et al. (2010)
clp$d13a = approx(md$age, md$d13C, xout = clp$age)$y
## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13c.stdev = 0.2
clp$d18c.stdev = 0.2
clp$d13o.stdev = 0.2
clp$d13a.stdev = 0.05
clp$MS.stdev = 10

## Make ages negative
clp$age = -clp$age
clp = clp %>% arrange(age)
D47$age = -D47$age
## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13c", "d13c.stdev")])
d18Oc = na.exclude(clp[c("age", "d18c", "d18c.stdev")])
d13Co = na.exclude(clp[c("age", "d13o", "d13o.stdev")])
D47c = na.exclude(D47[c("age", "D47", "D47.sd")])
d13Ca = na.exclude(clp[c("age", "d13a", "d13a.stdev")])
MS = na.exclude(clp[c("age", "MS", "MS.stdev")])

# ages = clp$age
# dt = abs(diff(ages, lag = 1))
dt = 0.1
ages = seq(-2.6, 0, by = dt)
d18Oc.ai = get.ind(d18Oc$age, ages)
d13Cc.ai = get.ind(d13Cc$age, ages)
d13Co.ai = get.ind(d13Co$age, ages)
D47c.ai = get.ind(D47c$age, ages)
MS.ai = get.ind(MS$age, ages)
d13Ca.ai = get.ind(d13Ca$age, ages)

d = list(ai = ages, dt = dt,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = d18Oc.ai,
         d18Oc.obs2 = d18Oc[, 2:3], d18Oc.ai2 = d18Oc.ai,
         d13Co.obs = d13Co[, 2:3], d13Co.ai = d13Co.ai,
         # D47c.obs = D47c[, 2:3], D47c.ai = D47c.ai,
         # MS.obs = MS[, 2:3], MS.ai = MS.ai,
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = d13Ca.ai)

parms = c("pCO2", "MAT", "PCQ_to", "Tsoil", "tsc", "MAP", "PCQ_pf", "PPCQ", "d18.p", 
          "d18O.s", "AET_PCQ", "z_m", "f_R", "S_z", "Ratio", "ETR", "d13Cs", "d13Cc", "d13Cc")

system.time({post.clp = jags.parallel(d, NULL, parms, "code/models/time_series_discrete.R", 
                        n.iter = 1e5, n.chains = 3, n.burnin = 3e4)})
sum_clp = post.clp$BUGSoutput$summary

save(post.clp, file = "out/ts_lc_1e5_30ppm_normal_MS.rda")
load("out/ts_fx_1e5_30ppm.rda")
# traceplot(post.clp, varname = "Ratio")