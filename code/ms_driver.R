library(R2jags)
library(tidyverse)
source("code/constructors.R")
source("code/helpers.R")

## Read data
clp = read.csv("data/data.csv") %>% filter(site == "Fuxian") %>% filter(age < 2.6)
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
         d13Co.obs = d13Co[, 2:3], d13Co.ai = tsi[[3]],
         # D47c.obs = D47c[, 2:3], D47c.ai = tsi[[4]],
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = tsi[[4]]
         )

parms = c("pCO2", "MAT", "Tsoil", "MAP", "PCQ_pf", "d18.p", 
          "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Co", "f_R", 
          "d13Cc", "d18Oc", "d13Ca", "ha", "ETR")

post.clp1 = jags.parallel(d, NULL, parms, "code/models/multi_sample.R",
                         n.iter = 1e3, n.chains = 3)

sum1 = post.clp1$BUGSoutput$summary
save(post.clp1, file = "out/clp1e3_ms.rda")
# traceplot(post.clp1, varname = "pCO2")

## plot ----
ages1 = ages$ts
clp.co2 = data.frame(cbind(ages1, 
                           t(apply(post.clp1$BUGSoutput$sims.list$pCO2, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.co2) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.co2, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.mat = data.frame(cbind(ages1,
                           t(apply(post.clp1$BUGSoutput$sims.list$MAT, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.mat) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.mat, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.st = data.frame(ages1,
                          t(apply(post.clp1$BUGSoutput$sims.list$Tsoil, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95))))
names(clp.st) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.st, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.map = data.frame(ages1,
                           t(apply(post.clp1$BUGSoutput$sims.list$MAP, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95))))
names(clp.map) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.map, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.ppcq = data.frame(cbind(ages1, 
                            t(apply(post.clp1$BUGSoutput$sims.list$PCQ_pf, 2, quantile, 
                                    c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.ppcq) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.ppcq, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))


clp.d18p = data.frame(ages1,
                            t(apply(post.clp1$BUGSoutput$sims.list$d18.p, 2, quantile, 
                                    c(0.05, 0.25, 0.5, 0.75, 0.95))))
names(clp.d18p) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.d18p, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.sz = data.frame(ages1,
                          t(apply(post.clp1$BUGSoutput$sims.list$S_z, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95))))
names(clp.sz) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.sz, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6)) +
  scale_y_continuous(limits = c(0, 5000))

clp.fr = data.frame(cbind(ages1,
                          t(apply(post.clp1$BUGSoutput$sims.list$f_R, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.fr) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.fr, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.ha = data.frame(cbind(ages1,
                          t(apply(post.clp1$BUGSoutput$sims.list$ha, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.ha) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.ha, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.d13Co = data.frame(cbind(ages1,
                          t(apply(post.clp1$BUGSoutput$sims.list$d13Co, 2, quantile, 
                                  c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.d13Co) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.d13Co, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

clp.d13Cc = data.frame(cbind(ages1,
                             t(apply(post.clp1$BUGSoutput$sims.list$d13Cc, 2, quantile, 
                                     c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(clp.d13Cc) = c("age", "x5", "x25", "median", "x75", "x95")
ggplot(clp.d13Cc, aes(x = age, y = median)) +
  geom_path() +
  scale_x_continuous(limits = c(0, 2.6))

