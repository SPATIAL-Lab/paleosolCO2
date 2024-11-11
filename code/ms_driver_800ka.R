library(R2jags)
library(tidyverse)
library(readxl)
source("code/constructors.R")
source("code/helpers.R")

## Read data
clp = read.csv("data/clp_800ka.csv") %>% drop_na()
clp = clp[order(clp$age),]
pCO2 = read_xlsx("/Users/jiawei/Documents/Work/2019- glacial CO2/data/related_data/global records/CO2/ice core pCO2.xlsx")
pCO2 = pCO2[,c(1,3)]
clp$co2 = approx(x = pCO2$`Age (Ma)`, y = pCO2$`CO2 (ppmv)`, xout = clp$age)$y
md = read.csv("data/d13Ca_tipple.csv") %>% filter(age <= 0.8)
clp$d13a = approx(x = md$age, y = md$d13C, xout = clp$age)$y

## Propagate measurement, within outcrop, and within age bin uncertainties
clp$d13c.sd = 0.2
clp$d18c.sd = 0.2
clp$d13o.sd = 0.2
clp$d13a.sd = 0.2
clp$co2.sd = 20
## Parse data into series
d13Cc = na.exclude(clp[c("age", "d13c", "d13c.sd")])
d18Oc = na.exclude(clp[c("age", "d18c", "d18c.sd")])
d13Co = na.exclude(clp[c("age", "d13o", "d13o.sd")])
d13Ca = na.exclude(clp[c("age", "d13a", "d13a.sd")])
ice = na.exclude(clp[c("age", "co2", "co2.sd")])

ages = ts(d13Cc$age, d18Oc$age, d13Co$age, d13Ca$age, ice$age) # D47c$age, 
tsi = ages$ts_ind
ai = ages$ts

d = list(ai = ages$ts, 
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = tsi[[1]],
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = tsi[[2]],
         d13Co.obs = d13Co[, 2:3], d13Co.ai = tsi[[3]],
         # D47c.obs = D47c[, 2:3], D47c.ai = tsi[[4]],
         d13Ca.obs = d13Ca[, 2:3], d13Ca.ai = tsi[[4]],
         ice.obs = ice[, 2:3], ice.ai = tsi[[5]]
)

parms = c("pCO2", "MAT", "PCQ_to", "Tsoil", "tsc", "MAP", "PCQ_pf", "PPCQ", "d18.p", 
          "d18O.s", "AET_PCQ", "z_m", "f_R", "S_z", "d13Cr", 
          "d13Cc", "d18Oc", "D47c", "d13Ca", "ha", "Ratio")

system.time({post.clp = jags.parallel(d, NULL, parms, "code/models/multi_sample_800ka.R",
                                      n.iter = 1e5, n.chains = 3, n.burnin = 3e4)})

sum.clp = post.clp$BUGSoutput$summary
save(post.clp, file = "out/ms_800ka_1e5.rda")
load("out/ms_800ka_1e5.rda")
# traceplot(post.clp, varname = "MAP")

Sz = as.data.frame(post.clp$BUGSoutput$mean$S_z)
Sz.sd = as.data.frame(post.clp$BUGSoutput$sd$S_z)
d18c = as.data.frame(post.clp$BUGSoutput$mean$d18Oc)
d18c.sd = as.data.frame(post.clp$BUGSoutput$sd$d18Oc)
post = as.data.frame(cbind(d18c, d18c.sd, Sz, Sz.sd))
names(post) = c("d18c", "d18c.sd", "Sz", "Sz.sd")
post = cbind(clp[c("age", "section", "MS")], post)
# post = post %>% filter(section != "Luochuan")
m1 = nls(data = post, d18c ~ a * log(Sz) + b, start = list(a = -2, b = 0.4))
summary(m1)
d18c.p = seq(min(post$d18c), max(post$d18c), 0.01)
Sz.p = predict(m1, newdata = data.frame(d18c = d18c.p))
fit = data.frame(d18c.p, Sz.p)

ggplot(post, aes(x = d18c, y = Sz, fill = section)) +
  geom_errorbar(aes(xmin = d18c - d18c.sd, xmax = d18c + d18c.sd), size = 0.2) +
  geom_errorbar(aes(ymin = Sz - Sz.sd, ymax = Sz + Sz.sd), size = 0.2) +
  geom_point(size = 4, shape = 21) +
  # geom_line(data = fit, aes(x = d18c.p, y = Sz.p), size = 1, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(size = 10)) +
  labs(x = expression(delta^"18"*"O"[carb]*" (\u2030)"),
       y = expression(S[(z)]*" (ppmv)")) +
  scale_y_continuous(breaks = seq(0, 2500, 500)) +
  scale_x_continuous(breaks = seq(-11, -8.5, 0.5))


ggplot(post, aes(x = MS, y = Sz, fill = section)) +
  geom_errorbar(aes(ymin = Sz - Sz.sd, ymax = Sz + Sz.sd), size = 0.2) +
  geom_point(size = 4, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(size = 10)) +
  labs(x = expression(delta^"18"*"O"[carb]*" (\u2030)"),
       y = expression(S[(z)]*" (ppmv)")) +
  scale_y_continuous(breaks = seq(0, 2500, 500))
post.lc = post %>% filter(section == "Luochuan")
m1 = lm(data = post.lc, MS ~ Sz)
summary(m1)
