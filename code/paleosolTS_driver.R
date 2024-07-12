library(R2jags)
library(openxlsx)
source("code/constructors.R")
source("code/helpers.R")

## Read data
dfs = list.files("data/TJ/", full.names = TRUE)
tds = list()
for(i in seq_along(dfs)){
  tds[[i]] = read.xlsx(dfs[i], sheet = "paleosol data", startRow = 3)
  tds[[i]]$si = rep(i)
}

td = tds[[1]]
for(i in 2:length(tds)){
  names(tds[[i]]) = names(td)
  td = rbind(td, tds[[i]])
}

td = td[td$`Age.(Ma)` < 251,]
td = td[td$`Age.(Ma)` > 140,]
td$`Age.(Ma)` = 0 - td$`Age.(Ma)`

## Parse data into series
d13Cc = na.exclude(td[c("Age.(Ma)", "d13Ccc.(‰.PDB)", "si")])
d13Cc$d13Cc.sd = rep(0.5)
names(d13Cc) = c("age", "d13Cc", "si", "d13Cc.sd")

d18Oc = na.exclude(td[c("Age.(Ma)", "d18Occ.(‰.PDB)", "si")])
d18Oc$d18Oc.sd = rep(0.5)
names(d18Oc) = c("age", "d18Oc", "si", "d18Oc.sd")

d13Cno = na.exclude(td[c("Age.(Ma)", "d13Coccluded.om.(‰.PDB)", "si")])
d13Cno$d13Cno.sd = rep(0.5)
names(d13Cno) = c("age", "d13Cno", "si", "d13Cno.sd")

d13Cbo = na.exclude(td[c("Age.(Ma)", "d13Cbulk.paleosol.om.(‰.PDB)", "si")])
d13Cbo$d13Cbo.sd = rep(1)
names(d13Cbo) = c("age", "d13Cbo", "si", "d13Cbo.sd")

st = na.exclude(td[c("Age.(Ma)", "Temperature.of.calcium.carbonate.formation.(°C)", "si")])
st$st.sd = rep(2)
names(st) = c("age", "st", "si", "st.sd")

map = na.exclude(td[c("Age.(Ma)", "Mean.Annual.Precipitation.(MAP).(mm)", "si")])
map$map.sd = map$`Mean.Annual.Precipitation.(MAP).(mm)` * 0.2
names(map) = c("age", "map", "si", "map.sd")

ages = data.frame("ages" = c(d13Cc$age, d18Oc$age, d13Cbo$age, 
                             d13Cno$age, st$age, map$age))
ages$si = c(d13Cc$si, d18Oc$si, d13Cbo$si, d13Cno$si, st$si, map$si)
ages = rbind(ages, data.frame("ages" = seq(-252, -140, by = 2),
                              "si" = rep(1)))

ages.v = sort(unique(ages$ages))
si = ages$si[match(ages.v, ages$ages)]
ages = ages.v
dt = ages[2:length(ages)] - ages[1:(length(ages) - 1)]

d18Oc.ai = get.ind(d18Oc$age, ages)
d13Cc.ai = get.ind(d13Cc$age, ages)
d13Cno.ai = get.ind(d13Cno$age, ages)
d13Cbo.ai = get.ind(d13Cbo$age, ages)
st.ai = get.ind(st$age, ages)
map.ai = get.ind(map$age, ages)

d = list(ai = ages, dt = dt, si = si, nsites = 8,
         d13Cc.obs = d13Cc[, c(2, 4)], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, c(2, 4)], d18Oc.ai = d18Oc.ai,
         d13Cno.obs = d13Cno[, c(2, 4)], d13Cno.ai = d13Cno.ai,
         d13Cbo.obs = d13Cbo[, c(2, 4)], d13Cbo.ai = d13Cbo.ai,
         st.obs = st[, c(2, 4)], st.ai = st.ai,
         map.obs = map[, c(2, 4)], map.ai = map.ai)

parms = c("pCO2", "MAT", "GMT", "MAP", "TmPCQ", "PPCQ", "d18.p", 
          "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Cr", "DIFC",
          "d13Cc", "d18Oc", "d13Ca", "ha")

system.time({post = jags.parallel(d, NULL, parms, "code/models/paleosolTS.R", 
                        n.iter = 1e4, n.chains = 3)})

save(post, file = "bigout/post.Rda")
