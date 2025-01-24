load("out/ms_fx_1e4.rda")
fx = post.clp
load("out/ms_zjc_1e4.rda")
zjc = post.clp

inv = read_xlsx("data/Dataset S1.xlsx")
inv = inv[1:40, c(1,3,6,8,10)]
names(inv) = c("site", "age", "d18c", "R","CO2")
inv$age = round(inv$age/1000, 3)

# Ratio 
param = "Ratio"
fx.r = fm("Fuxian", param, fx)
zjc.r = fm("Zhaojiachuan", param, zjc)
ms.r = rbind(fx.r, zjc.r)
ms.r = ms.r %>% filter(age <= 0.8)
dat.r = merge(ms.r, inv, by = "age")
ggplot(dat.r) +
  geom_errorbar()
  geom_point(aes(x = d18c, y= median), color = "blue", shape = 21, size = 3) +
  geom_point(aes(x = d18c, y = R), color = "red", shape = 21, size = 3)

# S_z
param = "S_z"
fx.sz = fm("Fuxian", param, fx)
zjc.sz = fm("Zhaojiachuan", param, zjc)
ms.sz = rbind(fx.sz, zjc.sz)
ms.sz = ms.sz %>% filter(age <= 0.8)
dat.sz = merge(ms.sz, inv, by = "age")
dat.r = dat.r %>%
  mutate(Sz.inv = CO2/median)
ggplot(dat.sz, aes(x = d18c, y = median)) +
  geom_point() +
  geom_point(data = dat.r, aes(x = d18c, y = Sz.inv), color = "red")

