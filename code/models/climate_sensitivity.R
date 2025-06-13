rm(list = ls())
pacman::p_load(tidyverse, readxl, ggpubr, viridisLite)
nsyth = 1e4
set.seed(42)

# load and groom data ----
GMST = read_xlsx("data/global_data/gmst_clark_2024.xlsx")[, c(1, 8)]
names(GMST) = c("age", "gmst")
GMST = GMST |>
  filter(age <= 2.58)

ice_forcing = read_xlsx("data/global_data/ice_forcing_stap_2018.xlsx", sheet = 2)
ice_forcing = ice_forcing[2:nrow(ice_forcing), 1:2]
names(ice_forcing) = c("age", "ice_forcing")
ice_forcing = ice_forcing |>
  mutate(across(everything(), as.numeric)) |>
  filter(age <= 2.58)
write.csv(ice_forcing, file = "data/global_data/ice_forcing.csv")

ice_co2 = read_xls("data/global_data/ice_core_bereiter_2015.xls", sheet = 3)
ice_co2 = ice_co2[15:nrow(ice_co2), 1:2]
names(ice_co2) = c("age", "co2")
ice_co2 = ice_co2 |>
  mutate(across(everything(), as.numeric)) |>
  mutate(age = age / 1e6) |>
  filter(age > 0)
write_csv(ice_co2, "data/global_data/ice_core_co2.csv")

boron_co2 = read_csv("data/global_data/co2_proxy_data.csv")[, c(4, 7:9, 10)] |>
  filter(age > 800) |>
  mutate(age = age / 1e3)
colnames(boron_co2)[3:4] = c("lower", "higher")
ref = unique(boron_co2$reference)
boron_co2 = boron_co2 |>
  mutate(paper = reference) |>
  mutate(paper = case_when(
    reference == ref[1] ~ "Boti et al. (2015)",
    reference == ref[2] ~ "Sosdian et al. (2018)",
    reference == ref[3] ~ "Dyez et al. (2018)",
    reference == ref[4] ~ "de la Vega et al. (2020)"
  )) |>
  mutate(co2_lower = co2 - lower,
         co2_higher = co2 + higher)
write_csv(boron_co2, "data/global_data/boron_co2.csv")

# time series plot ----
pal = mako(5)
png("figure/climate_sensitivity.png", width = 5, height = 5.7, units = "in", res = 500)
par(mar = c(4,4,4,4))
plot(-1, 0, xlim = c(0, 2.6), ylim = c(0,3), axes = FALSE, 
     xlab = "", ylab = "")

axis(3, mgp = c(2, 1.2, .5))
mtext("Age (Ma)", 3, line = 2.5)

yext = range(GMST$gmst)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 2)
gmst.rs = cbind(GMST$age,
                2 + (GMST$gmst - min(tix)) / diff(range(tix)))
lines(gmst.rs[, 1], gmst.rs[, 2], col = pal[3])
axis(2, 2 + (tix - min(tix)) / diff(range(tix)))
mtext(expression(paste(Delta*"GMST (", degree, "C)")), 2, line = 2.5, at = 2.5)

pal_proxy = pal[factor(boron_co2$paper)]
yext = range(boron_co2$co2_lower, boron_co2$co2_higher)
tix = seq(floor(min(yext)+1), ceiling(max(yext)+20), by = 100)
co2.rs = cbind(boron_co2$age,
               1 + (boron_co2$co2 - min(tix)) / diff(range(tix)),
               1 + (boron_co2$co2_lower - min(tix)) / diff(range(tix)),
               1 + (boron_co2$co2_higher - min(tix)) / diff(range(tix)))
arrows(co2.rs[,1], co2.rs[,3], co2.rs[,1], co2.rs[,4], col = "grey",
       angle = 90, length = 0, code = 0)
points(co2.rs[,1], co2.rs[,2], col = "black", bg = pal_proxy,
       pch = 21, cex = 1.3)
ice_co2.rs = cbind(ice_co2$age,
                   1 + (ice_co2$co2 - min(tix)) / diff(range(tix)))
lines(ice_co2.rs[,1], ice_co2.rs[,2], col = pal[2])
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("CO"[2]*" (ppmv)"), 4, line = 2.5, at = 1.5)
legend(x = .9, y = 2.3, legend = unique(boron_co2$paper),
       pt.bg = pal, pch = 21, cex = .8, pt.cex = 1, bty = "n")

yext = range(ice_forcing$ice_forcing)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 1)
R_ice.rs = cbind(ice_forcing$age,
                 0 + (ice_forcing$ice_forcing - min(tix)) / diff(range(tix)))
lines(R_ice.rs[,1], R_ice.rs[,2], col = pal[4])
axis(2, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta*"R"[ice]*" (W/K/m"^"2"*")"), 2, line = 2.5, at = .5)

axis(1)
mtext("Age (Ma)", 1, line = 2.5)
text(.2, 3, labels = "a", font = 2, cex = 1.2)
text(.2, 1.5, labels = "b", font = 2, cex = 1.2)
text(.2, 0, labels = "c", font = 2, cex = 1.2)

dev.off()

# ECS ----
ECS = data.frame(age = seq(0, 2.6, .1))
ECS_sum = data.frame(matrix(nrow = (nrow(ECS) - 1),
                            ncol = 5))
names(ECS_sum) = c("time", "R_sf", "R_sf_sd", "gmst", "gmst_sd")
co2_composite = rbind(ice_co2, boron_co2[, 1:2])
for (i in 1:(nrow(ECS)-1)) {
  age_min = ECS$age[i]
  age_max = ECS$age[i+1]
  ECS_sum$time[i] = mean(age_min, age_max)
  co2 = co2_composite |>
    filter(age > age_min & age < age_max)
  R_ice = ice_forcing |>
    filter(age > age_min & age < age_max)
  co2_s = sample(co2$co2, nsyth, replace = TRUE)
  ice_s = sample(R_ice$ice_forcing, nsyth, replace = TRUE)
  R_slow = 5.35 * log(co2_s / 278) + .45 * ice_s
  ECS_sum$R_sf[i] = mean(R_slow)
  ECS_sum$R_sf_sd[i] = sd(R_slow)
  temp = GMST |>
    filter(age > age_min & age < age_max) |>
    summarize(mean = mean(gmst),
              sd = sd(gmst))
  ECS_sum[i, 4:5] = temp
}
ECS_sum = ECS_sum |>
  filter(R_sf > -4)

m1 = lm(gmst ~ R_sf, data = ECS_sum)
summary(m1)
p2 = ggplot(ECS_sum, aes(x = R_sf, y = gmst)) +
  geom_errorbar(aes(xmin = R_sf - R_sf_sd, xmax = R_sf + R_sf_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = gmst - gmst_sd, ymax = gmst + gmst_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  geom_point(aes(fill = time), shape = 21, size = 4) +
  annotate("text", x = -2.5, y = 3, label = expression("R"^"2"*" = 0.67")) +
  annotate("text", x = -2.5, y = 2.2, label = expression(italic(p)*" < 0.001")) +
  annotate("text", x = 1.5, y = 3.9, label = "b",
           size = 8, face = "bold") +
  scale_fill_viridis_c(option = "mako") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.8, .3),
        legend.title = element_text(margin = margin(b = 10)),
        axis.text = element_text(size = 10, color = "black")) +
  labs(x = expression(Delta*"R"["CO2,LI"]*" (W/K/m"^"2"*")"),
       y = expression(paste(Delta*"GMST (", degree, "C)")),
       fill = "Age (Ma)")
# ggsave("figure/climate_sensitivity_2.png", width = 3.5, height = 3.8, dpi = 500)  

ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv")
ggsave("figure/empirical_relationship.png", width = 7.5, height = 4,
       dpi = 500, bg = "white")
