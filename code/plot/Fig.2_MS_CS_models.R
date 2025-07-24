rm(list = ls())
pacman::p_load(tidyverse, readxl, ggpubr, viridisLite)
nsyth = 1e4
set.seed(42)

# load data ----
modern_soil = read_csv("data/CLP_data/MS_porter2001.csv") |>
  drop_na(MAP, MAT, MS)

GMST = read_xlsx("data/global_data/gmst_clark_2024.xlsx")[, c(1, 8)]
names(GMST) = c("age", "gmst")
GMST = GMST |>
  filter(age <= 2.58)

ice_forcing = read_csv("data/global_data/ice_forcing.csv")
ice_co2 = read_csv("data/global_data/ice_core_co2.csv")
boron_co2 = read_csv("data/global_data/boron_co2.csv")

# plot ----
m2 = lm(log10(MS) ~ MAP, data = modern_soil)
summary(m2)
p1 = ggplot(modern_soil, aes(x = MAP, y = log10(MS))) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black") +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 500, y = 1.2,
           label = expression("log"[10]*"("*italic(chi)[lf]*") = 0.0018 * MAP + 0.945")) +
  annotate("text", x = 350, y = 1.9,
           label = expression("R"^"2"*"=0.59, "*italic(p)*" < 0.001")) +
  annotate("text", x = 250, y = 2.3,
           label = "a", size = 8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black")) +
  labs(x = "MAP (mm)",
       y = expression("log"[10]*"("*italic(chi)[lf]*")"))

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
ggsave("figure/Fig.2_MS_CS_models.png", width = 8, height = 4,
       dpi = 500, bg = "white")

