rm(list = ls())
pacman::p_load(tidyverse, readxl, ggpubr)
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(hjust = 0.1, vjust = -10))
post_ms = read_csv("out/ms_fuxian_D47_MS_ECS_2e5.csv")
# post_ts = read_csv("out/ts")

# load CO2 records ----
ice_co2 = read_csv("data/global_data/ice_core_co2.csv")
boron_co2 = read_csv("data/global_data/boron_co2.csv")
blue_co2 = read_csv("data/global_data/AllanHillsCO2CarbonIsotopes.csv")[, c(3,6,8)]
names(blue_co2) = c("co2", "age", "exclude")
blue_co2 = blue_co2 |>
  filter(is.na(exclude) & age > 800) |>
  mutate(age = age / 1e3)
D47 = read_csv("data/CLP_data/D47.csv")
D47$position = 5

# plot CO2 ----
p1 = ggplot(boron_co2, aes(x = age, y = co2)) +
  geom_ribbon(data = post_ms,
              aes(x = age, y = pCO2,
                  ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd),
              fill = "tomato", alpha = .2) +
  geom_errorbar(aes(ymin = co2_lower, ymax = co2_higher),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = paper),
             shape = 22, size = 3, alpha = .5) +
  geom_line(data = ice_co2, aes(x = age, y = co2), color = "grey50") +
  geom_point(data = blue_co2, aes(x = age, y = co2),
             shape = 23, size = 3, fill = "lightblue1") +
  geom_point(data = post_ms, aes(x = age, y = pCO2),
             shape = 21, size = 3, fill = "white", color = "tomato") +
  annotate("text", x = 2.5, y = 160, label = "a", size = 7) +
  scale_fill_viridis_d(option = "mako") +
  guides(fill = guide_legend(
    title = "",
    override.aes = list(alpha = 1)
  )) +
  theme_bw() + theme +
  theme(legend.position = c(.3, .8),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = seq(0, 2.6, .5)) +
  labs(x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

# load local climate data ----
gdgt = read_xlsx("data/CLP_data/lingtai_temperature_lu2022.xlsx")[, c(1, 18)]
names(gdgt) = c("age", "MAST")
p2 = ggplot(post_ms, aes(x = age, y = MAT)) +
  geom_line(data = gdgt, aes(x = age, y = MAST),
            color = "grey50") +
  geom_ribbon(aes(ymin = MAT - MAT_sd,
                  ymax = MAT + MAT_sd),
              fill = "tomato", alpha = .2) +
  geom_point(shape = 21, size = 3, fill = "white", color = "tomato") +
  geom_point(data = D47, aes(x = age, y = position),
             shape = 21, size = 3) +
  annotate("text", x = 2.9, y = 5, label = "b", size = 7) +
  theme_bw() + theme +
  scale_x_continuous(breaks = seq(0, 3, .5)) +
  labs(x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))
p2
# time series posteriors ----
ages = seq(-3, 0, .05)
load("out/ts_fuxian_D47_MS_1e4.rda")
post_ts = data.frame(ages, t(apply(post.ts$BUGSoutput$sims.list$pCO2, 2, quantile,
                                   c(.16, .5, .84))))
names(post_ts) = c("age", "lower", "median", "higher")
load("out/ts_fuxian_D47_MS_ECS_1e4.rda")
post_ts_ecs = data.frame(ages, t(apply(post.ts$BUGSoutput$sims.list$pCO2, 2, quantile,
                                   c(.16, .5, .84))))
names(post_ts_ecs) = c("age", "lower", "median", "higher")
post_ts = post_ts |>
  filter(age >= -2.6) |>
  mutate(age = -age)
post_ts_ecs = post_ts_ecs |>
  filter(age >= -2.6) |>
  mutate(age = -age)
p3 = ggplot(boron_co2, aes(x = age, y = co2)) +
  geom_ribbon(data = post_ts,
                aes(x = age, y = median,
                    ymin = lower, ymax = higher), 
                fill = "tomato", alpha = .2) +
  geom_ribbon(data = post_ts_ecs,
              aes(x = age, y = median,
                  ymin = lower, ymax = higher), 
              fill = "royalblue", alpha = .2) +
  geom_errorbar(aes(ymin = co2_lower, ymax = co2_higher),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = paper),
             shape = 22, size = 3, alpha = .5) +
  geom_line(data = ice_co2, aes(x = age, y = co2), color = "grey50") +
  geom_point(data = blue_co2, aes(x = age, y = co2),
             shape = 23, size = 3, fill = "lightblue1") +
  geom_line(data = post_ts,
            aes(x = age, y = median), 
            linewidth = 2, color = "tomato") +
  geom_line(data = post_ts_ecs,
            aes(x = age, y = median), 
            linewidth = 2, color = "royalblue") +
  annotate("text", x = 2.5, y = 120, label = "c", size = 7) +
  scale_fill_viridis_d(option = "mako") +
  guides(fill = guide_legend(
    title = "",
    override.aes = list(alpha = 1)
  )) +
  theme_bw() + theme +
  theme(legend.position = c(.3, .8),
        legend.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks = seq(0, 2.6, .5)) +
  labs(x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))
p3
load("out/ts_fuxian_D47_MS_1e4.rda")
post_ts = data.frame(ages, t(apply(post.ts$BUGSoutput$sims.list$MAT, 2, quantile,
                                   c(.16, .5, .84))))
names(post_ts) = c("age", "lower", "median", "higher")
load("out/ts_fuxian_D47_MS_ECS_1e4.rda")
post_ts_ecs = data.frame(ages, t(apply(post.ts$BUGSoutput$sims.list$MAT, 2, quantile,
                                       c(.16, .5, .84))))
names(post_ts_ecs) = c("age", "lower", "median", "higher")
post_ts = post_ts |>
  filter(age >= -2.6) |>
  mutate(age = -age)
post_ts_ecs = post_ts_ecs |>
  filter(age >= -2.6) |>
  mutate(age = -age)

p4 = ggplot(post_ts, aes(x = age, y = median)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = higher),
              fill = "tomato", alpha = .2) +
  geom_ribbon(data = post_ts_ecs,
              aes(ymin = lower,
                  ymax = higher),
              fill = "royalblue", alpha = .2) +
  geom_line(data = gdgt, aes(x = age, y = MAST),
            color = "grey50") +
  geom_line(shape = 21, linewidth = 2, fill = "white", color = "tomato") +
  geom_line(data = post_ts_ecs,
            aes(x = age, y = median),
            shape = 21, linewidth = 2, fill = "white", color = "royalblue") +
  geom_point(data = D47, aes(x = age, y = position),
             shape = 21, size = 3) +
  annotate("text", x = 2.9, y = 5, label = "d", size = 7) +
  theme_bw() + theme +
  scale_x_continuous(breaks = seq(0, 3, .5)) +
  labs(x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))
p4
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "v")
ggsave("figure/Fig.5_model_data_comparison.png", width = 9, height = 6, dpi = 500)
