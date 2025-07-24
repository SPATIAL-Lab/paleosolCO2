rm(list = ls())
pacman::p_load(tidyverse, readxl, ggpubr)
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(hjust = 0.1, vjust = -10))
post = read_csv("out/ms_fuxian_2e5.csv")
post_47 = read_csv("out/ms_fuxian_D47_2e5.csv")
post_47_ms = read_csv("out/ms_fuxian_D47_MS_2e5.csv")
post_47_ms_ecs = read_csv("out/ms_fuxian_D47_MS_ECS_2e5.csv")
# post_ts = read_csv("out/ts")
ice_co2 = read_csv("data/global_data/ice_core_co2.csv")

post_ice = function(post, label, title) {
  post = post |>
    select(age, pCO2, pCO2_sd) |>
    filter(age <= .8)

  for (i in 1:nrow(post)) {
    age_min = post$age[i] - 0.01
    age_max = post$age[i]
    ice_co2_sample = ice_co2 |>
      filter(age >= age_min & age <= age_max)
    post$ice_co2[i] = mean(ice_co2_sample$co2)
    post$ice_co2_sd[i] = sd(ice_co2_sample$co2)
  }
  
  post = post |>
    mutate(ice_co2_low = ice_co2 - ice_co2_sd,
           ice_co2_high = ice_co2 + ice_co2_sd,
           post_co2_low = pCO2 - pCO2_sd,
           post_co2_high = pCO2 + pCO2_sd)
  co2_range = range(post$ice_co2_low, post$ice_co2_high, post$post_co2_low, post$post_co2_high)
  ggplot(post, aes(x = ice_co2, y = pCO2)) +
    geom_errorbar(aes(xmin = ice_co2 - ice_co2_sd,
                      xmax = ice_co2 + ice_co2_sd),
                  linewidth = .2, width = 0, color = "grey80") +
    geom_errorbar(aes(ymin = pCO2 - pCO2_sd,
                      ymax = pCO2 + pCO2_sd),
                  linewidth = .2, width = 0, color = "grey80") +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(fill = age), shape = 21, size = 3) +
    annotate("text", x = floor(co2_range[1]) + 10, y = ceiling(co2_range[2]) - 10, label = label,
             size = 5, fontface = "bold") +
    annotate("text", x = mean(co2_range), y = ceiling(co2_range[2]) - 10, label = title,
             size = 5, fontface = "bold") +
    scale_fill_viridis_c(option = "mako", limits = c(0, .8)) +
    guides(fill = guide_colorbar(
      breaks = seq(0, .8, .2),
      labels = c("0", "0.2", "0.4", "0.6", "0.8")
    )) +
    theme_bw() + theme +
    theme(legend.position = c(.8, .3),
          legend.title = element_text(margin = margin(b = 15))) +
    scale_x_continuous(limits = c(floor(co2_range[1]), ceiling(co2_range[2]))) +
    scale_y_continuous(limits = c(floor(co2_range[1]), ceiling(co2_range[2]))) +
    labs(x = expression("ice-core CO"[2]*" (ppmv)"),
         y = expression("posterior CO"[2]*" (ppmv)"),
         fill = "Age (Ma)")
}


p1 = post_ice(post, "a", NA)
title = expression(Delta[47])
p2 = post_ice(post_47, "b", title)
title = expression(Delta[47]*" + MS")
p3 = post_ice(post_47_ms, "c", title)
title = expression(Delta[47]*" + MS + ECS")
p4 = post_ice(post_47_ms_ecs, "d", title)
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, align = "hv")
ggsave("figure/800ky_CO2.png", width = 6.7, height = 7, dpi = 500)
