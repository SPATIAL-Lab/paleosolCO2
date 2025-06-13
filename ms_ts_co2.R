rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/constructors.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              legend.position = c(.7, .9),
              legend.background = element_rect(fill = NA),
              legend.key = element_rect(fill = NA),
              legend.text = element_text(size = 10))

#### ms_time_series ----
post = read_csv("out/ms_fuxian_1e5.csv")
post_47 = read_csv("out/ms_fuxian_D47_1e5.csv")
post_47_ms = read_csv("out/ms_fuxian_D47_MS_1e5.csv")
post_47_ms_ecs = read_csv("out/ms_fuxian_D47_MS_ECS_1e5.csv")
mean(post_47_ms_ecs$pCO2_sd)

post1 = data.frame(post[, c(2:11)], type = "post")
post2 = data.frame(post_47[, c(2:11)], type = "post_47")
composite = rbind(post1, post2)
p_47 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd, ymax = pCO2 + pCO2_sd), 
                color = "grey80", linewidth = .2, width = 0,
                show.legend = FALSE) +
  scale_fill_discrete(labels = c(
    expression(delta^"13"*"C + "*delta^"18"*"O"),
    expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47])
  )) +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 0.2, y = 550, label = "a", fontface = "bold", size = 6) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

post1 = data.frame(post_47[, c(2:11)], type = "post_47")
post2 = data.frame(post_47_ms[, c(2:11)], type = "post_47_ms")
composite = rbind(post1, post2)
p_ms = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd, ymax = pCO2 + pCO2_sd), 
                color = "grey80", linewidth = .2, width = 0,
                show.legend = FALSE) +
  scale_fill_discrete(labels = c(
    expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]),
    expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]*" + MS")
  )) +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 0.2, y = 550, label = "b", fontface = "bold", size = 6) +
  theme_bw() + theme + theme(legend.position = c(.6, .9)) +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

post1 = data.frame(post_47_ms[, c(2:11)], type = "post_47_ms")
post2 = data.frame(post_47_ms_ecs[, c(2:11)], type = "post_47_ms_ecs")
composite = rbind(post1, post2)
p_ecs = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd, ymax = pCO2 + pCO2_sd), 
                color = "grey80", linewidth = .2, width = 0,
                show.legend = FALSE) +
  scale_fill_discrete(labels = c(
    expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]*" + MS"),
    expression(delta^"13"*"C+"*delta^"18"*"O+"*Delta[47]*"+MS+S"["CO2,LI"])
  )) +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 0.2, y = 550, label = "c", fontface = "bold", size = 6) +
  theme_bw() + theme + theme(legend.position = c(.6, .9)) +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

ggarrange(p_47, p_ms, p_ecs, nrow = 1, ncol = 3)
ggsave("figure/ms_CO2_time_series.png", width = 10.6, height = 3.5, dpi = 500)

# ts_time_series ----
ages = seq(-3, 0, .05)
load("out/ts_fuxian_1e4.rda")
post = post.ts
load("out/ts_fuxian_D47_1e4.rda")
post_47 = post.ts
load("out/ts_fuxian_D47_MS_1e4.rda")
post_47_ms = post.ts
load("out/ts_fuxian_D47_MS_ECS_1e4.rda")
post_47_ms_ecs = post.ts
cat("\014")
sample_res = data.frame(age = post1$age, position = 0)
plot_ts = function(num, dat, ages){
  for (i in 1:num) {
    post_df = data.frame(age = -ages, co2 = dat$BUGSoutput$sims.list$pCO2[i, ])
    if(i == 1){
      p1 = ggplot(post_df, aes(x = age, y = co2)) +
        geom_line(alpha = .1) +
        theme_bw() + theme +
        labs(x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)")) 
    } else {
      p1 = p1 + geom_line(data = post_df, aes(x = age, y = co2), alpha = .1)
    }
  }
  post_df_mean = data.frame(age = -ages, co2 = dat$BUGSoutput$mean$pCO2)
  p1 = p1 + 
    geom_line(data = post_df_mean, aes(x = age, y = co2), color = "tomato", size = 2) +
    geom_point(data = sample_res, aes(x = age, y = position), shape = 21, size = 2) +
    scale_x_continuous(breaks = seq(0, 3, .5))
  return(p1)
}
p_base = plot_ts(500, post, ages) +
  annotate("text", x = 2.8, y = 100, label = "a", size = 6, fontface = "bold") +
  annotate("text", x = 2.5, y = 1000, label = expression(delta^"13"*"C + "*delta^"18"*"O"))
p_47 = plot_ts(500, post_47, ages) +
  annotate("text", x = 2.8, y = 100, label = "b", size = 6, fontface = "bold") +
  annotate("text", x = 2.3, y = 1000, label = expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]))
p_ms = plot_ts(500, post_47_ms, ages) +
  annotate("text", x = 2.8, y = 50, label = "c", size = 6, fontface = "bold") +
  annotate("text", x = 1, y = 500, label = expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]*" + MS"))
p_ecs = plot_ts(500, post_47_ms_ecs, ages) +
  annotate("text", x = 2.8, y = 50, label = "d", size = 6, fontface = "bold") +
  annotate("text", x = 1, y = 500, label = expression(delta^"13"*"C + "*delta^"18"*"O + "*Delta[47]*" + MS + ECS"))
ggarrange(p_base, p_47, p_ms, p_ecs, nrow = 2, ncol = 2)
ggsave("figure/ts_CO2_time_series.png", width = 7.5, height = 6.4, dpi = 500)

