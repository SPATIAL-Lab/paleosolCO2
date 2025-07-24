rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/constructors.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              legend.background = element_blank(),
              legend.position = c(.8, .9))

#### time series ----
post = read_csv("out/ms_fuxian_2e5.csv")
post_47 = read_csv("out/ms_fuxian_D47_2e5.csv")
post_ms = read_csv("out/ms_fuxian_D47_MS_2e5.csv")
post_cs = read_csv("out/ms_fuxian_D47_MS_ECS_2e5.csv")
mean(post_cs$pCO2_sd)

# post vs post_47 ----
post1 = data.frame(post[, c(2:11)], type = "post")
post2 = data.frame(post_47[, c(2:11)], type = "post_47")
composite = rbind(post1, post2)
p1 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd,
                    ymax = pCO2 + pCO2_sd), alpha = .1) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post" = "JPI-base",
                               "post_47" = expression("JPI-"*Delta[47]))) +
  annotate("text", x = 0.1, y = 570, label = "a", fontface = "bold", size = 5) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

# post_47 vs post_47_MS ----
post1 = data.frame(post_47, type = "post_47")
post2 = data.frame(post_ms, type = "post_47_MS")
composite = rbind(post1, post2)
p2 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd), alpha = .1) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47" = expression("JPI-"*Delta[47]),
                               "post_47_MS" = "JPI-MS")) +
  annotate("text", x = 0.1, y = 570, label = "b", fontface = "bold", size = 5) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

# post_47_MS vs post_47_MS_ECS ----
post1 = data.frame(post_ms, type = "post_47_MS")
post2 = data.frame(post_cs, type = "post_47_MS_ECS")
composite = rbind(post1, post2)
p3 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_errorbar(aes(ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47_MS" = "JPI-MS",
                               "post_47_MS_ECS" = "JPI-CS")) +
  annotate("text", x = 0.1, y = 570, label = "c", fontface = "bold", size = 5) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(100, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

ggarrange(p1, p2, p3, nrow = 1, ncol = 3, align = "hv")
ggsave("figure/Fig.4_ms_co2_time_series.png", width = 10, height = 3.5, dpi = 500)
 
