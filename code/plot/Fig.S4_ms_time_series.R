rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/constructors.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              legend.position = c(.8, .9),
              legend.background = element_rect(fill = NA),
              legend.key = element_rect(fill = NA))

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
p1 = ggplot(composite, aes(x = age, y = MAT, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAT - MAT_sd,
                  ymax = MAT + MAT_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 0.1, y = 16.5, label = "a", fontface = "bold", size = 5) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post" = "JPI-base",
                               "post_47" = expression("JPI-"*Delta[47]))) +
  theme_bw() + theme + 
  scale_y_continuous(limits = c(4, 17)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))

p2 = ggplot(composite, aes(x = age, y = PCQ_to, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_to - PCQ_to_sd,
                  ymax = PCQ_to + PCQ_to_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 0.1, y = 14.8, label = "b", fontface = "bold", size = 5) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post" = "JPI-base",
                               "post_47" = expression("JPI-"*Delta[47]))) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(7, 15)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste(Delta*"T (", degree, "C)")))

# post_47 vs post_47_MS ----
post1 = data.frame(post_47, type = "post_47")
post2 = data.frame(post_ms, type = "post_47_MS")
composite = rbind(post1, post2)
p3 = ggplot(composite, aes(x = age, y = MAP, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAP - MAP_sd,
                  ymax = MAP + MAP_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47" = expression("JPI-"*Delta[47]),
                               "post_47_MS" = "JPI-MS")) +
  annotate("text", x = 0.1, y = 980, label = "c", fontface = "bold", size = 5) +
  theme_bw() + theme + 
  scale_y_continuous(limits = c(1e2, 1e3)) +
  labs(fill = "", x = "Age (Ma)", y = "MAP (mm)")

p4 = ggplot(composite, aes(x = age, y = PCQ_pf, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_pf - PCQ_pf_sd,
                  ymax = PCQ_pf + PCQ_pf_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47" = expression("JPI-"*Delta[47]),
                               "post_47_MS" = "JPI-MS")) +
  annotate("text", x = 0.1, y = .77, label = "d", fontface = "bold", size = 5) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(.3, 0.8)) +
  labs(fill = "", x = "Age (Ma)", y = expression(italic(f)[PPCQ]))

# post_47_MS vs post_47_MS_ECS ----
post1 = data.frame(post_ms, type = "post_47_MS")
post2 = data.frame(post_cs, type = "post_47_MS_ECS")
composite = rbind(post1, post2)
p5 = ggplot(composite, aes(x = age, y = MAT, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAT - MAT_sd,
                  ymax = MAT + MAT_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47_MS" = "JPI-MS",
                               "post_47_MS_ECS" = "JPI-CS")) +
  annotate("text", x = 0.1, y = 23, label = "e", fontface = "bold", size = 5) +
  theme_bw() + theme + 
  scale_y_continuous(limits = c(4, 25)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))

p6 = ggplot(composite, aes(x = age, y = PCQ_to, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_to - PCQ_to_sd,
                  ymax = PCQ_to + PCQ_to_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("post_47_MS" = "JPI-MS",
                               "post_47_MS_ECS" = "JPI-CS")) +
  annotate("text", x = 0.1, y = 14.5, label = "f", fontface = "bold", size = 5) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(7, 15)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste(Delta*"T (", degree, "C)")))


ggarrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2, align = "hv")
ggsave("figure/Fig.S4_ms_time_series.png", width = 7.5, height = 7.4, dpi = 500)

