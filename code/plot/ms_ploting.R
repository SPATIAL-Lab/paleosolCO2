rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/constructors.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(vjust = -10, hjust = .1),
              plot.margin = margin(2, 2, 2, 2))
parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf", "f_R", "spre")
prior_range = read_xlsx("data/input_params_range.xlsx")[, c(1, 4, 6, 7)]

# load Bayesian data ----
load("out/ms_fuxian_1e5.rda")
post = post.ms
load("out/ms_fuxian_D47_1e5.rda")
post_47 = post.ms
load("out/ms_fuxian_D47_MS_1e5.rda")
post_47_ms = post.ms
load("out/ms_fuxian_D47_MS_ECS_1e5.rda")
post_47_ms_ecs = post.ms
cat("\014")

# prior vs posterior distributions ----
prior_post = function(index, parms){
  plot_list = list()
  for (i in 1:length(parms)) {
    name = parms[i]
    prior_info = prior_range |>
      filter(variable == name)
    if (prior_info$distribution_type == "uniform") {
      prior = data.frame(type = "prior", value = runif(1e6, prior_info$value1, prior_info$value2))
    } else {
      prior = data.frame(type = "prior", value = rbeta(1e6, prior_info$value1, prior_info$value2))
    }
    post_pdf = data.frame(type = "post", value = post$BUGSoutput$sims.list[[name]][, index])
    post_47_pdf = data.frame(type = "post_47", value = post_47$BUGSoutput$sims.list[[name]][, index])
    post_47_ms_pdf = data.frame(type = "post_47_ms", value = post_47_ms$BUGSoutput$sims.list[[name]][, index])
    post_47_ms_ecs_pdf = data.frame(type = "post_47_ms_ecs", value = post_47_ms_ecs$BUGSoutput$sims.list[[name]][, index])
    composite = rbind(prior, post_pdf, post_47_pdf, post_47_ms_pdf, post_47_ms_ecs_pdf)
    composite$type = factor(composite$type, levels = c("prior", "post", "post_47", "post_47_ms", "post_47_ms_ecs"))
    ext = diff(range(composite$value))
    p = ggplot(data = composite) +
      geom_density(aes(x = value, group = type, fill = type), alpha = .5) +
      scale_fill_viridis_d(option = "mako") +
      scale_x_continuous(limits = c(floor((min(composite$value) - ext / 10) * 10) / 10, 
                                    ceiling((max(composite$value) + ext / 10) * 10) / 10)) +
      theme_bw() + theme +
      labs(x = name)
    plot_list[[i]] = p
  }
  ggarrange(plotlist = plot_list, nrow = 3, ncol = 3, 
            align = "hv", common.legend = TRUE, legend = "right")
}
prior_post(50, parms) 
ggsave("figure/prior_post_pdf.png", width = 8.3, height = 7.3, dpi = 500)

#### time series ----
post = read_csv("out/ms_fuxian_1e5.csv")
post_47 = read_csv("out/ms_fuxian_D47_1e5.csv")
post_47_ms = read_csv("out/ms_fuxian_D47_MS_1e5.csv")
post_47_ms_ecs = read_csv("out/ms_fuxian_D47_MS_ECS_1e5.csv")

# post vs post_47 ----
post1 = data.frame(post[, c(2:11)], type = "post")
post2 = data.frame(post_47[, c(2:11)], type = "post_47")
composite = rbind(post1, post2)
p1 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_ribbon(aes(ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  guides(fill = "none") +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

p2 = ggplot(composite, aes(x = age, y = MAT, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAT - MAT_sd,
                  ymax = MAT + MAT_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme + 
  guides(fill = "none") +
  scale_y_continuous(limits = c(4, 17)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))

p3 = ggplot(composite, aes(x = age, y = PCQ_to, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_to - PCQ_to_sd,
                  ymax = PCQ_to + PCQ_to_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_discrete(labels = c(expression("w/o "*Delta[47]), expression("w "*Delta[47]))) +
  theme_bw() + theme +
  theme(legend.position = c(.8, .2),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA)) +
  scale_y_continuous(limits = c(7, 15)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste(Delta*"T (", degree, "C)")))

ggarrange(p1, p2, p3, nrow = 3, ncol = 1, align = "hv")
ggsave("figure/time_series_D47.png", width = 4, height = 7, dpi = 500)

# post_47 vs post_47_MS ----
post1 = data.frame(post_47, type = "post_47")
post2 = data.frame(post_47_ms, type = "post_47_MS")
composite = rbind(post1, post2)
p1 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_ribbon(aes(ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  guides(fill = "none") +
  scale_y_continuous(limits = c(150, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

p2 = ggplot(composite, aes(x = age, y = MAP, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAP - MAP_sd,
                  ymax = MAP + MAP_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme + 
  guides(fill = "none") +
  scale_y_continuous(limits = c(1e2, 1e3)) +
  labs(fill = "", x = "Age (Ma)", y = "MAP (mm)")

p3 = ggplot(composite, aes(x = age, y = PCQ_pf, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_pf - PCQ_pf_sd,
                  ymax = PCQ_pf + PCQ_pf_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_discrete(labels = c("w/o MS", "w MS")) +
  theme_bw() + theme +
  theme(legend.position = c(.8, .8),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA)) +
  scale_y_continuous(limits = c(.3, 1)) +
  labs(fill = "", x = "Age (Ma)", y = expression("P"[PCQ]))

ggarrange(p1, p2, p3, nrow = 3, ncol = 1, align = "hv")
ggsave("figure/time_series_MS.png", width = 4, height = 7, dpi = 500)


# post_47_MS vs post_47_MS_ECS ----
post1 = data.frame(post_47_ms, type = "post_47_MS")
post2 = data.frame(post_47_ms_ecs[1:32], type = "post_47_MS_ECS")
composite = rbind(post1, post2)
p1 = ggplot(composite, aes(x = age, y = pCO2, group = type, fill = type)) +
  geom_ribbon(aes(ymin = pCO2 - pCO2_sd,
                  ymax = pCO2 + pCO2_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme +
  guides(fill = "none") +
  scale_y_continuous(limits = c(100, 600)) +
  labs(fill = "", x = "Age (Ma)", y = expression("CO"[2]*" (ppmv)"))

p2 = ggplot(composite, aes(x = age, y = MAT, group = type, fill = type)) +
  geom_ribbon(aes(ymin = MAT - MAT_sd,
                  ymax = MAT + MAT_sd), alpha = .2) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme + 
  guides(fill = "none") +
  scale_y_continuous(limits = c(4, 25)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste("MAT (", degree, "C)")))

p3 = ggplot(composite, aes(x = age, y = PCQ_to, group = type, fill = type)) +
  geom_ribbon(aes(ymin = PCQ_to - PCQ_to_sd,
                  ymax = PCQ_to + PCQ_to_sd), alpha = .2, show.legend = FALSE) +
  geom_point(shape = 21, size = 3) +
  scale_fill_discrete(labels = c("w/o ECS", "w ECS")) +
  theme_bw() + theme +
  theme(legend.position = c(.8, .2),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA)) +
  scale_y_continuous(limits = c(7, 15)) +
  labs(fill = "", x = "Age (Ma)", y = expression(paste(Delta*"T (", degree, "C)")))

ggarrange(p1, p2, p3, nrow = 3, ncol = 1, align = "hv")
ggsave("figure/time_series_ECS.png", width = 4, height = 7, dpi = 500)


# function to create data frame used to plot parameter curves ----
fm = function(site, var, param) {
  if(site == "Luochuan") {
    age = read.csv("data/loess_interglacial.csv")
  } else {
    age = read.csv("data/loess_glacial.csv") %>% filter(section == site)
  }
  post.var = data.frame(cbind(site, age$age,
                              t(apply(param$BUGSoutput$sims.list[[var]], 2, quantile, 
                                      c(0.05, 0.25, 0.5, 0.75, 0.95)))))
  names(post.var) = c("site", "age", "x5", "x25", "median", "x75", "x95")
  post.var[, 2:7] = lapply(post.var[, 2:7], as.numeric)
  results = post.var
}


