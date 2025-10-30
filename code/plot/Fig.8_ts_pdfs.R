rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/helpers.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(hjust = 0.1, vjust = -10))
prior_range = read_xlsx("data/input_params_range.xlsx")[, c(1, 4, 6, 7)]

# load Bayesian data ----
load("out/ts_fuxian_8e5.rda")
post = post.ts
load("out/ts_fuxian_D47_8e5.rda")
post_47 = post.ts
load("out/ts_fuxian_MS_8e5.rda")
post_ms = post.ts
load("out/ts_fuxian_CS_1e5.rda")
post_ecs = post.ts
cat("\014")
# mean(post_47_ms_ecs$BUGSoutput$sd$pCO2)

# prior vs posterior distributions ----
prior_post = function(index, name){
  prior_info = prior_range |>
    filter(variable == name)
  if (prior_info$distribution_type == "uniform") {
    prior = data.frame(type = "prior", value = runif(1e6, prior_info$value1, prior_info$value2))
  } else {
    prior = data.frame(type = "prior", value = rbeta(1e6, prior_info$value1, prior_info$value2))
  }
  post_pdf = data.frame(type = "post", value = post$BUGSoutput$sims.list[[name]][, index])
  post_47_pdf = data.frame(type = "post_47", value = post_47$BUGSoutput$sims.list[[name]][, index])
  post_ms_pdf = data.frame(type = "post_ms", value = post_ms$BUGSoutput$sims.list[[name]][, index])
  post_ecs_pdf = data.frame(type = "post_ecs", value = post_ecs$BUGSoutput$sims.list[[name]][, index])
  composite = rbind(prior, post_pdf, post_47_pdf, post_ms_pdf, post_ecs_pdf)
  composite$type = factor(composite$type, levels = c("prior", "post", "post_47", "post_ms", "post_ecs"))
  ext = diff(range(composite$value))
  p = ggplot(data = composite) +
    geom_density(aes(x = value, group = type, fill = type), alpha = .5) +
    scale_fill_viridis_d(option = "mako",
                         labels = c(
                           "prior", "JPI-base",
                           expression("JPI-"*Delta[47]),
                           "JPI-MS", "JPI-CS")) +
    scale_x_continuous(limits = c(floor((min(composite$value) - ext / 10) * 10) / 10, 
                                  ceiling((max(composite$value) + ext / 10) * 10) / 10)) +
    theme_bw() + theme
  return(p)
}
p1 = prior_post(50, "pCO2") + 
  labs(x = expression("CO"[2]*" (ppmv)"), y = "density", fill = "")
p2 = prior_post(50, "MAT") + 
  labs(x = expression(paste("MAT (", degree, "C)")), y = "density", fill = "")
p3 = prior_post(50, "PCQ_to") + 
  labs(x = expression(paste(Delta*"T (", degree, "C)")), y = "density", fill = "")
p4 = prior_post(50, "MAP") + 
  labs(x = "MAP (mm)", y = "density", fill = "")
p5 = prior_post(50, "PCQ_pf") + 
  labs(x = expression(italic(f)[PPCQ]), y = "density", fill = "")
p6 = prior_post(50, "f_R") + 
  labs(x = expression(italic(f)[R]), y = "density", fill = "") 

prior = data.frame(type = "prior", value = rnorm(1e6, 2.02, 0.29))
post_ecs_pdf = data.frame(type = "post_ecs", value = 0.42 + post_ecs$BUGSoutput$sims.list$ECS[, 50] / (5.35*log(2)*0.64))
composite = rbind(prior, post_ecs_pdf)
composite$type = factor(composite$type, levels = c("prior", "post_ecs"))
p7 = ggplot(data = composite) +
  geom_density(aes(x = value, group = type, fill = type), alpha = .5) +
  scale_fill_viridis_d(option = "mako",
                       labels = c("prior", "JPI-CS")) +
  theme_bw() + theme +
  labs(x = expression("S"["CO2,LI"]*" (K/W/m"^"2"*")"))
ggarrange(p1, p2, p3, p4, p5, p6, p7, nrow = 2, ncol = 4, align = "hv",
          labels = c("a", "b", "c", "d", "e", "f", "g"),
          common.legend = TRUE, legend = "right")
# ggsave("figure/prior_post_pdf_ts.png", width = 7.4, height = 5.5, dpi = 500, bg = "white")
ggsave("figure/Fig.8_ts_pdfs.pdf", width = 12, height = 5.5)
