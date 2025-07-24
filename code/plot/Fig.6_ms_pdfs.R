rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/constructors.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"))
prior_range = read_xlsx("data/input_params_range.xlsx")[, c(1, 4, 6, 7)]

# load Bayesian data ----
load("out/ms_fuxian_2e5.rda")
post = post.ms
load("out/ms_fuxian_D47_2e5.rda")
post_47 = post.ms
load("out/ms_fuxian_D47_MS_2e5.rda")
post_ms = post.ms
load("out/ms_fuxian_D47_MS_ECS_2e5.rda")
post_cs = post.ms
cat("\014")

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
  post_47_ms_pdf = data.frame(type = "post_ms", value = post_ms$BUGSoutput$sims.list[[name]][, index])
  post_47_ms_ecs_pdf = data.frame(type = "post_cs", value = post_cs$BUGSoutput$sims.list[[name]][, index])
  composite = rbind(prior, post_pdf, post_47_pdf, post_47_ms_pdf, post_47_ms_ecs_pdf)
  composite$type = factor(composite$type, levels = c("prior", "post", "post_47", "post_ms", "post_cs"))
  ext = diff(range(composite$value))
  p = ggplot(data = composite) +
    geom_density(aes(x = value, group = type, fill = type), alpha = .5) +
    scale_fill_viridis_d(option = "mako",
                         labels = c("prior", "JPI-base",
                                    expression("JPI-"*Delta[47]),
                                    "JPI-MS", "JPI-CS")) +
    scale_x_continuous(limits = c(floor((min(composite$value) - ext / 10) * 10) / 10, 
                                  ceiling((max(composite$value) + ext / 10) * 10) / 10)) +
    theme_bw() + theme
  return(p)
}
p1 = prior_post(50, "pCO2") + 
  # annotate("text", x = 600, y = .01, label = "a", fontface = "bold", size = 6) +
  labs(x = expression("CO"[2]*" (ppmv)"), y = "density", fill = "")
p2 = prior_post(50, "MAT") + 
  # annotate("text", x = 24, y = .12, label = "b", fontface = "bold", size = 6) +
  labs(x = expression(paste("MAT (", degree, "C)")), y = "density", fill = "")
p3 = prior_post(50, "PCQ_to") + 
  # annotate("text", x = 15, y = .16, label = "c", fontface = "bold", size = 6) +
  labs(x = expression(paste(Delta*"T (", degree, "C)")), y = "density", fill = "")
p4 = prior_post(50, "MAP") + 
  # annotate("text", x = 1e3, y = .0042, label = "d", fontface = "bold", size = 6) +
  labs(x = "MAP (mm)", y = "density", fill = "")
p5 = prior_post(50, "PCQ_pf") + 
  # annotate("text", x = 1, y = 4.1, label = "e", fontface = "bold", size = 6) +
  labs(x = expression(italic(f)[PPCQ]), y = "density", fill = "")
p6 = prior_post(50, "f_R") + 
  # annotate("text", x = .45, y = 9, label = "f", fontface = "bold", size = 6) +
  labs(x = expression(italic(f)[R]), y = "density", fill = "") +
  scale_x_continuous(limits = c(-.1, 0.5))
ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, align = "hv",
          labels = c("a", "b", "c", "d", "e", "f"),
          common.legend = TRUE, legend = "top")
ggsave("figure/Fig.6_ms_pdfs.png", width = 8.2, height = 5.9, dpi = 500, bg = "white")


