rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("code/helpers.R")
theme = theme(panel.grid = element_blank(),
              axis.text = element_text(size = 10, color = "black"),
              plot.title = element_text(hjust = 0.1, vjust = -10))
parms = c("pCO2", "MAT", "PCQ_to", "tsc", "MAP", "PCQ_pf", "f_R", "spre")
prior_range = read_xlsx("data/input_params_range.xlsx")[, c(1, 4, 6, 7)]

# load Bayesian data ----
load("out/ts_fuxian_1e4.rda")
post = post.ts
load("out/ts_fuxian_D47_1e4.rda")
post_47 = post.ts
load("out/ts_fuxian_D47_MS_1e4.rda")
post_47_ms = post.ts
load("out/ts_fuxian_D47_MS_ECS_1e4.rda")
post_47_ms_ecs = post.ts
cat("\014")
mean(post_47_ms_ecs$BUGSoutput$sd$pCO2)
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
  post_47_ms_pdf = data.frame(type = "post_47_ms", value = post_47_ms$BUGSoutput$sims.list[[name]][, index])
  post_47_ms_ecs_pdf = data.frame(type = "post_47_ms_ecs", value = post_47_ms_ecs$BUGSoutput$sims.list[[name]][, index])
  composite = rbind(prior, post_pdf, post_47_pdf, post_47_ms_pdf, post_47_ms_ecs_pdf)
  composite$type = factor(composite$type, levels = c("prior", "post", "post_47", "post_47_ms", "post_47_ms_ecs"))
  ext = diff(range(composite$value))
  p = ggplot(data = composite) +
    geom_density(aes(x = value, group = type, fill = type), alpha = .5) +
    scale_fill_viridis_d(option = "mako",
                         labels = c(
                           "prior", "post",
                           expression("post_"*Delta[47]),
                           expression("post_"*Delta[47]*"_MS"),
                           expression("post_"*Delta[47]*"_MS_S"["CO2,LI"])
                         )) +
    scale_x_continuous(limits = c(floor((min(composite$value) - ext / 10) * 10) / 10, 
                                  ceiling((max(composite$value) + ext / 10) * 10) / 10)) +
    theme_bw() + theme
  return(p)
}
p1 = prior_post(50, "pCO2") + 
  annotate("text", x = 1.1e3, y = .0083, label = "a", fontface = "bold", size = 6) +
  labs(x = expression("CO"[2]*" (ppmv)"), y = "density", fill = "")
p2 = prior_post(50, "MAT") + 
  annotate("text", x = -8, y = .13, label = "b", fontface = "bold", size = 6) +
  labs(x = expression(paste("MAT (", degree, "C)")), y = "density", fill = "")
p3 = prior_post(50, "PCQ_to") + 
  annotate("text", x = -22, y = .117, label = "c", fontface = "bold", size = 6) +
  labs(x = expression(paste(Delta*"T (", degree, "C)")), y = "density", fill = "")
p4 = prior_post(50, "MAP") + 
  annotate("text", x = 1.8e3, y = .012, label = "d", fontface = "bold", size = 6) +
  labs(x = "MAP (mm)", y = "density", fill = "")
p5 = prior_post(50, "PCQ_pf") + 
  annotate("text", x = 1, y = 6.8, label = "e", fontface = "bold", size = 6) +
  labs(x = expression(italic(f)[PPCQ]), y = "density", fill = "")
p6 = prior_post(50, "f_R") + 
  annotate("text", x = .45, y = 11.5, label = "f", fontface = "bold", size = 6) +
  labs(x = expression(italic(f)[R]), y = "density", fill = "") +
  scale_x_continuous(limits = c(-.1, 0.5))
ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, align = "hv",
          common.legend = TRUE, legend = "top")
ggsave("figure/prior_post_pdf_ts.png", width = 7.4, height = 5.5, dpi = 500, bg = "white")

# time series w/ iterations ----
ages = seq(-3, 0, .05)
# D47 ----
png("figure/ts_time_series_D47.png", width = 6.5, height = 7.2, units = "in", res = 500)
par(mfrow = c(3, 2), mar = margin(3, 3, 1, 1))
plot.jpi(ages, post$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         mgp = c(2, .8, 0))
text(-2.7, 1200, "a", cex = 1.5, font = 2)
plot.jpi(ages, post_47$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         mgp = c(2, .8, 0))
text(-2.7, 1000, "b", cex = 1.5, font = 2)
text(-1.5, 1000, labels = expression(Delta[47]), cex = 1.5)
plot.jpi(ages, post$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, -5, "c", cex = 1.5, font = 2)
plot.jpi(ages, post_47$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 22, "d", cex = 1.5, font = 2)
text(-1.5, 22, labels = expression(Delta[47]), cex = 1.5)
plot.jpi(ages, post$BUGSoutput$sims.list$PCQ_to, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste(Delta*"T (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 30, "e", cex = 1.5, font = 2)
plot.jpi(ages, post_47$BUGSoutput$sims.list$PCQ_to, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste(Delta*"T (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 32, "f", cex = 1.5, font = 2)
text(-1.5, 32, labels = expression(Delta[47]), cex = 1.5)
dev.off()

# MS ----
png("figure/ts_time_series_D47_MS.png", width = 6.5, height = 7.2, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(3, 3, 1, 1))
plot.jpi(ages, post_47$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         mgp = c(2, .8, 0))
text(-2.7, 1000, "a", cex = 1.5, font = 2)
text(-1.5, 1000, labels = expression(Delta[47]), cex = 1.5)
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$pCO2, n = 500, 
         xlab = "Age (Ma)", ylab = expression("CO"[2]*" (ppmv)"),
         mgp = c(2, .8, 0))
text(-2.7, 80, "b", cex = 1.5, font = 2)
text(-1, 500, labels = expression(Delta[47]*" + MS"), cex = 1.5)
plot.jpi(ages, post_47$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         mgp = c(2, .8, 0))
text(-2.7, 1900, "c", cex = 1.5, font = 2)
text(-1.5, 1900, labels = expression(Delta[47]), cex = 1.5)
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         mgp = c(2, .8, 0))
text(-.3, 180, "d", cex = 1.5, font = 2)
text(-1, 650, labels = expression(Delta[47]*" + MS"), cex = 1.5)
plot.jpi(ages, post_47$BUGSoutput$sims.list$PCQ_pf, n = 500, 
         xlab = "Age (Ma)", ylab = expression("P"[PCQ]),
         mgp = c(2, .8, 0))
text(-2.7, .15, "e", cex = 1.5, font = 2)
text(-1.5, .55, labels = expression(Delta[47]), cex = 1.5)
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$PCQ_pf, n = 500, 
         xlab = "Age (Ma)", ylab = expression("P"[PCQ]),
         mgp = c(2, .8, 0))
text(-2.7, .18, "f", cex = 1.5, font = 2)
text(-1, .9, labels = expression(Delta[47]*" + MS"), cex = 1.5)
dev.off()

# ECS ----
png("figure/ts_time_series_D47_MS_ECS.png", width = 6.5, height = 7.2, units = "in", res = 500)
par(mfrow = c(3, 2), mar = c(3, 3, 1, 1))
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         mgp = c(2, .8, 0))
text(-0.2, 150, "a", cex = 1.5, font = 2)
text(-1, 600, labels = expression(Delta[47]*" + MS"), cex = 1.5)
plot.jpi(ages, post_47_ms_ecs$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP",
         mgp = c(2, .8, 0))
text(-0.2, 1200, "b", cex = 1.5, font = 2)
text(-1.5, 1200, labels = expression(Delta[47]*" + MS + ECS"), cex = 1.5)
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 25, "c", cex = 1.5, font = 2)
text(-1, 25, labels = expression(Delta[47]*" + MS"), cex = 1.5)
plot.jpi(ages, post_47_ms_ecs$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, -10, "d", cex = 1.5, font = 2)
text(-1, -10, labels = expression(Delta[47]*" + MS + ECS"), cex = 1.5)
plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$Tsoil, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("T"[soil]*" (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 30, "e", cex = 1.5, font = 2)
text(-1, 30, labels = expression(Delta[47]*" + MS"), cex = 1.5)
plot.jpi(ages, post_47_ms_ecs$BUGSoutput$sims.list$Tsoil, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("T"[soil]*" (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, -10, "f", cex = 1.5, font = 2)
text(-1, -10, labels = expression(Delta[47]*" + MS + ECS"), cex = 1.5)
dev.off()


# key parameters for different versions ----
png("figure/ts_time_series.png", width = 6.5, height = 7.2, units = "in", res = 500)
par(mfrow = c(3, 2), mar = margin(3, 3, 1, 1))
plot.jpi(ages, post$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, -5, "a", cex = 1.5, font = 2)
text(-2, 20, labels = "base version", cex = 1.1)

plot.jpi(ages, post_47$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 0, "b", cex = 1.5, font = 2)
text(-2, 23, labels = expression(Delta[47]*" version"), cex = 1.1)

plot.jpi(ages, post_47$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         mgp = c(2, .8, 0))
text(-2.7, 1.9e3, "c", cex = 1.5, font = 2)
text(-1.5, 1.8e3, labels = expression(Delta[47]*" version"), cex = 1.1)

plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$MAP, n = 500, 
         xlab = "Age (Ma)", ylab = "MAP (mm)",
         mgp = c(2, .8, 0))
text(-.2, 150, "d", cex = 1.5, font = 2)
text(-1.5, 650, labels = "MS version", cex = 1.1)

plot.jpi(ages, post_47_ms$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, 25, "e", cex = 1.5, font = 2)
text(-1.5, 23, labels = expression(Delta[47]*" version"), cex = 1.1)

plot.jpi(ages, post_47_ms_ecs$BUGSoutput$sims.list$MAT, n = 500, 
         xlab = "Age (Ma)", ylab = expression(paste("MAT (", degree, "C)")),
         mgp = c(2, .8, 0))
text(-2.7, -10, "f", cex = 1.5, font = 2)
text(-1.5, -5, labels = expression("S"[CO2,LI]*" version"), cex = 1.1)
dev.off()

# function to create data frame used to plot parameter curves ----
fm.ts = function(site, var, param) {
  age = seq(-2.6, 0, by = 0.1)
  post.var = data.frame(cbind(site, age,
                              t(apply(param$BUGSoutput$sims.list[[var]], 2, quantile, 
                                      c(0.05, 0.25, 0.5, 0.75, 0.95)))))
  names(post.var) = c("site", "age", "x5", "x25", "median", "x75", "x95")
  post.var[, 2:7] = lapply(post.var[, 2:7], as.numeric)
  post.var$age = -post.var$age
  results = post.var
}

mean(post$BUGSoutput$sd$pCO2)
mean(post_47$BUGSoutput$sd$pCO2)
mean(post_47_ms$BUGSoutput$sd$pCO2)
mean(post_47_ms_ecs$BUGSoutput$sd$pCO2)
