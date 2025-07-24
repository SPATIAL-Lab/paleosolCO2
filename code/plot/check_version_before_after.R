load("out/ms_fx_1e4.rda")
before = post.clp
load("out/ms_fx_1e4_v2.rda")
after = post.clp
age = read_csv("data/loess_glacial.csv") %>% filter(section == "Fuxian")
params.before = data.frame(cbind("before", age$age, 
                           t(apply(before$BUGSoutput$sims.list$PCQ_to, 2, quantile, 
                                   c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(params.before) = c("version", "age", "x5", "x25", "median", "x75", "x95")
params.before[,2:7] = lapply(params.before[,2:7], as.numeric)
params.after = data.frame(cbind("after", age$age, 
                                 t(apply(after$BUGSoutput$sims.list$PCQ_to, 2, quantile, 
                                         c(0.05, 0.25, 0.5, 0.75, 0.95)))))
names(params.after) = c("version", "age", "x5", "x25", "median", "x75", "x95")
params.after[,2:7] = lapply(params.after[,2:7], as.numeric)
params = rbind(params.before, params.after)
ggplot(params, aes(x = age, y = median, fill = version)) +
  geom_ribbon(aes(ymin = x5, ymax = x95), alpha = 0.3) +
  geom_line(aes(color = version)) +
  scale_fill_manual(values = c("firebrick2", "royalblue")) +
  scale_color_manual(values = c("firebrick2", "royalblue")) +
  theme_bw() + theme +
  labs(x = "Age (Ma)", y = expression("CO2")) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5))
