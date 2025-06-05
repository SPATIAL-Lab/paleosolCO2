rm(list = ls())
pacman::p_load(tidyverse)
modern_soil = read_csv("data/CLP_data/MS_porter2001.csv") |>
  drop_na(MAP, MAT, MS)

m1 = nls(MS ~ a * exp(b * MAP), data = modern_soil, start = c(a = 8.8, b = 4e-3))
summary(m1)

base_text_size = 10
modern_soil$fit = predict(m1)
a = coef(m1)[1]
b = coef(m1)[2]
ggplot(modern_soil, aes(x = MAP, y = MS)) +
  # geom_line(aes(y = fit), color = "red", size = 1) +
  # stat_function(fun = function(x) a * exp(b * x), 
  #               color = "grey", size = 1.5, linetype = "dashed") +
  geom_point(shape = 21, size = 3) +
  # annotate("text", x = 300, y = 200, 
  #          label = expression(italic(chi)[lf]* "= 10.4 e"^"0.00399 * MAP")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        ) +
  labs(x = "MAP (mm)",
       y = expression(italic(chi)[lf]* " (10"^"-8"*"m"^"3"*"kg"^"-1"*")"))

m2 = lm(log10(MS) ~ MAP, data = modern_soil)
summary(m2)
ggplot(modern_soil, aes(x = MAP, y = log10(MS))) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black") +
  geom_point(shape = 21, size = 3) +
  annotate("text", x = 500, y = 1.2,
           label = expression("log"[10]*"("*italic(chi)[lf]*") = 0.0018 * MAP + 0.945")) +
  annotate("text", x = 350, y = 1.9,
           label = expression("R"^"2"*"=0.59, "*italic(p)*" < 0.001")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
  ) +
  labs(x = "MAP (mm)",
       y = expression("log"[10]*"("*italic(chi)[lf]*")"))
ggsave("figure/MAP_MS_model.png", width = 4, height = 4.5, dpi = 500)
