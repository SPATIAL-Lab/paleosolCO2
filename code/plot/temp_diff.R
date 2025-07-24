rm(list = ls())
pacman::p_load(tidyverse, readxl)
gmst = read_xlsx("data/global_data/gmst_clark_2024.xlsx")[, c(1, 8)]
lingtai = read_xlsx("data/CLP_data/lingtai_temperature_lu2022.xlsx")[, c(1, 18)]
names(lingtai) = c("age", "msat")
lingtai$gmst = approx(gmst$`Age (Ma)`, gmst$`Global mean surface temp (GMST) (oC) area weighted`, xout = lingtai$age)$y
lingtai = lingtai |>
  mutate(temp_diff = msat - gmst)
range(lingtai$temp_diff)

ggplot(lingtai, aes(x = age, y = msat)) +
  geom_line()
