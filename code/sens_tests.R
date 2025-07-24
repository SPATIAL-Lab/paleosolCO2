library(tidyverse)

source("code/models/forward_model.R")
vars = ctrl()
baseCase = fm(vars)

# MAT
vars = ctrl()
vars$MAT = (5:20)
sens.t = fm(vars)
sens.t$MAT = vars$MAT
sens.plot(sens.t, "MAT")

# MAP
vars = ctrl()
vars$MAP = (100:2000)
sens.t = fm(vars)
sens.t$MAP = vars$MAP
sens.plot(sens.t, "MAP")

# pCO2
vars = ctrl()
vars$pCO2 = (150:2000)
sens.t = fm(vars)
sens.t$pCO2 = vars$pCO2
sens.plot(sens.t, "pCO2")

# PCQ.pf
vars = ctrl()
vars$PCQ.pf = seq(0, 0.5, 0.01)
sens.t = fm(vars)
sens.t$PCQ.pf = vars$PCQ.pf
sens.plot(sens.t, "PCQ.pf")

# PCQ_to
vars = ctrl()
vars$PCQ_to = (-5:15)
sens.t = fm(vars)
sens.t$PCQ_to = vars$PCQ_to
sens.plot(sens.t, "PCQ_to")

# tsc
vars = ctrl()
vars$tsc = seq(0,1,0.01)
sens.t = fm(vars)
sens.t$tsc = vars$tsc
sens.plot(sens.t, "tsc")

# porosity
vars = ctrl()
vars$pore = seq(0.1, 0.5, 0.01)
sens.t = fm(vars)
sens.t$pore = vars$pore
sens.plot(sens.t, "pore")
