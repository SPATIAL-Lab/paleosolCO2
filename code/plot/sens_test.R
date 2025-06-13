rm(list = ls())
source('code/models/forward_model.R')

# CO2 ----
vars = ctrl()
for (i in 1:10) {
  vars$pCO2 = 1e2 * i + 1e2 # 100-1000 ppm
  sims = fm(vars)
  sims$pCO2 = vars$pCO2
  if (i == 1){
    results_CO2 = sims
  } else {
    results_CO2 = rbind(results_CO2, sims)
  }
}

# MAT ----
vars = ctrl()
for (i in 1:31) {
  vars$MAT = i-1 # 0-30 degree
  sims = fm(vars)
  sims$MAT = vars$MAT
  if (i == 1){
    results_MAT = sims
  } else {
    results_MAT = rbind(results_MAT, sims)
  }
}

# PCQ_to ----
vars = ctrl()
for (i in 1:31) {
  vars$PCQ_to = i-1 # 0-30 degree
  sims = fm(vars)
  sims$PCQ_to = vars$PCQ_to
  if (i == 1){
    results_PCQ_to = sims
  } else {
    results_PCQ_to = rbind(results_PCQ_to, sims)
  }
}

# MAP ----
vars = ctrl()
for (i in 1:100) {
  vars$MAP = 10 * i # 10-1000 ppm
  sims = fm(vars)
  sims$MAP = vars$MAP
  if (i == 1){
    results_MAP = sims
  } else {
    results_MAP = rbind(results_MAP, sims)
  }
}

# PCQ_pf ----
vars = ctrl()
for (i in 1:100) {
  vars$PCQ_pf = .01 * i # 0-1
  sims = fm(vars)
  sims$PCQ_pf = vars$PCQ_pf
  if (i == 1){
    results_PCQ_pf = sims
  } else {
    results_PCQ_pf = rbind(results_PCQ_pf, sims)
  }
}

# pore ----
vars = ctrl()
for (i in 1:60) {
  vars$pore = .01 * i + .2 # 0.2-0.8
  sims = fm(vars)
  sims$pore = vars$pore
  if (i == 1){
    results_pore = sims
  } else {
    results_pore = rbind(results_pore, sims)
  }
}

# f_R ----
vars = ctrl()
for (i in 1:40) {
  vars$f_R = .01 * i # 0.01-0.4
  sims = fm(vars)
  sims$f_R = vars$f_R
  if (i == 1){
    results_f_R = sims
  } else {
    results_f_R = rbind(results_f_R, sims)
  }
}

# spre ----
vars = ctrl()
for (i in 1:40) {
  vars$spre = .01 * i + .4 # 0.4-0.8
  sims = fm(vars)
  sims$spre = vars$spre
  if (i == 1){
    results_spre = sims
  } else {
    results_spre = rbind(results_spre, sims)
  }
}

# tsc ----
vars = ctrl()
for (i in 1:100) {
  vars$tsc = .005 * i # 0-0.5
  sims = fm(vars)
  sims$tsc = vars$tsc
  if (i == 1){
    results_tsc = sims
  } else {
    results_tsc = rbind(results_tsc, sims)
  }
}

# d13Cc plot ----
png("figure/sens_d13c.png", width = 6.2, height = 6.2, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(3,4,1,1))
plot(results_CO2$pCO2, results_CO2$d13Cc, type = "l",
     xlab = expression("CO"[2]*" (ppmv)"),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(1000, -9.2, "a", cex = 1.5, font = 2)

plot(results_MAP$MAP, results_MAP$d13Cc, type = "l",
     xlab = "MAP (mm)",
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(100, -9, "b", cex = 1.5, font = 2)

plot(results_PCQ_pf$PCQ_pf, results_PCQ_pf$d13Cc, type = "l",
     xlab = expression("P"[PCQ]),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.1, -10.3, "c", cex = 1.5, font = 2)

plot(results_MAT$MAT, results_MAT$d13Cc, type = "l",
     xlab = expression(paste("MAT (", degree, "C)")),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(2, -11.3, "d", cex = 1.5, font = 2)

plot(results_PCQ_to$PCQ_to, results_PCQ_to$d13Cc, type = "l",
     xlab = expression(paste(Delta*"T (", degree, "C)")),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(2, -11.5, "e", cex = 1.5, font = 2)

plot(results_pore$pore, results_pore$d13Cc, type = "l",
     xlab = expression(rho),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.75, -10, "f", cex = 1.5, font = 2)

plot(results_f_R$f_R, results_f_R$d13Cc, type = "l",
     xlab = expression("f"[R]),
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.35, -3, "g", cex = 1.5, font = 2)

plot(results_tsc$tsc, results_tsc$d13Cc, type = "l",
     xlab = "tsc",
     ylab = expression(delta^"13"*"C"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.05, -8.9, "h", cex = 1.5, font = 2)
dev.off()

# d18Oc ----
png("figure/sens_d18c.png", width = 6.2, height = 6.2, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(3,4,1,1))
plot(results_MAP$MAP, results_MAP$d18Oc, type = "l",
     xlab = "MAP (mm)",
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(100, -7, "a", cex = 1.5, font = 2)

plot(results_PCQ_pf$PCQ_pf, results_PCQ_pf$d18Oc, type = "l",
     xlab = expression("P"[PCQ]),
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.9, -7, "b", cex = 1.5, font = 2)

plot(results_MAT$MAT, results_MAT$d18Oc, type = "l",
     xlab = expression(paste("MAT (", degree, "C)")),
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(27, -11, "c", cex = 1.5, font = 2)

plot(results_PCQ_to$PCQ_to, results_PCQ_to$d18Oc, type = "l",
     xlab = expression(paste(Delta*"T (", degree, "C)")),
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(3, -8.4, "d", cex = 1.5, font = 2)

plot(results_pore$pore, results_pore$d18Oc, type = "l",
     xlab = expression(rho),
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.75, -9, "e", cex = 1.5, font = 2)

plot(results_spre$spre, results_spre$d18Oc, type = "l",
     xlab = expression("f"[np]),
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.44, -10.5, "f", cex = 1.5, font = 2)

plot(results_tsc$tsc, results_tsc$d18Oc, type = "l",
     xlab = "tsc",
     ylab = expression(delta^"18"*"O"[c]*" (\u2030)"),
     cex.lab = 1.2, cex.axis = 1.2, mgp = c(2,.8,0))
text(.05, -9.2, "g", cex = 1.5, font = 2)
dev.off()





