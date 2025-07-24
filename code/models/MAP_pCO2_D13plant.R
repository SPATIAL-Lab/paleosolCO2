ctrl = function(){
  vars = list(
    "MAP" = 500,
    "pCO2" = 250,
    "d13Ca" = -7
  )
}

# photosynthetic fractionation
pf = function(vars){
  list2env(vars, environment())
  DD13water = 25.09 - (1.2 * (MAP + 975) / (27.2 + 0.04 * (MAP + 975))) # Diefendorf et al (2010)
  D13plant = (28.25 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))) - DD13water
  d13Cr = d13Ca - D13plant
  results = data.frame("MAP" = rep(MAP), "pCO2" = rep(pCO2), "d13Ca" = rep(d13Ca),
                       "DD13water" = rep(DD13water), "D13plant" = rep(D13plant),
                       "d13Cr" = rep(d13Cr))
  return(results)
}

# sensitivity test 
vars = ctrl()
vars$MAP = seq(200, 700, 100)
sens.p = pf(vars)
plot(sens.p$MAP, sens.p$D13plant)

vars = ctrl()
vars$pCO2 = seq(150, 500, 50)
sens.c = pf(vars)
plot(sens.c$pCO2, sens.c$D13plant)

vars = ctrl()
vars$d13Ca = seq(-8, -6, 0.5)
sens.d = pf(vars)
plot(sens.d$d13Ca, sens.d$d13Cr)
