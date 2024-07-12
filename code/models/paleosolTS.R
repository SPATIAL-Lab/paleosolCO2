model{

  # Data model ----
  for(i in 1:length(d13Cc.ai)){
    d13Cc.obs[i, 1] ~ dnorm(d13Cc[d13Cc.ai[i]], 1 / d13Cc.obs[i, 2] ^ 2)
  }

  for(i in 1:length(d18Oc.ai)){
    d18Oc.obs[i, 1] ~ dnorm(d18Oc[d18Oc.ai[i]], 1 / d18Oc.obs[i, 2] ^ 2)
  }

  for(i in 1:length(d13Cno.ai)){
    d13Cno.obs[i, 1] ~ dnorm(d13Ca[d13Cno.ai[i]] - D13C_plant[d13Cno.ai[i]], 
                            1 / d13Cno.obs[i, 2] ^ 2)
  }

  for(i in 1:length(d13Cbo.ai)){
    d13Cbo.obs[i, 1] ~ dnorm(d13Ca[d13Cbo.ai[i]] - D13C_plant[d13Cbo.ai[i]], 
                            1 / d13Cbo.obs[i, 2] ^ 2)
  }
  
  for(i in 1:length(st.ai)){
    st.obs[i, 1] ~ dnorm(Tsoil[st.ai[i]], 1 / st.obs[i, 2] ^ 2)
  }
  
  for(i in 1:length(map.ai)){
    map.obs[i, 1] ~ dnorm(MAP[si[map.ai[i]], map.ai[i]], 
                          1 / map.obs[i, 2] ^2)
  }
  
  for(i in 1:length(ai)){  
    # Soil carbonate ----
    ## Depth to carbonate formation based on Retallack (2005) data, meters
    z.mean[i] = (0.093 * MAP[si[i], i] + 13.12)
    ### Gamma rate
    z.beta[i] = z.mean[i] / (22 ^ 2)
    ### Gamma shape
    z.alpha[i] = z.mean[i] * z.beta[i]
    z[i] ~ dgamma(z.alpha[i], z.beta[i])
    z_m[i] = z[i] / 100
    
    ## Soil temperatures at depth z
    Tsoil[i] = MAT[si[i], i] + (PCQ_to[si[i], i] * sin(2 * 3.1415 * tsc[si[i], i] - z[i] / d)) / 
      exp(z[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
    Tair_PCQ[i] = MAT[si[i], i] + PCQ_to[si[i], i]
    PET_PCQ_D.1[i] = ifelse(ha[si[i], i] < 0.5, 
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha[si[i], i]) / 0.7)),
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50))
    PET_PCQ_D[i] = max(PET_PCQ_D.1[i], 0.01)
    PET_PCQ[i] = PET_PCQ_D[i] * 90
    
    ## AET in mm/quarter from Budyko curve - Pike (1964)
    AET_var[i] ~ dgamma(1 / 0.2 ^ 2, 1 / 0.2 ^ 2) # noise parameter - Gentine (2012)
    AET_PCQ[i] = PPCQ[si[i], i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[si[i], i])) * AET_var[i])) ^ 2)))
    
    ## Carbon isotopes ----
    ### Free air porosity
    FAP.1[i] = min((pore - ((PPCQ[si[i], i] - AET_PCQ[i]) / (L * 10 * pore))), pore - 0.05)
    FAP[i] = max(FAP.1[i], 0.01)
    
    ### Soil respiration rate 
    R_PCQ_D_m1[i] = 1.25 * exp(0.05452 * TmPCQ[si[i], i]) * PPCQ[si[i], i] / (127.77 + PPCQ[si[i], i])
    R_PCQ_D_m[i] = R_PCQ_D_m1[i] * f_R[si[i], i] # (gC/m2/d)
    R_beta[i] = R_PCQ_D_m[i] / (R_PCQ_D_m[i] * 0.1) ^ 2
    R_alpha[i] = R_PCQ_D_m[i] * R_beta[i]
    R_PCQ_D[i] ~ dgamma(R_alpha[i], R_beta[i])
    
    ### Convert to molC/cm3/s
    R_PCQ_D.1[i] = R_PCQ_D[i] / (12.01 * 100 ^ 2)  # from gC/m2/d to molC/cm2/d
    R_PCQ_S[i] = R_PCQ_D.1[i] / (24 * 3600)  # molC/ cm2 / s
    R_PCQ_S_0[i]= R_PCQ_S[i] / (L * pore) # Quade et al. (2007)
    
    ### CO2 diffusion
    Dair[i] = 0.1369 * (Tsoil.K[i] / 273.15) ^ 1.958
    DIFC[i] = FAP[i] * tort * Dair[i]
    
    ### S(z)
    S_z_mol[i] = k ^ 2 * R_PCQ_S_0[i] / DIFC[i] * (1 - exp(-z[i] / k)) # (mol/cm3)
    S_z[i] = S_z_mol[i] * (0.08206 * Tsoil.K[i] * 10^9) # ppmv
    
    ### d13C of soil-respired CO2
    DD13_water[i] = 25.09 - 1.2 * (MAP[si[i], i] + 975) / (27.2 + 0.04 * (MAP[si[i], i] + 975))
    D13C_plant[i] = (28.26 * 0.22 * (pCO2[i] + 23.9)) / (28.26 + 0.22 * (pCO2[i] + 23.9)) - DD13_water[i] # schubert & Jahren (2015)
    D13C_off[i] ~ dnorm(0, 1 / 2 ^ 2) # Noise term
    d13Cr[i] = d13Ca[i] - D13C_plant[i] + D13C_off[i]
    
    ### d13C of pedogenic carbonate
    d13Cs[i] = (pCO2[i] * d13Ca[i] + S_z[i] * (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[i])
    d13Cc[i] = ((1 + (11.98 - 0.12 * Tsoil[i]) / 1000) * (d13Cs[i] + 1000)) - 1000
    
    ## Oxygen isotopes ----
    ### Rainfall isotopes
    R18.p[i] = (d18.p[si[i], i] / 1000 + 1) * R18.VSMOW
    
    ### Equilibrium fractionation (Horita and Wesolowski 1994)
    alpha18.eq[i] = 1 / exp(((1.137e6 / (Tsoil.K[i] ^ 2) - 0.4156e3/Tsoil.K[i] - 2.0667) /1000))
    
    ### Atmospheric water vapor isotopes
    R18.a[i] = R18.p[i] * alpha18.eq[i]
    
    ### Soil evaporation from AET
    E1[i] = ETR[si[i], i] * AET_PCQ[i]
    E[i] = max(E1[i], 1) # minimum of 1 mm
    E_s[i] = E[i] / (1000 * 90 * 24 * 3600) # soil evaporation rate in m/sec
    
    ### Water vapor diffusivity
    es[i] = (0.611 * exp(17.502 * Tsoil[i] / (Tsoil[i] + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
    N.sat[i] = 0.01802 * es[i] / (Rgas * Tsoil.K[i]) # saturated water vapor concentration at a given temperature
    z.bar[i] = N.sat[i] * Dv.soil / (E_s[i] * rho) # penetration depth (m)
    z.ef1[i] = (1 - ha[si[i], i]) * z.bar[i] # the thickness of the water vapor phase region (m)
    z.ef[i] = max(z.ef1[i], 1e-10)
    
    ### Liquid water diffusivity (m2/s) (Easteal 1984)
    Dl[i] = exp(1.6766 + 1.6817 * (1000 / Tsoil.K[i]) - 0.5773 * (1000 / Tsoil.K[i]) ^ 2) * 1e-9 
    Dl.soil[i] = Dl[i] * pore * tort # effective diffusivity of liquid water (m2/s)
    z.hat[i] = Dl.soil[i] / E_s[i] # the decay length (mean penetration depth)
    
    ### The evaporation front
    h.ef[i] = ha[si[i], i] + z.ef[i] / z.bar[i] # humidity at the evaporation front
    R18.ef[i] = (alpha18.diff * R18.p[i] * (z.ef[i] / z.bar[i]) + 
                   ha[si[i], i] * R18.a[i]) / (h.ef[i] * alpha18.eq[i]) # isotopic composition at the evaporation front
    
    ### Isotope composition of soil water at depth z
    hs[i] = min(ha[si[i], i] + z_m[i] / z.bar[i], 1)
    z.f[i] = (pore / a.theta) * log(z_m[i] / z.ef[i]) # the modified depth function
    R18.s[i] = ifelse(z_m[i] <= z.ef[i], 
                      (alpha18.diff * R18.p[i] * z_m[i] / z.bar[i] + ha[si[i], i] * R18.a[i]) / 
                        (hs[i] * alpha18.eq[i]),
                      (R18.ef[i] - R18.p[i]) * exp(-z.f[i] / z.hat[i]) + R18.p[i])
    d18O.s[i] = ((R18.s[i] / R18.VSMOW) - 1) * 1000
    
    ### Isotope composition of soil carbonate
    alpha18_c_w_eq[i] = exp((1.61e4 / Tsoil.K[i] - 24.6) / 1000) # Wostbrock (2020)
    R18.c[i] = R18.s[i] * alpha18_c_w_eq[i]
    d18Oc[i] = (R18.c[i] / R18.VPDB - 1) * 1000
    D47c[i] = 0.0417e6 / Tsoil.K[i] ^ 2 + 0.139
  }
  
  # Time dependent variables, time series ----

  for(i in 2:length(ai)){
    pCO2[i] = max(pCO2.p[i], 100)
    pCO2.p[i] = pCO2[1] * (1 + pCO2.eps[i])
    pCO2.eps[i] ~ dnorm(pCO2.eps[i - 1] * (pCO2.phi ^ dt[i - 1]), pCO2.pc[i])
    pCO2.pc[i] = pCO2.tau * ((1 - pCO2.phi ^ 2) / (1 - pCO2.phi ^ (2 * dt[i - 1])))
    
    d13Ca[i] = d13Ca[1] + d13Ca.eps[i]
    d13Ca.eps[i] ~ dnorm(d13Ca.eps[i - 1] * (d13Ca.phi ^ dt[i - 1]), d13Ca.pc[i])
    d13Ca.pc[i] = d13Ca.tau * ((1 - d13Ca.phi ^ 2) / (1 - d13Ca.phi ^ (2 * dt[i - 1])))

    GMT[i] = GMT[1] + GMT.eps[i]
    GMT.eps[i] ~ dnorm(GMT.eps[i - 1] * (GMT.phi ^ dt[i - 1]), GMT.pc[i])
    GMT.pc[i] = GMT.tau * ((1 - GMT.phi ^ 2) / (1 - GMT.phi ^ (2 * dt[i - 1])))
  }
  
  pCO2[1] ~ dunif(500, 2000) # atmospheric CO2 mixing ratio
  pCO2.eps[1] = 0
  
  d13Ca[1] ~ dnorm(-9, 1 / 1 ^ 2) # Atmospheric d13C, ppt
  d13Ca.eps[1] = 0
  
  GMT[1] ~ dnorm(25, 1 / 3 ^ 2)
  GMT.eps[1] = 0
  
  for(j in 1:nsites){
    for(i in 2:length(ai)){
      
      ## Derived values ----
      d18.p[j, i] ~ dnorm(-15 + 0.58 * (MAT[j, i] * (1 - PCQ_pf[j, i]) +
                                       TmPCQ[j, i] * PCQ_pf[j, i]), 1 / 1 ^ 2) # Precipitation d18O, ppt
      PPCQ[j, i] = MAP[j, i] * PCQ_pf[j, i] 
      TmPCQ[j, i] = MAT[j, i] + PCQ_to[j, i] 
      MAT[j, i] = GMT[i] + MATo[j, i]
      
      ## Primary environmental ----
      MATo[j, i] = MATo[j, 1] + MATo.eps[j, i]
      MATo.eps[j, i] ~ dnorm(MATo.eps[j, i - 1] * (MATo.phi ^ dt[i - 1]), MATo.pc[j, i])
      MATo.pc[j, i] = MATo.tau * ((1 - MATo.phi ^ 2) / (1 - MATo.phi ^ (2 * dt[i - 1])))
      
      PCQ_to[j, i] = PCQ_to[j, 1] + PCQ_to.eps[j, i]
      PCQ_to.eps[j, i] ~ dnorm(PCQ_to.eps[j, i - 1] * (PCQ_to.phi ^ dt[i - 1]), PCQ_to.pc[j, i])
      PCQ_to.pc[j, i] = PCQ_to.tau * ((1 - PCQ_to.phi ^ 2) / (1 - PCQ_to.phi ^ (2 * dt[i - 1])))
      
      MAP[j, i] = MAP[j, 1] * (1 + MAP.eps[j, i])
      MAP.eps[j, i] ~ dnorm(MAP.eps[j, i - 1] * (MAP.phi ^ dt[i - 1]), MAP.pc[j, i])T(-1,)
      MAP.pc[j, i] = MAP.tau * ((1 - MAP.phi ^ 2) / (1 - MAP.phi ^ (2 * dt[i - 1])))
      
      PCQ_pf[j, i] = max(min(PCQ_pf.p[j, i], 0.25), 0.01)
      PCQ_pf.p[j, i] = PCQ_pf[j, 1] * (1 + PCQ_pf.eps[j, i])
      PCQ_pf.eps[j, i] ~ dnorm(PCQ_pf.eps[j, i - 1] * (PCQ_pf.phi ^ dt[i - 1]), PCQ_pf.pc[j, i])
      PCQ_pf.pc[j, i] = PCQ_pf.tau * ((1 - PCQ_pf.phi ^ 2) / (1 - PCQ_pf.phi ^ (2 * dt[i - 1])))
      
      ## Secondary soil ----
      tsc[j, i] = tsc[j, 1] + tsc.eps[j, i]
      tsc.eps[j, i] ~ dnorm(tsc.eps[j, i - 1] * (tsc.phi ^ dt[i - 1]), tsc.pc[j, i])
      tsc.pc[j, i] = tsc.tau * ((1 - tsc.phi ^ 2) / (1 - tsc.phi ^ (2 * dt[i - 1])))
      
      ha[j, i] = max(min(ha.p[j, i], 0.95), 0.1)
      ha.p[j, i] = ha[j, 1] * (1 + ha.eps[j, i])
      ha.eps[j, i] ~ dnorm(ha.eps[j, i - 1] * (ha.phi ^ dt[i - 1]), ha.pc[j, i])
      ha.pc[j, i] = ha.tau * ((1 - ha.phi ^ 2) / (1 - ha.phi ^ (2 * dt[i - 1])))
      
      f_R[j, i] = max(min(f_R.p[j, i], 0.3), 0.02)
      f_R.p[j, i] = f_R[j, 1] * (1 + f_R.eps[j, i])
      f_R.eps[j, i] ~ dnorm(f_R.eps[j, i - 1] * (f_R.phi ^ dt[i - 1]), f_R.pc[j, i])
      f_R.pc[j, i] = f_R.tau * ((1 - f_R.phi ^ 2) / (1 - f_R.phi ^ (2 * dt[i - 1])))
      
      ETR[j, i] = max(min(ETR.p[j, i], 0.1), 0.01)
      ETR.p[j, i] = ETR[j, 1] * (1 + ETR.eps[j, i])
      ETR.eps[j, i] ~ dnorm(ETR.eps[j, i - 1] * (ETR.phi ^ dt[i - 1]), ETR.pc[j, i])
      ETR.pc[j, i] = ETR.tau * ((1 - ETR.phi ^ 2) / (1 - ETR.phi ^ (2 * dt[i - 1])))
    }
    
    # Time dependent variables, initial conditions ----
    ## Derived values ----
    d18.p[j, 1] = -17 + 0.58 * (MAT[j, 1] * (1 - PCQ_pf[j, 1]) +
                               TmPCQ[j, 1] * PCQ_pf[j, 1]) + MAP[j, 1] / 250 # Precipitation d18O, ppt
    PPCQ[j, 1] = MAP[j, 1] * PCQ_pf[j, 1] 
    TmPCQ[j, 1] = MAT[j, 1] + PCQ_to[j, 1] 
    MAT[j, 1] = GMT[1] + MATo[j, 1]
    
    ## Primary environmental ----
    MATo[j, 1] ~ dunif(-15, 15) # terrestrial temperature site offset, C
    MATo.eps[j, 1] = 0
    PCQ_to[j, 1] ~ dunif(-2, 14) # PCQ temperature offset, C
    PCQ_to.eps[j, 1] = 0
    MAP[j, 1] ~ dunif(200, 1400) # mean annual terrestrial site precipitation, mm
    MAP.eps[j, 1] = 0
    PCQ_pf[j, 1] ~ dunif(0.02, 0.25) # PCQ precipitation fraction
    PCQ_pf.eps[j, 1] = 0
    
    ## Secondary soil ----
    tsc[j, 1] ~ dbeta(0.25 * 1000 / 0.75, 1000) # seasonal offset of PCQ for thermal diffusion
    tsc.eps[j, 1] = 0
    ha[j, 1] ~ dunif(0.1, 0.7)
    ha.eps[j, 1] = 0
    f_R[j, 1] ~ dbeta(0.15 * 500 / 0.85, 500) # ratio of PCQ to mean annual respiration rate
    f_R.eps[j, 1] = 0
    ETR[j, 1] ~ dbeta(0.06 * 1000 / 0.94, 1000) # Soil evaporation / AET
    ETR.eps[j, 1] = 0
  }
  
  # Time dependent variables, ts parameters ----
  pCO2.tau ~ dgamma(10, 0.1)
  pCO2.phi ~ dbeta(5, 2)
  
  GMT.tau ~ dgamma(20, 5)
  GMT.phi ~ dbeta(5, 2)
  
  MATo.tau ~ dgamma(20, 5)
  MATo.phi ~ dbeta(5, 2)
  
  PCQ_to.tau ~ dgamma(100, 5)
  PCQ_to.phi ~ dbeta(5, 2)
  
  MAP.tau ~ dgamma(10, 0.1)
  MAP.phi ~ dbeta(5, 2)
  
  PCQ_pf.tau ~ dgamma(10, 0.01)
  PCQ_pf.phi ~ dbeta(5, 2)
  
  tsc.tau ~ dgamma(10, 0.01)
  tsc.phi ~ dbeta(5, 2)
  
  ha.tau ~ dgamma(10, 0.01)
  ha.phi ~ dbeta(5, 2)
  
  f_R.tau ~ dgamma(10, 0.001)
  f_R.phi ~ dbeta(5, 2)
  
  d13Ca.tau ~ dgamma(10, 1)
  d13Ca.phi ~ dbeta(5, 2)
  
  ETR.tau ~ dgamma(10, 0.001)
  ETR.phi ~ dbeta(5, 2)

  
  # Not time dependent ----
  lat = 30 # terrestrial site latitude
  Ra = 42.608 - 0.3538 * abs(lat) # total radiation at the top of the atmosphere
  Rs = Ra * 0.16 * sqrt(12) # daily temperature range assumed to be 12
  L ~ dgamma(50, 1) # mean rooting depth, cm
  k = L / 2 / log(2)   # Respiration characteristic production depth (cm) - Quade (2007)
  pore ~ dbeta(0.35 * 100 / 0.65, 100)T(0.06,) # soil porosity
  tort ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
  Dv.soil = Dv.air * tort * (pore - 0.05) # effective diffusivity of water vapor in soil (m2/s)
  
  ## Constants ----
  R13.VPDB = 0.011237
  R18.VSMOW = 0.0020052
  R18.VPDB = 0.0020672
  alpha18.diff = 1.028489
  a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
  Rgas = 8.314462 # gas constant
  rho = 1000 # liquid water density (kg/m3)
  Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  
}


