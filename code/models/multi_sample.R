model{

  # Data model ----
  for(i in 1:length(d13Cc.ai)){
    d13Cc.obs[i, 1] ~ dnorm(d13Cc[d13Cc.ai[i]], d13Cc.pre[i])
    d13Cc.pre[i] = 1 / d13Cc.obs[i, 2] ^ 2
  }
  
  for(i in 1:length(d18Oc.ai)){
    d18Oc.obs[i, 1] ~ dnorm(d18Oc[d18Oc.ai[i]], d18Oc.pre[i])
    d18Oc.pre[i] = 1 / d18Oc.obs[i, 2] ^ 2
  }

  # for(i in 1:length(d13Co.ai)){
  #   d13Co.obs[i, 1] ~ dnorm(d13Co[d13Co.ai[i]], d13Co.pre[i])
  #   d13Co.pre[i] = 1 / d13Co.obs[i, 2] ^ 2
  # }
    
  for (i in 1:length(d13Ca.ai)) {
    d13Ca.obs[i, 1] ~ dnorm(d13Ca[d13Ca.ai[i]], d13Ca.pre[i])
    d13Ca.pre[i] = 1 / d13Ca.obs[i, 2] ^ 2
  }
  
  # for(i in 1:length(D47c.ai)){
  #   D47c.obs[i, 1] ~ dnorm(D47c[D47c.ai[i]], D47c.pre[i])
  #   D47c.pre[i] = 1 / D47c.obs[i, 2] ^ 2
  # }
  
  for(i in 1:length(ai)){  
    # Soil carbonate ----
    ## Depth to carbonate formation based on Retallack (2005) data, meters
    z[i] = (0.093 * MAP[i] + 13.12)
    z_m[i] = z[i] / 100
    
    ## Soil temperatures at depth z
    Tsoil[i] = MAT[i] + (PCQ_to[i] * sin(2 * 3.1415 * tsc[i] - z[i] / d)) / 
      exp(z[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
    Tair_PCQ[i] = MAT[i] + PCQ_to[i]
    PET_PCQ_D.1[i] = ifelse(ha[i] < 0.5, 
                       0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                       0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50))
    PET_PCQ_D[i] = max(PET_PCQ_D.1[i], 0.01)
    PET_PCQ[i] = PET_PCQ_D[i] * 90
    
    ## AET in mm/quarter from Budyko curve - Pike (1964)
    PPCQ[i] = MAP[i] * PCQ_pf[i]
    AET_PCQ[i] = PPCQ[i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[i])))) ^ 2)))
    
    ## Carbon isotopes ----
    ### Free air porosity
    FAP.1[i] = min((pore - ((PPCQ[i] - AET_PCQ[i]) / (L * 10 * pore))), pore - 0.05)
    FAP[i] = max(FAP.1[i], 0.01)
    
    ### Soil respiration rate 
    R_PCQ_D_m1[i] = 1.25 * exp(0.05452 * Tair_PCQ[i]) * PPCQ[i] / (127.77 + PPCQ[i])
    R_PCQ_D[i] = R_PCQ_D_m1[i] * f_R[i] # (gC/m2/d)
    
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
    # DD13_water[i] = 25.09 - 1.2 * (MAP[i] + 975) / (27.2 + 0.04 * (MAP[i] + 975))
    # D13C_plant[i] = (28.26 * 0.22 * (pCO2[i] + 23.9)) / (28.26 + 0.22 * (pCO2[i] + 23.9)) - DD13_water[i] # schubert & Jahren (2015)
    # d13Cr[i] = d13Ca[i] - D13C_plant[i]
    d13Co[i] = d13Cr[i] + 1 + SOM.frac
    
    ### d13C of pedogenic carbonate
    d13Cs[i] = (pCO2[i] * d13Ca[i] + S_z[i] * (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[i])
    d13Cc[i] = ((1 + (11.98 - 0.12 * Tsoil[i]) / 1000) * (d13Cs[i] + 1000)) - 1000
    
    ## Oxygen isotopes ----
    ### Rainfall isotopes
    R18.p[i] = (d18.p[i] / 1000 + 1) * R18.VSMOW
    
    ### Equilibrium fractionation (Horita and Wesolowski 1994)
    alpha18.eq[i] = 1 / exp(((1.137e6 / (Tsoil.K[i] ^ 2) - 0.4156e3/Tsoil.K[i] - 2.0667) /1000))
    
    ### Atmospheric water vapor isotopes
    R18.a[i] = R18.p[i] * alpha18.eq[i]
    
    ### Soil evaporation from AET
    E1[i] = ETR[i] * AET_PCQ[i]
    E[i] = max(E1[i], 1) # minimum of 1 mm
    E_s[i] = E[i] / (1000 * 90 * 24 * 3600) # soil evaporation rate in m/sec
    
    ### Water vapor diffusivity
    es[i] = (0.611 * exp(17.502 * Tsoil[i] / (Tsoil[i] + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
    N.sat[i] = 0.01802 * es[i] / (Rgas * Tsoil.K[i]) # saturated water vapor concentration at a given temperature
    z.bar[i] = N.sat[i] * Dv.soil / (E_s[i] * rho) # penetration depth (m)
    z.ef1[i] = (1 - ha[i]) * z.bar[i] # the thickness of the water vapor phase region (m)
    z.ef[i] = max(z.ef1[i], 1e-10)
    
    ### Liquid water diffusivity (m2/s) (Easteal 1984)
    Dl[i] = exp(1.6766 + 1.6817 * (1000 / Tsoil.K[i]) - 0.5773 * (1000 / Tsoil.K[i]) ^ 2) * 1e-9 
    Dl.soil[i] = Dl[i] * pore * tort # effective diffusivity of liquid water (m2/s)
    z.hat[i] = Dl.soil[i] / E_s[i] # the decay length (mean penetration depth)
    
    ### The evaporation front
    h.ef[i] = ha[i] + z.ef[i] / z.bar[i] # humidity at the evaporation front
    R18.ef[i] = (alpha18.diff * R18.p[i] * (z.ef[i] / z.bar[i]) + 
                   ha[i] * R18.a[i]) / (h.ef[i] * alpha18.eq[i]) # isotopic composition at the evaporation front
    
    ### Isotope composition of soil water at depth z
    hs[i] = min(ha[i] + z_m[i] / z.bar[i], 1)
    z.f[i] = (pore / a.theta) * log(z_m[i] / z.ef[i]) # the modified depth function
    R18.s[i] = ifelse(z_m[i] <= z.ef[i], 
                      (alpha18.diff * R18.p[i] * z_m[i] / z.bar[i] + ha[i] * R18.a[i]) / 
                        (hs[i] * alpha18.eq[i]),
                      (R18.ef[i] - R18.p[i]) * exp(-z.f[i] / z.hat[i]) + R18.p[i])
    d18O.s[i] = ((R18.s[i] / R18.VSMOW) - 1) * 1000
    
    ### Isotope composition of soil carbonate
    alpha18_c_w_eq[i] = exp((1.61e4 / Tsoil.K[i] - 24.6) / 1000) # Wostbrock (2020)
    R18.c[i] = R18.s[i] * alpha18_c_w_eq[i]
    d18Oc[i] = (R18.c[i] / R18.VPDB - 1) * 1000
    # D47c[i] = 0.0391e6 / Tsoil.K[i] ^ 2 + 0.15
  }
  
  for(i in 1:length(ai)){
    # Time dependent variables ----
    ## Primary environmental ----
    pCO2[i] ~ dunif(100, 400) # atmospheric CO2 mixing ratio
    MAT[i] ~ dunif(0, 20)
    PCQ_to[i] ~ dunif(10, 16)
    MAP[i] ~ dunif(150, 750) # mean annual terrestrial site precipitation, mm
    PCQ_pf[i] ~ dnorm(0.55, 1 / 0.1 ^ 2)T(0.3, 0.8) # PCQ precipitation fraction
    Tair_OOS[i] = (4 * MAT[i] - Tair_PCQ[i]) / 3
    d18.p[i] ~ dnorm(-15 + 0.58 * (Tair_PCQ[i] * PCQ_pf[i] + Tair_OOS[i] * (1 - PCQ_pf[i])), 1 / 1 ^ 2)
    
    ## Secondary soil ----
    tsc[i] ~ dbeta(0.29 * 1000 / 0.71, 1000) # seasonal offset of PCQ for thermal diffusion
    h_m[i] = 0.25 + 0.7 * (PPCQ[i] / 900)
    ha[i] ~ dbeta(h_m[i] * 100 / (1 - h_m[i]), 100) # PCQ atmospheric humidity
    f_R[i] ~ dbeta(0.11 * 500 / 0.89, 500) # ratio of PCQ to mean annual respiration rate
    d13Ca[i] ~ dunif(-8, -5) # Atmospheric d13C, ppt
    ETR[i] ~ dbeta(0.06 * 1000 / 0.94, 1000) # Soil evaporation / AET
    d13Cr[i] ~ dunif(-27, -22)
  }
  
  # Not time dependent ----
  ## Terrestrial ----
  lat = 36 # terrestrial site latitude
  Ra = 42.608 - 0.3538 * abs(lat) # total radiation at the top of the atmosphere
  Rs = Ra * 0.16 * sqrt(12) # daily temperature range assumed to be 12
  L ~ dgamma(40, 1) # mean rooting depth, cm
  k = L / 2 / log(2)   # Respiration characteristic production depth (cm) - Quade (2007)
  pore ~ dbeta(0.45 * 100 / 0.55, 100)T(0.06,) # soil porosity
  tort ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
  Dv.soil = Dv.air * tort * (pore - 0.05) # effective diffusivity of water vapor in soil (m2/s)
  SOM.frac ~ dunif(-0.5, 0.5)

  ## Constants ----
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  ### Isotope ratio constants
  R13.VPDB = 0.011237
  R18.VSMOW = 0.0020052
  R18.VPDB = 0.0020672
  alpha18.diff = 1.028489
  a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
  Rgas = 8.314462 # gas constant
  rho = 1000 # liquid water density (kg/m3)
  Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)

}
  