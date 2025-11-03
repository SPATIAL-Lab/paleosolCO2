# paleosolCO2
PSM for paleosol-CO2 system


## code/ 
Scripts used to reproduce the study results of Da et al. (2025):

## code/models/ 
Scripts of the paleosol proxy system model (PSM) 
- **forward_model.R** The forward paleosol PSM 
- **multi_sample.R** The JAGS version of the paleosol PSM without the time-series model 
- **time_series.R** The JAGS version of the paleosol PSM with the time-series model 

## code/drivers/ 
Scripts used to perform Bayesian inversion 
- **ms_driver.R** Load data and run the inversion of the PSM without the time-series model 
- **ts_driver.R** Load data and run the inversion of the PSM with the time-series model 

## code/plot/ 
Scripts used to plot the figures in Da et al. (2025)

## data/
Houses any input data, including observation data (e.g., isotopic compositions of pedogenic carbonates) used for the JPI, and data used to plot figures.

## out/
Houses the output data of different JPI versions, including the mean and standard deviations of posterior distributions of various environmental parameters (e.g., CO<sub>2</sub>), as well as their Gelman-Rubin statistics (Rhat, n.eff).
