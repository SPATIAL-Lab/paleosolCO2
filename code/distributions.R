## For gamma distribution
# plot function - precision instead of sd
pg = function(shp, rt){
  pre = rgamma(1e6, shp, rt)
  sd = sqrt(1/pre)
  plot(density(sd))
}
pg(100, 100)
mean = 0.86
var = .14
rt = mean / var
shp = mean * rt
pg(shp, rt)

## For beta distribution
pb = function(shp, rt){
  pre = rbeta(1e6, shp, rt)
  plot(density(pre))
}
mean = .86
var_max = mean * (1 - mean)
var = .12
shp = mean * (mean * (1 - mean) / var - 1)
rt = (1 - mean) * (mean * (1 - mean) / var - 1)
pb(6.956, 1.139) 
pb(0.55 * 10 / 0.45, 10)
## For normal distribution
pn = function(mean, sd){
  pre = rnorm(1e6, mean, sd)
  plot(density(pre))
}
pn(300, 50)
