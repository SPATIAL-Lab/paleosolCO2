## For gamma distribution
# plot function - precision instead of sd
pg = function(shp, rt){
  pre = rgamma(1e6, shp, rt)
  sd = sqrt(1/pre)
  plot(density(sd))
}
pg(10, 5)

## For beta distribution
pb = function(shp, rt){
  pre = rbeta(1e6, shp, rt)
  plot(density(pre))
}
pb(0.06 * 1e3 / (1 - 0.06), 1e3) 

## For normal distribution
pn = function(mean, sd){
  pre = rnorm(1e6, mean, sd)
  plot(density(pre))
}
pn(0, .5)
