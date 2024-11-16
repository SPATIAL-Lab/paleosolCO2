## For gamma distribution
# plot function - precision instead of sd
pg = function(shp, rt){
  pre = rgamma(1e6, shp, rt)
  sd = sqrt(1/pre)
  plot(density(sd))
}
pg(1 / 0.2 ^ 2, 1 / 0.2 ^ 2)

## For beta distribution
pb = function(shp, rt){
  pre = rbeta(1e6, shp, rt)
  plot(density(pre))
}
pb(0.55 * 500 / 0.45, 500) 

## For normal distribution
pn = function(mean, sd){
  pre = rnorm(1e6, mean, sd)
  plot(density(pre))
}
pn(-8, 3)
