## For gamma distribution
# plot function - precision instead of sd
pg = function(shp, rt){
  pre = rgamma(1e6, shp, rt)
  sd = sqrt(1/pre)
  plot(density(sd))
}
pg(10, 100)

## For beta distribution
pb = function(shp, rt){
  pre = rbeta(1e6, shp, rt)
  plot(density(pre))
}
pb(0.11 * 1e4 / 0.89, 1e4) 

## For normal distribution
pn = function(mean, sd){
  pre = rnorm(1e6, mean, 1/sd^2)
  plot(density(pre))
}
pn(-6.5, 1/1^2)
