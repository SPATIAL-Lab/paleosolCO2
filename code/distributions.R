## For gamma distribution
# plot function - precision instead of sd
pg = function(shp, rt){
  pre = rgamma(1e6, shp, rt)
  sd = sqrt(1/pre)
  plot(density(sd))
}
pg(10, 1e-4)

## For beta distribution
pb = function(shp, rt){
  pre = rbeta(1e6, shp, rt)
  plot(density(pre))
}
pb(0.11 * 500 / 0.89, 500) 

## For normal distribution
pn = function(mean, sd){
  pre = rnorm(1e6, mean, 1/sd^2)
  plot(density(pre))
}
pn(0.6, 1/0.5^2)
