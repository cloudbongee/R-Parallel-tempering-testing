library(Rcpp)
library(bench)
library(SLHD)
library(here)
sourceCpp(here::here("maximin_geom_seq.cpp"))
sourceCpp(here::here("maxpro.cpp"))
sourceCpp(here::here("criteria.cpp"))

## bench_mean_sd_measures returns a list with the bench, the mean, and the standard
## deviation of a given function (fun) and its parameters, with alpha, on the standards of a phi with p 15
bench_mean_sd_measures <- function(R,fun,n,k,M,p,Nswap,Nmax,tolerance,alpha){

  ## 1 to R measures
  measures_geom <- lapply(1:R , function(i) { phi(fun(n,k,M,p,Nswap,Nmax,tolerance,alpha)$design,n,k,15) })
  measures_geom <- unlist(measures_geom)

  ## return mean, sd, benchmark an extra run.
  return(list(
    mean = mean(measures_geom),
    sd = sd(measures_geom),
    benchmark = bench::mark(fun(n,k,M,p,Nswap,Nmax,tolerance,alpha))
  ))

}

phi(maximinSLHD(t = 1, m = 27, k = 13, power= 20)$Design, 27,13,15)

## example with very good parameters
d <- bench_mean_sd_measures(100, maximinLHD_geom,27,13,8,20,10,1000000,5000,0.1) 
d$benchmark
d$mean
d$sd

# [1] 0.09044156
# [1] 0.09066566
