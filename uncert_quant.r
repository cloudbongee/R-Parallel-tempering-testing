rm(list = ls())

library(Rcpp)
library(bench)
library(here)
library(ggplot2)


sourceCpp(here::here("maximin_geom_seq.cpp")) ## load maximin with geometric sequence temperatures
sourceCpp(here::here("maxpro.cpp"))           ## load maxpro


## Simulation:

## By Jaime Meyer Beilis Michel
## At Bloomington, Indiana, 7.16.2025
## Optimization of parallel tempering for Latin Hypercube Designs

## Uncertainity quantification using PT
## Piston motion with a cylinder. Cycle time in seconds

## Experiment from SFDesign: An R package for Space-Filling Designs
## Shangkun Wang∗ , Weijun Xie and V. Roshan Joseph Georgia Institute of Technology

## Which references the experiment from:
## Kenett R, Zacks S (1998). 
## Modern industrial statistics: design and control of quality and reliability. 
## Duxbury Press, Pacific Grove, C



## normalization for an n row design
norm_design <- function(X){
  (apply(X, 2, rank) - 0.5) / nrow(X)
}

## maxpro temperature scheduling function
maxpro_temps <- function(n,k,M){
  crit1 = 1.0 / (n-1)
  crit2 = (1.0 / ((n - 1)^(k -1) * (n-2)))^(1.0/k)
  delta0 = crit2 - crit1;
  T0 = -delta0 * (1.0/log(0.99))
  L <- sapply(0:M-1, function(i) T0^(M - i))
}



## Takes a [0,1] normalized latin hypercube design of n = 50, k = 7.
## Sets the variables to their ranges
cycle_time <- function(LHD, norm = T, wt_upper =  60, wt_lower = 30,area_lower =  0.005,area_upper = 0.020,v_gas_lwr = 0.002,v_gas_upp = 0.01,k_lower = 1000,k_upper = 5000,p0_lower = 90000,p0_upper = 110000,temp_ambience_lower = 290,temp_ambience_upper = 296,filling_gas_lower = 340,filling_gas_upper = 360){
  if(norm){ ## normalize?
    X <- norm_design(LHD)
  }else{
    X <- LHD
  }

  X <- as.data.frame(X)
  names(X) = c("piston_wt", "surf_area", "init_gas_vol", "spring_coef", "atmospheric_press", "ambient_temp", "filling_gas_temp")


  ## Domain arrangement
  X[,"piston_wt"] = X[,"piston_wt"] * (wt_upper - wt_lower) + wt_lower                                            # [30:60]
  X[,"surf_area"] = X[,"surf_area"] * (area_upper - area_lower) + area_lower                                      # [0.005 : 0.020]
  X[,"init_gas_vol"] = X[,"init_gas_vol"] * (v_gas_upp - v_gas_lwr) + v_gas_lwr                                   # [0.002 : 0.010]
  X[,"spring_coef"] = X[,"spring_coef"] * (k_upper - k_lower) + k_lower                                           # [1000 : 5000]
  X[,"atmospheric_press"] = X[,"atmospheric_press"] * (p0_upper - p0_lower) + p0_lower                            # [90000 : 110000]
  X[,"ambient_temp"] = X[,"ambient_temp"] * (temp_ambience_upper - temp_ambience_lower) + temp_ambience_lower     # [290 : 296]
  X[,"filling_gas_temp"] = X[,"filling_gas_temp"] * (filling_gas_upper - filling_gas_lower) + filling_gas_lower   # [340 : 360]

  ## Computation
  A <- X[,"atmospheric_press"] * X[,"surf_area"] + 19.62 * X[,"piston_wt"] - X[,"spring_coef"] * X[,"init_gas_vol"]/X[,"surf_area"]
  V <- X[,"surf_area"]/(2 * X[,"spring_coef"]) * (sqrt(A^2 + 4 * X[,"spring_coef"]* X[,"atmospheric_press"] * X[,"init_gas_vol"]/ X[,"filling_gas_temp"] * X[,"ambient_temp"]) - A)
  C <- 2 * pi * sqrt(X[,"piston_wt"]/(X[,"spring_coef"] + X[,"surf_area"]^2 * (X[,"atmospheric_press"] * X[,"init_gas_vol"] * X[,"ambient_temp"]/ (X[,"filling_gas_temp"] * V^2)) ))

  return ( list(
    input = X,
    output = C,
    mean = mean(C),
    var = var(C),
    sd = sd(C)
  ) )}



## parallelize: (50 times) TODO 
## simulation with pre_filled values
sim_piston_cycles_time <- function(){
  temps = maxpro_temps(50,7,5)
  X1 <- maximinLHD_geom(50,7,8,20,10,1000000,1000,0.1)$design
  X2 <- pt_maxpro_lhd(50,7,5,100000,10,2500, temps)$design
  randX <- matrix(runif(50 * 7, min = 0, max = 1), nrow = 50, ncol = 7)
  maximin_sim <- cycle_time(X1)
  maxpro_sim <- cycle_time(X2, norm= F)
  random <- cycle_time( randX,norm = F)
  
  box <- ggplot(data = data.frame(
    second_cycle = c(maximin_sim$output, maxpro_sim$output, random$output),
    type = factor(c(rep("maximin", nrow(X1)), rep("maxpro", nrow(X2)), rep("random_unif", nrow(randX))))
  ), aes(x = second_cycle)) + geom_boxplot(fill = c("#f67e7d", "#843b62", "#621940")) + facet_wrap(~type) + coord_flip()

  ggsave(filename = "sim_cycle_time_boxplot.png", plot = box,width = 6,height = 4,units = "in",dpi = 300)
  
  res =  list(maximin = maximin_sim, maxpro = maxpro_sim)

  saveRDS( # save maxpro
    res,
    file = here::here("simulated_cycle_time.rds")
  )

  return(res)
  
}

result <- sim_piston_cycles_time()
