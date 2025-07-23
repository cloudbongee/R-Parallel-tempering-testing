rm(list = ls())

library(Rcpp)
library(bench)
library(here)
library(ggplot2)
library(spacefillr)

library(SLHD)                                 ## compare to SLHD library

sourceCpp(here::here("maximin_geom_seq.cpp")) ## load maximin with geometric sequence temperatures
sourceCpp(here::here("maxpro.cpp"))           ## load maxpro
sourceCpp(here::here("randomLHD.cpp"))        ## randomLHD function. Generate a random LHD.

reps <- 50

## Simulation:

## By Jaime Meyer Beilis Michel
## At Bloomington, Indiana, 7.16.2025
## Optimization of parallel tempering for Latin Hypercube Designs

## Uncertainity quantification using PT, EXPERIMENT:
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

  ## Computation of seconds per cycle in a piston
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



## simulation with pre_filled values
## It runs a 50 x 7 design cycle time, and boxplots its its means and variance
sim_piston_cycles_time <- function(reps){
  print("Running piston cycle per second simulation on a 50 x 7 matrix")
  if(file.exists(here::here("simulated_cycle_time.rds"))){
    response <- readline(prompt = "File simulated_cycle_time.rds exists already, do you want to overwrite it? [y/n]: ")
    if(response != "y"){ return(NULL) }
  }

  
  temps = maxpro_temps(50,7,5)


  maximin_sim = lapply(1:reps, function(i) cycle_time(maximinLHD_geom(50,7,9,20,10,1000000,1000,0.095)$design) )
  maxpro_sim  = lapply(1:reps, function(i) cycle_time(pt_maxpro_lhd(50,7,5,100000,10,2500, temps)$design)    )
  random_sim  = lapply(1:reps, function(i) cycle_time(randomLHD(50,7))                                       )
    SLHD_sim  = lapply(1:reps, function(i) cycle_time( maximinSLHD(t = 1, m = 50, k = 7)$Design )            )

  sobol_test =  spacefillr::generate_sobol_set(1e6, 7)
  sobol_sim  =  cycle_time(sobol_test, norm = F)
  mean.response  =  mean(sobol_sim$output)
  var.response   =  var(sobol_sim$output)

  df = data.frame(
    type =    factor(c(rep("maximin",reps), rep("maxpro",reps),rep("random_LHD",reps),rep("SLHD",reps))),
    mean =    c(sapply(1:reps, function(i) maximin_sim[[i]]$mean), sapply(1:reps, function(i) maxpro_sim[[i]]$mean), sapply(1:reps, function(i) random_sim[[i]]$mean), sapply(1:reps, function(i) SLHD_sim[[i]]$mean)),
    var =     c(sapply(1:reps, function(i) maximin_sim[[i]]$var), sapply(1:reps, function(i) maxpro_sim[[i]]$var), sapply(1:reps, function(i) random_sim[[i]]$var), sapply(1:reps, function(i) SLHD_sim[[i]]$var)),
    sd =      c(sapply(1:reps, function(i) maximin_sim[[i]]$sd), sapply(1:reps, function(i) maxpro_sim[[i]]$sd), sapply(1:reps, function(i) random_sim[[i]]$sd), sapply(1:reps, function(i) SLHD_sim[[i]]$sd))
   )

  
 
  ## plot means
  meanplot <- ggplot(data = df, aes(y = type, x = mean, fill = type) ) + 
    geom_boxplot(fill = c("#C7DAE2","#D6DFDF","#E0E3E0", "#EAEBE3")) + 
    theme_bw() + geom_vline(xintercept = mean.response, color = "red", size = 1) + 
    geom_point(data = df,aes(y = type, x = mean),color = "black",alpha = 0.4,size = 1.32) + 
    coord_flip()+
    labs(x = "seconds (mean)", y = "")+
    theme(legend.position = "none", legend.frame = element_blank())

  ggsave(filename = "sim_cycle_time_mean_boxplot.png", plot = meanplot,width = 16,height = 5,units = "in",dpi = 300)

  ## plot variances
  ## TODO: group instead of facet_wrap()
  varplot <- ggplot(data = df, aes(y = type, x = var, fill = type) ) +
    geom_boxplot(fill = c("#C7DAE3","#D6DFE0","#E0E3E1", "#EAEBE4")) + 
    theme_bw() + geom_vline(xintercept = var.response, color = "red", size = 1) + 
    geom_point(data = df,aes(y = type, x = var),color = "black",alpha = 0.4,size = 1.32) + 
    coord_flip() + 
    labs(x = "seconds (variance)", y = "")+
    theme(legend.position = "none", legend.frame = element_blank())

  ggsave(filename = "sim_cycle_time_variance_boxplot.png", plot = varplot,width = 16,height = 5,units = "in",dpi = 300)


  ## SAVE RESULTS
  res =  list(maximin = maximin_sim, maxpro = maxpro_sim, randLHD = random_sim, SLHD = SLHD_sim)

  saveRDS( # save maxpro
    res,
    file = here::here("simulated_cycle_time.rds")
  )

  saveRDS( # save maxpro
    res,
    file = here::here("simulated_cycle_time_data_frame.rds")
  )

  return(res)
  
}

result <- sim_piston_cycles_time(reps)
