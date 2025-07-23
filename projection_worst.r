
library(Rcpp)
library(here)
library(dplyr)
library(ggplot2)
library(SLHD)
library(MaxPro)

## Jaime Meyer Beilis Michel
## At Bloomington Indiana
## 7.17.2025

## Select the worst psi from all possible selections of an increasing n dimensions
## in a (normalized) latin hypercube design. Plots a line graph comparing the results
## between distinct procedures.


## combinatoric function. Out of an X design, given an n amount of dimensions
## extracts the worst measure across all possible combinations of the projection
## makes a line plot.

sourceCpp(here::here("maximin_geom_seq.cpp"))
sourceCpp(here::here("maxpro.cpp"))
sourceCpp(here::here("criteria.cpp"))



## scale the design
scale_design <- function(x) {

  (apply(x, 2, rank) - 0.5) / nrow(x)

}


## generates the temperature schedule for maxpro
maxpro_temps <- function(n,k,M){
  crit1 = 1.0 / (n-1)
  crit2 = (1.0 / ((n - 1)^(k -1) * (n-2)))^(1.0/k)
  delta0 = crit2 - crit1;
  T0 = -delta0 * (1.0/log(0.99))
  L <- sapply(0:M-1, function(i) T0^(M - i))
  return(L)
}

maxpro_temps_geom <- function(n,k,M, a){
  crit1 = 1.0 / (n-1)
  crit2 = (1.0 / ((n - 1)^(k -1) * (n-2)))^(1.0/k)
  delta0 = crit2 - crit1;
  T0 = -delta0 * (1.0/log(0.99))
  L <- sapply(0:M-1, function(i) T0 * a^(M - i))
  return(L)
}



## takes all possible combinations of a matrix for i from 3 to k, returns the worse out of each i projectio
comb_proj_wrst <- function(design){

  n <- nrow(design)
  k <- ncol(design)

  ## output all possible combinations
  powset <- sapply(3:k, function(i) combn(ncol(design), i, simplify = FALSE))
  result <- sapply(1:length(powset), function(i){
    
    max(sapply(powset[[i]], function(comb){ 
      psi(design[, comb], n, length(comb)) 
    }))
  })

  return(result)

}


## Call plot projections on k
## Get rid of pow T0 maxpro
## ensure consistency of alpha (more chains, hotter alpha)


## tests with a default of 18 , 19 and plots result of worse psi projection.
compare_projections <- function(n = 18, k = 19, M = 8, Nmax = 1000000, Nswap = 10){

  maxpro_sched <- maxpro_temps(n=n,k=k,M=M)
  maxpro_geom_sched <- maxpro_temps_geom(n=n,k=k,M=M,a=0.85)


  maxpro_X <- pt_maxpro_lhd(n = n,k = k,M = M,Nmax = Nmax,Nswap = Nswap,tolerance = 5000,temp_sched = maxpro_sched)$design
  maxpro_geom_X <- pt_maxpro_lhd(n = n,k = k,M = M,Nmax = Nmax,Nswap = Nswap,tolerance = 5000,temp_sched = maxpro_geom_sched)$design

  maximin_X <- scale_design(maximinLHD_geom(n=n,k=k,M=M,20,Nswap = Nswap,Nmax = Nmax,tolerance = 5000,alpha =0.1)$design)
  SLHD_X <- scale_design(maximinSLHD(t=1,m=n,k=k)$Design)
  maxpro_lib_X <- MaxProLHD(n=n,p=k)$Design

  maxpro_projs   = data.frame( psi = comb_proj_wrst(maxpro_X) ,    x = 3:k)
  maxpro_geom_proj= data.frame( psi = comb_proj_wrst(maxpro_geom_X),x = 3:k)
  maximin_projs  = data.frame( psi = comb_proj_wrst(maximin_X),    x = 3:k)
  SLHD_projs     = data.frame( psi = comb_proj_wrst(SLHD_X)   ,    x = 3:k)
  mxprolib_projs = data.frame( psi = comb_proj_wrst(maxpro_lib_X), x = 3:k)

  df <- bind_rows(
                  mutate(maxpro_projs,    design = "Maxpro"),
                  mutate(maxpro_geom_proj,design ="Maxpro with geom seq"),
                  mutate(maximin_projs,   design = "Maximin"),
                  mutate(SLHD_projs,      design = "SLHD"),
                  mutate(mxprolib_projs,  design = "MaxPro Library")
  )
    saveRDS( # save maxpro
    df,
    file = here::here("psi_projection_results.rds")
    )


  plot <- ggplot(df, aes(x = x, y = psi, color = design, linetype = design, shape = design))+
    geom_point() +
    geom_line(linewidth = 1) +
    labs(
      x = "k density",
      y = expression(psi ~ " measure"),
      title = ""
    ) +
    theme_bw()
  
  ggsave("psi_worse_proj.png", plot = plot , width = 10,height = 6,units = "in",dpi = 300)


  return(df)
}

