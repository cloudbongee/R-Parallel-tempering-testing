

maxpro_temps <- function(n,k,M, alpha){
  crit1 = 1.0 / (n-1)
  crit2 = (1.0 / ((n - 1)^(k -1) * (n-2)))^(1.0/k)
  delta0 = crit2 - crit1;
  T0 = -delta0 * (1.0/log(0.99))
  return(sapply(0:M-1, function(i) T0*alpha^(M - i)))
}


## Just a simple function to get a latex table of space-time performance of the functions.
## currently using the proposed maximin measures as based on hyperparameter tuning.
bench_out <- function(type = "maximin", alpha = 0.25, M = 9, n = 27, k = 13, p = 15){

  library(dplyr)
  library(Rcpp)
  library(here)
  library(tinytable)
  library(bench)

  if(type == "maximin"){
    sourceCpp(here("maximin_geom_seq.cpp"))
    marked <- mark(maximinLHD_geom(n = n, k = k, M = M, p = p, Nswap = 10, Nmax = 1000000, tolerance = 5000, alpha = alpha))

  }else if(type == "maxpro"){
    sourceCpp(here("maxpro.cpp"))
    temp <- maxpro_temps(n=n, k = k, M = M, alpha = alpha)
    marked <- mark(pt_maxpro_lhd(n = n, k = k, M = M, Nmax = 1000000, Nswap = 10, tolerance = 5000, temp_sched = temp))
    

  }else{
    print("No type specified , no results")
    return(NULL)
  }

  print( marked ) 

  print("WITH LATEX OUTPUT ~~~~~~~~~~~~~~~")
  print("")


  ## after the benchmark has been specified:
  ## return the latex list. Only the subset containing

  marked <- marked[c("min", "median", "itr/sec", "mem_alloc", "total_time")]
  benchmark <- as.data.frame(marked) |> 
    mutate(
      min = as.numeric(min),
      median = as.numeric(median),
      `itr/sec` = as.numeric(`itr/sec`),
      mem_alloc = format(mem_alloc) ,
      total_time = as.numeric(total_time)
    )

  return(benchmark |> 
    tinytable::tt(digits = 3) |> 
    tinytable::theme_tt("tabular") |>
    print("latex") )
  

}


## Ignore these, I am just testing the output and runtime
bench_out(alpha = 0.09, M = 9, n = 27, k = 3, p = 30)


bench_out("maxpro", alpha = 0.9, M = 11, n = 27, k = 13)
