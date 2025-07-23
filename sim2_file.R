
## This has to be ignored by git : not written by mes
## Jose Toledo

library(dplyr)
library(parallel)
library(Rcpp)
library(here)

sourceCpp(here::here("maxpro.cpp"))

reps <- 100

# simulation 1

n_rows <- 27
n_cols <- 2:13
alpha <- c(0.25, 0.5, 0.75, 0.9, 0.95)
n_chains <- c(7, 9, 11, 15)

single_grid <- expand.grid(n_rows = n_rows,n_cols = n_cols,alpha = alpha,n_chains = n_chains)

dat <- do.call(rbind, replicate(reps, single_grid, simplify = FALSE))
# 24000 x 4 # split into three people

start_row = 1
end_row = 24000                          # person 1   took the 24000
# start_row = 8001, end_row = 16000      # person 2
# start_row = 16001, end_row = 24000     # person 3

dat_subset <- dat[start_row:end_row, ]

maxpro_temps <- function(n,k,M, alpha){
  crit1 = 1.0 / (n-1)
  crit2 = (1.0 / ((n - 1)^(k -1) * (n-2)))^(1.0/k)
  delta0 = crit2 - crit1;
  T0 = -delta0 * (1.0/log(0.99))
  return(sapply(0:M-1, function(i) T0*alpha^(M - i)))
}

# Parallelize the loop within the job
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("dat", "dat_subset","single_grid","alpha","end_row","n_chains","n_cols","n_rows","num_cores","reps","start_row", "pt_maxpro_lhd", "maxpro_temps"))
clusterEvalQ(cl, {
  library(here)
  library(Rcpp)
  sourceCpp(here("maxpro.cpp"))
})
# do not change anything here
measures <- parLapply(
  cl,
  1:nrow(dat_subset),
  function(i) {
    return(pt_maxpro_lhd(
      n = dat_subset$n_rows[i],
      k = dat_subset$n_cols[i],
      M = dat_subset$n_chains[i],
      Nswap = 10,
      Nmax = 1000000,
      tolerance = 5000,
      temp_sched = maxpro_temps(
        n = dat_subset$n_rows[i],
        k= dat_subset$n_cols[i], 
        M = dat_subset$n_chains[i], 
        alpha = dat_subset$alpha[i]
      )
    )$measure)
  }
)

measures <- unlist(measures)
final_dat <- cbind(dat_subset, measure = measures)

saveRDS(final_dat, file = here::here("sim_results2.rds")) ## precious data

final_dat <- readRDS(here::here("sim_results2.rds"))
## separate by number of rows / number of columns
## then proceed to analyze means . . . analyze performance , separated in best 10

grouped <- final_dat |> group_by(n_cols, n_rows) 
grouped <- grouped |> summarise(mean_measure = mean(measure), sd_measure = sd(measure))
grouped$lower_confidence_bound <- grouped$mean_measure - grouped$sd_measure


  final_dat |>
    group_by(n_cols, n_rows, alpha, n_chains) |> 
    summarise(mean_measure = mean(measure), sd_measure = sd(measure),
  lower_bound = mean_measure - sd_measure) |> 
    ungroup() |> 
    group_by(n_cols, n_rows) |> 
    arrange( lower_bound , .by_group = TRUE) |>
    slice_head(n = 1) |> 
    tinytable::tt(digits = 3) |>
    tinytable::theme_tt("tabular") |>
    print("latex")



## show the minimum of the lwoer confidence bound in the data.
print(grouped[which.min( grouped$lower_confidence_bound),])
grouped <- grouped[order(grouped$lower_confidence_bound),]

## note that all these correspond to a similar number of columns [12,13] and perform on lower alphas
head(grouped, 10)
