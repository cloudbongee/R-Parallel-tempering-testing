
## graph sim  dot R

library(dplyr)
library(parallel)
library(Rcpp)
library(ggplot2)
library(here)
library(patchwork)


## groups by n_chains
# draws as many graphs as df$alpha 
graph_1 <- function(filename, measure_tag = "Measure"){

  df <- readRDS(here::here(filename))

  df$n_chains <- factor(df$n_chains)

  alphas <- sort(unique(df$alpha))

  ## filter to only have the number of chains 7,9,11
  df <- df |> 
    filter(n_chains %in% c(7,9,11))


  for(a in alphas){
    df_distinct <- df |> 
      distinct(n_rows, n_cols)

    plot_list <- lapply(seq_len(nrow(df_distinct)), function(i){
      dat <- df |>
        filter(n_rows == df_distinct$n_rows[i], n_cols == df_distinct$n_cols[i], alpha == a)

      p <- ggplot(data = dat, aes( x = measure, fill = n_chains, group = n_chains))  + ## faceted by number of chains 
        geom_density(alpha = 0.2) + 
        labs(
          title = paste0( "rows = ", df_distinct$n_rows[i], ", cols = ", df_distinct$n_cols[i]),
          x = measure_tag,
          y = "Density"
        ) + 
        theme_minimal()

        
        return(p)
    })


    plot <- wrap_plots(plot_list) + plot_layout(guides = 'collect')
    ggsave(filename = here(paste0(filename,"graph_by_alpha",a,".png")), plot, width = 10, height = 8, units = "in")
  }
}

graph_1("sim_results1.rds", expression(phi ~ "measure"))

graph_1("sim_results2.rds", expression(psi ~ "measure"))
