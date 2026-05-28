library(jsonlite)

n_v = c(50, 100, 200)
M_v = c(3, 9)
re_dist_v = c('Bernoulli', 'Poisson')
with_rd_slope_v = c(FALSE, TRUE)

i = 1
for (with_rd_slope in with_rd_slope_v){
  for (n in n_v){
    for (M in M_v){
      for (re_dist in re_dist_v){
        config = list(n = n, M = M, re_dist = re_dist, with_rd_slope = with_rd_slope)
        write_json(config, paste0('config/config_', i, '.json'), pretty = TRUE)
        i = i + 1
      }
    }
  }
}
