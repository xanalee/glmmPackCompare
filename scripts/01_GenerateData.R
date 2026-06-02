library(MASS)
library(jsonlite)

source('helpers/data_generate_tool.R')

set.seed(126) # for reproducibility

for (i in 1:24){
  config_i = fromJSON(sprintf('config/config_%s.json', i))
  gen_unbalanced_data(n = config_i$n, M = config_i$M, resp_dist = config_i$re_dist,
                      with_rd_slope = config_i$with_rd_slope)
}
