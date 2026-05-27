library(lme4)
library(GLMMadaptive)
library(glmmTMB)
library(MASS)
library(hglm)
library(jsonlite)

source(file.path('helpers/estimation_tool.R'))

# Read config file
c = 10
config = fromJSON(sprintf('config/config_%s.json', c))
comb = sprintf('n%s_M%s_%s_%s', config$n, config$M, config$re_dist,
               ifelse(config$with_rd_slope, 'rdis', 'rdi'))

# Load datasets
load(sprintf('data/%s.RData', comb))
n_sim = length(simulate_data_ls)

# Create/Load results
res_path = sprintf('results/%s.RData', comb)
pack_v = c('lme4_LA', 'lme4_AGQ', 'GLMMadaptive', 'glmmTMB', 'MASS', 'hglm')
result_v = c('conv_status', 'compute_time', 'beta0_hat', 'beta1_hat',
             'beta2_hat', 'beta3_hat', 'tau0_hat',
             if (config$with_rd_slope) c('tau1_hat', 'rho01_hat'),
             'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni', 'power_hat_mult')
if (!file.exists(res_path)) {
  for (result in result_v){
    NA_df = data.frame(matrix(NA, nrow = n_sim, ncol = length(pack_v)))
    colnames(NA_df) = pack_v
    assign(sprintf('%s_df', result), NA_df)
  }
  save(list = paste0(result_v, '_df'), file = res_path)
}

load(res_path)

# Estimate

for (pack in pack_v){
  for (i in 1:n_sim){
    cat(sprintf('Fitting GLMM on dataset %s by %s...\n', i, pack))
    pack_result = match.fun(sprintf("est_%s", pack))(
      data = simulate_data_ls[[i]], resp_dist = config$re_dist, with_rd_slope = config$with_rd_slope)
    for (result in result_v){
      result_df = get(sprintf('%s_df', result))
      result_df[i, pack] = pack_result[result]
      assign(sprintf('%s_df', result), result_df)
    }
    if (i %% 10 == 0){
      save(list = paste0(result_v, '_df'), file = res_path)
    }
  }
}

