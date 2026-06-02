source('helpers/results_gathering_tool.R')

alpha = 0.05
true_beta_v = c(0.1, 0.3, -0.2, 0.1)
true_variance_com_v = c(0.8, 0.4, 0.2)
n_beta = length(true_beta_v)

result_file_v = list.files('results', pattern = "^n\\d+_M\\d+_[A-Za-z]+_rdi?s?\\.RData$")
n_result_file = length(result_file_v)

result_df_ls = vector('list', length(n_result_file))
for (i in 1:n_result_file){
  result_file = result_file_v[i]
  n = as.integer(sub('^n(\\d+)_.*', '\\1', result_file))
  M = as.integer(sub('.*_M(\\d+)_.*', '\\1', result_file))
  re_dist = sub('.*_M\\d+_([^_]+)_.*', '\\1', result_file)
  with_rd_slope = ifelse(sub('.*_([^_]+)\\.[^.]*$', '\\1', result_file) == 'rdis', TRUE, FALSE)
  freq_pack_v = c('lme4_LA', if (! with_rd_slope) 'lme4_AGQ', 'GLMMadaptive',
                  'glmmTMB', 'MASS', 'hglm')
  bayes_pack_v = c(if (! with_rd_slope) 'brms', 'rstanarm')
  pack_v = c(freq_pack_v, bayes_pack_v)
  n_pack = length(pack_v)
  load_obj_v = load(sprintf('results/%s', result_file))
  fail_pos_df = conv_status_df != 'good'
  for (load_obj in load_obj_v) {
    obj = get(load_obj)
    if (! load_obj %in% c('conv_status_df', 'compute_time_df')){
      obj[fail_pos_df] = NA
    }
    assign(load_obj, obj[, pack_v], envir = .GlobalEnv)
}

  result_df = data.frame(
    'N.subjects' = rep(n, n_pack),
    'Max.m' = rep(M, n_pack),
    'Response.distribution' = rep(re_dist, n_pack),
    'Random.effect' = rep(ifelse(with_rd_slope, 'Intercepts+Slopes', 'Intercepts'), n_pack),
    'Packages' = pack_v)

  result_df[,'Convergency.ratio.official'] = apply(conv_status_df, 2, official_conv_ratio)
  result_df[,'Convergency.ratio.real'] = apply(conv_status_df, 2, real_conv_ratio)
  result_df[,'Computational.time'] = apply(compute_time_df, 2, mean)
  result_df[,'Mean.beta0_hat'] = apply(beta0_hat_df, 2, mean, na.rm = TRUE)
  result_df[,'RMSE.beta0_hat'] = apply(beta0_hat_df, 2, rmse, true_theta = true_beta_v[1])
  result_df[,'Mean.beta1_hat'] = apply(beta1_hat_df, 2, mean, na.rm = TRUE)
  result_df[,'RMSE.beta1_hat'] = apply(beta1_hat_df, 2, rmse, true_theta = true_beta_v[2])
  result_df[,'Mean.beta2_hat'] = apply(beta2_hat_df, 2, mean, na.rm = TRUE)
  result_df[,'RMSE.beta2_hat'] = apply(beta2_hat_df, 2, rmse, true_theta = true_beta_v[3])
  result_df[,'Mean.beta3_hat'] = apply(beta3_hat_df, 2, mean, na.rm = TRUE)
  result_df[,'RMSE.beta3_hat'] = apply(beta3_hat_df, 2, rmse, true_theta = true_beta_v[4])

  # Testing
  alpha_hat_uni_df[, freq_pack_v] = alpha_hat_uni_df[, freq_pack_v] < alpha
  result_df[, 'Alpha_hat.univariate'] = apply(alpha_hat_uni_df, 2, mean, na.rm = TRUE)
  power_hat_uni_df[, freq_pack_v] = power_hat_uni_df[, freq_pack_v] < alpha
  result_df[, 'Power_hat.univariate'] = apply(power_hat_uni_df, 2, mean, na.rm = TRUE)
  alpha_hat_mult_df[, freq_pack_v] = alpha_hat_mult_df[, freq_pack_v] < alpha
  result_df[, 'Alpha_hat.multivariate'] = apply(alpha_hat_mult_df, 2, mean, na.rm = TRUE)
  power_hat_mult_df[, freq_pack_v] = power_hat_mult_df[, freq_pack_v] < alpha
  result_df[, 'Power_hat.multivariate'] = apply(power_hat_mult_df, 2, mean, na.rm = TRUE)

  # Random effects
  result_df[, 'Mean.tau0_hat'] = apply(tau0_hat_df, 2, mean, na.rm = TRUE)
  result_df[, 'RMSE.tau0_hat'] = apply(tau0_hat_df, 2, rmse, true_theta = true_variance_com_v[1])
  if (with_rd_slope){
    result_df[, 'Mean.tau1_hat'] = apply(tau1_hat_df, 2, mean, na.rm = TRUE)
    result_df[, 'RMSE.tau1_hat'] = apply(tau1_hat_df, 2, rmse, true_theta = true_variance_com_v[2])
    result_df[, 'Mean.rho01_hat'] = apply(rho01_hat_df, 2, mean, na.rm = TRUE)
    result_df[, 'RMSE.rho01_hat'] = apply(rho01_hat_df, 2, rmse, true_theta = true_variance_com_v[3])
  } else{
    result_df[, 'Mean.tau1_hat'] = rep(NA, n_pack)
    result_df[, 'RMSE.tau1_hat'] = rep(NA, n_pack)
    result_df[, 'Mean.rho01_hat'] = rep(NA, n_pack)
    result_df[, 'RMSE.rho01_hat'] = rep(NA, n_pack)
  }
  result_df_ls[[i]] = result_df
}
final_result_df = do.call('rbind', result_df_ls)

# Mean absolute bias
mean_df = final_result_df[, c('Packages', 'Mean.beta0_hat', 'Mean.beta1_hat',
                              'Mean.beta2_hat', 'Mean.beta3_hat', 'Mean.tau0_hat',
                              'Mean.tau1_hat', 'Mean.rho01_hat')]
mean_df[, 2:8] = abs(sweep(mean_df[, 2:8], 2, c(true_beta_v, true_variance_com_v), '-'))
mean_abs_bias_df = aggregate(mean_df[-1], by = list(Group = mean_df[[1]]), FUN = mean)

# Mean RMSE
rmse_df = final_result_df[, c('Packages', 'RMSE.beta0_hat', 'RMSE.beta1_hat',
                              'RMSE.beta2_hat', 'RMSE.beta3_hat', 'RMSE.tau0_hat',
                              'RMSE.tau1_hat', 'RMSE.rho01_hat')]
mean_rmse_df = aggregate(rmse_df[-1], by = list(Group = rmse_df[[1]]), FUN = mean)

save(list = c('final_result_df', 'mean_abs_bias_df', 'mean_rmse_df'), file = 'results/final.RData')
