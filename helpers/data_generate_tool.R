gen_random_effects = function(with_rd_slope, n, subject_id){
  if (!with_rd_slope) {
    # Random intercepts for each subject
    tau = 0.8
    random_effects = rnorm(n, mean = 0, sd = tau)
    random_intercept = random_effects[subject_id]
    random_slope = 0
  } else {
    # Random intercepts and slopes for each subject
    tau0 = 0.8
    tau1 = 0.4
    rho01 = 0.2
    sigma_mat = matrix(c(tau0^2, rho01*tau0*tau1, rho01*tau0*tau1, tau1^2),
                       nrow = 2, ncol = 2, byrow = T)
    random_effects = mvrnorm(n, mu = c(0, 0), Sigma = sigma_mat)
    random_intercept = random_effects[subject_id, 1]
    random_slope = random_effects[subject_id, 2]
  }
  return(list('random_intercept' = random_intercept,
              'random_slope' = random_slope))
}

gen_response = function(z, num_obs, resp_dist){
  if (resp_dist == 'Bernoulli') {
    pi_full = 1/(1+exp(-z[['z_full']]))
    pi_0 = 1/(1+exp(-z[['z_0']]))
    pi_00 = 1/(1+exp(-z[['z_00']]))

    y_full = rbinom(num_obs, size = 1, prob = pi_full)
    y_0 = rbinom(num_obs, size = 1, prob = pi_0)
    y_00 = rbinom(num_obs, size = 1, prob = pi_00)
  } else if (resp_dist == 'Poisson'){
    mu_full = exp(z[['z_full']])
    mu_0 = exp(z[['z_0']])
    mu_00 = exp(z[['z_00']])

    y_full = rpois(num_obs, lambda = mu_full)
    y_0 = rpois(num_obs, lambda = mu_0)
    y_00 = rpois(num_obs, lambda = mu_00)
} else {
    warning('The given response distribution is not supported.')
  }

  return(list('y_full' = y_full, 'y_0' = y_0, 'y_00' = y_00))
}

gen_unbalanced_data = function(n, M, resp_dist, with_rd_slope){

  num_sim = 1000

  # Fixed effects
  beta_v = c(0.1, 0.3, -0.2, 0.1)
  beta_0_v = c(0.1, 0.3, -0.2, 0)
  beta_00_v = c(0.1, 0.3, 0, 0)

  # Possible set of time points
  possi_time_set = (0:(M-1))/(M-1)

  simulate_data_ls = list()
  for (i in 1:num_sim){
    # Time variable
    num_time_v = sample(1:M, n, replace = TRUE)
    time_variable = unlist(lapply(
      num_time_v,
      FUN = function(x) sort(sample(possi_time_set, size = x, replace = FALSE))))

    # Sample size
    num_obs = sum(num_time_v)

    # Subject ID
    subject_id = factor(rep(1:n, times = num_time_v))

    # Group variable
    group = rbinom(n, size = 1, prob = 0.4)
    group_variable = group[subject_id]

    # Design matrix of fixed effects
    X_fixed = cbind(rep(1, num_obs), time_variable, group_variable,
                    time_variable * group_variable)

    # Random effects
    temp = gen_random_effects(with_rd_slope = with_rd_slope, n = n,
                              subject_id = subject_id)
    random_intercept = temp[['random_intercept']]
    random_slope = temp[['random_slope']]

    # Linear predictors
    z_full = X_fixed %*% beta_v + random_intercept + random_slope * time_variable
    z_0 = X_fixed %*% beta_0_v + random_intercept + random_slope * time_variable
    z_00 = X_fixed %*% beta_00_v + random_intercept + random_slope * time_variable

    # Response
    temp = gen_response(z = list('z_full' = z_full, 'z_0' = z_0, 'z_00' = z_00),
                        num_obs = num_obs, resp_dist = resp_dist)

    # Dataset
    simulate_data = data.frame(
      y_full = temp[['y_full']],
      y_0 = temp[['y_0']],
      y_00 = temp[['y_00']],
      id = subject_id,
      time = time_variable,
      group = group_variable,
      random_intercept = random_intercept,
      random_slope = random_slope
    )
    simulate_data_ls[[i]] = simulate_data
  }
  comb = paste0('n', n, '_M', M, '_', resp_dist,
                ifelse(with_rd_slope, '_rdis', '_rdi'))
  save(simulate_data_ls, file = sprintf('data/%s.RData', comb))
  cat(sprintf('%s.RData done.\n', comb))
}