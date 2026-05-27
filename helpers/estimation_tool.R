build_formula = function(response, fixed_terms, random_terms) {
  form = paste(response, '~', fixed_terms, '+', random_terms)
  return(as.formula(form))
}

get_p_value = function(model, term) {
  # lme4
  if (inherits(model, 'glmerMod')) {
    return(summary(model)$coefficients[term, 'Pr(>|z|)'])
  }
  # GLMMadaptive
  if (inherits(model, 'MixMod')) {
    return(summary(model)$coef_table[term, 'p-value'])
  }
  # glmmTMB
  if (inherits(model, 'glmmTMB')) {
    return(summary(model)$coefficients$cond[term, 'Pr(>|z|)'])
  }
  # MASS
  if (inherits(model, 'glmmPQL')) {
    return(summary(model)$tTable[term, 'p-value'])
  }
  # hglm
  if (inherits(model, 'hglm')) {
    return(summary(model)$FixCoefMat[term, 'Pr(>|t|)'])
  }
  stop('Unknown model class – cannot extract p-value.')
}

check_CI = function(model, term) {
    CI = posterior_interval(model, prob = 0.95)
    !((CI[term, '2.5%'] < 0) & (CI[term, '97.5%'] > 0))
  }

check_convergence = function(model, ...) {
  if (inherits(model, 'try-error')) return('error')

  # lme4
  if (inherits(model, 'glmerMod')) {
    k = list(...)$k
    if (is.null(k)) stop('For lme4 models, provide k = number of fixed effects.')
    if (model@optinfo$conv$opt != 0) return('official_fail')
    if (length(model@optinfo$conv$lme4) > 0 || length(model@beta) != k) return('extra_fail')
    return('good')
  }

  # GLMMadaptive
  if (inherits(model, 'MixMod')) {
    if (!model$converged) return('official_fail')
    return('good')
  }

  # glmmTMB
  if (inherits(model, 'glmmTMB')) {
    if (model$fit$convergence != 0) return('official_fail')
    return('good')
  }

  # MASS
  if (inherits(model, 'glmmPQL')) {
    if (!is.matrix(model$apVar)) return('extra_fail')
    return('good')
  }

  # hglm
  if (inherits(model, 'hglm')) {
    if (model$Converge != 'converged') return('official_fail')
    return('good')
  }

  # brms
  if (inherits(model, 'brmsfit')){
    n_term = list(...)$n_term
    if (is.null(n_term)) stop('For brms models, provide n_term = number of parameters need to be estimated.')
    Rhat = rhat(model)
    if (any(abs(Rhat[1:n_term] - 1) >= 0.1)) return('official_fail')
    return('good')
  }

  # rstanarm
  if (inherits(model, 'stanreg')){
    fixed_names = list(...)$fixed_names
    if (is.null(fixed_names)) stop("For rstanarm models, provide fixed_names = c('time', 'group', ...).")
    summ = summary(model)
    if (any(! fixed_names %in% row.names(summ))) return('error')
    with_rd_slope = list(...)$with_rd_slope
    if (is.null(with_rd_slope)) stop("For rstanarm models, provide with_rd_slope = if the model includes random slope.")
    random_names = c('Sigma[id:(Intercept),(Intercept)]',
                     if (config$with_rd_slope) c('Sigma[id:time,(Intercept)]', 'Sigma[id:time,time]'))
    if (any(abs(summ[c(fixed_names, random_names), 'Rhat'] - 1) >= 0.1)) return('official_fail')
    return('good')
  }

  stop('Unknown model class – cannot check convergence.')
}

est_lme4 = function(data, resp_dist, with_rd_slope, nAGQ){

  # Specify random effects
  rd_term = if (with_rd_slope) '(1 + time|id)' else '(1|id)'

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial
  }
  if (resp_dist == 'Poisson'){
    fami = poisson
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(
      glmer(build_formula('y_full', 'time + group + time:group', rd_term),
            family = fami, data = data, nAGQ = nAGQ))
  })[['elapsed']]
  conv_status = check_convergence(model_full, k = 4)

  if (conv_status %in% c('good', 'extra_fail')){

    # Record beta and theta after fitting full model on full dataset
    beta0_hat = model_full@beta[1]
    beta1_hat = model_full@beta[2]
    beta2_hat = model_full@beta[3]
    beta3_hat = model_full@beta[4]

    vc = VarCorr(model_full)$id
    if (with_rd_slope){
      tau0_hat = attr(vc,'stddev')[['(Intercept)']]
      tau1_hat = attr(vc,'stddev')[['time']]
      rho01_hat = attr(vc,'correlation')['(Intercept)', 'time']
    } else {
      tau0_hat = attr(vc,'stddev')[['(Intercept)']]
    }
    if (conv_status == 'good'){
      # Calculate empirical power of univariate test (Wald)
      power_hat_uni = get_p_value(model_full, 'time:group')
      # Fit null model on full dataset
      null_model_full = try(
        glmer(build_formula('y_full', 'time', rd_term),
              family = fami, data = data, nAGQ = nAGQ))

      if (check_convergence(null_model_full, k = 2) == 'good'){
        # Calculate empirical power of multivariate test (LRT)
        ano = anova(null_model_full, model_full)
        power_hat_mult = ano['model_full', 'Pr(>Chisq)']
      }
    }

    # Fit full model on dataset with beta3 = 0
    model_uni = try(
      glmer(build_formula('y_0', 'time + group + time:group', rd_term),
            family = fami, data = data, nAGQ = nAGQ))
    if (check_convergence(model_uni, k = 4) == 'good'){
      # Calculate empirical alpha of univariate test (Wald)
      alpha_hat_uni = get_p_value(model_uni, 'time:group')
    }

    # Fit full model on dataset with beta2 = beta3 = 0
    model_mult = try(
      glmer(build_formula('y_00', 'time + group + time:group', rd_term),
            family = fami, data = data, nAGQ = nAGQ))
    if (check_convergence(model_mult, k = 4) == 'good'){
      # Fit null model on dataset with beta2 = beta3 = 0
      null_model_mult = try(
        glmer(build_formula('y_00', 'time', rd_term),
              family = fami, data = data, nAGQ = nAGQ))
      if (check_convergence(null_model_mult, k = 2) == 'good'){
        # Calculate empirical alpha of multivariate test (LRT)
        ano = anova(null_model_mult, model_mult)
        alpha_hat_mult = ano['model_mult', 'Pr(>Chisq)']
      }
    }
  }

  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_lme4_LA = pryr::partial(est_lme4, nAGQ = 1)
est_lme4_AGQ = pryr::partial(est_lme4, nAGQ = 11)

est_GLMMadaptive = function(data, resp_dist, with_rd_slope){

  # Specify random effects
  rd_term = as.formula(ifelse(with_rd_slope, '~ 1 + time | id', '~ 1 | id'))

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial()
  }
  if (resp_dist == 'Poisson'){
    fami = poisson()
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(mixed_model(fixed = y_full ~ time + group + time:group,
                                 random = rd_term,
                                 family = fami,
                                 data = data
    ))})[['elapsed']]
  conv_status = check_convergence(model_full)

  if (conv_status == 'good'){

    # Record beta and theta after fitting full model on full dataset
    beta0_hat = model_full$coefficients[['(Intercept)']]
    beta1_hat = model_full$coefficients[['time']]
    beta2_hat = model_full$coefficients[['group']]
    beta3_hat = model_full$coefficients[['time:group']]
    if (with_rd_slope){
      tau0_hat = model_full$D['(Intercept)', '(Intercept)']^0.5
      tau0_hat = tau0_hat
      tau1_hat = model_full$D['time', 'time']^0.5
      rho01_hat = model_full$D['(Intercept)', 'time']/(tau0_hat*tau1_hat)
    } else {
      tau_hat = model_full$D['(Intercept)', '(Intercept)']^0.5
    }

    # Calculate empirical power of univariate test (Wald)
    power_hat_uni = get_p_value(model_full, 'time:group')

    # Fit null model on full dataset
    null_model_full = try(mixed_model(fixed = y_full ~ time,
                                      random = rd_term,
                                      family = fami,
                                      data = data
    ))
    if (check_convergence(null_model_full) == 'good') {
      # Calculate empirical power of multivariate test (LRT)
      ano = anova(null_model_full, model_full)
      power_hat_mult = ano$p.value
    }

    # Fit full model on dataset with beta3 = 0
    model_uni = try(mixed_model(fixed = y_0 ~ time + group + time:group,
                                random = rd_term,
                                family = fami,
                                data = data
    ))
    if (check_convergence(model_uni) == 'good'){
      # Calculate empirical alpha of univariate test (Wald)
      alpha_hat_uni = get_p_value(model_uni, 'time:group')
    }

    # Fit full model on dataset with beta2 = beta3 = 0
    model_mult = try(mixed_model(fixed = y_00 ~ time + group + time:group,
                                 random = rd_term,
                                 family = fami,
                                 data = data
    ))
    if (check_convergence(model_mult) == 'good') {
      # Fit null model on dataset with beta2 = beta3 = 0
      null_model_mult = try(mixed_model(fixed = y_00 ~ time,
                                        random = rd_term,
                                        family = fami,
                                        data = data
      ))
      if (check_convergence(null_model_mult) == 'good'){
        # Calculate empirical alpha of multivariate test (LRT)
        ano = anova(null_model_mult, model_mult)
        alpha_hat_mult = ano$p.value
      }
    }
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_glmmTMB = function(data, resp_dist, with_rd_slope){

  # Specify random effects
  rd_term = if (with_rd_slope) '(1 + time|id)' else '(1|id)'

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial
  }
  if (resp_dist == 'Poisson'){
    fami = poisson
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(glmmTMB(build_formula('y_full', 'time + group + time:group', rd_term),
                             family = fami,
                             data = data
    ))})[['elapsed']]
  conv_status = check_convergence(model_full)

  if (conv_status == 'good'){

    # Record beta and theta after fitting full model on full dataset
    beta0_hat = model_full$fit$par[[1]]
    beta1_hat = model_full$fit$par[[2]]
    beta2_hat = model_full$fit$par[[3]]
    beta3_hat = model_full$fit$par[[4]]
    summ = summary(model_full)
    if (with_rd_slope){
      tau0_hat = as.numeric(attr(summ$varcor$cond$id, 'stddev')['(Intercept)'])
      tau1_hat = as.numeric(attr(summ$varcor$cond$id, 'stddev')['time'])
      rho01_hat = attr(summ$varcor$cond$id, 'correlation')['(Intercept)', 'time']
    } else {
      tau_hat = as.numeric(attr(summ$varcor$cond$id, 'stddev'))
    }

    # Calculate empirical power of univariate test (Wald)
    power_hat_uni = get_p_value(model_full, 'time:group')

    # Fit null model on full dataset
    null_model_full = try(glmmTMB(build_formula('y_full', 'time', rd_term),
                                  family = fami,
                                  data = data
    ))
    if (check_convergence(null_model_full) == 'good'){
      # Calculate empirical power of multivariate test (LRT)
      ano = anova(null_model_full, model_full)
      power_hat_mult = ano['model_full', 'Pr(>Chisq)']
    }

    # Fit full model on dataset with beta3 = 0
    model_uni = try(glmmTMB(build_formula('y_0', 'time + group + time:group', rd_term),,
                            family = fami,
                            data = data
    ))
    if (check_convergence(model_uni) == 'good'){
      # Calculate empirical alpha of univariate test (Wald)
      alpha_hat_uni = get_p_value(model_uni, 'time:group')
    }

    # Fit full model on dataset with beta2 = beta3 = 0
    model_mult = try(glmmTMB(build_formula('y_00', 'time + group + time:group', rd_term),
                             family = fami,
                             data = data
    ))
    if (check_convergence(model_mult) == 'good'){
      # Fit null model on dataset with beta2 = beta3 = 0
      null_model_mult = try(glmmTMB(build_formula('y_00', 'time', rd_term),
                                    family = fami,
                                    data = data
      ))
      if (check_convergence(null_model_mult) == 'good'){
        # Calculate empirical alpha of multivariate test (LRT)
        ano = anova(null_model_mult, model_mult)
        alpha_hat_mult = ano['model_mult', 'Pr(>Chisq)']
      }
    }
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_MASS = function(data, resp_dist, with_rd_slope){
  # Specify random effects
  rd_term = as.formula(ifelse(with_rd_slope, '~ 1 + time | id', '~ 1 | id'))

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial
  }
  if (resp_dist == 'Poisson'){
    fami = poisson
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(glmmPQL(y_full ~ time + group + time:group,
                             random = rd_term,
                             family = fami,
                             data = data,
                             control = list(returnObject = T)
    ))})[['elapsed']]
  conv_status = check_convergence(model_full)

  if (conv_status %in% c('good', 'extra_fail')) {

    # Record beta and theta after fitting full model on full dataset
    beta0_hat = model_full$coefficients$fixed[['(Intercept)']]
    beta1_hat = model_full$coefficients$fixed[['time']]
    beta2_hat = model_full$coefficients$fixed[['group']]
    beta3_hat = model_full$coefficients$fixed[['time:group']]
    temp = VarCorr(model_full)
    if (with_rd_slope){
      tau0_hat = as.numeric(temp['(Intercept)', 'StdDev'])
      tau1_hat = as.numeric(temp['time', 'StdDev'])
      rho01_hat = as.numeric(temp['time', 'Corr'])
    } else {
      tau_hat = as.numeric(temp['(Intercept)', 'StdDev'])
    }

    if (conv_status == 'good'){
      # Calculate empirical power of univariate test (Wald)
      power_hat_uni = get_p_value(model_full, 'time:group')
    }

    # Fit full model on dataset with beta3 = 0
      model_uni = try(glmmPQL(y_0 ~ time + group + time:group,
                              random = rd_term,
                              family = fami,
                              data = data,
                              control = list(returnObject = T)
      ))
      if (check_convergence(model_uni) == 'good'){
        # Calculate empirical alpha of univariate test (Wald)
        alpha_hat_uni = get_p_value(model_uni, 'time:group')
      }
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_hglm = function(data, resp_dist, with_rd_slope){

  # Specify random effects
  rd_term = if (with_rd_slope) '(1 + time | id)' else '(1 | id)'

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial(link = 'logit')
  }
  if (resp_dist == 'Poisson'){
    fami = poisson(link = 'log')
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(
      hglm2(meanmodel = build_formula('y_full', 'time + group + time:group', rd_term),
            family = fami,
            rand.family = gaussian(link = identity),
            data = data)
    )})[['elapsed']]
  conv_status = check_convergence(model_full)

  if (conv_status == 'good'){

    # Record beta and theta after fitting full model on full dataset
    beta0_hat = model_full$fixef[['(Intercept)']]
    beta1_hat = model_full$fixef[['time']]
    beta2_hat = model_full$fixef[['group']]
    beta3_hat = model_full$fixef[['time:group']]
    if (with_rd_slope){
      tau0_hat = model_full$varRanef[1]^0.5
      tau1_hat = model_full$varRanef[2]^0.5
      rho01_hat = 0
    } else {
      tau_hat = model_full$varRanef^0.5
    }

    # Calculate empirical power of univariate test (Wald)
    power_hat_uni = get_p_value(model_full, 'time:group')

    # Fit full model on dataset with beta3 = 0
    model_uni = try(
      hglm2(meanmodel = build_formula('y_0', 'time + group + time:group', rd_term),
            family = fami,
            rand.family = gaussian(link = identity),
            data = data
        )
    )
    if (check_convergence(model_uni, warning_list) == 'good'){
      # Calculate empirical alpha of univariate test (Wald)
      alpha_hat_uni = get_p_value(model_uni, 'time:group')
    }
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_brms = function(data, resp_dist, with_rd_slope, n_cores = 4, n_threads = 6){

  # Specify random effects
  rd_term = if (with_rd_slope) '(1 + time | id)' else '(1|id)'

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = bernoulli(link = 'logit')
  }
  if (resp_dist == 'Poisson'){
    fami = poisson(link = 'log')
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(brm(build_formula('y_full', 'time + group + time:group', rd_term),
                         family = fami,
                         data = data,
                         threads = threading(n_threads),
                         cores = n_cores,
                         seed = 123
    ))})[['elapsed']]

  conv_status = check_convergence(model_full, n_term = 5)

  if (conv_status == 'good'){

    # Record betas and tau after fitting full model on full dataset
    summ = summary(model_full)
    beta0_hat = summ$fixed['Intercept', 'Estimate']
    beta1_hat = summ$fixed['time', 'Estimate']
    beta2_hat = summ$fixed['group', 'Estimate']
    beta3_hat = summ$fixed['time:group', 'Estimate']
    if (with_rd_slope){
      tau0_hat = summ$random$id['sd(Intercept)', 'Estimate']
      tau1_hat = summ$random$id['sd(time)', 'Estimate']
      rho01_hat = summ$random$id['cor(Intercept,time)', 'Estimate']
    } else {
      tau_hat = summ$random$id['sd(Intercept)', 'Estimate']
    }
    # Calculate empirical power of univariate test (CI)
    power_hat_uni = check_CI(model_full, 'b_time:group')

    # Fit null model on full dataset
    null_model_full = try(brm(build_formula('y_full', 'time', rd_term),
                              family = fami,
                              data = data,
                              threads = threading(n_threads),
                              cores = n_cores,
                              seed = 123
    ))
    if (check_convergence(null_model_full, n_term = 3) == 'good'){

      # Calculate empirical power of multivariate test (LOOIC)
      null_loo = loo(null_model_full)
      alt_loo = loo(model_full)
      power_hat_mult = null_loo$estimates['looic', 'Estimate'] > alt_loo$estimates['looic', 'Estimate']
    }
    rm(model_full, null_model_full)

    # Fit full model on dataset with beta3 = 0
    model_uni = try(brm(build_formula('y_0', 'time + group + time:group', rd_term),
                        family = fami,
                        data = data,
                        threads = threading(n_threads),
                        cores = n_cores,
                        seed = 123
    ))
    if (check_convergence(model_uni, n_term = 5) == 'good'){

      # Calculate empirical alpha of univariate test (CI)
      alpha_hat_uni = check_CI(model_uni, 'b_time:group')
    }
    rm(model_uni)

    # Fit full model on dataset with beta2 = beta3 = 0
    model_mult = try(brm(build_formula('y_00', 'time + group + time:group', rd_term),
                         family = fami,
                         data = data,
                         threads = threading(n_threads),
                         cores = n_cores,
                         seed = 123
    ))
    if (check_convergence(model_mult, n_term = 5) == 'good'){

      # Fit null model on dataset with beta2 = beta3 = 0
      null_model_mult = try(brm(build_formula('y_00', 'time', rd_term),
                                family = fami,
                                data = data,
                                threads = threading(n_threads),
                                cores = n_cores,
                                seed = 123
      ))
      if (check_convergence(null_model_mult, n_term = 3) == 'good') {

        # Calculate empirical alpha of multivariate test (LOOIC)
        null_loo = loo(null_model_mult)
        alt_loo = loo(model_mult)
        alpha_hat_mult = null_loo$estimates['looic', 'Estimate'] > alt_loo$estimates['looic', 'Estimate']
      }
      rm(null_model_mult)
    }
    rm(model_mult)
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}

est_rstanarm = function(data, resp_dist, with_rd_slope, n_cores = 4){

  # Specify random effects
  rd_term = if (with_rd_slope) '(1 + time | id)' else '(1 | id)'

  # Specify response distribution
  if (resp_dist == 'Bernoulli'){
    fami = binomial(link = 'logit')
  }
  if (resp_dist == 'Poisson'){
    fami = poisson(link = 'log')
  }

  # Fit full model on full dataset
  compute_time = system.time({
    model_full = try(stan_glmer(build_formula('y_full', 'time + group + time:group', rd_term),
                                family = fami,
                                data = data,
                                cores = n_cores,
                                seed = 123
    ))})[['elapsed']]

  conv_status = check_convergence(model_full, fixed_names = c('time', 'group', 'time:group'), with_rd_slope = with_rd_slope)

  if (conv_status == 'good'){

    # Record betas and tau after fitting full model on full dataset
    beta0_hat = model_full$coefficients[['(Intercept)']]
    beta1_hat = model_full$coefficients[['time']]
    beta2_hat = model_full$coefficients[['group']]
    beta3_hat = model_full$coefficients[['time:group']]
    if (with_rd_slope){
      tau0_hat = attr(VarCorr(model_full)$id,'stddev')[['(Intercept)']]
      tau1_hat = attr(VarCorr(model_full)$id,'stddev')[['time']]
      rho01_hat = attr(VarCorr(model_full)$id,'correlation')['(Intercept)', 'time']
    } else {
      tau_hat = as.numeric(attr(VarCorr(model_full)$id,'stddev'))
    }
    # Calculate empirical power of univariate test (CI)
    power_hat_uni = check_CI(model_full, 'time:group')

    # Fit null model on full dataset
    null_model_full = try(stan_glmer(build_formula('y_full', 'time', rd_term),
                                     family = fami,
                                     data = data,
                                     cores = n_cores,
                                     seed = 123
    ))
    if (check_convergence(null_model_full, fixed_names = c('time'), with_rd_slope = with_rd_slope) == 'good'){
      # Calculate empirical power of multivariate test (LOOIC)
      null_loo = loo(null_model_full)
      alt_loo = loo(model_full)
      power_hat_mult = null_loo$estimates['looic', 'Estimate'] > alt_loo$estimates['looic', 'Estimate']
    }
    rm(model_full, null_model_full)

    # Fit full model on dataset with beta3 = 0
    model_uni = try(stan_glmer(build_formula('y_0', 'time + group + time:group', rd_term),
                               family = fami,
                               data = data,
                               cores = n_cores,
                               seed = 123
    ))
    if (check_convergence(model_uni, fixed_names = c('time', 'group', 'time:group'), with_rd_slope = with_rd_slope) == 'good') {
      # Calculate empirical alpha of univariate test (CI)
      alpha_hat_uni = check_CI(model_uni, 'time:group')
    }
    rm(model_uni)
    # Fit full model on dataset with beta2 = beta3 = 0
    model_mult = try(stan_glmer(build_formula('y_00', 'time + group + time:group', rd_term),
                                family = fami,
                                data = data,
                                cores = n_cores,
                                seed = 123
    ))
    if (check_convergence(model_mult, fixed_names = c('time', 'group', 'time:group'), with_rd_slope = with_rd_slope) == 'good'){
      # Fit null model on dataset with beta2 = beta3 = 0
      null_model_mult = try(stan_glmer(build_formula('y_00', 'time', rd_term),
                                       family = fami,
                                       data = data,
                                       cores = n_cores,
                                       seed = 123
      ))
      if (check_convergence(null_model_mult, fixed_names = c('time'), with_rd_slope = with_rd_slope) == 'good'){
        # Calculate empirical alpha of multivariate test (LRT)
        null_loo = loo(null_model_mult)
        alt_loo = loo(model_mult)
        alpha_hat_mult = null_loo$estimates['looic', 'Estimate'] > alt_loo$estimates['looic', 'Estimate']
      }
    }
    rm(model_mult, null_model_mult)
  }
  # Save results
  result_names = c('compute_time', 'conv_status', 'beta0_hat', 'beta1_hat',
                   'beta2_hat', 'beta3_hat', 'tau0_hat', 'tau1_hat', 'rho01_hat',
                   'alpha_hat_uni', 'alpha_hat_mult', 'power_hat_uni',
                   'power_hat_mult')
  result = unlist(mget(result_names, ifnotfound = list(NA), inherits = FALSE))
  return(result)
}