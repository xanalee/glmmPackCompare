official_conv_ratio = function(converge_v){
  conv_ratio = sum(converge_v %in% c('extra_fail', 'good'))/length(converge_v)
  return(conv_ratio)
}

real_conv_ratio = function(converge_v){
  conv_ratio = sum(converge_v == 'good')/length(converge_v)
  return(conv_ratio)
}

rmse = function(theta_v, true_theta){
  v_theta_v = na.omit(theta_v)
  return(sqrt(mean((v_theta_v - true_theta)^2)))
}