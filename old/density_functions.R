library(purrr)
library(flexsurv)
library(mvtnorm) 

#gaussian copula density for 4-variate vector
copula_density = function(Sigma, x_vec){
  #transform uniform variables to standard normal variables
  q_vec = qnorm(p = x_vec)
  #put normal transformed variables into appropriate vector
  q_vec = matrix(q_vec, nrow = 4)
  #inverse of correlation matrix
  inv_Sigma = tryCatch(exp = solve(Sigma), finally = NaN)
  #determinant of correlation matrix
  det_Sigma = det(Sigma)
  #required product in the exponent
  exp_product = -0.5*t(q_vec)%*%(inv_Sigma - diag(1, 4))%*%q_vec
  #return copula density
  return(as.numeric((det_Sigma**-0.5)*exp(exp_product)))
}

#gaussian copula density for 2-variate vector 
#(marginalised over two other two variables)
copula_density_marg = function(rho, x_vec){
  #correlation matrix Sigma
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  #transform uniform variables to standard normal variables
  q_vec = qnorm(p = x_vec)
  #put normal transformed variables into appropriate vector
  q_vec = matrix(q_vec, nrow = 2)
  #inverse of correlation matrix
  inv_Sigma = tryCatch(expr = solve(Sigma), finally = NaN) 
  #determinant of correlation matrix
  det_Sigma = det(Sigma)
  #required product in the exponent
  exp_product = -0.5*t(q_vec)%*%(inv_Sigma - diag(1, 2))%*%q_vec
  #return copula density
  return(as.numeric((det_Sigma**-0.5)*exp(exp_product)))
}


#SPECIFIC MARGINALS
#compute density for a censored observation for copula model
model_density_copula_marg_single_cens = function(rho, x_vec, 
                                             cens, prec = 1000){
  #transform uniform variables to standard normal variables
  q_vec = qnorm(p = x_vec)
  #density for uncensored observation
  dens_uncens = dnorm(q_vec[-cens])
  #probability for censored observation, given censored observation
  prob_cens = 1 - pnorm(q = q_vec[cens], mean = rho*q_vec[-cens], sd = 1 - rho**2)
  return(dens_uncens*prob_cens)
}

#compute density of both observations are censored copula model
model_density_copula_marg_double_cens = function(rho, x_vec){
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  q_vec = qnorm(p = x_vec)
  return(pmvnorm(mean = c(0, 0), corr = Sigma, lower = q_vec))
}



#four variate gaussian copula model density for weibull marginals
model_density_wb = function(Sigma, t_vec, lambda_vec, k_vec){
  x_vec = pmap_dbl(.l = list(shape = k_vec, scale = lambda_vec, q = t_vec), 
                   .f = pweibull)
  t_dens = pmap_dbl(.l = list(shape = k_vec, scale = lambda_vec, x = t_vec),
                    .f = dweibull)
  return(copula_density(Sigma, x_vec)*t_dens[1]*t_dens[2]*t_dens[3]*t_dens[4])
}

#compute density for a censored observation for weibull
model_density_wb_marg_single_cens = function(rho, t_vec, lambda_vec, k_vec, 
                                             cens, prec = 1000){
  upper = qweibull(p = 0.95, shape = k_vec[cens], scale = lambda_vec[cens])
  grid = seq(from = t_vec[cens], to = upper, length.out = prec)
  h = grid[2] - grid[1]
  grid = grid + h/2
  dens = 0
  for(i in 1:prec){
    dens = dens + h*model_density_wb_marg(rho, c(t_vec[1], grid[i]), 
                                            lambda_vec, k_vec)
  }
  
  upper_2 = qweibull(p = 1 - (0.05/prec), shape = k_vec[cens], scale = lambda_vec[cens])
  grid_2 = seq(from = upper + h, to = upper_2, length.out = prec)
  h_2 = grid_2[2] - grid_2[1]
  grid_2 = grid_2 + h_2/2
  for(i in 1:prec){
    dens = dens + h_2*model_density_wb_marg(rho, c(t_vec[1], grid_2[i]), 
                                              lambda_vec, k_vec)
  }  
  return(dens)
}

#compute density of both observations are censored weibull
model_density_wb_marg_double_cens = function(rho, t_vec, lambda_vec, k_vec){
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  x_vec = pmap_dbl(.l = list(x = t_vec, scale = lambda_vec, shape = k_vec), 
                   .f = pweibull)
  q_vec = qnorm(p = x_vec)
  return(pmvnorm(mean = c(0, 0), corr = Sigma, lower = q_vec))
}



#four variate gaussian copula model density for loglogistic marginals
model_density_llog = function(Sigma, t_vec, alpha_vec, beta_vec){
  x_vec = pmap_dbl(.l = list(shape = beta_vec, scale = alpha_vec, q = t_vec), 
                   .f = pllogis)
  t_dens = pmap_dbl(.l = list(shape = beta_vec, scale = alpha_vec, x = t_vec),
                    .f = dllogis)
  return(copula_density(Sigma, x_vec)*t_dens[1]*t_dens[2]*t_dens[3]*t_dens[4])
}

#compute density for a censored observation for loglogistic
model_density_llog_marg_single_cens = function(rho, t_vec, alpha_vec, beta_vec, 
                                             cens, prec = 50){
  upper = qllogis(p = 0.95, shape = beta_vec[cens], scale = alpha_vec[cens])
  grid = seq(from = t_vec[cens], to = upper, length.out = prec)
  h = grid[2] - grid[1]
  grid = grid + h/2
  dens = 0
  for(i in 1:prec){
    dens = dens + h*model_density_llog_marg(rho, c(t_vec[1], grid[i]), 
                                            alpha_vec, beta_vec)
  }
  
  upper_2 = qllogis(p = 1 - (0.05/prec), shape = beta_vec[cens], scale = alpha_vec[cens])
  grid_2 = seq(from = upper + h, to = upper_2, length.out = prec)
  h_2 = grid_2[2] - grid_2[1]
  grid_2 = grid_2 + h_2/2
  for(i in 1:prec){
    dens = dens + h_2*model_density_llog_marg(rho, c(t_vec[1], grid_2[i]), 
                                            alpha_vec, beta_vec)
  }  
  return(dens)
}

#compute density of both observations are censored loglogistic
model_density_llog_marg_double_cens = function(rho, t_vec, alpha_vec, beta_vec){
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  x_vec = pmap_dbl(.l = list(q = t_vec, scale = alpha_vec, shape = beta_vec), 
                   .f = pllogis)
  q_vec = qnorm(p = x_vec)
  return(pmvnorm(mean = c(0, 0), corr = Sigma, lower = q_vec))
}



#four variate gaussian copula model density for generalized gamma marginals
model_density_gg = function(Sigma, t_vec, mu_vec, sigma_vec, Q_vec){
  x_vec = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, Q = Q_vec,
                             q = t_vec), 
                   .f = pgengamma)
  t_dens = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, Q = Q_vec,
                              x = t_vec),
                    .f = dgengamma)
  return(copula_density(Sigma, x_vec)*t_dens[1]*t_dens[2]*t_dens[3]*t_dens[4])
}

#compute density for a censored observation for generalized gamma
model_density_gg_marg_single_cens = function(rho, t_vec, mu_vec, sigma_vec, Q_vec,
                                             cens, prec = 50){
  upper = qgengamma(p = 0.95, mu = mu_vec[cens], sigma = sigma_vec[cens],
                    Q = Q_vec[cens])
  grid = seq(from = t_vec[cens], to = upper, length.out = prec)
  h = grid[2] - grid[1]
  grid = grid + h/2
  dens = 0
  for(i in 1:prec){
    dens = dens + h*model_density_gg_marg(rho, c(t_vec[1], grid[i]), 
                                          mu_vec, sigma_vec, Q_vec)
  }
  
  upper_2 = qgengamma(p = 1 - (0.05/prec), mu = mu_vec[cens], sigma = sigma_vec[cens],
                      Q = Q_vec[cens])
  grid_2 = seq(from = upper + h, to = upper_2, length.out = prec)
  h_2 = grid_2[2] - grid_2[1]
  grid_2 = grid_2 + h_2/2
  for(i in 1:prec){
    dens = dens + h_2*model_density_gg_marg(rho, c(t_vec[1], grid_2[i]), 
                                            mu_vec, sigma_vec, Q_vec)
  }  
  return(dens)
}

#compute density of both observations are censored generalized gamma
model_density_gg_marg_double_cens = function(rho, t_vec, 
                                               mu_vec, sigma_vec, Q_vec){
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  x_vec = pmap_dbl(.l = list(q = t_vec, mu = mu_vec, 
                             sigma = sigma_vec, Q = Q_vec), 
                   .f = pgengamma)
  q_vec = qnorm(p = x_vec)
  return(pmvnorm(mean = c(0, 0), corr = Sigma, lower = q_vec))
}




#four variate gaussian copula model density for generalized F
model_density_genf = function(Sigma, t_vec, mu_vec, sigma_vec, Q_vec, P_vec){
  x_vec = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, Q = Q_vec, P = P_vec,
                             q = t_vec), 
                   .f = pgenf)
  t_dens = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, Q = Q_vec, P = P_vec,
                              x = t_vec),
                    .f = dgenf)
  return(copula_density(Sigma, x_vec)*t_dens[1]*t_dens[2]*t_dens[3]*t_dens[4])
}

#compute density for a censored observation for generalized F
model_density_genf_marg_single_cens = function(rho, t_vec, mu_vec, sigma_vec, Q_vec, P_vec,
                                             cens, prec = 50){
  upper = qgenf(p = 0.95, mu = mu_vec[cens], sigma = sigma_vec[cens],
                    Q = Q_vec[cens], P = P_vec[cens])
  grid = seq(from = t_vec[cens], to = upper, length.out = prec)
  h = grid[2] - grid[1]
  grid = grid + h/2
  dens = 0
  for(i in 1:prec){
    dens = dens + h*model_density_genf_marg(rho, c(t_vec[1], grid[i]), 
                                          mu_vec, sigma_vec, Q_vec, P_vec)
  }
  
  upper_2 = qgenf(p = 1 - (0.05/prec), mu = mu_vec[cens], sigma = sigma_vec[cens],
                      Q = Q_vec[cens], P = P_vec[cens])
  grid_2 = seq(from = upper + h, to = upper_2, length.out = prec)
  h_2 = grid_2[2] - grid_2[1]
  grid_2 = grid_2 + h_2/2
  for(i in 1:prec){
    dens = dens + h_2*model_density_genf_marg(rho, c(t_vec[1], grid_2[i]), 
                                            mu_vec, sigma_vec, Q_vec, P_vec)
  }  
  return(dens)
}

#compute density of both observations are censored generalized F
model_density_genf_marg_double_cens = function(rho, t_vec, 
                                             mu_vec, sigma_vec, Q_vec, P_vec){
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  x_vec = pmap_dbl(.l = list(q = t_vec, mu = mu_vec, 
                             sigma = sigma_vec, Q = Q_vec, P = P_vec), 
                   .f = pgenf)
  q_vec = qnorm(p = x_vec)
  return(pmvnorm(mean = c(0, 0), corr = Sigma, lower = q_vec))
}



#two variate gaussian copula model density for weibull marginals
model_density_wb_marg = function(rho, t_vec, lambda_vec, k_vec){
  x_vec = pmap_dbl(.l = list(shape = k_vec, scale = lambda_vec, q = t_vec), 
                   .f = pweibull)
  t_dens = pmap_dbl(.l = list(shape = k_vec, scale = lambda_vec, x = t_vec),
                    .f = dweibull)
  return(copula_density_marg(rho, x_vec)*t_dens[1]*t_dens[2])  
}

#two variate gaussian copula model density for loglogistic marginals
model_density_llog_marg = function(rho, t_vec, alpha_vec, beta_vec){
  x_vec = pmap_dbl(.l = list(shape = beta_vec, scale = alpha_vec, q = t_vec), 
                   .f = pllogis)
  t_dens = pmap_dbl(.l = list(shape = beta_vec, scale = alpha_vec, x = t_vec),
                    .f = dllogis)
  return(copula_density_marg(rho, x_vec)*t_dens[1]*t_dens[2])  
}

#two variate gaussian copula model density for generalized gamma marginals
model_density_gg_marg = function(rho, t_vec, mu_vec, sigma_vec, Q_vec){
  x_vec = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, 
                             Q = Q_vec, q = t_vec), 
                   .f = pgengamma)
  t_dens = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, 
                              Q = Q_vec, x = t_vec),
                    .f = dgengamma)
  return(copula_density_marg(rho, x_vec)*t_dens[1]*t_dens[2])  
}

#two variate gaussian copula model density for generalized F marginals
model_density_genf_marg = function(rho, t_vec, mu_vec, sigma_vec, Q_vec, P_vec){
  x_vec = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, 
                             Q = Q_vec, P = P_vec, q = t_vec), 
                   .f = pgenf)
  t_dens = pmap_dbl(.l = list(mu = mu_vec, sigma = sigma_vec, 
                              Q = Q_vec, P = P_vec, x = t_vec),
                    .f = dgenf)
  return(copula_density_marg(rho, x_vec)*t_dens[1]*t_dens[2])  
}

#FUNCTIONS FOR FULL DATA SETS
full_loglik_wb = function(data, Sigma, lambda_vec, k_vec){
  #loglikelihood contributions for individual observations
  loglik_obs = sapply(X = data, MARGIN = 1, FUN = model_density_wb,
                      lambda_vec = lambda_vec, k_vec = k_vec)
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}

full_loglik_llog = function(data, Sigma, alpha_vec, beta_vec){
  #loglikelihood contributions for individual observations
  loglik_obs = sapply(X = data, MARGIN = 1, FUN = model_density_llog,
                      alpha_vec = alpha_vec, beta_vec = beta_vec)
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs)) 
}




#FUNCTIONS FOR REAL DATA SETS
marg_loglik_wb = function(data, rho_vec, lambda_vec, k_vec){
  n = nrow(data)
  trt = data[, 3]
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    trt_ind = trt[i] + 1
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_wb_marg(rho = rho_vec[trt_ind], t_vec = t_vec,
                                            lambda_vec = lambda_vec[c(trt_ind, trt_ind + 2)],
                                            k_vec = k_vec[c(trt_ind, trt_ind + 2)])      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_wb_marg_double_cens(rho = rho_vec[trt_ind],
                                                        t_vec = t_vec, 
                                                        lambda_vec = lambda_vec[c(trt_ind, trt_ind + 2)],
                                                        k_vec = k_vec[c(trt_ind, trt_ind + 2)]
                                                        )
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_wb_marg_single_cens(rho = rho_vec[trt_ind],
                                                        t_vec = t_vec, 
                                                        lambda_vec = lambda_vec[c(trt_ind, trt_ind + 2)],
                                                        k_vec = k_vec[c(trt_ind, trt_ind + 2)],
                                                        cens = cens)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_wb_help = function(parms, data){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho_vec = (exp(parms[1:2]) - 1)/(exp(parms[1:2]) + 1)
  lambda_vec = exp(parms[3:6])
  k_vec = exp(parms[7:10])
  return(-2*marg_loglik_wb(data, rho_vec, lambda_vec, k_vec))
}

marg_loglik_llog = function(data, rho_vec, alpha_vec, beta_vec){
  n = nrow(data)
  trt = data[, 3]
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    trt_ind = trt[i] + 1
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_llog_marg(rho = rho_vec[trt_ind], t_vec = t_vec,
                                            alpha_vec = alpha_vec[c(trt_ind, trt_ind + 2)],
                                            beta_vec = beta_vec[c(trt_ind, trt_ind + 2)])      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_llog_marg_double_cens(rho = rho_vec[trt_ind],
                                                        t_vec = t_vec, 
                                                        alpha_vec = alpha_vec[c(trt_ind, trt_ind + 2)],
                                                        beta_vec = beta_vec[c(trt_ind, trt_ind + 2)])
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_llog_marg_single_cens(rho = rho_vec[trt_ind],
                                                        t_vec = t_vec, 
                                                        alpha_vec = alpha_vec[c(trt_ind, trt_ind + 2)],
                                                        beta_vec = beta_vec[c(trt_ind, trt_ind + 2)],
                                                        cens = cens)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_llog_help = function(parms, data){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho_vec = (exp(parms[1:2]) - 1)/(exp(parms[1:2]) + 1)
  alpha_vec = exp(parms[3:6])
  beta_vec = exp(parms[7:10])
  return(-2*marg_loglik_llog(data, rho_vec, alpha_vec, beta_vec))
}

fit_wb = function(data, inits = c(2, 2.1, log(1), log(5), log(12), log(15), 
                                  log(1), log(1), log(2), log(2)),
                  maxit = 100){
  optim_out = optim(par = inits, fn = marg_min2loglik_wb_help, method = "L-BFGS-B",
                    data = data, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)
}

fit_llog = function(data, inits = c(2, 2.1, log(1), log(1.2), log(2), log(3), 
                                    log(1), log(1), log(2), log(2)),
                    maxit = 100){
  optim_out = optim(par = inits, fn = marg_min2loglik_llog_help, method = "L-BFGS-B", 
                    data = data, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)  
}





#FUNCTIONS FOR REAL DATA SETS, all in parts
#COPULA
marg_loglik_copula_part = function(data, rho, prec){
  n = nrow(data)
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    x_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = copula_density_marg(rho = rho, x_vec = x_vec)      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_copula_marg_double_cens(rho = rho,
                                                        x_vec = x_vec)
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_copula_marg_single_cens(rho = rho,
                                                        x_vec = x_vec,
                                                        cens = cens, prec = prec)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_copula_help = function(data, parms, prec){
  rho = (exp(parms[1]) - 1)/(exp(parms[1]) + 1)
  return(-2*marg_loglik_copula_part(data, rho, prec))
}
#WEIBULL
marg_loglik_wb_part = function(data, rho, lambda_vec, k_vec, prec){
  n = nrow(data)
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_wb_marg(rho = rho, t_vec = t_vec,
                                            lambda_vec = lambda_vec,
                                            k_vec = k_vec)      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_wb_marg_double_cens(rho = rho,
                                                        t_vec = t_vec, 
                                                        lambda_vec = lambda_vec,
                                                        k_vec = k_vec)
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_wb_marg_single_cens(rho = rho,
                                                        t_vec = t_vec, 
                                                        lambda_vec = lambda_vec,
                                                        k_vec = k_vec,
                                                        cens = cens, prec = prec)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_wb_help_part = function(parms, data, prec){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho = (exp(parms[1]) - 1)/(exp(parms[1]) + 1)
  lambda_vec = exp(parms[2:3])
  k_vec = exp(parms[4:5])
  return(-2*marg_loglik_wb_part(data, rho, lambda_vec, k_vec, prec))
}

#LOGLOGISTIC
marg_loglik_llog_part = function(data, rho, alpha_vec, beta_vec, prec){
  n = nrow(data)
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_llog_marg(rho = rho, t_vec = t_vec,
                                              alpha_vec = alpha_vec,
                                              beta_vec = beta_vec)      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_llog_marg_double_cens(rho = rho,
                                                          t_vec = t_vec, 
                                                          alpha_vec = alpha_vec,
                                                          beta_vec = beta_vec)
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_llog_marg_single_cens(rho = rho,
                                                          t_vec = t_vec, 
                                                          alpha_vec = alpha_vec,
                                                          beta_vec = beta_vec,
                                                          cens = cens, prec = prec)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_llog_help_part = function(parms, data, prec){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho = (exp(parms[1]) - 1)/(exp(parms[1]) + 1)
  alpha_vec = exp(parms[2:3])
  beta_vec = exp(parms[4:5])
  return(-2*marg_loglik_llog_part(data, rho, alpha_vec, beta_vec, prec))
}

#GENERALIZED GAMMA
marg_loglik_gg_part = function(data, rho, mu_vec, sigma_vec, Q_vec, prec){
  n = nrow(data)
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_gg_marg(rho = rho, t_vec = t_vec,
                                            mu_vec = mu_vec, sigma_vec = sigma_vec,
                                            Q_vec = Q_vec)      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_gg_marg_double_cens(rho = rho,
                                                          t_vec = t_vec, 
                                                        mu_vec = mu_vec, sigma_vec = sigma_vec,
                                                        Q_vec = Q_vec)
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_gg_marg_single_cens(rho = rho,
                                                          t_vec = t_vec, 
                                                        mu_vec = mu_vec, sigma_vec = sigma_vec,
                                                        Q_vec = Q_vec,
                                                          cens = cens, prec = prec)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_gg_help_part = function(parms, data, prec){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho = (exp(parms[1]) - 1)/(exp(parms[1]) + 1)
  mu_vec = parms[2:3]
  sigma_vec = exp(parms[4:5])
  Q_vec = exp(parms[6:7])
  
  return(-2*marg_loglik_gg_part(data, rho, mu_vec, sigma_vec, Q_vec, prec))
}

#GENERALIZED F
marg_loglik_genf_part = function(data, rho, mu_vec, sigma_vec, Q_vec, P_vec, prec){
  n = nrow(data)
  s_cens = data[, 4]
  t_cens = data[, 5]
  #loglikelihood contributions for individual observations
  loglik_obs = rep(0, n)
  for(i in 1:n){
    t_vec = as.numeric(data[i, 1:2])
    if((s_cens[i] == 1) & (t_cens[i] = 1)){
      loglik_obs[i] = model_density_genf_marg(rho = rho, t_vec = t_vec,
                                            mu_vec = mu_vec, sigma_vec = sigma_vec,
                                            Q_vec = Q_vec, P_vec = P_vec)      
    }
    else if((s_cens[i] == 0) & (t_cens[i] = 0)){
      loglik_obs[i] = model_density_genf_marg_double_cens(rho = rho,
                                                        t_vec = t_vec, 
                                                        mu_vec = mu_vec, sigma_vec = sigma_vec,
                                                        Q_vec = Q_vec, P_vec = P_vec)
    }
    else{
      if(s_cens[i] == 0) cens = 1
      else cens = 2
      loglik_obs[i] = model_density_genf_marg_single_cens(rho = rho,
                                                        t_vec = t_vec, 
                                                        mu_vec = mu_vec, sigma_vec = sigma_vec,
                                                        Q_vec = Q_vec, P_vec = P_vec,
                                                        cens = cens, prec = prec)
    }
  }
  loglik_obs = log(loglik_obs)
  return(sum(loglik_obs))
}
#wrapper to fit model
marg_min2loglik_genf_help_part = function(parms, data, prec){
  #parms contains transformed values for the parameters
  #such that these parameters are defined on the entire real line
  rho = (exp(parms[1]) - 1)/(exp(parms[1]) + 1)
  mu_vec = parms[2:3]
  sigma_vec = exp(parms[4:5])
  Q_vec = exp(parms[6:7])
  P_vec = exp(parms[8:9])
  
  return(-2*marg_loglik_genf_part(data, rho, mu_vec, sigma_vec, Q_vec, P_vec, prec))
}



#FITTING FUNCTIONS
fit_wb_part = function(data, inits = c(2, 2.1, log(1), log(5), log(1), log(1)),
                       maxit = 100, prec = 50){
  optim_out = optim(par = inits, fn = marg_min2loglik_wb_help_part, method = "L-BFGS-B",
                    data = data, prec = prec, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)
}

fit_llog_part = function(data, inits = c(2, 2.1, log(1), log(1.2), log(1), log(1)),
                         maxit = 100, prec = 50){
  optim_out = optim(par = inits, fn = marg_min2loglik_llog_help_part, method = "L-BFGS-B", 
                    data = data, prec = prec, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)  
}

fit_gg_part = function(data, inits = c(2, 1, 1, 1, 1, 1, 1),
                       maxit = 100, prec = 50){
  optim_out = optim(par = inits, fn = marg_min2loglik_gg_help_part, method = "L-BFGS-B", 
                    data = data, prec = prec, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)    
}

fit_genf_part = function(data, inits = c(2, 1, 1, 1, 1, 1, 1, 1, 1),
                       maxit = 100, prec = 50){
  optim_out = optim(par = inits, fn = marg_min2loglik_genf_help_part, method = "CG", 
                    data = data, prec = prec, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)    
}

fit_copula = function(data, inits = c(1), maxit = 100, prec = 50){
  optim_out = optim(par = inits, fn = marg_min2loglik_copula_help, method = "L-BFGS-B", 
                    data = data, prec = prec, control = list(maxit = maxit, REPORT = 2, trace = TRUE), 
                    hessian = TRUE)
  return(optim_out)    
}


