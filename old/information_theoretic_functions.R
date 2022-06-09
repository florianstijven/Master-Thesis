density_fast_s0s1_rh = function(Sigma_det_root_inv, Sigma_inv, 
                                u0, u1, s0, s1, f_s){
  qs0 = qnorm(u0)
  qs1 = qnorm(u1)
  z = matrix(c(qs0, qs1), nrow = 1)
  cop_dens = (Sigma_det_root_inv*
                exp(-0.5*z%*%(Sigma_inv - diag(1, 2))%*%t(z)))[1,1]
  return(cop_dens*f_s)
}

density_fast_t0t1_rh = function(Sigma_det_root_inv, Sigma_inv, 
                                u0, u1, t0, t1, f_t){
  qt0 = qnorm(u0)
  qt1 = qnorm(u1)
  z = matrix(c(qt0, qt1), nrow = 1)
  cop_dens = (Sigma_det_root_inv*
                exp(-0.5*z%*%(Sigma_inv - diag(1, 2))%*%t(z)))[1,1]
  return(cop_dens*f_t)
}

density_fast_st_rh = function(Sigma_det_root_inv, Sigma_inv, 
                              us0, us1, s0, s1, ut0, ut1, t0, t1,
                              f_t, f_s){

  qs0 = qnorm(us0)
  qs1 = qnorm(us1)
  qt0 = qnorm(ut0)
  qt1 = qnorm(ut1)
  z = matrix(c(qs0, qs1, qt0, qt1), nrow = 1)
  cop_dens = (Sigma_det_root_inv*
                exp(-0.5*z%*%(Sigma_inv - diag(1, 4))%*%t(z)))[1,1]
  return(cop_dens*f_t*f_s)
}


sample_rh = function(delta_s, delta_t, Sigma, N = 1000,
                     Sigma_det_root_inv, Sigma_inv,
                     SigmaS_det_root_inv, SigmaS_inv,
                     SigmaT_det_root_inv, SigmaT_inv){
  #compute density for delta s
  if(delta_s < 0){
    Us1 = runif(n = N)
    S1 = predict(fit_s1, type = "quantile", newdata = data[1,], p = Us1)$.pred[[1]]$.pred
    S0 = S1 - delta_s
    Us0 = 1 - predict(fit_s0, type = "survival", newdata = data[1,], times = S0)$.pred[[1]]$.pred
    f_s = (predict(fit_s0, type = "hazard", newdata = data[1,], times = S0)$.pred[[1]]$.pred)*
      (predict(fit_s0, type = "survival", newdata = data[1,], times = S0)$.pred[[1]]$.pred)
  }
  else{
    Us0 = runif(n = N)
    S0 = predict(fit_s0, type = "quantile", newdata = data[1,], p = Us0)$.pred[[1]]$.pred
    S1 = S0 + delta_s
    Us1 = 1 - predict(fit_s1, type = "survival", newdata = data[1,], times = S1)$.pred[[1]]$.pred
    f_s = (predict(fit_s1, type = "hazard", newdata = data[1,], times = S1)$.pred[[1]]$.pred)*
      (predict(fit_s1, type = "survival", newdata = data[1,], times = S1)$.pred[[1]]$.pred)
  }
  f_deltaS = sum(pmap_dbl(.l = list(u0 = Us0, u1 = Us1, s0 = S0, s1 = S1, f_s = f_s), 
                           .f = density_fast_s0s1_rh,
                           Sigma_det_root_inv = SigmaS_det_root_inv, 
                           Sigma_inv = SigmaS_inv))/N
  
  #compute density for delta T
  if(delta_t < 0){
    Ut1 = runif(n = N)
    T1 = predict(fit_t1, type = "quantile", newdata = data[1,], p = Ut1)$.pred[[1]]$.pred
    T0 = T1 - delta_t
    Ut0 = 1 - predict(fit_t0, type = "survival", newdata = data[1,], times = T0)$.pred[[1]]$.pred
    f_t = (predict(fit_t0, type = "hazard", newdata = data[1,], times = T0)$.pred[[1]]$.pred)*
      (predict(fit_t0, type = "survival", newdata = data[1,], times = T0)$.pred[[1]]$.pred) 
  }
  else{
    Ut0 = runif(n = N)
    T0 = predict(fit_t0, type = "quantile", newdata = data[4,], p = Ut0)$.pred[[1]]$.pred
    T1 = T0 + delta_t
    Ut1 = 1 - predict(fit_t1, type = "survival", newdata = data[1,], times = T1)$.pred[[1]]$.pred
    f_t = (predict(fit_t1, type = "hazard", newdata = data[1,], times = T1)$.pred[[1]]$.pred)*
      (predict(fit_t1, type = "survival", newdata = data[1,], times = T1)$.pred[[1]]$.pred) 
  }
  f_deltaT = sum(pmap_dbl(.l = list(u0 = Ut0, u1 = Ut1, t0 = T0, t1 = T1, f_t = f_t), 
                           .f = density_fast_t0t1_rh,
                           Sigma_det_root_inv = SigmaT_det_root_inv, 
                           Sigma_inv = SigmaT_inv))/N
  #compute density for delta S and T
  f_deltaST = sum(pmap_dbl(.l = list(us0 = Us0, us1 = Us1, s0 = S0, s1 = S1, 
                                      ut0 = Ut0, ut1 = Ut1, t0 = T0, t1 = T1,
                                      f_s = f_s, f_t = f_t), 
                            .f = density_fast_st_rh,
                            Sigma_det_root_inv = Sigma_det_root_inv, 
                            Sigma_inv = Sigma_inv))/N
  return(log(f_deltaST/(f_deltaS*f_deltaT)))
}

r_h = function(Sigma, n = 1000, N = 1000){
  #compute different matrices once
  rho13 = Sigma[1, 3]
  rho24 = Sigma[2, 4]
  Sigma_det_root_inv = det(Sigma)**-0.5
  Sigma_inv = solve(Sigma)
  SigmaS = matrix(c(1, rho13, rho13, 1), nrow = 2)
  SigmaS_det_root_inv = det(SigmaS)**-0.5
  SigmaS_inv = solve(SigmaS)
  SigmaT = matrix(c(1, rho24, rho24, 1), nrow = 2)
  SigmaT_det_root_inv = det(SigmaT)**-0.5
  SigmaT_inv = solve(SigmaT)
  
  #sample from the full joint distribution to get a sample
  #of delta S and delta T
  q_vec = rmvnorm(n = n, sigma = Sigma)
  u_vec = pnorm(q = q_vec)
  S_delta = predict(fit_s0, type = "quantile", newdata = data[4,], p = u_vec[,1])$.pred[[1]]$.pred -
    predict(fit_s1, type = "quantile", newdata = data[1,], p = u_vec[,2])$.pred[[1]]$.pred
  T_delta = predict(fit_t0, type = "quantile", newdata = data[4,], p = u_vec[,3])$.pred[[1]]$.pred -
    predict(fit_t1, type = "quantile", newdata = data[1,], p = u_vec[,4])$.pred[[1]]$.pred
  mut_inf = sum(pmap_dbl(.l = list(delta_s = S_delta, delta_t = T_delta), 
                          .f = sample_rh, Sigma = Sigma, N = N,
                          Sigma_det_root_inv = Sigma_det_root_inv,
                          SigmaS_det_root_inv = SigmaS_det_root_inv,
                          SigmaT_det_root_inv = SigmaT_det_root_inv,
                          Sigma_inv = Sigma_inv,
                          SigmaS_inv = SigmaS_inv,
                          SigmaT_inv = SigmaT_inv))/n
  return(mut_inf)
}

