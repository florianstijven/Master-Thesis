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
                              us0, us1, ut0, ut1,
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
                     SigmaT_det_root_inv, SigmaT_inv,
                     gammas0, gammas1, gammat0, gammat1,
                     knots0, knots1, knott0, knott1){
  #compute density for delta s
  if(delta_s < 0){
    Us1 = runif(n = N)
    S1 = qsurvspline(p = Us1, gamma = gammas1, knots = knots1, timescale = "log")
    S0 = S1 - delta_s
    Us0 = psurvspline(q = S0, gamma = gammas0, knots = knots0, timescale = "log")
    f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0, timescale = "log")
  }
  else{
    Us0 = runif(n = N)
    S0 = qsurvspline(p = Us0, gamma = gammas0, knots = knots0, timescale = "log")
    S1 = S0 + delta_s
    Us1 = psurvspline(q = S1, gamma = gammas1, knots = knots1, timescale = "log")
    f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1, timescale = "log")
  }
  f_deltaS = sum(pmap_dbl(.l = list(u0 = Us0, u1 = Us1, s0 = S0, s1 = S1, f_s = f_s), 
                          .f = density_fast_s0s1_rh,
                          Sigma_det_root_inv = SigmaS_det_root_inv, 
                          Sigma_inv = SigmaS_inv))/N
  
  #compute density for delta T
  if(delta_t < 0){
    Ut1 = runif(n = N)
    T1 = qsurvspline(p = Ut1, gamma = gammat1, knots = knott1, timescale = "log")
    T0 = T1 - delta_t
    Ut0 = psurvspline(q = T0, gamma = gammat0, knots = knott0, timescale = "log")
    f_t = dsurvspline(x = T0, gamma = gammat0, knots = knott0, timescale = "log")
  }
  else{
    Ut0 = runif(n = N)
    T0 = qsurvspline(p = Ut0, gamma = gammat0, knots = knott0, timescale = "log")
    T1 = T0 + delta_t
    Ut1 = psurvspline(q = T1, gamma = gammat1, knots = knott1, timescale = "log")
    f_t = dsurvspline(x = T1, gamma = gammat1, knots = knott1, timescale = "log")
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

r_h = function(Sigma, n = 1000, N = 2000,
               gammas0, gammas1, gammat0, gammat1,
               knots0, knots1, knott0, knott1){
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
  S_delta = -1*qsurvspline(p = u_vec[,1], gamma = gammas0, knots = knots0, timescale = "log") + 
    qsurvspline(p = u_vec[,2], gamma = gammas1, knots = knots1, timescale = "log")
  T_delta = -1*qsurvspline(p = u_vec[,3], gamma = gammat0, knots = knott0, timescale = "log") + 
    qsurvspline(p = u_vec[,4], gamma = gammat1, knots = knott1, timescale = "log")
  mut_inf = pmap_dbl(.l = list(delta_s = S_delta, delta_t = T_delta), 
                         .f = sample_rh, Sigma = Sigma, N = N,
                         Sigma_det_root_inv = Sigma_det_root_inv,
                         SigmaS_det_root_inv = SigmaS_det_root_inv,
                         SigmaT_det_root_inv = SigmaT_det_root_inv,
                         Sigma_inv = Sigma_inv,
                         SigmaS_inv = SigmaS_inv,
                         SigmaT_inv = SigmaT_inv,
                         gammas0 = gammas0, gammas1 = gammas1, 
                         gammat0 = gammat0, gammat1 = gammat1,
                         knots0 = knots0, knots1 = knots1, 
                         knott0 = knott0, knott1 = knott1)
  return(mean(mut_inf))
}


#FUNCTIONS SPECIFICALLY FOR NEW APPROACH
density_fast_s0s1_new = function(Sigma11, Sigma12, Sigma21,
                                 q00, q01, q10, q11,
                                 Sigma22_det_root_inv, Sigma22_inv, 
                                 f_s){
  q = matrix(c(q00, q01), nrow = 1)
  cop_dens = (Sigma22_det_root_inv*
                exp(-0.5*q%*%(Sigma22_inv - diag(1, 2))%*%t(q)))[1,1]
  P = pmvnorm(lower = c(q10, q11), upper = +Inf, 
              mean = as.vector(Sigma12%*%Sigma22_inv%*%t(q)), 
              sigma = Sigma11 - Sigma12%*%Sigma22_inv%*%Sigma21)
  return(P*cop_dens*f_s)
}

# density_deltaS_new = function(Sigma, delta_s, N){
#   if(delta_s < 0){
#     U = runif(n = N)
#     #PART 1
#     S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
#     S0 = S1 - delta_s
#     f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
#     q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
#     q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
#     q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
#     q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part2
#     S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
#     T0 = S1 - delta_s
#     f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
#     q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
#     q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
#     q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
#     q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part3
#     T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
#     S0 = T1 - delta_s
#     f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
#     q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
#     q01 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
#     q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
#     q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part4
#     T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
#     T0 = T1 - delta_s
#     f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott1)
#     q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
#     q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
#     q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
#     q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,4,1,2), c(3,4,1,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     return(part1 + part2 + part3 + part4)
#   }
#   else{
#     U = runif(n = N)
#     #PART 1
#     S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
#     S1 = S0 + delta_s
#     f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
#     q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
#     q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
#     q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
#     q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                             .f = density_fast_s0s1_new,
#                             Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                             Sigma22_inv = solve(Sigma22), 
#                             Sigma11 = Sigma11,
#                             Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part2
#     T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
#     S1 = T0 + delta_s
#     f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
#     q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
#     q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
#     q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
#     q10 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part3
#     S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
#     T1 = S0 + delta_s
#     f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott0)
#     q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
#     q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
#     q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
#     q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     #part4
#     T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
#     T1 = T0 + delta_s
#     f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott0)
#     q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
#     q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
#     q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
#     q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
#     finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
#     Sigma_new = Sigma[c(3,4,1,2), c(3,4,1,2)]
#     Sigma11 = Sigma_new[1:2, 1:2]
#     Sigma22 = Sigma_new[3:4, 3:4]
#     Sigma12 = Sigma_new[1:2, 3:4]
#     part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
#                                    f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]), 
#                          .f = density_fast_s0s1_new,
#                          Sigma22_det_root_inv = det(Sigma22)**-0.5, 
#                          Sigma22_inv = solve(Sigma22), 
#                          Sigma11 = Sigma11,
#                          Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
#     return(part1 + part2 + part3 + part4)
#   }
# 
# }



density_deltaS_new = function(Sigma, delta_s, N,
                                gammas0, gammas1, gammat0, gammat1,
                                knots0, knots1, knott0, knott1){
  if(delta_s < 0){
    U = runif(n = N)
    #PART 1
    S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
    S0 = S1 - delta_s
    f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part2
    S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
    T0 = S1 - delta_s
    f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part3
    T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
    S0 = T1 - delta_s
    f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part4
    T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
    T0 = T1 - delta_s
    f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,2,3,4), c(1,2,3,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    return(part1 + part2 + part3 + part4)
  }
  else{
    U = runif(n = N)
    #PART 1
    S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
    S1 = S0 + delta_s
    f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part2
    T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
    S1 = T0 + delta_s
    f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part3
    S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
    T1 = S0 + delta_s
    f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott1)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    # part4
    T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
    T1 = T0 + delta_s
    f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott1)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,2,3,4), c(1,2,3,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    return(part1 + part2 + part3 + part4)
  }
}

density_deltaST_new = function(Sigma, delta_s, delta_t, N,
                              gammas0, gammas1, gammat0, gammat1,
                              knots0, knots1, knott0, knott1){
  if(delta_s < 0){
    U = runif(n = N)
    V = runif(n = N)
    #PART 1
    S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
    S0 = S1 - delta_s
    T1 = qsurvspline(p = V, gamma = gammat1, knots = knott1)
    T0 = T1 - delta_t
    f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
    f_t = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
    
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part2
    S1 = qsurvspline(p = U, gamma = gammas1, knots = knots1)
    T0 = S1 - delta_s
    f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part3
    T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
    S0 = T1 - delta_s
    f_s = dsurvspline(x = S0, gamma = gammas0, knots = knots0)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part4
    T1 = qsurvspline(p = U, gamma = gammat1, knots = knott1)
    T0 = T1 - delta_s
    f_s = dsurvspline(x = T0, gamma = gammat0, knots = knott0)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,2,3,4), c(1,2,3,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    return(part1 + part2 + part3 + part4)
  }
  else{
    U = runif(n = N)
    #PART 1
    S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
    S1 = S0 + delta_s
    f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,4,1,2),c(3,4,1,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part1 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part2
    T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
    S1 = T0 + delta_s
    f_s = dsurvspline(x = S1, gamma = gammas1, knots = knots1)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = S1, gamma = gammas1, knots = knots1))
    q11 = qnorm(p = psurvspline(q = S1, gamma = gammat1, knots = knott1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,4,3,2), c(1,4,3,2)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part2 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    #part3
    S0 = qsurvspline(p = U, gamma = gammas0, knots = knots0)
    T1 = S0 + delta_s
    f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott1)
    q00 = qnorm(p = psurvspline(q = S0, gamma = gammas0, knots = knots0))
    q10 = qnorm(p = psurvspline(q = S0, gamma = gammat0, knots = knott0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(3,2,1,4), c(3,2,1,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part3 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    # part4
    T0 = qsurvspline(p = U, gamma = gammat0, knots = knott0)
    T1 = T0 + delta_s
    f_s = dsurvspline(x = T1, gamma = gammat1, knots = knott1)
    q00 = qnorm(p = psurvspline(q = T0, gamma = gammat0, knots = knott0))
    q10 = qnorm(p = psurvspline(q = T0, gamma = gammas0, knots = knots0))
    q01 = qnorm(p = psurvspline(q = T1, gamma = gammat1, knots = knott1))
    q11 = qnorm(p = psurvspline(q = T1, gamma = gammas1, knots = knots1))
    finite = is.finite(q00) & is.finite(q01) & is.finite(f_s) & is.finite(q10) & is.finite(q11)
    Sigma_new = Sigma[c(1,2,3,4), c(1,2,3,4)]
    Sigma11 = Sigma_new[1:2, 1:2]
    Sigma22 = Sigma_new[3:4, 3:4]
    Sigma12 = Sigma_new[1:2, 3:4]
    part4 = sum(pmap_dbl(.l = list(q00 = q00[finite], q01 = q01[finite],
                                   f_s = f_s[finite], q10 = q10[finite], q11 = q11[finite]),
                         .f = density_fast_s0s1_new,
                         Sigma22_det_root_inv = det(Sigma22)**-0.5,
                         Sigma22_inv = solve(Sigma22),
                         Sigma11 = Sigma11,
                         Sigma12 = Sigma12, Sigma21 = t(Sigma12)))/N
    return(part1 + part2 + part3 + part4)
  }
}
