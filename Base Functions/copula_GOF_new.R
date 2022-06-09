marginal_gof_scr = function(fit_0_par, fit_1_par, 
                            knots0, knots1, knott0, knott1,
                            data, cop_type,
                            grid){
  par(mfrow = c(2,2))
  #goodness of fit of marginal survival function of OS
  Surv_t0 = 1 - psurvspline(q = grid, gamma = fit_0_par[5:8], knots = knott0)
  KM_est = survfit(Surv(Surv, SurvInd)~1, data = data, subset = data$Treat == 0)
  ggsurvplot(fit = KM_est, data = data) 
  plot(survfit(Surv(Surv, SurvInd)~1, data = data, subset = data$Treat == 0),
       xlim = c(min(grid), max(grid)), main = "OS (0)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, Surv_t0, col = "red")

  Surv_t1 = 1 - psurvspline(q = grid, gamma = fit_1_par[5:8], knots = knott1)
  plot(survfit(Surv(Surv, SurvInd)~1, data = data, subset = data$Treat == 1),
       xlim = c(min(grid), max(grid)), main = "OS (1)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, Surv_t1, col = "red")

  #goodness of fit for marginal distribution of PFS
  pfs_surv = function(s, gammas, gammat, knots, knott, theta){
    u = 1 - psurvspline(q = s, gamma = gammas, knots = knots)
    v = 1 - psurvspline(q = s, gamma = gammat, knots = knott)
    if(cop_type == "frank"){
      C = (-1/theta)*log(((1-exp(-theta)-(1-exp(-theta*u))*(1-exp(-theta*v))))/(1-exp(-theta)))
    }
    else if(cop_type == "gaussian"){
      rho = (exp(theta) - 1)/(exp(theta) + 1)
      Sigma = matrix(data = c(1, rho, rho, 1), ncol = 2)
      V = qnorm(c(u, v))
      C = pmvnorm(lower = -Inf,upper=V, sigma=Sigma, mean=c(0,0))[1]
    }
    else if(cop_type == "clayton"){
      C = (u^(-theta) + v^(-theta) - 1)^(-1/theta)
    }
    else if(cop_type == "gumbel"){
      C = exp(-((-log(u))^(theta)+(-log(v))^(theta))^(1/theta))
    }
    return(C)
  }

  probs0 = sapply(grid, pfs_surv, gammas = fit_0_par[1:4], gammat = fit_0_par[5:8],
                  knots = knots0, knott = knott0, theta = fit_0_par[9])
  KM_est = survfit(Surv(data$Pfs, pmax(data$PfsInd, data$SurvInd))~1, data = data,
                   subset = data$Treat == 0)
  ggsurvplot(KM_est) +
    scale_x_continuous(name = "time (months)",
                       limits = c(min(grid), max(grid))) +
    scale_y_continuous(name = "S(t)") +
    geom_line(data = data.frame(grid = grid, probs0 = probs0),
              aes(x = grid, y = probs0), color = "red")
    ggtitle("PFS (0)")

  probs1 = sapply(grid, pfs_surv, gammas = fit_1_par[1:4], gammat = fit_1_par[5:8],
                  knots = knots1, knott = knott1, theta = fit_1_par[9])
  plot(survfit(Surv(data$Pfs, pmax(data$PfsInd, data$SurvInd))~1, data = data,
               subset = data$Treat == 1),
       xlim = c(min(grid), max(grid)), main = "PFS (1)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, probs1, col = "red")
}

marginal_gof_no = function(fit_0_par, fit_1_par, 
                           knots0, knots1, knott0, knott1,
                           data, grid){
  par(mfrow = c(2,2))
  #goodness of fit of marginal survival function of OS
  Surv_t0 = 1 - psurvspline(q = grid, gamma = fit_0_par[5:8], knots = knott0)
  plot(survfit(Surv(Surv, SurvInd)~1, data = data, subset = data$Treat == 0),
       xlim = c(min(grid), max(grid)), main = "OS (0)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, Surv_t0, col = "red")
  
  Surv_t1 = 1 - psurvspline(q = grid, gamma = fit_1_par[5:8], knots = knott1)
  plot(survfit(Surv(Surv, SurvInd)~1, data = data, subset = data$Treat == 1),
       xlim = c(min(grid), max(grid)), main = "OS (1)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, Surv_t1, col = "red")
  
  probs0 = 1 - psurvspline(q = grid, gamma = fit_0_par[1:4], knots = knots0)
  plot(survfit(Surv(data$Pfs, pmax(data_pfs$PfsInd, data_pfs$SurvInd))~1, data = data_pfs, 
               subset = data_pfs$Treat == 0),
       xlim = c(min(grid), max(grid)), main = "PFS (0)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, probs0, col = "red")
  
  probs1 = 1 - psurvspline(q = grid, gamma = fit_1_par[1:4], knots = knots1)
  plot(survfit(Surv(data_pfs$Pfs, pmax(data_pfs$PfsInd, data_pfs$SurvInd))~1, data = data_pfs, 
               subset = data_pfs$Treat == 1),
       xlim = c(min(grid), max(grid)), main = "PFS (1)", xlab = "time (months)", ylab = "S(t)")
  lines(grid, probs1, col = "red")
}

association_gof_scr = function(fit_0_par, fit_1_par, 
                               knots0, knots1, knott0, knott1,
                               data, cop_type){
  #estimate survival function for censoring
  km_censor = flexsurvspline(Surv(pmax(Pfs, Surv), 1 - SurvInd)~1, 
                             data = data, k = 3)
  #sample from KM-curve
  U = runif(n = nrow(data))
  C = qsurvspline(p = U, gamma = km_censor$coefficients, knots = km_censor$knots)
  
  #sample from model for control
  data_temp = data[data$Treat == 0,]
  data_control = data_temp
  n_temp = nrow(data_temp)
  if(cop_type == "frank"){
    fitted_copula = bicop_dist(family = "frank", parameters = fit_0_par[9])
  }
  else if(cop_type == "gaussian"){
    rho = (exp(fit_0_par[9]) - 1)/(exp(fit_0_par[9]) + 1) 
    fitted_copula = bicop_dist(family = "gaussian", parameters = rho)
  }
  else if(cop_type == "clayton"){
    fitted_copula = bicop_dist(family = "clayton", parameters = fit_0_par[9],
                               rotation = 180)
  }
  else if(cop_type == "gumbel"){
    fitted_copula = bicop_dist(family = "gumbel", parameters = fit_0_par[9],
                               rotation = 180)
  }
  copula_structure = dvine_structure(order = 1:2)
  vine_cop = vinecop_dist(pair_copulas = list(list(fitted_copula)), structure = copula_structure)
  X = rvinecop(n_temp, vine_cop)
  data_temp$Surv = pmin(qsurvspline(p = X[,2], gamma = fit_0_par[5:8], 
                                    knots = knott0),
                        C[1:n_temp])
  data_temp$Pfs = pmin(qsurvspline(p = X[,1], gamma = fit_0_par[1:4], 
                                   knots = knots0), 
                       C[1:n_temp], data_temp$Surv)
  cens = pmax(data_temp$Pfs, data_temp$Surv) == C[1:n_temp]
  #plot sample and true data
  par(mfrow = c(2,2))
  plot(data_temp$Pfs[!cens], data_temp$Surv[!cens], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Fitted Model", xlab = "S_0", ylab = "T_0")
  points(data_temp$Pfs[cens], data_temp$Surv[cens])
  plot(data_control$Pfs[(data_control$SurvInd == 1)], 
       data_control$Surv[(data_control$SurvInd == 1)], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Observed Data", xlab = "S_0", ylab = "T_0")
  points(data_control$Pfs[!(data_control$SurvInd == 1)], data_control$Surv[!(data_control$SurvInd == 1)])
  
  #sample from model for treated
  data_temp = data[data$Treat == 1,]
  data_control = data_temp
  n_temp = nrow(data_temp)
  if(cop_type == "frank"){
    fitted_copula = bicop_dist(family = "frank", parameters = fit_1_par[9])
  }
  else if(cop_type == "gaussian"){
    rho = (exp(fit_1_par[9]) - 1)/(exp(fit_1_par[9]) + 1) 
    fitted_copula = bicop_dist(family = "gaussian", parameters = rho)
  }
  else if(cop_type == "clayton"){
    fitted_copula = bicop_dist(family = "clayton", parameters = fit_1_par[9],
                               rotation = 180)
  }
  else if(cop_type == "gumbel"){
    fitted_copula = bicop_dist(family = "gumbel", parameters = fit_1_par[9],
                               rotation = 180)
  }
  copula_structure = dvine_structure(order = 1:2)
  vine_cop = vinecop_dist(pair_copulas = list(list(fitted_copula)), structure = copula_structure)
  X = rvinecop(n_temp, vine_cop)
  data_temp$Surv = pmin(qsurvspline(p = X[,2], gamma = fit_1_par[5:8], 
                                    knots = knott1),
                        C[1:n_temp])
  data_temp$Pfs = pmin(qsurvspline(p = X[,1], gamma = fit_1_par[1:4], 
                                   knots = knots1), 
                       C[1:n_temp], data_temp$Surv)
  
  cens = data_temp$Surv == C[1:n_temp]
  #plot sample and true data
  plot(data_temp$Pfs[!cens], data_temp$Surv[!cens], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Fitted Model", xlab = "S_1", ylab = "T_1")
  points(data_temp$Pfs[cens], data_temp$Surv[cens])
  plot(data_control$Pfs[(data_control$SurvInd == 1)], 
       data_control$Surv[(data_control$SurvInd == 1)], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Observed Data", xlab = "S_1", ylab = "T_1")
  points(data_control$Pfs[!(data_control$SurvInd == 1)], data_control$Surv[!(data_control$SurvInd == 1)])
}

association_no = function(fit_0_par, fit_1_par, 
                          knots0, knots1, knott0, knott1,
                          data, cop_type){
  #estimate survival function for censoring
  km_censor = flexsurvspline(Surv(pmax(Pfs, Surv), 1 - SurvInd)~1, 
                             data = data, k = 3)
  #sample from KM-curve
  U = runif(n = nrow(data))
  C = qsurvspline(p = U, gamma = km_censor$coefficients, knots = km_censor$knots)
  
  #sample from model for control
  data_temp = data_pfs[data_pfs$Treat == 0,]
  data_control = data_temp
  n_temp = nrow(data_temp)
  if(cop_type == "frank"){
    fitted_copula = bicop_dist(family = "frank", parameters = fit_0_par[9])
  }
  else if(cop_type == "gaussian"){
    rho = (exp(fit_0_par[9]) - 1)/(exp(fit_0_par[9]) + 1) 
    fitted_copula = bicop_dist(family = "gaussian", parameters = rho)
  }
  else if(cop_type == "clayton"){
    fitted_copula = bicop_dist(family = "clayton", parameters = fit_0_par[9],
                               rotation = 180)
  }
  else if(cop_type == "gumbel"){
    fitted_copula = bicop_dist(family = "gumbel", parameters = fit_0_par[9],
                               rotation = 180)
  }
  copula_structure = dvine_structure(order = 1:2)
  vine_cop = vinecop_dist(pair_copulas = list(list(fitted_copula)), structure = copula_structure)
  X = rvinecop(n_temp, vine_cop)
  data_temp$Surv = pmin(qsurvspline(p = X[,2], gamma = fit_0_par[5:8], 
                                    knots = knott0),
                        C[1:n_temp])
  data_temp$Pfs = pmin(qsurvspline(p = X[,1], gamma = fit_0_par[1:4], 
                                   knots = knots0), 
                       C[1:n_temp])
  cens = pmax(data_temp$Pfs, data_temp$Surv) == C[1:n_temp]
  #plot sample and true data
  par(mfrow = c(2,2))
  plot(data_temp$Pfs[!cens], data_temp$Surv[!cens], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Fitted Model", xlab = "S_0", ylab = "T_0")
  points(data_temp$Pfs[cens], data_temp$Surv[cens])
  plot(data_control$Pfs[(data_control$SurvInd == 1)], 
       data_control$Surv[(data_control$SurvInd == 1)], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Observed Data", xlab = "S_0", ylab = "T_0")
  points(data_control$Pfs[!(data_control$SurvInd == 1)], data_control$Surv[!(data_control$SurvInd == 1)])
  
  #sample from model for treated
  data_temp = data_pfs[data_pfs$Treat == 1,]
  data_control = data_temp
  n_temp = nrow(data_temp)
  if(cop_type == "frank"){
    fitted_copula = bicop_dist(family = "frank", parameters = fit_1_par[9])
  }
  else if(cop_type == "gaussian"){
    rho = (exp(fit_1_par[9]) - 1)/(exp(fit_1_par[9]) + 1) 
    fitted_copula = bicop_dist(family = "gaussian", parameters = rho)
  }
  else if(cop_type == "clayton"){
    fitted_copula = bicop_dist(family = "clayton", parameters = fit_1_par[9],
                               rotation = 180)
  }
  else if(cop_type == "gumbel"){
    fitted_copula = bicop_dist(family = "gumbel", parameters = fit_1_par[9],
                               rotation = 180)
  }
  copula_structure = dvine_structure(order = 1:2)
  vine_cop = vinecop_dist(pair_copulas = list(list(fitted_copula)), structure = copula_structure)
  X = rvinecop(n_temp, vine_cop)
  data_temp$Surv = pmin(qsurvspline(p = X[,2], gamma = fit_1_par[5:8], 
                                    knots = knott1),
                        C[1:n_temp])
  data_temp$Pfs = pmin(qsurvspline(p = X[,1], gamma = fit_1_par[1:4], 
                                   knots = knots1), 
                       C[1:n_temp])
  
  cens = data_temp$Surv == C[1:n_temp]
  #plot sample and true data
  plot(data_temp$Pfs[!cens], data_temp$Surv[!cens], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Fitted Model", xlab = "S_1", ylab = "T_1")
  points(data_temp$Pfs[cens], data_temp$Surv[cens])
  plot(data_control$Pfs[(data_control$SurvInd == 1)], 
       data_control$Surv[(data_control$SurvInd == 1)], 
       col = "red", xlim = c(0,24), ylim = c(0,24), 
       main = "Observed Data", xlab = "S_1", ylab = "T_1")
  points(data_control$Pfs[!(data_control$SurvInd == 1)], data_control$Surv[!(data_control$SurvInd == 1)])
}



