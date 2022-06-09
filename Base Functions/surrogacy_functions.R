surrogacy_measures_sens = function(cop_type, fit_0_par, fit_1_par, 
                              n_sim, n_prec,
                              knots0, knots1, knott0, knott1, 
                              minfo_prec = 0, restr = TRUE,
                              get_unid = FALSE, cop_type2, option = 1, 
                              ncores = 1, get_marg_tau = FALSE){
  #number of parameters
  k = length(fit_0_par)
  
  #initialize vectors to store measures of surrogacy in
  kendall <- sp_rho <- minfo <- 1:n_sim
  
  #pull association parameters from estimated parameter vectors
  #the association parameter is always the last one in the corresponding vector
  c12 = fit_0_par[k]
  c34 = fit_1_par[k]
  #"strongest" association parameter
  max_grid = max(c12, c34)
  if(option == 1){
    #sample from restricted space determined by max observed kendall's tau
    c = unid_cop_sample_restricted(pm1 = max_grid, cop_type1 = cop_type, 
                                   cop_type2 = cop_type2, n_sim = n_sim)
  }
  else if(option == 2){
    c = unid_cop_sample_unrestricted(cop_type2 = cop_type2, n_sim = n_sim)
  }
  else if(option == 3){
    c = unid_cop_sample_restricted_tau(pm1 = max_grid, cop_type1 = cop_type, 
                                       cop_type2 = cop_type2, n_sim = n_sim)
  }
  else if(option == 4){
    c = unid_cop_sample_hybrid(pm1 = max_grid, cop_type1 = cop_type, 
                               cop_type2 = cop_type2, n_sim = n_sim)
  }
  #sample rotation parameters, these are only used if the copula
  #cannot model positive and negative associations
  if(cop_type2 %in% c("gumbel", "clayton")){
    r = sample(x = c(0, 90, 180, 270), size = 4*n_sim, replace = TRUE)
  }
  else{
    r = rep(0, 4*n_sim)
  }
  #put the sampled unidentifiable parameters into a matrix and list
  r_matrix = matrix(data = r, ncol = 4)
  c_matrix = matrix(data = c, ncol = 4)
  c_list = as.list(data.frame(t(c_matrix)))
  r_list = as.list(data.frame(t(r_matrix)))
  #put all other arguments in a list for the apply function
  MoreArgs = list(fit_0_par = fit_0_par, fit_1_par = fit_1_par,
                  n_prec = n_prec,
                  knots0 = knots0, knots1 = knots1, 
                  knott0 = knott0, knott1 = knott1, 
                  cop_type = cop_type, restr = restr, 
                  cop_type2 = cop_type2, minfo_prec = minfo_prec, 
                  get_marg_tau = get_marg_tau)
  if(ncores > 1){
    cl  <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl = cl, "surrogacy_sample_sens")
    parallel::clusterEvalQ(cl = cl, expr = library(flexsurv))
    parallel::clusterEvalQ(cl = cl, expr = library(rvinecopulib))
    temp = parallel::clusterMap(cl = cl, fun = compute_mci, 
               c = c_list, r = r_list, seed = 1:n_sim, MoreArgs = MoreArgs)
    on.exit(expr = {parallel::stopCluster(cl)
      rm("cl")})
  }
  else if (ncores == 1){
    temp = mapply(FUN = compute_mci, 
                  c = c_list, r = r_list, MoreArgs = MoreArgs, 
                  SIMPLIFY = FALSE)
  }
  
  measures_df = t(as.data.frame(temp))
  if(get_marg_tau){
    colnames(measures_df) = c("kendall", "sp_rho", "minfo", 
                              "tau_s0s1", "tau_s0t0", "tau_s0t1",
                              "tau_s1t0", "tau_s1t1", 
                              "tau_t0t1")
  }
  else{
    colnames(measures_df) = c("kendall", "sp_rho", "minfo")
  }
  rownames(measures_df) = NULL
  if(!get_unid){
    return(as.data.frame(measures_df))
  }
  else{
    c = as.data.frame(x = c_matrix)
    colnames(c) = c("c23", "c13_2", "c24_3", "c14_23")
    r = as.data.frame(x = r_matrix)
    colnames(r) = c("r23", "r13_2", "r24_3", "r14_23")
    return(dplyr::bind_cols(as.data.frame(measures_df), c, r))
  }
}


compute_mci = function(fit_0_par, fit_1_par, n_prec,
                       knots0, knots1, knott0, knott1, 
                       c, r,
                       cop_type, restr, cop_type2, minfo_prec, 
                       get_marg_tau = FALSE, seed = 1){
  set.seed(seed)
  #sample data, which depends on the copula type
  data_sample = surrogacy_sample_sens(fit_0_par, fit_1_par, n_prec,
                                      knots0, knots1, knott0, knott1, 
                                      c, r,
                                      cop_type, restr = restr, cop_type2 = cop_type2,
                                      get_marg_tau = get_marg_tau)
  
  deltaS = data_sample$deltaS
  deltaT = data_sample$deltaT
  kendall = cor(deltaS, deltaT, method = "kendall")
  sp_rho = cor(deltaS, deltaT, method = "spearman")
  minfo = NA
  if(minfo_prec != 0){
    a = rank(deltaS)/(length(deltaS)+1)
    b = rank(deltaT)/(length(deltaT)+1)
    temp_matrix = matrix(data = c(a,b), ncol = 2)
    tryCatch(expr = {
      t_kde = kdecopula::kdecop(udata = temp_matrix, method = "TLL2nn")
      minfo = kdecopula::dep_measures(object = t_kde, n_qmc = minfo_prec + 1, 
                                      measures = "minfo", seed = seed + 1)
    }
    )
  }
  if(get_marg_tau){
    marginal_tau = cor(data_sample[, c("s0", "s1", "t0", "t1")], method = "kendall")
    return(c(kendall, sp_rho, minfo,
             marginal_tau[1, 2:4], marginal_tau[2, 3:4], marginal_tau[3, 4]))
  }
  else{
    return(c(kendall, sp_rho, minfo))
  }
}

surrogacy_sample_sens = function(fit_0_par, fit_1_par, n_prec,
                            knots0, knots1, knott0, knott1, c_unid, r_unid,
                            cop_type, restr = TRUE, cop_type2,
                            get_marg_tau = FALSE){
  #number of knots k
  n_k = (length(fit_0_par) - 1)/2
  if (cop_type == "gaussian"){
    c12 = (exp(fit_0_par[n_k*2 + 1]) - 1)/(exp(fit_0_par[n_k*2 + 1]) + 1)
    c34 = (exp(fit_1_par[n_k*2 + 1]) - 1)/(exp(fit_1_par[n_k*2 + 1]) + 1)
  }
  else{
    c12 = fit_0_par[n_k*2 + 1]
    c34 = fit_1_par[n_k*2 + 1]
  }
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  #survival copula rotation
  if(cop_type %in% c("clayton", "gumbel")){
    rotation = 180
  }
  else{
    rotation = 0
  }
  
  pair_copulas = list(list(bicop_dist(family = cop_type, rotation = rotation, parameters = c12),
                           bicop_dist(family = cop_type2, rotation = r_unid[1], parameters = c23),
                           bicop_dist(family = cop_type, rotation = rotation, parameters = c34)),
                      list(bicop_dist(family = cop_type2, rotation = r_unid[2], parameters = c13_2),
                           bicop_dist(family = cop_type2, rotation = r_unid[3], parameters = c24_3)),
                      list(bicop_dist(family = cop_type2, rotation = r_unid[4], parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vinecop = vine_cop, cores = 1)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_0_par[1:n_k], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_0_par[(n_k+1):(2*n_k)], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_1_par[1:n_k], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_1_par[(n_k+1):(2*n_k)], knots = knott1)
  if(restr){
    pfs0 = pmin(s0, t0)
    pfs1 = pmin(s1, t1)
  }
  else{
    pfs0 = s0
    pfs1 = s1
  }
  deltaS = pfs1 - pfs0
  deltaT = t1 - t0
  if(get_marg_tau){
    return(data.frame(deltaS = deltaS, deltaT = deltaT,
                      s0 = pfs0, s1 = pfs1, t0 = t0, t1 = t1))
  }
  else{
    return(data.frame(deltaS = deltaS, deltaT = deltaT))
  }
  
}

unid_cop_sample_restricted = function(pm1, cop_type1, cop_type2, n_sim){
  if(cop_type1 == "frank"){
    cop = frankCopula(param = pm1, dim = 2)
  }
  else if(cop_type1 == "gaussian"){
    cop = ellipCopula(family = "normal", param = pm1)
  }
  else if(cop_type1 == "clayton"){
    cop = claytonCopula(param = pm1)
  }
  else if(cop_type1 == "gumbel"){
    cop = gumbelCopula(param = pm1)
  }
  rho_max = rho(cop)
  
  #sample unidentifiable parameters
  #the ranges depend on the type of copula
  if(cop_type2 == "frank"){
    u = runif(n = 4*n_sim, min = -rho_max, max = rho_max)
    c = sapply(u, iRho, copula = frankCopula())
  }
  else if(cop_type2 == "gaussian"){
    u = runif(n = 4*n_sim, min = -rho_max, max = rho_max)
    c = sapply(u, iRho, copula = ellipCopula(family = "normal"))
  }
  else if(cop_type2 == "clayton"){
    u = runif(n = 4*n_sim, min = 0, max = rho_max)
    c = sapply(u, iRho, copula = claytonCopula())
  }
  else if(cop_type2 == "gumbel"){
    u = runif(n = 4*n_sim, min = 0, max = rho_max)
    c = sapply(u, iRho, copula = gumbelCopula())
  }
  return(c)
}

unid_cop_sample_restricted_tau = function(pm1, cop_type1, cop_type2, n_sim){
  if(cop_type1 == "frank"){
    cop = frankCopula(param = pm1, dim = 2)
  }
  else if(cop_type1 == "gaussian"){
    cop = ellipCopula(family = "normal", param = pm1)
  }
  else if(cop_type1 == "clayton"){
    cop = claytonCopula(param = pm1)
  }
  else if(cop_type1 == "gumbel"){
    cop = gumbelCopula(param = pm1)
  }
  tau_max = tau(cop)
  #sample unidentifiable parameters
  #the ranges depend on the type of copula
  #sample uniformly on tau scale
  if(cop_type2 == "frank"){
    u = runif(n = 4*n_sim, min = -tau_max, max = tau_max)
    c = sapply(u, iTau, copula = frankCopula())
  }
  else if(cop_type2 == "gaussian"){
    u = runif(n = 4*n_sim, min = -tau_max, max = tau_max)
    c = sapply(u, iTau, copula = ellipCopula(family = "normal"))
  }
  else if(cop_type2 == "clayton"){
    u = runif(n = 4*n_sim, min = 0, max = tau_max)
    c = sapply(u, iTau, copula = claytonCopula())
  }
  else if(cop_type2 == "gumbel"){
    u = runif(n = 4*n_sim, min = 0, max = tau_max)
    c = sapply(u, iTau, copula = gumbelCopula())
  }
  return(c)
}

unid_cop_sample_unrestricted = function(cop_type2, n_sim){
  #sample uniformly parameters from spearman's correlation
  if(cop_type2 == "frank"){
    u = runif(n = 4*n_sim, min = -1, max = 1)
    c = sapply(X = u, FUN = iRho, copula = frankCopula())
    #functions for frank copula cannot handle parameter larger than abs(35)
    c = ifelse(abs(c) > 35, sign(c)*35, c)
  }
  else if(cop_type2 == "gaussian"){
    u = runif(n = 4*n_sim, min = -1, max = 1)
    c = sapply(X = u, FUN = iRho, copula = ellipCopula(family = "normal"))
  }
  else if(cop_type2 == "clayton"){
    u = runif(n = 4*n_sim, min = 0, max = 1)
    c = sapply(X = u, FUN = iRho, copula = claytonCopula())
    #functions for clayton copula cannot handle parameter larger than 28
    c = ifelse(c > 28, 28, c)
  }
  else if(cop_type2 == "gumbel"){
    u = runif(n = 4*n_sim, min = 0, max = 1)
    c = sapply(X = u, FUN = iRho, copula = gumbelCopula())
    #functions for gumbel copula cannot handle parameter values larger than 50
    c = ifelse(c > 50, 50, c)
  }
  return(c)
}

unid_cop_sample_hybrid = function(pm1, cop_type1, cop_type2, n_sim){
  if(cop_type1 == "frank"){
    cop = frankCopula(param = pm1, dim = 2)
  }
  else if(cop_type1 == "gaussian"){
    cop = ellipCopula(family = "normal", param = pm1)
  }
  else if(cop_type1 == "clayton"){
    cop = claytonCopula(param = pm1)
  }
  else if(cop_type1 == "gumbel"){
    cop = gumbelCopula(param = pm1)
  }
  rho_max = rho(cop)
  #sample unidentifiable parameters
  #the ranges depend on the type of copula
  if(cop_type2 == "frank"){
    u1 = runif(n = 2*n_sim, min = -rho_max, max = rho_max)
    u2 = runif(n = 2*n_sim, min = -1, max = 1)
    u = 1:(4*n_sim)
    u[1:(4*n_sim)%%4 %in% 0:1] = u1
    u[1:(4*n_sim)%%4 %in% 2:3] = u2
    c = sapply(u, iRho, copula = frankCopula())
    c = ifelse(abs(c) > 35, sign(c)*35, c)
  }
  else if(cop_type2 == "gaussian"){
    u1 = runif(n = 2*n_sim, min = -rho_max, max = rho_max)
    u2 = runif(n = 2*n_sim, min = -1, max = 1)
    u = 1:(4*n_sim)
    u[1:(4*n_sim)%%4 %in% 0:1] = u1
    u[1:(4*n_sim)%%4 %in% 2:3] = u2
    c = sapply(u, iRho, copula = ellipCopula(family = "normal"))
  }
  else if(cop_type2 == "clayton"){
    u1 = runif(n = 2*n_sim, min = 0, max = rho_max)
    u2 = runif(n = 2*n_sim, min = 0, max = 1)
    u = 1:(4*n_sim)
    u[1:(4*n_sim)%%4 %in% 0:1] = u1
    u[1:(4*n_sim)%%4 %in% 2:3] = u2
    c = sapply(u, iRho, copula = claytonCopula())
    c = ifelse(c > 28, 28, c)
  }
  else if(cop_type2 == "gumbel"){
    u1 = runif(n = 2*n_sim, min = 0, max = rho_max)
    u2 = runif(n = 2*n_sim, min = 0, max = 1)
    u = 1:(4*n_sim)
    u[1:(4*n_sim)%%4 %in% 0:1] = u1
    u[1:(4*n_sim)%%4 %in% 2:3] = u2
    c = sapply(u, iRho, copula = gumbelCopula())
    c = ifelse(c > 50, 50, c)
  }
  return(c)
}




