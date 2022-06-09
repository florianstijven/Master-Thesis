surrogacy_sample = function(fit_0_par, fit_1_par, n_prec,
                            knots0, knots1, knott0, knott1, c_unid, r_unid,
                            cop_type, restr = TRUE){
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
                           bicop_dist(family = cop_type, rotation = r_unid[1], parameters = c23),
                           bicop_dist(family = cop_type, rotation = rotation, parameters = c34)),
                      list(bicop_dist(family = cop_type, rotation = r_unid[2], parameters = c13_2),
                           bicop_dist(family = cop_type, rotation = r_unid[3], parameters = c24_3)),
                      list(bicop_dist(family = cop_type, rotation = r_unid[4], parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
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
  deltaS_log = log(pfs1/pfs0)
  deltaT_log = log(t1/t0)
  deltaS = pfs1 - pfs0
  deltaT = t1 - t0
  return(data.frame(deltaS = deltaS, deltaT = deltaT,
                    deltaS_log = deltaS_log, deltaT_log = deltaT_log))
}


surrogacy_measures = function(cop_type, fit_0_par, fit_1_par, 
                                n_sim, n_prec,
                                knots0, knots1, knott0, knott1, 
                                minfo_prec = 0, restr = TRUE,
                                get_unid = FALSE){
  #number of parameters
  k = length(fit_0_par)
  
  #initialize vectors to store measures of surrogacy in
  kendall <- sp_rho <- kendall_log <- sp_rho_log <- minfo <- minfo_log <- 1:n_sim
  
  #pull association parameters from estimated parameter vectors
  #the association parameter is always the last one in the corresponding vector
  c12 = fit_0_par[k]
  c34 = fit_1_par[k]
  #"strongest" association parameter
  max_grid = max(c12, c34)
  #sample unidentifiable parameters
  #the ranges depend on the type of copula
  if(cop_type == "frank"){
    c = runif(n = 4*n_sim, min = -1*max_grid, max = max_grid)
  }
  else if(cop_type == "gaussian"){
    c = runif(n = 4*n_sim, min = -1, max = 1)
  }
  else if(cop_type == "clayton"){
    c = runif(n = 4*n_sim, min = 0, max = max_grid)
  }
  else if(cop_type == "gumbel"){
    c = runif(n = 4*n_sim, min = 1, max = max_grid)
  }
  #sample rotation parameters, these are only used if the copula
  #cannot model positive and negative associations
  if(cop_type %in% c("gumbel", "clayton")){
    r = sample(x = c(0, 90, 180, 270), size = 4*n_sim, replace = TRUE)
  }
  else{
    r = rep(0, 4*n_sim)
  }
  r = matrix(data = r, ncol = 4)
  #put the sampled unidentifiable parameters into a matrix
  c = matrix(data = c, ncol = 4)
  for(i in 1:n_sim){
    #sample data, which depends on the copula type
    data_sample = surrogacy_sample(fit_0_par, fit_1_par, n_prec,
                                   knots0, knots1, knott0, knott1, 
                                   c[i,], r[i,],
                                   cop_type, restr = restr)

    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix, method = "TLL2")
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    # if(minfo_prec != 0){
    #   a = rank(deltaS)/(n_prec+1)
    #   b = rank(deltaT)/(n_prec+1)
    #   temp_matrix = matrix(data = c(a,b), ncol = 2)
    #   t_kde = kdecop(udata = temp_matrix, method = "TLL2")
    #   minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    # }
  }
  measures_df = data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                           kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log)
  if(!get_unid){
    return(measures_df)
  }
  else{
    c = as.data.frame(x = c)
    colnames(c) = c("c23", "c13_2", "c24_3", "c14_23")
    r = as.data.frame(x = r)
    colnames(r) = c("r23", "r13_2", "r24_3", "r14_23")
    return(bind_cols(measures_df, c, r))
  }
}



surrogacy_measures_sens = function(cop_type, fit_0_par, fit_1_par, 
                              n_sim, n_prec,
                              knots0, knots1, knott0, knott1, 
                              minfo_prec = 0, restr = TRUE,
                              get_unid = FALSE, cop_type2, option = 1){
  #number of parameters
  k = length(fit_0_par)
  
  #initialize vectors to store measures of surrogacy in
  kendall <- sp_rho <- kendall_log <- sp_rho_log <- minfo <- minfo_log <- 1:n_sim
  
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
  #sample rotation parameters, these are only used if the copula
  #cannot model positive and negative associations
  if(cop_type2 %in% c("gumbel", "clayton")){
    r = sample(x = c(0, 90, 180, 270), size = 4*n_sim, replace = TRUE)
  }
  else{
    r = rep(0, 4*n_sim)
  }
  r = matrix(data = r, ncol = 4)
  #put the sampled unidentifiable parameters into a matrix
  c = matrix(data = c, ncol = 4)
  for(i in 1:n_sim){
    #sample data, which depends on the copula type
    data_sample = surrogacy_sample_sens(fit_0_par, fit_1_par, n_prec,
                                   knots0, knots1, knott0, knott1, 
                                   c[i,], r[i,],
                                   cop_type, restr = restr, cop_type2 = cop_type2)
    
    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix, method = "TLL2")
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    # if(minfo_prec != 0){
    #   a = rank(deltaS)/(n_prec+1)
    #   b = rank(deltaT)/(n_prec+1)
    #   temp_matrix = matrix(data = c(a,b), ncol = 2)
    #   t_kde = kdecop(udata = temp_matrix, method = "TLL2")
    #   minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    # }
  }
  measures_df = data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                           kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log)
  if(!get_unid){
    return(measures_df)
  }
  else{
    c = as.data.frame(x = c)
    colnames(c) = c("c23", "c13_2", "c24_3", "c14_23")
    r = as.data.frame(x = r)
    colnames(r) = c("r23", "r13_2", "r24_3", "r14_23")
    return(bind_cols(measures_df, c, r))
  }
}

surrogacy_sample_sens = function(fit_0_par, fit_1_par, n_prec,
                            knots0, knots1, knott0, knott1, c_unid, r_unid,
                            cop_type, restr = TRUE, cop_type2){
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
  
  u_vec = rvinecop(n = n_prec, vine_cop)
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
  deltaS_log = log(pfs1/pfs0)
  deltaT_log = log(t1/t0)
  deltaS = pfs1 - pfs0
  deltaT = t1 - t0
  return(data.frame(deltaS = deltaS, deltaT = deltaT,
                    deltaS_log = deltaS_log, deltaT_log = deltaT_log))
}

unid_cop_sample_restricted = function(pm1, cop_type1, cop_type2, n_sim){
  if(cop_type1 == "frank"){
    cop = frankCopula(param = pm1, dim = 2)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "gaussian"){
    cop = ellipCopula(family = "normal", param = pm1)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "clayton"){
    cop = claytonCopula(param = pm1)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "gumbel"){
    cop = gumbelCopula(param = pm1)
    tau_max = tau(cop)
  }
  #sample unidentifiable parameters
  #the ranges depend on the type of copula
  if(cop_type2 == "frank"){
    tau_max_inv = iTau(copula = frankCopula(), tau = tau_max)
    c = runif(n = 4*n_sim, min = -tau_max_inv, max = tau_max_inv)
  }
  else if(cop_type2 == "gaussian"){
    tau_max_inv = iTau(copula = ellipCopula(family = "normal"), tau = tau_max)
    c = runif(n = 4*n_sim, min = -tau_max_inv, max = tau_max_inv)
  }
  else if(cop_type2 == "clayton"){
    tau_max_inv = iTau(copula = claytonCopula(), tau = tau_max)
    c = runif(n = 4*n_sim, min = 0, max = tau_max_inv)
  }
  else if(cop_type2 == "gumbel"){
    tau_max_inv = iTau(copula = gumbelCopula(), tau = tau_max)
    c = runif(n = 4*n_sim, min = 1, max = tau_max_inv)
  }
  return(c)
}

unid_cop_sample_restricted_tau = function(pm1, cop_type1, cop_type2, n_sim){
  if(cop_type1 == "frank"){
    cop = frankCopula(param = pm1, dim = 2)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "gaussian"){
    cop = ellipCopula(family = "normal", param = pm1)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "clayton"){
    cop = claytonCopula(param = pm1)
    tau_max = tau(cop)
  }
  else if(cop_type1 == "gumbel"){
    cop = gumbelCopula(param = pm1)
    tau_max = tau(cop)
  }
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
  u = runif(n = 4*n_sim, min = 0, max = 1)
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

