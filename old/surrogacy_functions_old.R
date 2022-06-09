surrogacy_measures = function(cop_type, fit_0_par, fit_1_par, 
                              n_sim, n_prec,
                              knots0, knots1, knott0, knott1, 
                              minfo_prec = 0, restr = TRUE){
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
  r = sample(x = c(0, 90, 180, 270), size = 4*n_sim, replace = TRUE)
  r = matrix(data = r, ncol = 4)
  #put the sampled unidentifiable parameters into a matrix
  c = matrix(data = c, ncol = 4)
  for(i in 1:n_sim){
    #sample data, which depends on the copula type
    if(cop_type == "frank"){
      data_sample = frank_surrogacy_sample(fit_0_par, fit_1_par, n_prec,
                                           knots0, knots1, knott0, knott1, c[i,],
                                           restr = restr)
    }
    else if(cop_type == "gaussian"){
      data_sample = gaussian_surrogacy_sample(fit_0_par, fit_1_par, n_prec,
                                           knots0, knots1, knott0, knott1, c[i,],
                                           restr = restr)
    }
    else if(cop_type == "clayton"){
      data_sample = clayton_surrogacy_sample(fit_0_par, fit_1_par, n_prec,
                                           knots0, knots1, knott0, knott1, 
                                           c[i,], r[i,],
                                           restr = restr)
    }
    else if(cop_type == "gumbel"){
      data_sample = gumbel_surrogacy_sample(fit_0_par, fit_1_par, n_prec,
                                             knots0, knots1, knott0, knott1, 
                                             c[i,], r[i,],
                                             restr = restr)
    }
    
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
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix, method = "TLL2")
      minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
  }
  return(data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                    kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log))
}

surrogacy_sample = function(fit_0_par, fit_1_par, n_prec,
                            knots0, knots1, knott0, knott1, c_unid, r_unid,
                            cop_type, restr = TRUE){
  if (cop_type == "gaussian"){
    c12 = (exp(fit_gaussian_0_par[9]) - 1)/(exp(fit_gaussian_0_par[9]) + 1)
    c34 = (exp(fit_gaussian_1_par[9]) - 1)/(exp(fit_gaussian_1_par[9]) + 1)
  }
  else{
    c12 = fit_0_par[9]
    c34 = fit_1_par[9]
  }
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  
  pair_copulas = list(list(bicop_dist(family = cop_type, rotation = 180, parameters = c12),
                           bicop_dist(family = cop_type, rotation = r_unid[1], parameters = c23),
                           bicop_dist(family = cop_type, rotation = 180, parameters = c34)),
                      list(bicop_dist(family = cop_type, rotation = r_unid[2], parameters = c13_2),
                           bicop_dist(family = cop_type, rotation = r_unid[3], parameters = c24_3)),
                      list(bicop_dist(family = cop_type, rotation = r_unid[4], parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_0_par[1:4], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_0_par[5:8], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_1_par[1:4], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_1_par[5:8], knots = knott1)
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

surrogacy_measures_new = function(cop_type, fit_0_par, fit_1_par, 
                                n_sim, n_prec,
                                knots0, knots1, knott0, knott1, 
                                minfo_prec = 0, restr = TRUE){
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
                                   knots0, knots1, knott0, knott1, c[i,],
                                   cop_type, restr = restr)

    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
  }
  return(data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                    kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log))
}


frank_surrogacy = function(fit_frank_0_par, fit_frank_1_par, n_sim, n_prec,
                           knots0, knots1, knott0, knott1, minfo_prec = 0,
                           restr = TRUE){
  kendall = 1:n_sim
  sp_rho = 1:n_sim
  kendall_log = 1:n_sim
  sp_rho_log = 1:n_sim
  minfo = 1:n_sim
  minfo_log = 1:n_sim
  c12 = fit_frank_0_par[9]
  c34 = fit_frank_1_par[9]
  max_grid = max(c12, c34)
  c = runif(n = 4*n_sim, min = -35, max = 35) #sample unidentifiable parameters
  c = matrix(data = c, ncol = 4)
  for(i in 1:n_sim){
    data_sample = frank_surrogacy_sample(fit_frank_0_par, fit_frank_1_par, n_prec,
                                         knots0, knots1, knott0, knott1, c[i,],
                                         restr = restr)
    
    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
  }
  return(data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                    kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log))
}

frank_surrogacy_sample = function(fit_frank_0_par, fit_frank_1_par, n_prec,
                           knots0, knots1, knott0, knott1, c_unid,
                           restr = TRUE){
  c12 = fit_frank_0_par[9]
  c34 = fit_frank_1_par[9]
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  pair_copulas = list(list(bicop_dist(family = "frank", rotation = 0, parameters = c12),
                           bicop_dist(family = "frank", rotation = 0, parameters = c23),
                           bicop_dist(family = "frank", rotation = 0, parameters = c34)),
                      list(bicop_dist(family = "frank", rotation = 0, parameters = c13_2),
                           bicop_dist(family = "frank", rotation = 0, parameters = c24_3)),
                      list(bicop_dist(family = "frank", rotation = 0, parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_frank_0_par[1:4], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_frank_0_par[5:8], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_frank_1_par[1:4], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_frank_1_par[5:8], knots = knott1)
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


gaussian_surrogacy = function(fit_gaussian_0_par, fit_gaussian_1_par, n_sim, n_prec,
                             knots0, knots1, knott0, knott1, minfo_prec = 0,
                             restr = TRUE){
  kendall = 1:n_sim
  sp_rho = 1:n_sim
  kendall_log = 1:n_sim
  sp_rho_log = 1:n_sim
  minfo = 1:n_sim
  minfo_log = 1:n_sim
  c12 = (exp(fit_gaussian_0_par[9]) - 1)/(exp(fit_gaussian_0_par[9]) + 1)
  c34 = (exp(fit_gaussian_1_par[9]) - 1)/(exp(fit_gaussian_1_par[9]) + 1)
  max_grid = max(c12, c34)
  c = runif(n = 4*n_sim, min = -1, max = 1) #sample unidentifiable parameters
  c = matrix(data = c, ncol = 4)
  for(i in 1:n_sim){
    data_sample = gaussian_surrogacy_sample(fit_gaussian_0_par, fit_gaussian_1_par, n_prec,
                                         knots0, knots1, knott0, knott1, c[i,],
                                         restr = restr)
    
    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
  }
  return(data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                    kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log))
}

gaussian_surrogacy_sample = function(fit_gaussian_0_par, fit_gaussian_1_par, n_prec,
                                    knots0, knots1, knott0, knott1, c_unid,
                                    restr = TRUE){
  c12 = (exp(fit_gaussian_0_par[9]) - 1)/(exp(fit_gaussian_0_par[9]) + 1)
  c34 = (exp(fit_gaussian_1_par[9]) - 1)/(exp(fit_gaussian_1_par[9]) + 1)
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  
  pair_copulas = list(list(bicop_dist(family = "gaussian", rotation = 0, parameters = c12),
                           bicop_dist(family = "gaussian", rotation = 0, parameters = c23),
                           bicop_dist(family = "gaussian", rotation = 0, parameters = c34)),
                      list(bicop_dist(family = "gaussian", rotation = 0, parameters = c13_2),
                           bicop_dist(family = "gaussian", rotation = 0, parameters = c24_3)),
                      list(bicop_dist(family = "gaussian", rotation = 0, parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_gaussian_0_par[1:4], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_gaussian_0_par[5:8], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_gaussian_1_par[1:4], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_gaussian_1_par[5:8], knots = knott1)
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


clayton_surrogacy = function(fit_clayton_0_par, fit_clayton_1_par, n_sim, n_prec,
                            knots0, knots1, knott0, knott1, minfo_prec = 0,
                            restr = TRUE){
  kendall = 1:n_sim
  sp_rho = 1:n_sim
  kendall_log = 1:n_sim
  sp_rho_log = 1:n_sim
  minfo = 1:n_sim
  minfo_log = 1:n_sim
  c12 = fit_clayton_0_par[9]
  c34 = fit_clayton_1_par[9]
  max_grid = max(c12, c34)
  c = runif(n = 4*n_sim, min = 0, max = max_grid) #sample unidentifiable parameters
  c = matrix(data = c, ncol = 4)
  r = sample(x = c(0, 90, 180, 270), size = 4*n_sim, replace = TRUE)
  r = matrix(data = r, ncol = 4)
  for(i in 1:n_sim){
    data_sample = clayton_surrogacy_sample(fit_clayton_0_par, fit_clayton_1_par, n_prec,
                                          knots0, knots1, knott0, knott1, c[i,], r[i,],
                                          restr = restr)
    
    deltaS = data_sample$deltaS
    deltaT = data_sample$deltaT
    kendall[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo_log[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
    
    deltaS = data_sample$deltaS_log
    deltaT = data_sample$deltaT_log
    kendall_log[i] = cor.test(deltaS, deltaT, method = "kendall")$estimate
    sp_rho_log[i] = cor.test(deltaS, deltaT, method = "spearman")$estimate
    if(minfo_prec != 0){
      a = rank(deltaS)/(n_prec+1)
      b = rank(deltaT)/(n_prec+1)
      temp_matrix = matrix(data = c(a,b), ncol = 2)
      t_kde = kdecop(udata = temp_matrix)
      minfo[i] = dep_measures(object = t_kde, n_qmc = minfo_prec, measures = "minfo")
    }
  }
  return(data.frame(kendall = kendall, sp_rho = sp_rho, minfo = minfo,
                    kendall_log = kendall_log, sp_rho_log = sp_rho_log, minfo_log = minfo_log))
}

clayton_surrogacy_sample = function(fit_clayton_0_par, fit_clayton_1_par, n_prec,
                                    knots0, knots1, knott0, knott1, c_unid, r_unid,
                                    restr = TRUE){
  c12 = fit_clayton_0_par[9]
  c34 = fit_clayton_1_par[9]
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  
  pair_copulas = list(list(bicop_dist(family = "clayton", rotation = 180, parameters = c12),
                           bicop_dist(family = "clayton", rotation = r_unid[1], parameters = c23),
                           bicop_dist(family = "clayton", rotation = 180, parameters = c34)),
                      list(bicop_dist(family = "clayton", rotation = r_unid[2], parameters = c13_2),
                           bicop_dist(family = "clayton", rotation = r_unid[3], parameters = c24_3)),
                      list(bicop_dist(family = "clayton", rotation = r_unid[4], parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_clayton_0_par[1:4], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_clayton_0_par[5:8], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_clayton_1_par[1:4], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_clayton_1_par[5:8], knots = knott1)
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

gumbel_surrogacy_sample = function(fit_gumbel_0_par, fit_gumbel_1_par, n_prec,
                                    knots0, knots1, knott0, knott1, c_unid, r_unid,
                                    restr = TRUE){
  c12 = fit_gumbel_0_par[9]
  c34 = fit_gumbel_1_par[9]
  c23 = c_unid[1]
  c13_2 = c_unid[2]
  c24_3 = c_unid[3]
  c14_23 = c_unid[4]
  
  pair_copulas = list(list(bicop_dist(family = "gumbel", rotation = 180, parameters = c12),
                           bicop_dist(family = "gumbel", rotation = r_unid[1], parameters = c23),
                           bicop_dist(family = "gumbel", rotation = 180, parameters = c34)),
                      list(bicop_dist(family = "gumbel", rotation = r_unid[2], parameters = c13_2),
                           bicop_dist(family = "gumbel", rotation = r_unid[3], parameters = c24_3)),
                      list(bicop_dist(family = "gumbel", rotation = r_unid[4], parameters = c14_23)))
  copula_structure = dvine_structure(order = 1:4)
  vine_cop = vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
  
  u_vec = rvinecop(n = n_prec, vine_cop)
  s0 = qsurvspline(p = u_vec[,1], gamma = fit_gumbel_0_par[1:4], knots = knots0)
  t0 = qsurvspline(p = u_vec[,2], gamma = fit_gumbel_0_par[5:8], knots = knott0)
  s1 = qsurvspline(p = u_vec[,3], gamma = fit_gumbel_1_par[1:4], knots = knots1)
  t1 = qsurvspline(p = u_vec[,4], gamma = fit_gumbel_1_par[5:8], knots = knott1)
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

