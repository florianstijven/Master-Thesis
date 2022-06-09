source("surrogacy_functions_new.R")
library(flexsurv)
library(survival)
library(copula)
library(rvinecopulib)
library(kdecopula)
library(dplyr)


surrogacy_sens = function(restr, option = 2){
  temp_f = surrogacy_measures_sens(cop_type = cop_type, fit_0_par = ctrl_par, fit_1_par = treat_par,
                                   n_sim = N, n_prec = n_prec,
                                   knots0 = knot, knots1 = knot, knott0 = knot, knott1 = knot, 
                                   minfo_prec = n_prec, restr = restr, get_unid = TRUE, cop_type2 = "clayton",
                                   option = option, ncores = 8, get_marg_tau = TRUE)
  collect_df = data.frame(unid = "clayton", temp_f)
  
  temp_f = surrogacy_measures_sens(cop_type = cop_type, fit_0_par = ctrl_par, fit_1_par = treat_par,
                                   n_sim = N, n_prec = n_prec,
                                   knots0 = knot, knots1 = knot, knott0 = knot, knott1 = knot,
                                   minfo_prec = n_prec, restr = restr, get_unid = TRUE, cop_type2 = "gaussian",
                                   option = option, ncores = 8, get_marg_tau = TRUE)
  collect_df = dplyr::bind_rows(collect_df,
                         data.frame(unid = "gaussian", temp_f))
  
  temp_f = surrogacy_measures_sens(cop_type = cop_type, fit_0_par = ctrl_par, fit_1_par = treat_par,
                                   n_sim = N, n_prec = n_prec,
                                   knots0 = knot, knots1 = knot, knott0 = knot, knott1 = knot, 
                                   minfo_prec = n_prec, restr = restr, get_unid = TRUE, cop_type2 = "frank",
                                   option = option, ncores = 8, get_marg_tau = TRUE)
  collect_df = dplyr::bind_rows(collect_df,
                         data.frame(unid = "frank", temp_f))
  
  temp_f = surrogacy_measures_sens(cop_type = cop_type, fit_0_par = ctrl_par, fit_1_par = treat_par,
                                   n_sim = N, n_prec = n_prec,
                                   knots0 = knot, knots1 = knot, knott0 = knot, knott1 = knot, 
                                   minfo_prec = n_prec, restr = restr, get_unid = TRUE, cop_type2 = "gumbel",
                                   option = option, ncores = 8, get_marg_tau = TRUE)
  collect_df = dplyr::bind_rows(collect_df,
                         data.frame(unid = "gumbel", temp_f))
  return(collect_df)
}

set.seed(1)
n_prec = 2000
N = 5000
cop_type = "clayton"
knot = c(1, 5)
restr = FALSE

#very strong observable association parameters
theta_c12 = 11.33
theta_c34 = 11.33
ctrl_par = c(2, 2, 5, 2, theta_c12)
treat_par = c(1.5, 2, 4, 2, theta_c34)

temp = surrogacy_sens(restr = FALSE)
data_sens_all = data.frame(temp, restr = FALSE, setting = "very strong")

#strong observable association parameters
theta_c12 = 6
theta_c34 = 6
ctrl_par = c(2, 2, 5, 2, theta_c12)
treat_par = c(1.5, 2, 4, 2, theta_c34)


temp = surrogacy_sens(restr = FALSE)
data_sens_all = dplyr::bind_rows(data_sens_all,
                                 data.frame(temp, restr = FALSE, setting = "strong"))


#moderate observable association parameters
theta_c12 = 3.715
theta_c34 = 3.715
ctrl_par = c(2, 2, 5, 2, theta_c12)
treat_par = c(1.5, 2, 4, 2, theta_c34)

temp = surrogacy_sens(restr = FALSE)
data_sens_all = dplyr::bind_rows(data_sens_all,
                                 data.frame(temp, restr = FALSE, setting = "moderate"))


set.seed(1)
#WITH RESTRICTION
#construct data set to record the amount of dependent censoring
dep_censoring_data = data.frame(setting = character(), prop = numeric(), group = numeric())
#very strong observable association parameters
theta_c12 = 11.3
theta_c34 = 11.3
ctrl_par = c(2, 2, 1.5, 1.85, theta_c12)
treat_par = c(1.5, 2, 1, 1.8, theta_c34)

clayton_copula = copula::claytonCopula(param = 11.3, dim = 2)
U = rCopula(copula = clayton_copula, n = 100000)
s0 = qsurvspline0(p = 1 - U[,1], gamma0 = 2, gamma1 = 2, knots = knot)
t0 = qsurvspline0(p = 1 - U[,2], gamma0 = 1.5, gamma1 = 1.85, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "very strong", prop = mean(t0 < s0), group = 0)
  )
U = rCopula(copula = clayton_copula, n = 100000)
s1 = qsurvspline0(p = 1 - U[,1], gamma0 = 1.5, gamma1 = 2, knots = knot)
t1 = qsurvspline0(p = 1 - U[,2], gamma0 = 1, gamma1 = 1.8, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "very strong", prop = mean(t1 < s1), group = 1)
  )

temp = surrogacy_sens(restr = FALSE)
data_sens_all = dplyr::bind_rows(data_sens_all,
                                 data.frame(temp, restr = TRUE, setting = "very strong"))


#strong observable association parameters
theta_c12 = 6
theta_c34 = 6
ctrl_par = c(2, 2, 1.5, 2, theta_c12)
treat_par = c(1.5, 2, 1, 2, theta_c34)

clayton_copula = copula::claytonCopula(param = 6, dim = 2)
U = rCopula(copula = clayton_copula, n = 100000)
s0 = qsurvspline0(p = 1 - U[,1], gamma0 = 2, gamma1 = 2, knots = knot)
t0 = qsurvspline0(p = 1 - U[,2], gamma0 = 1.5, gamma1 = 2, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "strong", prop = mean(t0 < s0), group = 0)
  )
U = rCopula(copula = clayton_copula, n = 100000)
s1 = qsurvspline0(p = 1 - U[,1], gamma0 = 1.5, gamma1 = 2, knots = knot)
t1 = qsurvspline0(p = 1 - U[,2], gamma0 = 1, gamma1 = 2, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "strong", prop = mean(t1 < s1), group = 1)
  )

temp = surrogacy_sens(restr = FALSE)
data_sens_all = dplyr::bind_rows(data_sens_all,
                                 data.frame(temp, restr = TRUE, setting = "strong"))


#moderate observable association parameters
theta_c12 = 3.715
theta_c34 = 3.715
ctrl_par = c(1.75, 1.65, 1.5, 2, theta_c12)
treat_par = c(1.25, 1.6, 1, 2, theta_c34)

clayton_copula = copula::claytonCopula(param = 3.715, dim = 2)
U = rCopula(copula = clayton_copula, n = 100000)
s0 = qsurvspline0(p = 1 - U[,1], gamma0 = 1.75, gamma1 = 1.65, knots = knot)
t0 = qsurvspline0(p = 1 - U[,2], gamma0 = 1.5, gamma1 = 2, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "moderate", prop = mean(t0 < s0), group = 0)
  )
U = rCopula(copula = clayton_copula, n = 100000)
s1 = qsurvspline0(p = 1 - U[,1], gamma0 = 1.25, gamma1 = 1.6, knots = knot)
t1 = qsurvspline0(p = 1 - U[,2], gamma0 = 1, gamma1 = 2, knots = knot)
dep_censoring_data =  dep_censoring_data %>%
  bind_rows(
    data.frame(setting = "moderate", prop = mean(t1 < s1), group = 1)
  )

temp = surrogacy_sens(restr = FALSE)
data_sens_all = dplyr::bind_rows(data_sens_all,
                                 data.frame(temp, restr = TRUE, setting = "moderate"))



save(data_sens_all, dep_censoring_data, file = "simulations_results_data.RData")

