#Load files with necessary functions
source(file = "copula_SCR_fitting.R")
source(file = "copula_GOF.R")
source(file = "surrogacy_functions_new.R")
#load libraries
library(flexsurv)
library(Surrogate)
library(tidyverse)
library(survival)
library(survminer)
library(mvtnorm)
library(copula)
library(rvinecopulib)
library(kdecopula)
library(latex2exp)

#put data in correct format
data = read.csv("ovarian.csv")
data = data[,-1]
data$Pfs = data$Pfs*12
data$Surv = data$Surv*12

#number of knots
nknots = 3

#fit models with time ordering
fit_normal = fit_model(data, "gaussian", nknots)
fit_normal_0 = fit_normal$fit_0
fit_normal_1 = fit_normal$fit_1

fit_clayton = fit_model(data, "clayton", nknots)
fit_clayton_0 = fit_clayton$fit_0
fit_clayton_1 = fit_clayton$fit_1

fit_frank = fit_model(data, "frank", nknots)
fit_frank_0 = fit_frank$fit_0
fit_frank_1 = fit_frank$fit_1

fit_gumbel = fit_model(data, "gumbel", nknots)
fit_gumbel_0 = fit_gumbel$fit_0
fit_gumbel_1 = fit_gumbel$fit_1

knots0 = fit_normal$knots0
knots1 = fit_normal$knots1
knott0 = fit_normal$knott0
knott1 = fit_normal$knott1

#add results to a dataframe
model_comparison = data.frame(model = character(0),
                              loglik = numeric(0),
                              kendall0 = numeric(0),
                              kendall1 = numeric(0))
#normal model
loglik = fit_normal$log_lik
copula_par0 = (exp(fit_normal_0$par[2*(nknots + 2) + 1]) - 1)/(exp(fit_normal_0$par[2*(nknots + 2) + 1]) + 1)
copula_par1 = (exp(fit_normal_1$par[2*(nknots + 2) + 1]) - 1)/(exp(fit_normal_1$par[2*(nknots + 2) + 1]) + 1)
normal_copula = normalCopula(param = copula_par0, dim = 2, dispstr = "un")
kendall0 = tau(copula = normal_copula)
normal_copula = normalCopula(param = copula_par1, dim = 2, dispstr = "un")
kendall1 = tau(copula = normal_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "normal (ord)",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))
#clayton copula
loglik = fit_clayton$log_lik
copula_par0 = fit_clayton_0$par[2*(nknots + 2) + 1]
copula_par1 = fit_clayton_1$par[2*(nknots + 2) + 1]
clayton_copula = claytonCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = clayton_copula)
clayton_copula = claytonCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = clayton_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "Clayton (ord)",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))

#Frank Copula
loglik = fit_frank$log_lik
copula_par0 = fit_frank_0$par[2*(nknots + 2) + 1]
copula_par1 = fit_frank_1$par[2*(nknots + 2) + 1]
frank_copula = frankCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = frank_copula)
frank_copula = frankCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = frank_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "Frank (ord)",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))

#Gumbel Copula
loglik = fit_gumbel$log_lik
copula_par0 = fit_gumbel_0$par[2*(nknots + 2) + 1]
copula_par1 = fit_gumbel_1$par[2*(nknots + 2) + 1]
gumbel_copula = gumbelCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = gumbel_copula)
gumbel_copula = gumbelCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = gumbel_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "gumbel (ord)",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))

#load Ovarian data with PFS
data("Ovarian")
data_pfs = Ovarian %>% select(Pfs, Surv, Treat, PfsInd, SurvInd)
#filter observations that have impossible values, i.e. a PFS larger than OS
data_pfs = data_pfs %>% filter(Pfs <= Surv)
data_pfs = data_pfs %>% mutate(Pfs = 12*Pfs, Surv = 12*Surv)

#fit models without time ordering
fit_normal_no = fit_model(data_pfs, "gaussian", nknots)
fit_normal_0_no = fit_normal_no$fit_0
fit_normal_1_no = fit_normal_no$fit_1

fit_clayton_no = fit_model(data_pfs, "clayton", nknots)
fit_clayton_0_no = fit_clayton_no$fit_0
fit_clayton_1_no = fit_clayton_no$fit_1

fit_frank_no = fit_model(data_pfs, "frank", nknots)
fit_frank_0_no = fit_frank_no$fit_0
fit_frank_1_no = fit_frank_no$fit_1

fit_gumbel_no = fit_model(data_pfs, "gumbel", nknots)
fit_gumbel_0_no = fit_gumbel_no$fit_0
fit_gumbel_1_no = fit_gumbel_no$fit_1

knots0_no = fit_normal_no$knots0
knots1_no = fit_normal_no$knots1
knott0_no = fit_normal_no$knott0
knott1_no = fit_normal_no$knott1

#collect results in a dataframe
#normal model
loglik = fit_normal_no$log_lik
copula_par0 = (exp(fit_normal_0_no$par[2*(nknots + 2) + 1]) - 1)/(exp(fit_normal_0_no$par[2*(nknots + 2) + 1]) + 1)
copula_par1 = (exp(fit_normal_1_no$par[2*(nknots + 2) + 1]) - 1)/(exp(fit_normal_1_no$par[2*(nknots + 2) + 1]) + 1)
normal_copula = normalCopula(param = copula_par0, dim = 2, dispstr = "un")
kendall0 = tau(copula = normal_copula)
normal_copula = normalCopula(param = copula_par1, dim = 2, dispstr = "un")
kendall1 = tau(copula = normal_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "normal",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))
#clayton copula
loglik = fit_clayton_no$log_lik
copula_par0 = fit_clayton_0_no$par[2*(nknots + 2) + 1]
copula_par1 = fit_clayton_1_no$par[2*(nknots + 2) + 1]
clayton_copula = claytonCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = clayton_copula)
clayton_copula = claytonCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = clayton_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "Clayton",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))

#Frank Copula
loglik = fit_frank_no$log_lik
copula_par0 = fit_frank_0_no$par[2*(nknots + 2) + 1]
copula_par1 = fit_frank_1_no$par[2*(nknots + 2) + 1]
frank_copula = frankCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = frank_copula)
frank_copula = frankCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = frank_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "Frank",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))

#Gumbel Copula
loglik = fit_gumbel_no$log_lik
copula_par0 = fit_gumbel_0_no$par[2*(nknots + 2) + 1]
copula_par1 = fit_gumbel_1_no$par[2*(nknots + 2) + 1]
gumbel_copula = gumbelCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = gumbel_copula)
gumbel_copula = gumbelCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = gumbel_copula)
model_comparison = bind_rows(model_comparison,
                             data.frame(model = "gumbel",
                                        loglik = loglik,
                                        kendall0 = kendall0,
                                        kendall1 = kendall1))



grid = seq(0.01, 36, 0.1)
png(filename = "Figures/GOF_ovarian.png", width = 650, height = 520)
marginal_gof_scr(fit_clayton_0$par, fit_clayton_1$par, knots0, knots1, 
                 knott0, knott1,
                 data, "clayton", grid)
dev.off()
png(filename = "Figures/GOF_ovarian_no.png", width = 650, height = 520)
marginal_gof_no(fit_clayton_0_no$par, fit_clayton_1_no$par, 
                knots0_no, knots1_no,
                knott0_no, knott1_no, data_pfs, grid)
dev.off()


# SENSITIVITY ANALYSIS

sens_data_marg_i = surrogacy_measures_sens(cop_type = "clayton",
                               fit_0_par = fit_clayton_0$par, 
                               fit_1_par = fit_clayton_1$par,
                               n_sim = 5000, n_prec = 2000, 
                               knots0 = knots0, knots1 = knots1, 
                               knott0 = knott0, knott1 = knott1,
                               minfo_prec = 2000, restr = TRUE,
                               get_unid = FALSE, cop_type2 = "clayton", 
                               option = 2, ncores = 8, 
                               get_marg_tau = TRUE)
sens_data_marg = sens_data_marg_i %>%
  mutate(ordering = "Ordering")
sens_data_marg_i = surrogacy_measures_sens(cop_type = "clayton",
                                           fit_0_par = fit_clayton_0_no$par, 
                                           fit_1_par = fit_clayton_1_no$par,
                                           n_sim = 5000, n_prec = 2000, 
                                           knots0 = knots0_no, knots1 = knots1_no, 
                                           knott0 = knott0_no, knott1 = knott1_no,
                                           minfo_prec = 2000, restr = FALSE,
                                           get_unid = FALSE, cop_type2 = "clayton", 
                                           option = 2, ncores = 8, 
                                           get_marg_tau = TRUE)
sens_data_marg = sens_data_marg %>% 
  bind_rows(
    sens_data_marg_i %>%
      mutate(ordering = "No Ordering")
  )
                           


#save results of sensitivity analysis to file
save(model_comparison, sens_data_marg, file = "ovarian_sens_results.RData")
