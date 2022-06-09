#Load files with necessary functions
source(file = "Code/Base Functions/copula_SCR_fitting.R")
source(file = "Code/Base Functions/copula_GOF.R")
source(file = "Code/Base Functions/surrogacy_functions_new.R")
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
data = read.csv("DATA/Ovarian Cancer/ovarian.csv")
data = data[,-1]
data$Pfs = data$Pfs*12
data$Surv = data$Surv*12

#number of knots
nknots = 3

#fit models with time ordering
fit_clayton = fit_model(data, "clayton", nknots)
fit_clayton_0 = fit_clayton$fit_0
fit_clayton_1 = fit_clayton$fit_1

knots0 = fit_clayton$knots0
knots1 = fit_clayton$knots1
knott0 = fit_clayton$knott0
knott1 = fit_clayton$knott1


#clayton copula
loglik = fit_clayton$log_lik
copula_par0 = fit_clayton_0$par[2*(nknots + 2) + 1]
copula_par1 = fit_clayton_1$par[2*(nknots + 2) + 1]
clayton_copula = claytonCopula(param = copula_par0, dim = 2)
kendall0 = tau(copula = clayton_copula)
clayton_copula = claytonCopula(param = copula_par1, dim = 2)
kendall1 = tau(copula = clayton_copula)



sens_data_marg_i = surrogacy_measures_sens(cop_type = "clayton",
                                           fit_0_par = fit_clayton_0$par, 
                                           fit_1_par = fit_clayton_1$par,
                                           n_sim = 100000, n_prec = 2000, 
                                           knots0 = knots0, knots1 = knots1, 
                                           knott0 = knott0, knott1 = knott1,
                                           minfo_prec = 0, restr = TRUE,
                                           get_unid = FALSE, cop_type2 = "clayton", 
                                           option = 2, ncores = 8, 
                                           get_marg_tau = TRUE)
sens_data_marg = sens_data_marg_i %>%
  mutate(ordering = "Ordering")

save(object = sens_data_marg, file = "results_ovarian_100k.RData")

