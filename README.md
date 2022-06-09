# Master-Thesis
This repository contains the R-code used in the Master thesis. Three directories with R-files are provided: "Base Functions", "Ovarian Cancer Case Study", and "Simulations".
The content of these files is described next.

## Base Functions

This directory contains the R-code to conduct the surrogacy evaluation analysis with the methods described in the master thesis.
The methods are based on the causal-inference framework, and are developed for time-to-event surrogate and true endpoints. 
More detail on the methods can be found in the master thesis text.

### copula_SCR_fitting

This R-file contains functions to fit the proposed models. 
The following functions are helper-functions that compute the log-likelihood, given the model parameters:
- For normal identifiable copulas: normal_loglik
- For Clayton identifiable copulas: clayton_loglik
- For Frank identifiable copulas: frank_loglik
- For Gumbel identifiablek copulas: gumbel_loglik
Above functions do not need to be used directly by the end-user.

The fit_model function fits the model with following arguments:
- data: R dataframe with the following rows (S, S event indicator, T, T event indicator, treatment indicator)
- cop_type: type of identifiable copula, one of the following:
  - "gaussian"
  - "clayton"
  - "frank"
  - "gumbel"
- nknots: the number of internal knots for the Royston-Parmar model
This functions returns a list with the following elements:
- fit_0: parameter estimates for the control group
- fit_1: parameter estimates for the treated group
- log_lik: maximized log-likelihood
- knots0: knots for the survival model of S_0
- knots1: knots for the survival model of S_1
- knott0: knots for the survival model of T_0
- knott1: knots for the survival model of T_1
Note that the content of the fitted parameter vector for the respective treatment group is as follows, where k is the number of internal knots. 
- element 1 through k + 2: gamma_0 through gamma_(k+1) for S
- element k + 3 through 2k + 4: gamma_0 through gamma_(k+1) for T
- last element: copula parameter

### copula_GOF



### surrogacy_functions



