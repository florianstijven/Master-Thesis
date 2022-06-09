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

This R-file contains some basic functions to check the goodness of fit.
The marginal_gof_scr and marginal_gof_no plot the model based fitted survival functions against the KM estimates for the model with and without time orderings, respectively. These functions only check the goodness of fit marginally. The association_gof_scr and association_no functions can be used in addition to check the appropriateness of the association structure, for the model with and without time orderings, respectively. These functions work by sampling the from the fitted model, and by comparing this sample with the actually observed data.

The following arguments need to be provided for these functions:
- fit_0_par and fit_1_par: vectors of parameter estimates, as returned by the fit_model function
- knots0, knots1, knott0, knott1: vector of knots for the marginal survival models, as returned by the fit_model function
- data: dataframe with observed data, should have the same structure as the data argument in fit_model
- cop_type: type of identifiable copula. This argument is not required for the marginal_gof_no function.

### surrogacy_functions

The main function in this file is the surrogacy_measures_sens function. The other functions are merely helper functions, and should not be used by the end-user.
The following arguments need to be provided to the surrogacy_measures_sens function:
- copt_type: type of the identifiable copula
- fit_0_par and fit_1_par: vectors of parameter estimates, as returned by the fit_model function
- knots0, knots1, knott0, knott1: vector of knots for the marginal survival models, as returned by the fit_model function
- n_sim: number of replications in the sensitivity analysis
- n_prec: number of MC samples to compute the measures of surrogacy
- minfo_prec: number of quasi MC samples to compute the mutual information. If this is set to zero, the mutual information is not computed. Only Kendall's tau and Spearman's rho are then computed. This can save quite a lot of computation time.
- restr: TRUE if time orderings are taken into account, else FALSE
- get_unid: TRUE of the sampled unidentifiable copula parameters should be returned as well, else FALSE
- cop_type2: type of the unidentifiable copulas
- option: determines method of sampling for the unidentifiable copula parameters, should not be changed
- ncores: default is 1, then no parallel computing is performed. If this value is larger than 1, parallel computing is used. This is recomended to speed up the computations, especially when the mutual information is computed as well.
- get_marg_tau: TRUE if the marginal kendall's tau's should be computed and returned, else FALSE
 

## Ovarian Cancer Case Study

This directory contains all code for the application of the methods to the ovarian cancer data. 

### ovarian_final_analysis

This file contains the analysis of the ovarian cancer data. Note that the results of this analysis are saved to the ovarian_sens_results.RData file. This are further "processed" in the sens_analysis_ovarian_processing file. 

### sens_analysis_ovarian_processing

The results of the sensitivity analysis are processed in this file. The results are presented in graphs and tables. The processing is separated from the sensitivity analysis itself, because the latter takes quite some computational time.

### ovarian_analysis_full

This file contains supplementary sensitivity analysis with 100k replications, instead of 5k. The results of this analysis are saved to the results_ovarian_100k.RData file and processed in the sens_analysis_ovarian_full_proccessing file. 

### sens_analysis_ovarian_full_proccessing

The results of the sensitivity analysis with 100k replications are processde in this file

## Simulations

This directory contains the code for the small simulation study regarding the accuracy of the monte carlo integration: MCI.R. In addition, the other files contain an additional simulation study in which mock-up settings are analyzed. These are not presented in the master thesis. 





