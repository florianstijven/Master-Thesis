#Normal loglikelihood with Royston-Parmar marginal distributions
#This function is an adaptation of the R-functions provided in Sorell_2022
#This function was orginally written for exponential marginals, but the extension 
#to other distributions is relatively straightforward
normal_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is rho
  rho <- (exp(para[2*(k + 2) + 1]) - 1)/(exp(para[2*(k + 2) + 1]) + 1)
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  df = data.frame(X, Y)
  
  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]
  
  part1 <- ifelse((sum(df.1)>0),sum(-0.5*log(1-rho^2)+(((2*rho*qnorm(psurvspline(q = X.1, gamma = gammax, knots = knotsx))*qnorm(psurvspline(q = Y.1, gamma = gammay, knots = knotsy))-
                                                           rho^2*(qnorm(psurvspline(q = X.1, gamma = gammax, knots = knotsx))^2 + qnorm(psurvspline(q = Y.1, gamma = gammay, knots = knotsy))^2)))/((2*(1-rho^2))))
                                      + log(dsurvspline(x = X.1, gamma = gammax, knots = knotsx)) + log(dsurvspline(x = Y.1, gamma = gammay, knots = knotsy))),0)
  
  part2 <- ifelse((sum(df.2)>0),sum(log(pnorm(qnorm(psurvspline(q = Y.2, gamma = gammay, knots = knotsy, timescale = "log")), mean=rho*qnorm(psurvspline(q = X.2, gamma = gammax, knots = knotsx, timescale = "log")),
                                              sd=sqrt(1-rho^2), lower.tail=F)) + log(dsurvspline(x = X.2, gamma = gammax, knots = knotsx, timescale = "log"))),0)
  
  part3 <- ifelse((sum(df.3)>0),sum(log(pnorm(qnorm(psurvspline(q = X.3, gamma = gammax, knots = knotsx, timescale = "log")), mean=rho*qnorm(psurvspline(q = Y.3, gamma = gammay, knots = knotsy, timescale = "log")),
                                              sd=sqrt(1-rho^2), lower.tail=F)) + log(dsurvspline(x = Y.3, gamma = gammay, knots = knotsy, timescale = "log"))),0)
  
  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }
  
  part4 <- ifelse((sum(df.4)>0),sum(log(apply(qnorm(cbind(psurvspline(q = X.4, gamma = gammax, knots = knotsx, timescale = "log"), psurvspline(q = Y.4, gamma = gammay, knots = knotsy, timescale = "log"))),1,normal_cdf,cov_matrix))),0)
  
  loglik <- (part1+part2+part3+part4) 
  return(loglik)
}


clayton_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]
  
  #survival probabilities
  u = psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = dsurvspline(x = Y, gamma = gammay, knots = knotsy)
  
  C <- (u^(-theta) + v^(-theta) - 1)^(-1/theta) #copula
  
  part1 <- ifelse(d1*d2==1,log(1+theta)+(1+2*theta)*log(C) - (theta+1)*log(u)-(theta+1)*log(v) + log(du) + log(dv),0) #both events
  part2 <- ifelse(d1*(1-d2)==1, (theta+1)*log(C) - (theta+1)*log(u) + log(du),0) #non-terminal event only
  part3 <- ifelse((1-d1)*d2==1, (theta+1)*log(C) - (theta+1)*log(v) + log(dv),0) #terminal event only
  part4 <- ifelse((1-d1)*(1-d2)==1,log(C),0) #both events censored
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

frank_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]
  
  #survival probabilities
  u = psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = dsurvspline(x = Y, gamma = gammay, knots = knotsy)
  
  #copula
  C <- (-1/theta)*log(((1-exp(-theta)-(1-exp(-theta*u))*(1-exp(-theta*v))))/(1-exp(-theta)))
  
  part1 <- ifelse(d1*d2==1,(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*u)-1)-log(exp(theta*v)-1)+log(du)+log(dv)),0)
  part2 <- ifelse(d1*(1-d2)==1,log((1-exp(theta*C))/(1-exp(theta*u)))+log(du),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log((1-exp(theta*C))/(1-exp(theta*v)))+log(dv)),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

gumbel_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]
  
  #survival probabilities
  u = psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = dsurvspline(x = Y, gamma = gammay, knots = knotsy)
  
  #copula
  C <- exp(-((-log(u))^(theta)+(-log(v))^(theta))^(1/theta))
  
  part1 <- ifelse(d1*d2==1,log(C)+(theta-1)*log(-log(u))+(theta-1)*log(-log(v))+
                              log(theta-1-log(C))-log(u)-log(v)-(2*theta-1)*log(-log(C))+
                              log(dv)+log(du),0)
  part2 <- ifelse(d1*(1-d2)==1,log(C)+(theta-1)*log(-log(u))-log(u)-(theta-1)*log(-log(C))+log(du),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log(C)+(theta-1)*log(-log(v))-log(v)-(theta-1)*log(-log(C))+log(dv)),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
  
  loglik <- sum(part1+part2+part3+part4) 
  
  return(loglik)
}

#data should be in the prespecified format
#(S, S event indicator, T, T event indicator, treatment indicator)
fit_model = function(data, cop_type, nknots = 2){
  maxit = 500
  #colnames are added to make the intrepretation of the furthe code easier
  #Pfs refers to the surrogate, Surv refers to the true endpoint 
  colnames(data) = c("Pfs", "Surv", "Treat", "PfsInd", "SurvInd")
  #choose correct log-likelihood function
  #starting value for the association parameter is obtained by 
  #estimating the copula parameter through kendall's tau, ignoring censoring
  tau_0 = cor(data$Pfs[data$Treat == 0], data$Surv[data$Treat == 0], 
              method = "kendall")
  tau_1 = cor(data$Pfs[data$Treat == 1], data$Surv[data$Treat == 1],
              method = "kendall")
  if(cop_type == "gaussian"){
    log_lik_fn = normal_loglik
    inv_tau_0 = iTau(copula = ellipCopula(family = "normal"),
                     tau = tau_0)
    inv_tau_0 = log(1 + inv_tau_0) - log(1 - inv_tau_0)
    inv_tau_1 = iTau(copula = ellipCopula(family = "normal"),
                     tau = tau_1)
    inv_tau_1 = log(1 + inv_tau_1) - log(1 - inv_tau_1)
  }
  else if(cop_type == "clayton"){
    log_lik_fn = clayton_loglik
    inv_tau_0 = iTau(copula = claytonCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = claytonCopula(),
                     tau = tau_1)
  }
  else if(cop_type == "frank"){
    log_lik_fn = frank_loglik
    inv_tau_0 = iTau(copula = frankCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = frankCopula(),
                     tau = tau_1)
  }
  else if(cop_type == "gumbel"){
    log_lik_fn = gumbel_loglik
    inv_tau_0 = iTau(copula = gumbelCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = gumbelCopula(),
                     tau = tau_1)
  }
  #fit univariate models to provide starting values
  fit_s0 = flexsurvspline(formula = Surv(Pfs, PfsInd)~1, data = data, 
                          subset = data$Treat == 0, k = nknots, scale = "hazard")
  fit_t0 = flexsurvspline(formula = Surv(Surv, SurvInd)~1, data = data, 
                          subset = data$Treat == 0, k = nknots, scale = "hazard")
  fit_s1 = flexsurvspline(formula = Surv(Pfs, PfsInd)~1, data = data, 
                          subset = data$Treat == 1, k = nknots, scale = "hazard")
  fit_t1 = flexsurvspline(formula = Surv(Surv, SurvInd)~1, data = data, 
                          subset = data$Treat == 1, k = nknots, scale = "hazard")
  
  inits_0 = c(fit_s0$coefficients, fit_t0$coefficients, inv_tau_0)
  inits_1 = c(fit_s1$coefficients, fit_t1$coefficients, inv_tau_1)
  
  fit_0 = optim(par = inits_0, fn = log_lik_fn, method = "BFGS",
                X = data$Pfs[data$Treat == 0], Y = data$Surv[data$Treat == 0], 
                d1 = data$PfsInd[data$Treat == 0], d2 = data$SurvInd[data$Treat == 0],
                k = nknots, knotsx = fit_s0$knots, knotsy = fit_t0$knots, 
                control = list(maxit = maxit, fnscale = -1, reltol = 1e-8, 
                               ndeps = rep(1e-5, 2*(nknots + 2) + 1)))
  fit_1 = optim(par = inits_1, fn = log_lik_fn, method = "BFGS",
                X = data$Pfs[data$Treat == 1], Y = data$Surv[data$Treat == 1], 
                d1 = data$PfsInd[data$Treat == 1], d2 = data$SurvInd[data$Treat == 1], 
                k = nknots, knotsx = fit_s1$knots, knotsy = fit_t1$knots, 
                control = list(maxit = maxit, fnscale = -1, reltol = 1e-8, 
                               ndeps = rep(1e-5, 2*(nknots + 2) + 1)))
  
  return(list(fit_0 = fit_0, fit_1 = fit_1, 
              log_lik = fit_0$value + fit_1$value,
              knots0 = fit_s0$knots, knots1 = fit_s1$knots,
              knott0 = fit_t0$knots, knott1 = fit_t1$knots))
}
