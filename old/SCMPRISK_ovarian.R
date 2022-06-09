library(SemiCompRisks)
data = read.csv("ovarian.csv")
#FREQUENTIST
form = Formula(Pfs + PfsInd | Surv + SurvInd ~ 1 | 1 | 1)
fit = FreqID_HReg(Formula = form, data = data, model = "Markov", frailty = TRUE, 
            subset = data$Treat == 0)
pred <- predict(fit, time=seq(0,2,0.01), tseq=seq(from=0,to=2,by=0.2))
summary(fit)
plot(pred, plot.est = "Surv")
plot(pred, plot.est = "Haz")
fit = FreqID_HReg(Formula = form, data = data, model = "semi-Markov", frailty = TRUE, 
                  subset = data$Treat == 0)
summary(fit)
pred <- predict(fit, time=seq(0,2,0.01), tseq=seq(from=0,to=2,by=0.2))
summary(fit)
plot(pred, plot.est = "Surv")
plot(pred, plot.est = "Haz")

#BAYESIAN
data = data[data$Treat == 0,]
startValues <- initiate.startValues_HReg(form, data=data,
                                         model=c("semi-Markov","PEM"), nChain=3)
hyperParams <- list(theta=c(0.5,0.05), 
                    PEM=list(PEM.ab1=c(0.5,0.05),
                             PEM.ab2=c(0.5,0.05), 
                             PEM.ab3=c(0.5,0.05), 
                             PEM.alpha1=5, PEM.alpha2=5, PEM.alpha3=5))
sg_max <- c(2,2,2)
mcmcParams <- list(run=list(numReps=5e4, thin=1e2, burninPerc=0.5),
                   storage=list(nGam_save=0, storeV=rep(FALSE,3)),
                   tuning=list(mhProp_theta_var=0.05, Cg=rep(0.2,3), delPertg=rep(0.5,3),
                   rj.scheme=1, Kg_max=rep(30,3), sg_max=sg_max, time_lambda1=seq(1,sg_max[1],1),
                   time_lambda2=seq(1,sg_max[2],1), time_lambda3=seq(1,sg_max[3],1)))
fit = BayesID_HReg(form, data=data, model = c("semi-Markov", "PEM"),
                   startValues = startValues, hyperParams = hyperParams, mcmcParams = mcmcParams)


