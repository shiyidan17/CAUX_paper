library(survival)
library(SurvTrunc)

# load functions
f.files<-list.files("source_funs", full.names = TRUE)
for(ff in f.files) source(ff)

# load test data
load("testData.rda")

# Naive approach
fit.vs.cox<-coxph(Surv(test.data.CAUX$L,test.data.CAUX$X,test.data.CAUX$delta)~test.data.CAUX$Z)
summary(fit.vs.cox)

# CAUX
system.time(fit.SurroCOX<-SurroCOX(data=test.data.CAUX[1:10],
                                   phi.methods="Cox",             # model for phi
                                   bootstrap.se = TRUE,           # calculate bootstrap standard errors
                                   n.boot = 100,                  # number of bootstraps
                                   seed.boot = 12345,             # seed for bootstrap, default is NA
                                   LC.as.exact.t_star=TRUE,       # if there are too few left-censored individuals, combine them with events
                                   min.nLC.t_star=15))            # benchmark for doing so

# point estimate
fit.SurroCOX$fitcox$coefficients
# bootstrap se
fit.SurroCOX$bootstrap.se
