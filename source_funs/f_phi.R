calculate.phi.binary.cov<-function(data.obj,
                                   phi.methods.ds1=c("lm", "splines", "Cox"),
                                   ind.X_star=FALSE,
                                   bs.df=3, no.t_star.ind=FALSE,
                                   lt.nvs=FALSE,
                                   weight.type=0){

  list2env(data.obj, envir = environment())
  events.only.nvs<-min(delta_star.nvs)

  fit.ds0<-NULL
  fit.ds2<-NULL

  if(!lt.nvs) phi.methods.ds1<-"Cox"

  W<-rep(1, length(X))
  if(weight.type){
    data.missing<-data.frame(VS=rep(c(1,0), c(length(X), length(X_star.nvs))), Z=c(Z, Z.nvs))
    fit.missing<-glm(VS~Z, data = data.missing, family = binomial(link = "logit"))
    W<-predict(fit.missing, newdata = data.missing[1:length(X), ], type="response")
    if(weight.type==1) W<-1/W
    else W<-(1-W)/W
  }

  #--------------------------------------------------
  # For delta_star=1
  #--------------------------------------------------
  if(ind.X_star){
    data.fit.ds1<-data.frame(X=X,
                             X_star=X_star,
                             delta=delta,
                             Z=Z, L=L, W=W)[delta_star==1, ]
    fit.ds1<-coxph(Surv(L, X, delta)~Z, data = data.fit.ds1, weights = W)
    # summary(fit.ds1)

    haz.ds1<-basehaz(fit.ds1, centered = FALSE)

    # Z=0
    haz.ds1$surv<-exp(-haz.ds1$hazard)
    cdf.ds1<-stepfun(x=c(haz.ds1$time,Inf), y=1-c(1, haz.ds1$surv, 0))
    cdf.ds1<-cdf.ds1(t_dv.w)
    phi.ds1.Z0<-cdf.ds1-c(0,cdf.ds1[-length(cdf.ds1)])
    # Z=1
    haz.ds1$surv<-exp(-haz.ds1$hazard*exp(fit.ds1$coefficients))
    cdf.ds1<-stepfun(x=c(haz.ds1$time,Inf), y=1-c(1, haz.ds1$surv, 0))
    cdf.ds1<-cdf.ds1(t_dv.w)
    phi.ds1.Z1<-cdf.ds1-c(0,cdf.ds1[-length(cdf.ds1)])

    pred.id<-delta_star.nvs==1

    if(lt.nvs) phi<-sapply(1:sum(pred.id),
                           f.phi.adj.for.lt.rc,
                           phi.Z0=phi.ds1.Z0, phi.Z1=phi.ds1.Z1,
                           l=L.nvs[pred.id], x.values=t_dv.w,
                           z=Z.nvs[pred.id])
    else phi<-sapply(1:sum(pred.id),
                     f.phi.adj.for.lt.rc,
                     phi.Z0=phi.ds1.Z0, phi.Z1=phi.ds1.Z1,
                     l=rep(0, sum(pred.id)), x.values=t_dv.w,
                     z=Z.nvs[pred.id])
  }
  else{
    if(phi.methods.ds1[1]=="Cox"){
      # data for model and prediction
      data.fit.ds1<-data.frame(X=X,
                               X_star=X_star,
                               delta=delta,
                               Z=Z, L=L, W=W)[delta_star==1, ]
      if(lt.nvs) data.pred.ds1<-data.frame(X_star=X_star.nvs, L=L.nvs, Z=Z.nvs)[delta_star.nvs==1, ]
      else data.pred.ds1<-data.frame(X_star=X_star.nvs, L=rep(0, length(Z.nvs)), Z=Z.nvs)[delta_star.nvs==1, ]
      # fit model for phi
      fit.ds1<-coxph(Surv(L, X, delta)~Z*X_star, data = data.fit.ds1, weights = W)

      if(no.t_star.ind){
        if(summary(fit.ds1)$coef["X_star", 5]>0.05) fit.ds1<-coxph(Surv(L, X, delta)~Z, data = data.fit.ds1)
      }
      # summary(fit.ds1)
      haz.ds1<-basehaz(fit.ds1, centered = FALSE)
      pred.ds1<-predict(fit.ds1, newdata=data.pred.ds1, type="risk", reference="zero")

      #row - t_dv.w, col - subject
      phi<-sapply(1:nrow(data.pred.ds1), function(i, bh, pred.risk, lt, t_dv.w){
        pred.surv<-exp(-bh$hazard*pred.risk[i])
        cdf.ds1<-stepfun(x=c(bh$time,Inf), y=1-c(1, pred.surv, 0))
        cdf.ds1<-cdf.ds1(t_dv.w)
        phi.ds1.Z1<-cdf.ds1-c(0,cdf.ds1[-length(cdf.ds1)])
        phi.ds1.Z1[t_dv.w<lt[i]]<-0
        return(phi.ds1.Z1/sum(phi.ds1.Z1))
      }, bh=haz.ds1, pred.risk=pred.ds1, lt=data.pred.ds1$L, t_dv.w=t_dv.w)
    }
    else{
      # data for model and prediction
      data.fit.ds1<-data.frame(log.X.fit=log(X),
                               log.X_star.fit=log(X_star),
                               delta=delta,
                               Z=Z, L=L, log.L=log(L), W=W)[delta_star==1, ]
      if(lt.nvs) data.pred.ds1<-data.frame(log.X_star.fit=log(X_star.nvs), L=L.nvs, Z=Z.nvs, log.L=log(L.nvs))[delta_star.nvs==1, ]
      else data.pred.ds1<-data.frame(log.X_star.fit=log(X_star.nvs), L=rep(0,length(Z.nvs)), Z=Z.nvs, log.L=rep(-Inf,length(Z.nvs)))[delta_star.nvs==1, ]

      # fit model for phi
      if(phi.methods.ds1[1]=="lm") fit.ds1<-survreg(Surv(log.X.fit, delta)~log.X_star.fit, dist = "gaussian", data = data.fit.ds1, weights = W)
      if(phi.methods.ds1[1]=="splines"){ # general case
        if(any(L!=0)) fit.ds1<-survreg(Surv(log.X.fit, delta)~bs(log.X_star.fit, df=bs.df)*bs(log.L, df=bs.df)*Z, dist = "gaussian", data = data.fit.ds1, weights = W)
        else fit.ds1<-survreg(Surv(log.X.fit, delta)~bs(log.X_star.fit, df=bs.df)*Z, dist = "gaussian", data = data.fit.ds1, weights = W)
      }
      if(phi.methods.ds1[1]=="poly") fit.ds1<-survreg(Surv(log.X.fit, delta)~log.X_star.fit+log.X_star.fit^2, dist = "gaussian", data = data.fit.ds1, weights = W)

      # predict
      pred.phi<-predict(fit.ds1, newdata = data.pred.ds1, se.fit = TRUE)
      phi<-sapply(log(t_dv.w), pnorm,
                  mean=pred.phi$fit, sd=(pred.phi$se.fit^2+fit.ds1$scale^2)^(1/2)) # row - subject, col - t_dv
      phi<-phi-cbind(rep(0,sum(delta_star.nvs==1)), phi[, -J.w])
      #phi<-phi/apply(phi, sum, MARGIN=1)

      # CHECK: apply(phi, sum, MARGIN=1) #check
      # adjust for left-truncation: row - t_dv.w, col - subject
      phi<-sapply(1:nrow(data.pred.ds1),
                  f.phi.adj.for.lt,
                  phi=phi, l=data.pred.ds1$L, x.values=t_dv.w)
    }
  }


  # #CHECK:
  # apply(phi, sum, MARGIN=2)
  # dim(phi)
  # par(mfrow=c(1,2))
  # hist(phi)
  # hist(phi)
  # par(mfrow=c(1,1))


  #--------------------------------------------------
  # For delta_star=0
  #--------------------------------------------------
  if(any(delta_star.nvs==0)){
    if(all(delta_star==1)) stop("Validation set does not contain censored surrogate outcome!")

    fit.id<-delta_star==0
    fit.ds0<-coxph(Surv(L[fit.id], X[fit.id], delta[fit.id])~Z[fit.id], weights=W[fit.id])
    # summary(fit.ds0)

    haz.ds0<-basehaz(fit.ds0, centered = FALSE)

    # Z=0
    haz.ds0$surv<-exp(-haz.ds0$hazard)
    cdf.ds0<-stepfun(x=c(haz.ds0$time,Inf), y=1-c(1, haz.ds0$surv, 0))
    cdf.ds0<-cdf.ds0(t_dv.w)
    phi.rc.Z0<-cdf.ds0-c(0,cdf.ds0[-length(cdf.ds0)])
    # Z=1
    haz.ds0$surv<-exp(-haz.ds0$hazard*exp(fit.ds0$coefficients))
    cdf.ds0<-stepfun(x=c(haz.ds0$time,Inf), y=1-c(1, haz.ds0$surv, 0))
    cdf.ds0<-cdf.ds0(t_dv.w)
    phi.rc.Z1<-cdf.ds0-c(0,cdf.ds0[-length(cdf.ds0)])

    pred.id<-delta_star.nvs==0

    if(lt.nvs) phi.rc.adj<-sapply(1:sum(pred.id),
                                  f.phi.adj.for.lt.rc,
                                  phi.Z0=phi.rc.Z0, phi.Z1=phi.rc.Z1,
                                  l=L.nvs[pred.id], x.values=t_dv.w,
                                  z=Z.nvs[pred.id])
    else phi.rc.adj<-sapply(1:sum(pred.id),
                            f.phi.adj.for.lt.rc,
                            phi.Z0=phi.rc.Z0, phi.Z1=phi.rc.Z1,
                            l=rep(0, sum(pred.id)), x.values=t_dv.w,
                            z=Z.nvs[pred.id])
    #phi.rc.adj<-t(phi.rc.adj)
    # phi<-t(cbind(phi, phi.rc.adj))
  }
  #apply(phi.rc.adj, sum, MARGIN=2) #check
  #dim(phi.rc.adj)

  else phi.rc.adj<-NULL
  # phi: row - subject; column - t_dv.w

  #--------------------------------------------------
  # For delta_star=2
  #--------------------------------------------------
  if(any(delta_star.nvs==2)){
    if(all(delta_star==1)) stop("Validation set does not contain censored surrogate outcome!")

    fit.id<-delta_star==2
    fit.ds2<-coxph(Surv(L[fit.id], X[fit.id], delta[fit.id])~Z[fit.id], weights=W[fit.id])
    # summary(fit.ds0)

    haz.ds2<-basehaz(fit.ds2, centered = FALSE)

    # Z=0
    haz.ds2$surv<-exp(-haz.ds2$hazard)
    cdf.ds2<-stepfun(x=c(haz.ds2$time,Inf), y=1-c(1, haz.ds2$surv, 0))
    cdf.ds2<-cdf.ds2(t_dv.w)
    phi.lc.Z0<-cdf.ds2-c(0,cdf.ds2[-length(cdf.ds2)])
    # Z=1
    haz.ds2$surv<-exp(-haz.ds2$hazard*exp(fit.ds2$coefficients))
    cdf.ds2<-stepfun(x=c(haz.ds2$time,Inf), y=1-c(1, haz.ds2$surv, 0))
    cdf.ds2<-cdf.ds2(t_dv.w)
    phi.lc.Z1<-cdf.ds2-c(0,cdf.ds2[-length(cdf.ds2)])

    pred.id<-delta_star.nvs==2

    if(lt.nvs) phi.lc.adj<-sapply(1:sum(pred.id),
                                  f.phi.adj.for.lt.rc,
                                  phi.Z0=phi.lc.Z0, phi.Z1=phi.lc.Z1,
                                  l=L.nvs[pred.id], x.values=t_dv.w,
                                  z=Z.nvs[pred.id])
    else phi.lc.adj<-sapply(1:sum(pred.id),
                            f.phi.adj.for.lt.rc,
                            phi.Z0=phi.lc.Z0, phi.Z1=phi.lc.Z1,
                            l=rep(0, sum(pred.id)), x.values=t_dv.w,
                            z=Z.nvs[pred.id])
    #phi.rc.adj<-t(phi.rc.adj)
    # phi<-t(cbind(phi, phi.rc.adj))
  }
  #apply(phi.rc.adj, sum, MARGIN=2) #check
  #dim(phi.rc.adj)

  else phi.lc.adj<-NULL


  phi<-t(cbind(phi, phi.rc.adj, phi.lc.adj)) # phi: row - subject; column - t_dv.w

  Z.nvs.ds1<-Z.nvs[delta_star.nvs==1]
  Z.nvs.ds0<-Z.nvs[delta_star.nvs==0]
  Z.nvs.ds2<-Z.nvs[delta_star.nvs==2]
  which.Z0<-c(Z.nvs.ds1==0, Z.nvs.ds0==0, Z.nvs.ds2==0)
  which.Z1<-!which.Z0

  phi.Z0<-apply(phi[which.Z0, ], sum, MARGIN=2)
  phi.Z1<-apply(phi[which.Z1, ], sum, MARGIN=2)

  return(list(phi=cbind(phi.Z0, phi.Z1),
              fit.surro.res.events=fit.ds1, fit.surro.res.rc=fit.ds0,
              fit.surro.res.lc=fit.ds2))

}

f.smooth.ecdf<-function(ef){
  jump.value<-unique(ef)[2]-unique(ef)[1]
  which.jump<-c(1:length(ef))[!duplicated(ef)]

  value<-0
  for(i in 1:length(which.jump[-1])){
    int.start<-which.jump[i]
    int.end<-which.jump[i+1]
    value<-c(value, ef[which.jump[i]]+seq(0,jump.value,length.out=int.end-int.start+1)[-1])
  }
  if(int.end<length(ef)) value<-c(value, rep(1,length(ef)-int.end))
  return(value)
}

f.pdf.to.cdf<-function(pdf.value){
  value<-sapply(1:length(pdf.value), function(i,v) sum(v[1:i]),
                v=pdf.value)
  return(value)
}


f.phi.adj.for.lt<-function(r, phi, l, x.values){
  which.obs<-which.max(x.values>l[r])
  if(which.obs==1) return(phi[r, ])
  s.v<-sum(phi[r, -(1:(which.obs-1))])
  return(c(rep(0, (which.obs-1)), phi[r, -(1:(which.obs-1))]/s.v))
}

f.phi.adj.for.lt.rc<-function(r, l, x.values, phi.Z0, phi.Z1, z){

  if(z[r]==1) phi<-phi.Z1
  else phi<-phi.Z0
  #print(r)
  which.obs<-which.max(x.values>=l[r])
  if(which.obs==1) return(phi)
  else{
    # surv.l<-sum(phi[1:(which.obs-1)])+(l[r]-x.values[which.obs-1])/(x.values[which.obs]-x.values[which.obs-1])*phi[which.obs]
    # surv.l<-1-surv.l
    surv.l<-sum(phi[-(1:(which.obs-1))])
    return(c(rep(0, (which.obs-1)), phi[-(1:(which.obs-1))]/surv.l))
  }
}

f.pred.surv<-function(r, fitobj, pred.data, pred.points){
  cdf.t<-predict(fitobj, newdata=pred.data[r,], type='quantile',p=1:9998/10000)
  fun.cdf<-stepfun(x=exp(c(cdf.t,Inf)), y=c(0,1:9998/10000,1))
  return(fun.cdf(pred.points))
}

