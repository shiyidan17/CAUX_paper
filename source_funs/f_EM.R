# Fit COX model with surrogate outcome and nonvalidation set

# X:
# delta:
# L:
# Z:
# X_star:
# delta_star:
# bs.df: degree-of-freedom for splines model, default is 3
# tol: tolerance of convergence, default is 1e-7
# lt.value:
# phi.lm:
# check: if TRUE, iteration information will be printed when fitting model
# haz.method:
# value.at.failure.times:
# adj.pred: TBDeleted
# left.cont: TBDeleted

f.E1.rc<-function(rc.value, pcox.m, J.w=nrow(pcox.m)){
  #rc.value: col-1: rc.pos, col-2: k
  pcox.m<-pcox.m[, rc.value[2]]
  if(rc.value[1]>1){
    pcox.m[1:(rc.value[1]-1)]<-0
    pcox.m<-pcox.m/sum(pcox.m)
  }
  return(pcox.m)
}

SurroCOX<-function(data=NULL,
                   X=NULL, Z=NULL, delta=rep(1,length(X)), L=rep(0,length(X)), X_star=NULL, delta_star=rep(1,length(X_star)),
                   L.nvs=rep(0,length(X_star.nvs)), Z.nvs=NULL, X_star.nvs=NULL, delta_star.nvs=rep(1,length(X_star.nvs)),
                   tol=1e-5,
                   dist.lt=c("unknown", "uniform"),
                   phi.methods=c("Cox","splines","lm","poly"),
                   haz.method=c("EM", "coxph"),
                   left.cont=FALSE,
                   check=FALSE,
                   ind.X_star=FALSE,
                   no.t_star.ind=FALSE,
                   lt.nvs=FALSE,
                   bootstrap.se=FALSE,
                   n.boot=100,
                   seed.boot=NA,
                   bootout=FALSE,
                   CI.level=0.95, weight.type=0,
                   perturb.se=FALSE,
                   LC.as.exact.t_star=FALSE,
                   min.nLC.t_star=15){


  #--------------------------------------------------
  # data objects
  #--------------------------------------------------
  if(length(data)>0){
    list2env(data, envir = environment())
    rm(X.nvs, delta.nvs,
       X.all, delta.all,
       L.all, Z.all,
       X_star.all, delta_star.all,
       T.vs, T.nvs, T.all,
       T_star.vs, T_star.nvs, T_star.all,
       data)
  }
  if(LC.as.exact.t_star)
    if(sum(delta_star==2)<min.nLC.t_star){
      delta_star[delta_star==2]<-1
      delta_star.nvs[delta_star.nvs==2]<-1
    }
  #--------------------------------------------------
  # warnings and errors
  # scenarios
  #--------------------------------------------------
  if(max(delta)==0) stop("At least one failure time needed!")
  left.truncation<-1-(max(L)==0) # T - no left truncation
  surrogate.response<-(length(X_star.nvs)>0) # nvs
  right.censoring<-1-min(delta) # T - at least one right censoring, F - no right censoring

  #--------------------------------------------------
  # useful variables
  #--------------------------------------------------
  n<-length(X)
  n.nvs<-length(X_star.nvs)
  t_dv<-sort(unique(X[delta==1]))
  z_dv<-sort(unique(c(Z, Z.nvs)))
  J<-length(t_dv)
  K<-length(z_dv)
  t_dv.w<-c(t_dv,Inf) # With right-censoring, add p_{J+1, k}
  J.w<-J+1

  # data format for EM, (J+1)*K rows
  data.EM<-data.frame(X=rep(t_dv.w, K), Z=rep(z_dv, each=J.w), delta=rep(1, J.w*K))
  data.EM[data.EM$X==Inf, ]$delta<-0
  data.EM[data.EM$X==Inf, ]$X<-t_dv[J]
  data.EM$fixed.w<-apply(X=data.EM, MARGIN=1, function(value, data.X, data.Z, data.delta){
    if(value[3]==0) return(0)
    else return(sum((data.X==value[1])&(data.Z==value[2])&(data.delta==1)))
  }, data.X=X, data.Z=Z, data.delta=delta) # n.obs.events

  #--------------------------------------------------
  # Right censoring
  #--------------------------------------------------
  # sum(delta==0)*2 matrix: col-1: min j s.t. t_j>x, col-2: Z
  if(right.censoring) rc.value<-cbind(sapply(X[delta==0], function(x, t_dv.w) return(which.max(t_dv.w>x)), t_dv.w=t_dv.w),
                                      sapply(Z[delta==0], function(z, z_dv) return(which(z_dv==z)), z_dv=z_dv))

  #--------------------------------------------------
  # Left-truncation
  #--------------------------------------------------
  if(left.truncation){
    if(lt.nvs) {
      l_dv<-sort(unique(c(L, L.nvs)))
      M<-length(l_dv)
      n.k.lt<-sapply(z_dv, function(z, data.z) sum(data.z==z), data.z=c(Z, Z.nvs))

      # Left-truncation time distribution: qm
      if(dist.lt[1]=="unknown"){
        p.lt<-cdfDT(y=c(L, L.nvs), l=rep(-Inf, n+n.nvs), r=c(X, X_star.nvs), display=FALSE)
        p.lt<-p.lt$F-c(0,p.lt$F[-M])
      }
      if(dist.lt[1]=="uniform") p.lt=rep(1/M, M)
    }
    else{
      l_dv<-sort(unique(L))
      M<-length(l_dv)
      n.k.lt<-sapply(z_dv, function(z, data.z) sum(data.z==z), data.z=Z)

      # Left-truncation time distribution: qm
      if(dist.lt[1]=="unknown"){
        p.lt<-cdfDT(y=L, l=rep(-Inf, n), r=X, display=FALSE)
        p.lt<-p.lt$F-c(0,p.lt$F[-M])
      }
      if(dist.lt[1]=="uniform") p.lt=rep(1/M, M)
    }

    # J*(K+1) matrix: col-1~K: sum lm>t_j q_m*(n_k+n.nvs.k), col-K+1: sum lm<=t_j q_m
    lt.value<-sapply(t_dv.w, function(t, l_dv, p.lt) return(sum(p.lt[l_dv>t])), l_dv=l_dv, p.lt=p.lt)
    lt.value<-cbind(lt.value%*%t(n.k.lt), 1-lt.value)
  }

  #--------------------------------------------------
  # With validation set: Certain given uncertain
  #--------------------------------------------------
  if(surrogate.response){
    phi<-calculate.phi.binary.cov(data.obj=list(X=X, Z=Z, L=L, delta=delta, X_star=X_star, delta_star=delta_star,
                                                X_star.nvs=X_star.nvs, L.nvs=L.nvs, delta_star.nvs=delta_star.nvs, Z.nvs=Z.nvs,
                                                t_dv.w=t_dv.w, J.w=J.w),
                                  phi.methods.ds1=phi.methods[1], ind.X_star=ind.X_star, no.t_star.ind=no.t_star.ind,
                                  lt.nvs=lt.nvs, weight.type=weight.type)
    fit.surro.res.events<-phi$fit.surro.res.events
    fit.surro.res.rc<-phi$fit.surro.res.rc
    data.EM$fixed.w<-data.EM$fixed.w+as.numeric(phi$phi)
  }
  else{
    fit.surro.res.events<-NULL
    fit.surro.res.rc<-NULL
  }

  #--------------------------------------------------
  # EM algorithm
  #--------------------------------------------------
  # Initial fit
  fitcox<-coxph(Surv(L,X,delta)~Z)
  coef.0<-fitcox$coefficients
  pcox.m<-survfit(fitcox, newdata=data.frame(Z=0))
  pcox.m<-pcox.m$surv[pcox.m$time%in%t_dv.w[-J.w]]
  pcox.m<-cbind(c(1, pcox.m)-c(pcox.m, 0), c(1, pcox.m^exp(coef.0))-c(pcox.m^exp(coef.0), 0))

  # delete redundant variables
  # rm(list=setdiff(ls(), c("fitcox", "data.EM", "coef.0", "pcox.m",
  #                         "t_dv.w", "z_dv",
  #                         "rc.value", "lt.value", "J.w", "K",
  #                         "fit.surro.res.events",
  #                         "fit.surro.res.rc",
  #                         "haz.method",
  #                         "phi.methods",
  #                         "left.cont",
  #                         "tol",
  #                         "left.truncation",
  #                         "surrogate.response",
  #                         "right.censoring",
  #                         "check",
  #                         "bootstrap.se", "n.boot","bootout")))

  # CHECK
  n.step<-0
  iter.info<-c(n.step, fitcox$coefficients, NA)
  if(check) print(iter.info)

  # iterations
  repeat{

    #--------------------------------------------
    # calculate weights
    #--------------------------------------------
    data.EM$w<-data.EM$fixed.w

    # right censoring
    if(right.censoring){
      update.w<-apply(rc.value, MARGIN=1, f.E1.rc, pcox.m=pcox.m, J.w=J.w)
      update.w<-as.numeric(sapply(1:K, function(k, raw.w, rc.z) return(apply(raw.w[, rc.z==k], sum, MARGIN=1)), raw.w=update.w, rc.z=rc.value[,2]))
      data.EM$w<-data.EM$w+update.w
    }
    if(left.truncation){
      update.w<-t(lt.value[,1:K]*pcox.m)
      update.w<-as.numeric(t(update.w/apply(lt.value[,K+1]*pcox.m, sum, MARGIN=2)))
      data.EM$w<-data.EM$w+update.w
    }

    #--------------------------------------------
    # fit with new weight
    #--------------------------------------------
    fitcox<-coxph(Surv(X, delta)~Z, weights = w, data=data.EM[data.EM$w>0, ])

    #CHECK
    n.step<-n.step+1
    iter.info<-rbind(iter.info, c(n.step, fitcox$coefficients, sum(abs(fitcox$coefficients-coef.0))))
    if(check) print(c(n.step, fitcox$coefficients, sum(abs(fitcox$coefficients-coef.0))))

    if(sum(abs(fitcox$coefficients-coef.0))<=tol) break

    #--------------------------------------------
    # p_jk calculated using last iteration
    #--------------------------------------------
    coef.0<-fitcox$coefficients
    # if(haz.method[1]=="EM") fit.haz<-f.CoxHaz.EM(coef.0=coef.0, z_dv=z_dv, t_dv=t_dv.w, weights=data.EM$w)
    # else fit.haz<-f.CoxHaz.coxph(fitcox=fitcox, t_dv=t_dv)

    # return a (J+1)*K matrix where
    # p.cox[j,k]=p_jk=S(t_j|z_k)*lambda(t_j,z_k)
    # p.cox[J+1,k]=p_jk=S(t_J|z_k)
    # pcox.m<-pcox.mat(haz=fit.haz$haz, cumhaz=fit.haz$cumhaz, coef=coef.0, z_dv=z_dv, left.cont=left.cont, events.only=!right.censoring)

    pcox.m<-survfit(fitcox, newdata=data.frame(Z=0))
    pcox.m<-pcox.m$surv[pcox.m$time%in%t_dv.w[-J.w]]
    pcox.m<-cbind(c(1, pcox.m)-c(pcox.m, 0), c(1, pcox.m^exp(coef.0))-c(pcox.m^exp(coef.0), 0))
  }

  if(perturb.se){

    # coef.0<-fitcox$coefficients
    # pcox.m<-survfit(fitcox, newdata=data.frame(Z=0))
    # pcox.m<-pcox.m$surv[pcox.m$time%in%t_dv.w[-J.w]]
    # pcox.m<-cbind(c(1, pcox.m)-c(pcox.m, 0), c(1, pcox.m^exp(coef.0))-c(pcox.m^exp(coef.0), 0))
    # w.eps<-data.EM$fixed.w
    # # right censoring
    # if(right.censoring){
    #   update.w<-apply(rc.value, MARGIN=1, f.E1.rc, pcox.m=pcox.m, J.w=J.w)
    #   update.w<-as.numeric(sapply(1:K, function(k, raw.w, rc.z) return(apply(raw.w[, rc.z==k], sum, MARGIN=1)), raw.w=update.w, rc.z=rc.value[,2]))
    #   w.eps<-w.eps+update.w
    # }
    # if(left.truncation){
    #   update.w<-t(lt.value[,1:K]*pcox.m)
    #   update.w<-as.numeric(t(update.w/apply(lt.value[,K+1]*pcox.m, sum, MARGIN=2)))
    #   w.eps<-w.eps+update.w
    # }

    w.eps<-data.EM$w
    c.jk<-w.eps*data.EM$Z
    deriv.eps<-numeric(0)
    for(eps in c(-1, 1)/(n+n.nvs)){
      # coef.eps<-fitcox$coefficients+eps
      # a.h.eps<-matrix(c.jk*exp(coef.eps*data.EM$Z), ncol=K)
      # a.h.eps<-apply(a.h.eps, sum, MARGIN=1)
      # haz.eps<-f.CoxHaz.EM(coef.0=coef.eps, weights=w.eps,
      #                      z_dv=z_dv, t_dv=t_dv,
      #                      mat.by.row=FALSE)$haz
      # b.j.eps<-sapply(1:J, function(h, value) sum(value[h:length(value)]), value=a.h.eps)
      # deriv.eps<-c(deriv.eps, (sum(c.jk[data.EM$delta==1])+sum(b.j.eps*haz.eps))/eps)

      coef.eps<-fitcox$coefficients+eps
      cumhaz.eps<-f.CoxHaz.EM(coef.0=coef.eps, weights=w.eps,
                              z_dv=z_dv, t_dv=t_dv,
                              mat.by.row=FALSE)$cumhaz
      score.eps<-rep(c(cumhaz.eps, cumhaz.eps[J]), 2)*exp(data.EM$Z*coef.eps)
      score.eps<-c.jk*(data.EM$delta-score.eps)
      deriv.eps<-c(deriv.eps, -sum(score.eps)/eps)
    }
    A<-1/mean(deriv.eps)

    eps<-0
    coef.eps<-fitcox$coefficients+eps
    cumhaz.eps<-f.CoxHaz.EM(coef.0=coef.eps, weights=w.eps,
                            z_dv=z_dv, t_dv=t_dv,
                            mat.by.row=FALSE)$cumhaz
    score.eps<-rep(c(cumhaz.eps, cumhaz.eps[J]), 2)*exp(data.EM$Z*coef.eps)
    score.eps<-c.jk*(data.EM$delta-score.eps)
    B<-sum(score.eps^2)


    perturb.se<-data.frame(se=sqrt(A), robust.se=sqrt(A*B*A))
  }

  if(bootstrap.se){
    if(!is.na(seed.boot)) set.seed(seed.boot)
    fit.boot<-replicate(n.boot, boot.SurroCOX(X=X, Z=Z, delta=delta, L=L, X_star=X_star, delta_star=delta_star,
                                              L.nvs=L.nvs, Z.nvs=Z.nvs, X_star.nvs=X_star.nvs, delta_star.nvs=delta_star.nvs,
                                              tol=tol,
                                              dist.lt=dist.lt,
                                              phi.methods=phi.methods,
                                              haz.method=haz.method,
                                              left.cont=left.cont,
                                              ind.X_star=ind.X_star,
                                              no.t_star.ind=no.t_star.ind,
                                              lt.nvs=lt.nvs))
    fit.boot<-sort(fit.boot)
    bootstrap.se<-data.frame(se=sd(fit.boot), CI.L=fit.boot[round(n.boot*(1-CI.level)/2)], CI.U=fit.boot[round(n.boot*(1+CI.level)/2)])
  }

  if(bootout) return(fitcox$coefficient)
  else return(list(fitcox=fitcox, iter.info=iter.info, data.EM=data.EM,
                   haz.method=haz.method,
                   phi.model=phi.methods,
                   surv.cumhaz=c("current", "previous")[left.cont+1],
                   fit.surro.res.events=fit.surro.res.events,
                   fit.surro.res.rc=fit.surro.res.rc,
                   bootstrap.se=bootstrap.se,
                   perturb.se=perturb.se,
                   LC.as.exact.t_star=LC.as.exact.t_star))

}





