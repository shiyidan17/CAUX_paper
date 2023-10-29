boot.SurroCOX<-function(X=NULL, Z=NULL, delta=rep(1,length(X)), L=rep(0,length(X)), X_star=NULL, delta_star=rep(1,length(X_star)),
                        L.nvs=rep(0,length(X_star.nvs)), Z.nvs=NULL, X_star.nvs=NULL, delta_star.nvs=rep(1,length(X_star.nvs)),
                        tol=1e-5,
                        dist.lt=c("unknown", "uniform"),
                        phi.methods=c("Cox","splines","lm","poly"),
                        haz.method=c("EM", "coxph"),
                        left.cont=TRUE,
                        ind.X_star=FALSE,
                        no.t_star.ind=FALSE,
                        lt.nvs=FALSE){
  est.boot<-NULL

  repeat{
    sample.vs<-sample.int(length(X),replace = TRUE)
    sample.nvs<-sample.int(length(X_star.nvs),replace = TRUE)

    datamat.vs<-cbind(X,Z,delta,L,X_star,delta_star)[sample.vs, ]
    datamat.nvs<-cbind(L.nvs,Z.nvs,X_star.nvs,delta_star.nvs)[sample.nvs, ]
    try(est.boot<-SurroCOX(X=datamat.vs[,1], Z=datamat.vs[,2], delta=datamat.vs[,3],
                           L=datamat.vs[,4], X_star=datamat.vs[,5], delta_star=datamat.vs[,6],
                           L.nvs=datamat.nvs[,1], Z.nvs=datamat.nvs[,2],
                           X_star.nvs=datamat.nvs[,3], delta_star.nvs=datamat.nvs[,4],
                           tol=tol,
                           dist.lt=dist.lt,
                           phi.methods=phi.methods,
                           haz.method=haz.method,
                           left.cont=left.cont,
                           check=FALSE,
                           ind.X_star=ind.X_star,
                           no.t_star.ind=no.t_star.ind,
                           lt.nvs=lt.nvs,
                           bootstrap.se=FALSE,
                           bootout=TRUE, perturb.se = FALSE))
    if(!is.null(est.boot)) break
  }

  return(est.boot)
}



