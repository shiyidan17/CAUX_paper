#=====================================
# Getting estimates from coxph
#=====================================
# f.estCox<-function(fitcox, J, z_dv, t_dv){
#   coef.est<-fitcox$coefficients
#   cumhaz.est<-basehaz(fitcox, centered = FALSE)$hazard
#   haz.est<-cumhaz.est-c(0, cumhaz.est[-J])
#
#   return(list(coef.0=coef.est,
#               cumhaz.0=cumhaz.est,
#               haz.0=haz.est))
# }

# pcox<-function(j, k, haz=haz[j], cumhaz=cumhaz[j], coef, z.value=z_dv[k]){
#   # return only p.cox[j,k]
#   return(exp(log(haz)+coef*z.value-cumhaz*exp(coef*z.value)))
# }

# Initial weights, count validation set
f.ini.weights<-function(compare.value, data.X, data.Z, data.delta){

  data.X<-data.X[(data.Z==compare.value[2])&(data.delta==compare.value[3])]

  # for delta==1
  if(compare.value[3]==1) return(sum((data.X==compare.value[1])))
  # for delta==0
  else return(sum((data.X>=compare.value[1])))
}


pcox.mat<-function(haz, cumhaz, coef, z_dv, left.cont=FALSE, events.only=FALSE){
  # return a (J+1)*K matrix where
  # p.cox[j,k]=p_jk=S(t_j|z_k)*lambda(t_j,z_k)
  # p.cox[J+1,k]=p_jk=S(t_J|z_k)
  cov.value<-t(exp(coef*z_dv))
  v1<-log((haz%*%cov.value))
  v2<--cumhaz%*%cov.value

  if(left.cont) v2<-rbind(c(0,0), v2[-nrow(v2), ])

  if(events.only) return(exp(v1+v2))
  value<-rbind(v1+v2, v2[nrow(v2), ])
  return(exp(value))
}

#=====================================
# Calculate weights
#=====================================
f.E1<-function(j, k, pcox.k=pcox.k, t_dv.w=t_dv.w, z_dv=z_dv, X=X, delta=delta, Z=Z){
  # j<-which(t_dv.w==max(X[delta==0]))

  t.j<-t_dv.w[j]
  pcox.jk<-pcox.k[j]

  # Keep only Z=z.k and X=t.j
  V1<-(Z==z_dv[k])
  V2<-delta&(X==t.j)
  V3.nmr<-rep(0,length(X))
  if(sum((1-delta)&(X<t.j))){
    tmp.id<-(1-delta)&(X<t.j)
    #tmp.id<-(1-delta)&(X<t.j)&V1
    V3.nmr[tmp.id]<-pcox.jk #!!!!Browse[2]> pcox.jk [1] 0
  }

  # Browse[2]> sum(pcox.k==0)
  # [1] 49
  # Browse[2]> sum(delta==0)
  # [1] 49

  pos<-sapply(t_dv.w, function(t.u, X.vec){(X.vec<t.u)}, X.vec=X) #n*J+1, pos[i, j]=I(Xi<t.j)*I(Z_i=z.k)
  pos<-pos*V1
  V3.dnm<-pos%*%pcox.k

  value<-rep(0, length(V3.dnm))
  keep.id<-(V3.dnm!=0)
  value[keep.id]<-exp(log(V3.nmr[keep.id])-log(V3.dnm[keep.id]))
  value[!is.finite(value)]<-0

  #value<-exp(log(V3.nmr)-log(V3.dnm))

  value<-sum((V2+value)[V1])

  #print(c(j,k,value)) #CHECK





  return(value)
}

#=========================================================

#=========================================================
# Entire weights
# This function has redundant calculation, so it's for checking purpose.
f.weights<-function(t.value, z.value,
                    inputObjs=inputObjs,
                    pcox.m=pcox.m,
                    check=TRUE, new.E2=FALSE){
  # pass environment variables
  # inputObjs<-list(X=X, delta=delta, L=L, Z=Z, X_star=X_star,delta_star=delta_star,
  #                 L.nvs=L.nvs, Z.nvs=Z.nvs, X_star.nvs=X_star.nvs,delta_star.nvs=delta_star.nvs,
  #                 n=n, n.nvs=n.nvs, J=J, K=K, M=M,
  #                 t_dv=t_dv, l_dv=l_dv, z_dv=z_dv, t_dv.w=t_dv.w,
  #                 phi_delta0=phi_delta0, phi_delta1=phi_delta1,
  #                 p.lt=p.lt)
  list2env(inputObjs, envir = environment())
  rm(inputObjs)


  # j, k
  j<-which(t_dv.w==t.value)
  k<-which(z_dv==z.value)
  pcox.k<-pcox.m[, k]
  rm(pcox.m)

  # validation set
  E1<-f.E1(j=j, k=k, pcox.k=pcox.k, t_dv.w=t_dv.w, z_dv=z_dv, X=X, delta=delta, Z=Z)

  # nonvalidation set
  if(!n.nvs) E2<-0
  else{
    if(new.E2) E2<-phi[j,k]*pcox.k[j]
    else E2<-phi[j,k]
  }

  # truncation sample
  if(max(l_dv)){
    n.k<-sum(Z==z.value)
    n.nvs.k<-sum(Z.nvs==z.value)
    pos<-sapply(t_dv.w, function(t, l_dv){(l_dv>t)}, l_dv=l_dv) #M*J pos[m, j]=I(l_m>t.j)
    nmr<-pcox.k[j]*t(pos[,j])%*%p.lt
    dnm<-t(pcox.k)%*%t(!pos)%*%p.lt
    E3<-(n.k+n.nvs.k)*nmr/dnm
  }
  else E3<-0


  if(check) return(c(E1,E2,E3))
  else return(E1+E2+E3)

}
#=========================================================







