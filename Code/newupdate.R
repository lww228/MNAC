

############################ Update {g1,...,gN} ##########################################

update_g=function(Fc,Lamc,Fj,X,beta,Y,gmem,A2,Am12,S,alpha,T,N){
  
  ghatk=gmem

  for(i in 1:N){
    giobj=rep(0,S)
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,]
    for(j in 1:S){
      Lamji=t(Fj[[j]])%*%prei/T
      eps=prei-Fj[[j]]%*%Lamji
      giobj[j]=t(eps)%*%eps
    }
    
    for(l in (1:N)[-i]){
      inl=ghatk[l]
      giobj[inl]=giobj[inl]+2*alpha*Am12[i,l]
      giobj[-inl]=giobj[-inl]+2*alpha*A2[i,l]
    }
    
    ghatk[i]=which.min(giobj)
    
  }
  

  return(ghatk)
}

########################### End for Update {g1,...,gN} ############################


############################ Update {beta1,...,betaN} ##########################################

update_beta=function(Fc,Lamc,Fj,X,beta,Y,gmem,T,N,kapN){
  betahatk=vector("list",N)
  Ystar=matrix(0,T,N)
  for(i in 1:N){
    gik=gmem[i]
    Lamj_gik=t(Fj[[gik]])%*%(Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,])/T
    Ystar[,i]=Y[,i]-Fc%*%Lamc[i,]-Fj[[gik]]%*%Lamj_gik
    fit=ncvreg(X[[i]][,-1],Ystar[,i],penalty = "SCAD",lambda = kapN[i], warn = FALSE)
    betahatk[[i]]=fit$beta
  }
  return(betahatk)
}

######################### End for Update {beta1,...,betaN} ##############################



####################### Update Fc and Lamc #######################
update_FcLamc=function(Fc,Lamc,Fj,X,beta,Y,gmem,T,N,r){
  Zc=matrix(0,T,N)
  for(i in 1:N){
    gik=gmem[i]
    Lamj_gik=t(Fj[[gik]])%*%(Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,])/T
    Zc[,i]=Y[,i]-X[[i]]%*%beta[[i]]-Fj[[gik]]%*%Lamj_gik
  }
  Fchatk=sqrt(T)*eigen(Zc%*%t(Zc))$vectors[,1:r]
  Lamchatk=t(Zc)%*%Fchatk/T
  
  out=list()
  out$Fchatk=Fchatk
  out$Lamchatk=Lamchatk
  return(out)
}

######################### End for Update Fc and Lamc ####################


###################### Update Fj and Lamj #####################
update_FjLamj=function(Fc,Lamc,X,beta,Y,gmem,S,T,rS){
  Fjhatk=vector("list",S)
  Lamjhatk=vector("list",S)
  for(j in 1:S){
    iIndex=which(gmem==j)
    Njhat=length(iIndex)
    Zj=matrix(0,T,Njhat)
    ll=1
    for(l in iIndex){
      Zj[,ll]=Y[,l]-X[[l]]%*%beta[[l]]-Fc%*%Lamc[l,]
      ll=ll+1
    }
    Fjhatk[[j]]=sqrt(T)*eigen(Zj%*%t(Zj))$vectors[,1:rS[j]]
    Lamjhatk[[j]]=t(Zj)%*%Fjhatk[[j]]/T
  }
  out=list()
  out$Fjhatk=Fjhatk
  out$Lamjhatk=Lamjhatk
  return(out)
}


######################## End for Update Fj and Lamj ########################




###################### Implement Ando and Bai (2017) JASA #########################

############################ Update {g1,...,gN} ##########################################

Baiupdate_g=function(Fc,Lamc,Fj,X,beta,Y,S,T,N){
  

  ghatk=rep(0,N)
  
  for(i in 1:N){
    giobj=rep(0,S)
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,]
    for(j in 1:S){
      Lamji=t(Fj[[j]])%*%prei/T
      eps=prei-Fj[[j]]%*%Lamji
      giobj[j]=t(eps)%*%eps
    }
    ghatk[i]=which.min(giobj)
    
  }
  

  return(ghatk)
}

########################### End for Update {g1,...,gN} ############################












