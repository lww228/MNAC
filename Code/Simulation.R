

## DGP 1: Errors with homoscedasticity and cross-sectional independence
DGP1E=function(T,N){
  E=mvrnorm(T,rep(0,N),diag(N))
  return(E)
}



## DGP 2: Errors with heteroscedasticity and cross-sectional correlation
DGP2E=function(T,N,S){
  E=matrix(0,T,N)
  
  
  for(t in 1:T){
    if(t%%2==0){
      E[t,]=0.2*mvrnorm(1,rep(0,N),S)
    }
    else{
      E[t,]=0.2*mvrnorm(1,rep(0,N),S)+mvrnorm(1,rep(0,N),S)
    }
  }
  
  return(E)
}



## DGP 3: Errors with heteroscedasticity and cross-sectional correlation
DGP3E=function(T,N,S){
  E=matrix(0,T,N)
  
  for(t in 2:T){
    E[t,]=mvrnorm(1,rep(0,N),S)+E[t-1,]
  }
  return(E)
}



sim_generate=function(DGP,T,N,NS,r,rS,p,q,pre_beta,G,S){
  
  SS=matrix(0,N,N)
  for(i in 1:N){
    for(j in i:N){
      SS[i,j]=0.3^{abs(i-j)}
      SS[j,i]=SS[i,j]
    }
  }
  
  if(DGP==1){
    E=DGP1E(T,N)
  }
  else if(DGP==2){
    E=DGP2E(T,N,SS)
  }
  else{
    E=DGP3E(T,N,SS)
  }
  
  
  Fc=matrix(runif(T*r),T,r)
  Lamc=matrix(runif(N*r,min=-2,max=2),N,r)
  
  Fj=vector("list",S)
  Lamj=vector("list",S)
  for(j in 1:S){
    Fj[[j]]=matrix(rnorm(T*rS[j]),T,rS[j])
    Lamj[[j]]=matrix(rnorm(NS[j]*rS[j]),NS[j],rS[j])
  }
  
  
  Y=matrix(0,T,N)
  X=vector("list",N)
  beta=vector("list",N)
  
  for(i in 1:N){
    Xi=matrix(runif(T*p[i],min=-1,max=1),T,p[i])
    Xi=scale(Xi)
    Xi[,1]=1
    
    X[[i]]=Xi
    
    
    gi=G[i] 
    betai=runif(3,min=1.5,max=2.5)
   
    
    betai=c(betai,rep(0,p[i]-q[i]))
    
    beta[[i]]=betai
    
    rani=which(i==Group[[gi]]) 
    Y[,i]=Fc%*%Lamc[i,]+Fj[[gi]]%*%Lamj[[gi]][rani,]+E[,i]
    
  }
  
  
  out=list()
  out$E=E
  out$Fc=Fc
  out$Lamc=Lamc
  out$Fj=Fj
  out$Lamj=Lamj
  out$X=X
  out$beta=beta
  out$Y=Y
  
  return(out)
  
  
}






################### Performance measures ###################
betaMSE=function(newbeta,beta0){
  N=length(newbeta)
  preMSE=rep(0,N)
  for(i in 1:N){
    dif2=(as.vector(newbeta[[i]])-beta0[[i]])^2
    preMSE[i]=mean(dif2)
  }
  finMSE=mean(preMSE)
  return(finMSE)
}



betaabs=function(newbeta,beta0){
  N=length(newbeta)
  preabs=rep(0,N)
  for(i in 1:N){
    dif=abs(as.vector(newbeta[[i]])-beta0[[i]])
    preabs[i]=mean(dif)
  }
  finabs=mean(preabs)
  return(finabs)
}


betaMAX=function(newbeta,beta0){
  N=length(newbeta)
  preMAX=rep(0,N)
  for(i in 1:N){
    dif=abs(as.vector(newbeta[[i]])-beta0[[i]])
    preMAX[i]=max(dif)
  }
  finMAX=max(preMAX)
  return(finMAX)
}




## Some Indicator functions ##
## Indicator function 1
Indic=function(a,b){
  if(a==b){
    return(1)}
  else{
    return(0)}
}



## Indicator function 2

Indic2=function(a,b){
  if(a!=b){
    return(1)
  }
  else{
    return(0)
  }
}


## Indicator function 3
Indic3=function(a,b,c,d){
  if((a==b)&(c==d)){
    return(1)
  }
  else{
    return(0)
  }
}

## Indicator function 4
Indic4=function(a,b,c,d){
  if((a!=b)&(c!=d)){
    return(1)
  }
  else{
    return(0)
  }
}




betaTPR=function(newbeta,beta0){
  denom=0 
  numer=0 
  for(i in 1:length(newbeta)){
    for(k in 1:length(newbeta[[i]])){
      denom=denom+Indic2(beta0[[i]][k],0)
      numer=numer+Indic4(newbeta[[i]][k],0,beta0[[i]][k],0)
    }
  }
  finTPR=numer/denom
  return(finTPR)
}




betaTNR=function(newbeta,beta0){
  denom=0 
  numer=0 
  for(i in 1:length(newbeta)){
    for(k in 1:length(newbeta[[i]])){
      denom=denom+Indic(beta0[[i]][k],0)
      numer=numer+Indic3(newbeta[[i]][k],0,beta0[[i]][k],0)
    }
  }
  finTNR=numer/denom
  return(finTNR)
}


gmempro=function(newG,G0){
  finpro=sum(newG==G0)/length(G0)
  return(finpro)
}





maxbf=function(Bhat,Fhat,B,F){
  Nj=nrow(B)
  T=nrow(F)
  k=1
  reBF=rep(0,Nj*T)
  
  for(u in 1:Nj){
    for(v in 1:T){
      reBF[k]=abs(sum(Bhat[u,]*Fhat[v,])-sum(B[u,]*F[v,]))
      k=k+1
    }
  }
  out=max(reBF)
  return(out)
}


normF=function(C){
  resF=sqrt(sum(C^2))
  return(resF)
}






