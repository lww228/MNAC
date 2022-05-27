
####################### Objective function ######################
SCADfun=function(the,kap,gam){
  if(the<=kap){
    resSCAD=kap*the
  }
  else if((kap<the)&(the<=(gam*kap))){
    resSCAD=(gam*kap*the-0.5*(the*the+kap*kap))/(gam-1)
  }
  else{
    resSCAD=kap*kap*(gam*gam-1)/(2*(gam-1))
  }
  return(resSCAD)
}



obj=function(Y,X,beta,gmem,Fc,Lamc,Fj,T,N,alpha,A,kapN,gamma){
  Lobj=0
  for(i in 1:N){
    gi=gmem[i]
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,]
    Lgii=t(Fj[[gi]])%*%prei/T
    pre_eps1=prei-Fj[[gi]]%*%Lgii
    eps1=t(pre_eps1)%*%pre_eps1
    
    Lobj=Lobj+eps1
  }
  
  
  A0hat=matrix(0,N,N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(gmem[i]==gmem[j]){A0hat[i,j]=A0hat[j,i]=1}
    }
  }
  
  Lobj=Lobj+alpha*sum((A0hat-A)^2)
  
  return(Lobj)
}



Baiobj=function(Y,X,beta,gmem,Fc,Lamc,Fj,T,N,kapN,gamma){
  Lobj=0
  for(i in 1:N){
    gi=gmem[i]
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,]
    Lgii=t(Fj[[gi]])%*%prei/T
    pre_eps1=prei-Fj[[gi]]%*%Lgii
    eps1=t(pre_eps1)%*%pre_eps1
    
    Lobj=Lobj+eps1
  }
  
  
  return(Lobj)
}

####################### End for Objective function ######################

###### How much variation is left in the error term ######

R12=function(Y,X,beta,gmem,Fc,Lamc,Fj,T,N){
  Lobj=0
  
  YXb=0
  for(i in 1:N){
    gi=gmem[i]
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%Lamc[i,]
    Lgii=t(Fj[[gi]])%*%prei/T
    pre_eps1=prei-Fj[[gi]]%*%Lgii
    eps1=t(pre_eps1)%*%pre_eps1
  
    Lobj=Lobj+eps1
    
    addprei=Y[,i]-X[[i]]%*%beta[[i]]
    YXb=YXb+t(addprei)%*%addprei
  }
  
  R1=Lobj/sum(Y^2)
  R2=Lobj/YXb
  
  out=list()
  out$R1=R1
  out$R2=R2
  
  return(out)
}


############################### Add Price of Risk ##############################

PoR=function(esg,esbeta,esLamc,esLamj,Y,X,N,S){
  frp=matrix(0,S,10)
  profrp=matrix(0,S,10)
  
  lr=ncol(esLamc)
  
  rhat=rep(0,N)
  for(i in 1:N){
    rhat[i]=mean(Y[,i]-X[[i]]%*%esbeta[[i]])
  }
  
  for(j in 1:S){
    ind=which(esg==j)
    
    rjhat=rhat[ind]
    Lamcjhat=esLamc[ind,]
    Lamjhat=esLamj[[j]]
    lrj=ncol(Lamjhat)
    
    res=lm(rjhat~Lamcjhat+Lamjhat)
    coesig=summary(res)$coefficients[,"Pr(>|t|)"][-1]
    
    coeLamcj=coesig[1:lr]
    coeLamj=tail(coesig,lrj)
    
    frp[j,]=c(sum(coeLamcj<0.001),sum(coeLamcj<0.01),sum(coeLamcj<0.05),sum(coeLamcj<0.1),lr,
              sum(coeLamj<0.001),sum(coeLamj<0.01),sum(coeLamj<0.05),sum(coeLamj<0.1),lrj)
    profrp[j,]=c(sum(coeLamcj<0.001)/lr,sum(coeLamcj<0.01)/lr,sum(coeLamcj<0.05)/lr,sum(coeLamcj<0.1)/lr,lr/lr,
                 sum(coeLamj<0.001)/lrj,sum(coeLamj<0.01)/lrj,sum(coeLamj<0.05)/lrj,sum(coeLamj<0.1)/lrj,lrj/lrj)
    
  }
  
  out=list()
  out$frp=frp
  out$profrp=profrp
  return(out)
  
}


########################### End for Add Price of Risk ##########################







Outloop=function(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,A,A2,Am12,S,alpha,T,N,kapN,r,rS,gamma,max_step,show_step,tol0){
  
  obj_fun=NULL
  iter_gmem=NULL
  
  oldg=Ingmem
  oldbeta=Inbeta
  oldFc=InFc
  oldLamc=InLamc
  oldFj=InFj
  oldLamj=InLamj
  
  k=1
  
  newg=update_g(oldFc,oldLamc,oldFj,X,oldbeta,Y,oldg,A2,Am12,S,alpha,T,N)
  #newbeta=update_beta(oldFc,oldLamc,oldFj,X,oldbeta,Y,newg,T,N,kapN)
  newbeta=oldbeta
  newFcLamc=update_FcLamc(oldFc,oldLamc,oldFj,X,newbeta,Y,newg,T,N,r)
  newFc=newFcLamc$Fchatk
  newLamc=newFcLamc$Lamchatk
  newFjLamj=update_FjLamj(newFc,newLamc,X,newbeta,Y,newg,S,T,rS)
  newFj=newFjLamj$Fjhatk
  newLamj=newFjLamj$Lamjhatk
  
  
  obj_fun[k]=obj(Y,X,newbeta,newg,newFc,newLamc,newFj,T,N,alpha,A,kapN,gamma)
  iter_gmem=rbind(iter_gmem,newg)
  elapsed0 = proc.time()[3]
  
  tol=max(abs(newFc%*%t(newLamc)-oldFc%*%t(oldLamc)))
  
  while((tol>tol0)&(k<max_step)){
    k=k+1
    
    oldg=newg
    oldbeta=newbeta
    oldFc=newFc
    oldLamc=newLamc
    oldFj=newFj
    oldLamj=newLamj
    
    newg=update_g(oldFc,oldLamc,oldFj,X,oldbeta,Y,oldg,A2,Am12,S,alpha,T,N)
    #newbeta=update_beta(oldFc,oldLamc,oldFj,X,oldbeta,Y,newg,T,N,kapN)
    newbeta=oldbeta
    newFcLamc=update_FcLamc(oldFc,oldLamc,oldFj,X,newbeta,Y,newg,T,N,r)
    newFc=newFcLamc$Fchatk
    newLamc=newFcLamc$Lamchatk
    newFjLamj=update_FjLamj(newFc,newLamc,X,newbeta,Y,newg,S,T,rS)
    newFj=newFjLamj$Fjhatk
    newLamj=newFjLamj$Lamjhatk
    
    
    obj_fun[k]=obj(Y,X,newbeta,newg,newFc,newLamc,newFj,T,N,alpha,A,kapN,gamma)
    iter_gmem=rbind(iter_gmem,newg)
    tol=max(abs(newFc%*%t(newLamc)-oldFc%*%t(oldLamc)))
    
    #if(k%%show_step==0){
    #   print(paste0("step:",k,"time elapsed:",proc.time()[3] - elapsed0))
    #}
    
  }
  
  out=list()
  out$newg=newg
  out$newbeta=newbeta
  out$newFc=newFc
  out$newLamc=newLamc
  out$newFj=newFj
  out$newLamj=newLamj
  out$obj_fun=obj_fun
  out$iter_gmem=iter_gmem
  
  return(out)
  
}







BaiOutloop=function(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,S,T,N,kapN,r,rS,gamma,max_step,show_step,tol0){
  
  obj_fun=NULL
  iter_gmem=NULL
  
  oldg=Ingmem
  oldbeta=Inbeta
  oldFc=InFc
  oldLamc=InLamc
  oldFj=InFj
  oldLamj=InLamj
  
  k=1
  
  newg=Baiupdate_g(oldFc,oldLamc,oldFj,X,oldbeta,Y,S,T,N)
  #newbeta=update_beta(oldFc,oldLamc,oldFj,X,oldbeta,Y,newg,T,N,kapN)
  newbeta=oldbeta
  newFcLamc=update_FcLamc(oldFc,oldLamc,oldFj,X,newbeta,Y,newg,T,N,r)
  newFc=newFcLamc$Fchatk
  newLamc=newFcLamc$Lamchatk
  newFjLamj=update_FjLamj(newFc,newLamc,X,newbeta,Y,newg,S,T,rS)
  newFj=newFjLamj$Fjhatk
  newLamj=newFjLamj$Lamjhatk
  
  
  obj_fun[k]=Baiobj(Y,X,newbeta,newg,newFc,newLamc,newFj,T,N,kapN,gamma)
  iter_gmem=rbind(iter_gmem,newg)
  elapsed0 = proc.time()[3]
  
  tol=max(abs(newFc%*%t(newLamc)-oldFc%*%t(oldLamc)))
  
  while((tol>tol0)&(k<max_step)){
    k=k+1
    
    oldg=newg
    oldbeta=newbeta
    oldFc=newFc
    oldLamc=newLamc
    oldFj=newFj
    oldLamj=newLamj
    
    newg=Baiupdate_g(oldFc,oldLamc,oldFj,X,oldbeta,Y,S,T,N)
    #newbeta=update_beta(oldFc,oldLamc,oldFj,X,oldbeta,Y,newg,T,N,kapN)
    newbeta=oldbeta
    newFcLamc=update_FcLamc(oldFc,oldLamc,oldFj,X,newbeta,Y,newg,T,N,r)
    newFc=newFcLamc$Fchatk
    newLamc=newFcLamc$Lamchatk
    newFjLamj=update_FjLamj(newFc,newLamc,X,newbeta,Y,newg,S,T,rS)
    newFj=newFjLamj$Fjhatk
    newLamj=newFjLamj$Lamjhatk
    
    
    
    obj_fun[k]=Baiobj(Y,X,newbeta,newg,newFc,newLamc,newFj,T,N,kapN,gamma)
    iter_gmem=rbind(iter_gmem,newg)
    tol=max(abs(newFc%*%t(newLamc)-oldFc%*%t(oldLamc)))
    
    #if(k%%show_step==0){
    #  print(paste0("step:",k,"time elapsed:",proc.time()[3] - elapsed0))
    #}
    
  }
  
  out=list()
  out$newg=newg
  out$newbeta=newbeta
  out$newFc=newFc
  out$newLamc=newLamc
  out$newFj=newFj
  out$newLamj=newLamj
  out$obj_fun=obj_fun
  out$iter_gmem=iter_gmem
  
  return(out)
  
}



comparison=function(DGP,n,T,N,NS,r,rS,p,q,pre_beta,A,A2,Am12,S,gamma,max_step,show_step,tol0,
                    Gperm,nG,Gcom,G,cankapi,canS,cank,cankj,canalp,C,tolMS,maxMS,ep,nuA,obsA,zA,zA2,zAm12){
  
  
  ##################### beta coefficients estimation #########################
  reMSE=matrix(0,n,3)
  reTPR=matrix(0,n,3)
  reTNR=matrix(0,n,3)
  
  sumMSE=matrix(0,2,3)
  sumTPR=matrix(0,2,3)
  sumTNR=matrix(0,2,3)
  ##################### Group membership estimation ##########################
  rePRO=NULL
  
  sumPRO=matrix(0,2,3)
  
  
  ##################### Factor model estimation ####################
  FcLcMSE=matrix(0,n,3)
  FcLcUE=matrix(0,n,3)
  
  sumFcLcMSE=matrix(0,2,6)
  sumFcLcUE=matrix(0,2,6)
  
  
  FjLjMSE=NULL
  sumFjLjMSE=matrix(0,2,3)
 
  
  
  #################### R1 and R2 #####################################
  eR1=matrix(0,n,3)
  eR2=matrix(0,n,3)
  
  sumeR1=matrix(0,2,6)
  sumeR2=matrix(0,2,6)
  
  #################### Compute the Price of Risk ###################
  ProPoR=list()
  Pro1PoR=list()
  BaiPoR=list()
  
  sumProfrp=matrix(0,S,10)
  sumPro1frp=matrix(0,S,10)
  sumBaifrp=matrix(0,S,10)
  
  TProfrp=matrix(0,S,10)
  TPro1frp=matrix(0,S,10)
  TBaifrp=matrix(0,S,10)
  
  indfrp1=0
  indfrp2=0
  
  
  ########################## Record time ###################################
  Pro_time=matrix(0,n,3)
  Pro1_time=matrix(0,n,3)
  Bai_time=matrix(0,n,3)
  
  Pro_timems=matrix(0,2,3)
  Pro1_timems=matrix(0,2,3)
  Bai_timems=matrix(0,2,3)
  ######################### Save results ###############################
  
  ## Save weighting coefficients estimation ##    
  saveprea=matrix(0,n,nuA)
  saveweia=matrix(0,n,nuA)
  sumweia=matrix(0,2,nuA)
  saveA=vector("list",n)
  
  
  
  ## Proposed ## 
  saveprog=matrix(0,n,N)
  saveprobeta=vector("list",n)
  saveproFc=vector("list",n)
  saveproLamc=vector("list",n)
  saveproFj=vector("list",n)
  saveproLamj=vector("list",n)
  saveproobj=vector("list",n)
  saveprogmem=vector("list",n)
  


  ## Proposed 1 ## 
  savepro1g=matrix(0,n,N)
  savepro1beta=vector("list",n)
  savepro1Fc=vector("list",n)
  savepro1Lamc=vector("list",n)
  savepro1Fj=vector("list",n)
  savepro1Lamj=vector("list",n)
  savepro1obj=vector("list",n)
  savepro1gmem=vector("list",n)
  
  #########
  
  ## Bai ##
  
  saveBaig=matrix(0,n,N)
  saveBaibeta=vector("list",n)
  saveBaiFc=vector("list",n)
  saveBaiLamc=vector("list",n)
  saveBaiFj=vector("list",n)
  saveBaiLamj=vector("list",n)
  saveBaiobj=vector("list",n)
  saveBaigmem=vector("list",n)
  
  
  ## Save Initial Values ##
  
  saveIngmem=matrix(0,n,N)
  saveInbeta=vector("list",n)
  saveInFc=vector("list",n)
  saveInLamc=vector("list",n)
  saveInFj=vector("list",n)
  saveInLamj=vector("list",n)
  
  
  ## Save optimal tuning parameters ##
  
  saveopS=rep(0,n)
  saveopits=rep(0,n)
  saveopkapN=matrix(0,n,N)
  saveopk=rep(0,n)
  saveopkS=vector("list",n)
  saveopalp=rep(0,n)
  saveopalp1=rep(0,n)
  savefiPIC=rep(0,n)
  savefiPIC1=rep(0,n)

  ## Compute the percentages of under- (U), correct (C), and overidentification (O) of S, r, r1,...,rS ## 
  SUCO=rep(0,3)
  rUCO=rep(0,3)
  rjUCO=matrix(0,S,3)
  torjUCO=rep(0,3)


  for(s in 1:n){
    
    
    sim_res=sim_generate(DGP,T,N,NS,r,rS,p,q,pre_beta,G,S)
    E=sim_res$E
    Fc=sim_res$Fc
    Lamc=sim_res$Lamc
    Fj=sim_res$Fj
    Lamj=sim_res$Lamj
    X=sim_res$X
    beta=sim_res$beta
    Y=sim_res$Y
    
    

    ###################### Ando and Bai (2017) ######################
    Baisighat2=Baiscasig2(X,Y,canS,T,N,cankapi,cank,cankj,gamma,max_step,show_step,tol0,NS)
    
    
    opS=S
    opits=maxMS[1]
    opkapN=rep(0,N)
    opk=r
    opkS=rS
    

    ############################## Initial Parameter Values ##########################
    
    Init=Ini(X,Y,N,T,opS,opk,opkS,opkapN)
    Ingmem=Init$Ingmem
    Inbeta=Init$Inbeta
    InFc=Init$InFc
    InLamc=Init$InLamc
    InFj=Init$InFj
    InLamj=Init$InLamj
    
    
    ## Save Initial Values ##
    saveIngmem[s,]=Ingmem
    saveInbeta[[s]]=Inbeta
    saveInFc[[s]]=InFc
    saveInLamc[[s]]=InLamc
    saveInFj[[s]]=InFj
    saveInLamj[[s]]=InLamj
    
    
    ################# End for Initial Parameter Values ##########################
    
    ####################### Estimation ########################## 
    
    
    Baitime=proc.time()
    Baires=BaiOutloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,opS,T,N,opkapN,opk,opkS,gamma,max_step,show_step,tol0)
    Baiendtime=proc.time()
    
    Bai_time[s,]=c((Baiendtime-Baitime)[1],
                   (Baiendtime-Baitime)[2],
                   (Baiendtime-Baitime)[3])
    
    
    

    GA=Baires$newg
    
    BaiA=matrix(0,N,N)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        if(GA[i]==GA[j]){BaiA[i,j]=BaiA[j,i]=1}
      }
    }
    
    
    ## Calculate weighting coefficients ##
    
    prea=rep(0,nuA)
    for(nu in 1:nuA){
      prea[nu]=1/(sum((BaiA-obsA[[nu]])^2)+ep)
    }
    
    weia=prea/sum(prea)
    
    A=matrix(0,N,N)
    for(nu in 1:nuA){
      A=A+weia[nu]*obsA[[nu]]
    }
    
    
    A2=A^2
    Am12=(1-A)^2
    
    
    
    ## Save weighting coefficients estimation ##  
    saveprea[s,]=prea
    saveweia[s,]=weia
    saveA[[s]]=A
    
    #################
    
  
    
    upalp=updatealp(X,Y,A,A2,Am12,opS,T,N,opkapN,opk,opkS,gamma,max_step,show_step,tol0,
                    C,Baisighat2,canalp)
    opalp=upalp$opt
    fiPIC=upalp$minPIC
    

    
    
    upalp1=updatealp(X,Y,zA,zA2,zAm12,opS,T,N,opkapN,opk,opkS,gamma,max_step,show_step,tol0,
                    C,Baisighat2,canalp)
    opalp1=upalp1$opt
    fiPIC1=upalp1$minPIC
    
    
    
    #########
    
    #end_time=proc.time()-start_time
    #print(end_time)
    #########
    
  
    
    ## Save optimal tuning parameters ##
    
    saveopS[s]=opS
    saveopits[s]=opits
    saveopkapN[s,]=opkapN
    saveopk[s]=opk
    saveopkS[[s]]=opkS
    saveopalp[s]=opalp
    saveopalp1[s]=opalp1
    savefiPIC[s]=fiPIC
    savefiPIC1[s]=fiPIC1
    
    ################################## End for Model Selection ############################################
    

    
    ####################### Estimation ########################## 
    protime=proc.time()
    res=Outloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,A,A2,Am12,opS,opalp,T,N,opkapN,opk,opkS,gamma,max_step,show_step,tol0)
    proendtime=proc.time()
    
    Pro_time[s,]=c((proendtime-protime)[1],
                   (proendtime-protime)[2],
                   (proendtime-protime)[3])
    
    
    

    pro1time=proc.time()
    res1=Outloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,zA,zA2,zAm12,opS,opalp1,T,N,opkapN,opk,opkS,gamma,max_step,show_step,tol0)
    pro1endtime=proc.time()
    
    Pro1_time[s,]=c((pro1endtime-pro1time)[1],
                   (pro1endtime-pro1time)[2],
                   (pro1endtime-pro1time)[3])
    
    
    
    ########################################## Result ############################################
    ## Proposed ##
    pronewg=res$newg
    pronewbeta=res$newbeta
    pronewFc=res$newFc
    pronewLamc=res$newLamc
    pronewFj=res$newFj
    pronewLamj=res$newLamj
    proobj_fun=res$obj_fun
    proiter_gmem=res$iter_gmem
    
  
    
    ## Proposed 1 ##
    pro1newg=res1$newg
    pro1newbeta=res1$newbeta
    pro1newFc=res1$newFc
    pro1newLamc=res1$newLamc
    pro1newFj=res1$newFj
    pro1newLamj=res1$newLamj
    pro1obj_fun=res1$obj_fun
    pro1iter_gmem=res1$iter_gmem
    
  
    #########
    
    
    ## Bai ##
    Bainewg=Baires$newg
    Bainewbeta=Baires$newbeta
    BainewFc=Baires$newFc
    BainewLamc=Baires$newLamc
    BainewFj=Baires$newFj
    BainewLamj=Baires$newLamj
    Baiobj_fun=Baires$obj_fun
    Baiiter_gmem=Baires$iter_gmem
    
  
    
    ##################### beta coefficients estimation  #########################
    reMSE[s,]=c(betaMSE(pronewbeta,beta),betaMSE(pro1newbeta,beta),betaMSE(Bainewbeta,beta))
    reTPR[s,]=c(betaTPR(pronewbeta,beta),betaTPR(pro1newbeta,beta),betaTPR(Bainewbeta,beta))
    reTNR[s,]=c(betaTNR(pronewbeta,beta),betaTNR(pro1newbeta,beta),betaTNR(Bainewbeta,beta))
    

    ##################### Group membership estimation ##########################
  
    if(opS==S){
      preproG=rep(0,nG)
      prepro1G=rep(0,nG)
      preBaiG=rep(0,nG)
      for(nn in 1:nG){
        preproG[nn]=gmempro(pronewg,Gcom[nn,])
        prepro1G[nn]=gmempro(pro1newg,Gcom[nn,])
        preBaiG[nn]=gmempro(Bainewg,Gcom[nn,])
      }
      
      rePRO=rbind(rePRO,c(max(preproG),max(prepro1G),max(preBaiG)))
    
    

 
    
    
    
    ##################### Factor model estimation ####################
    

 
      FjLj0=NULL
      proFjLj=vector("list",S)
      pro1FjLj=vector("list",S)
      BaiFjLj=vector("list",S)
      
      
      for(j in 1:S){
        FjLjj=Fj[[j]]%*%t(Lamj[[j]])
        FjLj0=cbind(FjLj0,FjLjj)
        
        proFjLjj=pronewFj[[j]]%*%t(pronewLamj[[j]])
        proFjLj[[j]]=proFjLjj
        
       
        pro1FjLjj=pro1newFj[[j]]%*%t(pro1newLamj[[j]])
        pro1FjLj[[j]]=pro1FjLjj
       
        
        BaiFjLjj=BainewFj[[j]]%*%t(BainewLamj[[j]])
        BaiFjLj[[j]]=BaiFjLjj
      }
      
      
      
      proFjLjMSE=rep(0,nG)
      pro1FjLjMSE=rep(0,nG)
      BaiFjLjMSE=rep(0,nG)
      
      for(nn in 1:nG){
        proM=NULL
        pro1M=NULL
        BaiM=NULL
        for(j in 1:S){
          ind=Gperm[nn,j]
          proM=cbind(proM,proFjLj[[ind]])
          pro1M=cbind(pro1M,pro1FjLj[[ind]])
          BaiM=cbind(BaiM,BaiFjLj[[ind]])
        }
        proFjLjMSE[nn]=sum((FjLj0-proM)^2)/(N*T)
        pro1FjLjMSE[nn]=sum((FjLj0-pro1M)^2)/(N*T)
        BaiFjLjMSE[nn]=sum((FjLj0-BaiM)^2)/(N*T)
      }
      
      FjLjMSE=rbind(FjLjMSE,c(min(proFjLjMSE),min(pro1FjLjMSE),min(BaiFjLjMSE)))
      
      
      
      
      ############################### Add Price of Risk ##############################
      ProPoRs=PoR(pronewg,pronewbeta,pronewLamc,pronewLamj,Y,X,N,S)
      Profrp=ProPoRs$frp
      Proprofrp=ProPoRs$profrp
      
      
      Pro1PoRs=PoR(pro1newg,pro1newbeta,pro1newLamc,pro1newLamj,Y,X,N,S)
      Pro1frp=Pro1PoRs$frp
      Pro1profrp=Pro1PoRs$profrp
      
      
      BaiPoRs=PoR(Bainewg,Bainewbeta,BainewLamc,BainewLamj,Y,X,N,S)
      Baifrp=BaiPoRs$frp
      Baiprofrp=BaiPoRs$profrp
      
      indfrp1=indfrp1+1
      ProPoR[[indfrp1]]=ProPoRs
      Pro1PoR[[indfrp1]]=Pro1PoRs
      BaiPoR[[indfrp1]]=BaiPoRs
      
      sumProfrp=sumProfrp+Proprofrp
      sumPro1frp=sumPro1frp+Pro1profrp
      sumBaifrp=sumBaifrp+Baiprofrp
      
      if((opk==r)&(sum(opkS==rS)==S)){
        indfrp2=indfrp2+1
        TProfrp=TProfrp+Proprofrp
        TPro1frp=TPro1frp+Pro1profrp
        TBaifrp=TBaifrp+Baiprofrp
      }
      
    }
    
    FcLc0=Fc%*%t(Lamc)
    proFcLc=pronewFc%*%t(pronewLamc)
    pro1FcLc=pro1newFc%*%t(pro1newLamc)
    BaiFcLc=BainewFc%*%t(BainewLamc)
    
    
    
    FcLcMSE[s,]=c(sum((proFcLc-FcLc0)^2)/(N*T),sum((pro1FcLc-FcLc0)^2)/(N*T),sum((BaiFcLc-FcLc0)^2)/(N*T))
    
    
    
    FcLcUE[s,]=c(maxbf(pronewLamc,pronewFc,Lamc,Fc),maxbf(pro1newLamc,pro1newFc,Lamc,Fc),maxbf(BainewLamc,BainewFc,Lamc,Fc))
    
    
  
    ######################
    
    ## R1 and R2 ##
    proR12=R12(Y,X,pronewbeta,pronewg,pronewFc,pronewLamc,pronewFj,T,N)
    pro1R12=R12(Y,X,pro1newbeta,pro1newg,pro1newFc,pro1newLamc,pro1newFj,T,N)
    BaiR12=R12(Y,X,Bainewbeta,Bainewg,BainewFc,BainewLamc,BainewFj,T,N)
    
    eR1[s,]=c(proR12$R1,pro1R12$R1,BaiR12$R1)
    eR2[s,]=c(proR12$R2,pro1R12$R2,BaiR12$R2)
    
    ###############
    

    ######################### Save results ###############################
    ## Proposed ## 
    saveprog[s,]=pronewg
    saveprobeta[[s]]=pronewbeta
    saveproFc[[s]]=pronewFc
    saveproLamc[[s]]=pronewLamc
    saveproFj[[s]]=pronewFj
    saveproLamj[[s]]= pronewLamj
    saveproobj[[s]]=proobj_fun
    saveprogmem[[s]]=proiter_gmem
    
  
    ## Proposed 1 ## 
    savepro1g[s,]=pro1newg
    savepro1beta[[s]]=pro1newbeta
    savepro1Fc[[s]]=pro1newFc
    savepro1Lamc[[s]]=pro1newLamc
    savepro1Fj[[s]]=pro1newFj
    savepro1Lamj[[s]]= pro1newLamj
    savepro1obj[[s]]=pro1obj_fun
    savepro1gmem[[s]]=pro1iter_gmem
    
    #########

    ## Bai ##
    saveBaig[s,]=Bainewg
    saveBaibeta[[s]]=Bainewbeta
    saveBaiFc[[s]]=BainewFc
    saveBaiLamc[[s]]=BainewLamc
    saveBaiFj[[s]]=BainewFj
    saveBaiLamj[[s]]=BainewLamj
    saveBaiobj[[s]]=Baiobj_fun
    saveBaigmem[[s]]=Baiiter_gmem
    

    
  }
  
  ####################### Model Selection summary #######################
  
  ## Compute the percentages of under- (U), correct (C), and overidentification (O) of S, r, r1,...,rS ## 
  
  CSind=which(saveopS==S)
  lCS=length(CSind)
  SUCO[1]=sum(saveopS<S)/n; SUCO[2]=lCS/n; SUCO[3]=sum(saveopS>S)/n
  
  CSopk=saveopk[CSind]
  rUCO[1]=sum(CSopk<r)/lCS; rUCO[2]=sum(CSopk==r)/lCS; rUCO[3]=sum(CSopk>r)/lCS
  
  CSopkS=NULL
  for(l in CSind){
    CSopkS=cbind(CSopkS,saveopkS[[l]])
  }
  
  for(j in 1:S){
    rjUCO[j,1]=sum(CSopkS[j,]<rS[j])/lCS
    rjUCO[j,2]=sum(CSopkS[j,]==rS[j])/lCS
    rjUCO[j,3]=sum(CSopkS[j,]>rS[j])/lCS
  }
  
  torjUCO[1]=sum(CSopkS<rS[1])/(S*lCS); torjUCO[2]=sum(CSopkS==rS[1])/(S*lCS); torjUCO[3]=sum(CSopkS>rS[1])/(S*lCS)
  

  
  ##################### beta coefficients estimation summary #########################
  
  colnames(reMSE)=c("Proposed","Proposed1","Bai")
  sumMSE=rbind(colMeans(reMSE),apply(reMSE,2,sd))
  colnames(sumMSE)=c("Proposed","Proposed1","Bai")
  rownames(sumMSE)=c("mean","sd")
  
  colnames(reTPR)=c("Proposed","Proposed1","Bai")
  sumTPR=rbind(colMeans(reTPR),apply(reTPR,2,sd))
  colnames(sumTPR)=c("Proposed","Proposed1","Bai")
  rownames(sumTPR)=c("mean","sd")
  
  colnames(reTNR)=c("Proposed","Proposed1","Bai")
  sumTNR=rbind(colMeans(reTNR),apply(reTNR,2,sd))
  colnames(sumTNR)=c("Proposed","Proposed1","Bai")
  rownames(sumTNR)=c("mean","sd")
  
  ##################### Group membership estimation summary ##########################
  if(lCS>0){
    colnames(rePRO)=c("Proposed","Proposed1","Bai")
    sumPRO=rbind(colMeans(rePRO),apply(rePRO,2,sd))
    colnames(sumPRO)=c("Proposed","Proposed1","Bai")
    rownames(sumPRO)=c("mean","sd")
    
  
  ##################### Factor model estimation summary ####################
  
  
    colnames(FjLjMSE)=c("Proposed","Proposed1","Bai")
    sumFjLjMSE=rbind(colMeans(FjLjMSE),apply(FjLjMSE,2,sd))
    colnames(sumFjLjMSE)=c("Proposed","Proposed1","Bai")
    rownames(sumFjLjMSE)=c("mean","sd")
    
  }
  
 
  
  
  if(lCS!=1){
    CSFcLcMSE=FcLcMSE[CSind,]
    
    colnames(FcLcMSE)=c("Proposed","Proposed1","Bai")
    colnames(CSFcLcMSE)=c("Proposed","Proposed1","Bai")
    sumFcLcMSE=cbind(rbind(colMeans(FcLcMSE),apply(FcLcMSE,2,sd)),
                     rbind(colMeans(CSFcLcMSE),apply(CSFcLcMSE,2,sd)))
    colnames(sumFcLcMSE)=c("Proposed","Proposed1","Bai","CSProposed","CSProposed1","CSBai")
    rownames(sumFcLcMSE)=c("mean","sd")
    
    
    CSFcLcUE=FcLcUE[CSind,]
    
    colnames(FcLcUE)=c("Proposed","Proposed1","Bai")
    colnames(CSFcLcUE)=c("Proposed","Proposed1","Bai")
    sumFcLcUE=cbind(rbind(colMeans(FcLcUE),apply(FcLcUE,2,sd)),
                    rbind(colMeans(CSFcLcUE),apply(CSFcLcUE,2,sd)))
    colnames(sumFcLcUE)=c("Proposed","Proposed1","Bai","CSProposed","CSProposed1","CSBai")
    rownames(sumFcLcUE)=c("mean","sd")
    
    
    ###### R1 and R2 ######
    CSeR1=eR1[CSind,]
    
    colnames(eR1)=c("Proposed","Proposed1","Bai")
    colnames(CSeR1)=c("Proposed","Proposed1","Bai")
    sumeR1=cbind(rbind(colMeans(eR1),apply(eR1,2,sd)),
                     rbind(colMeans(CSeR1),apply(CSeR1,2,sd)))
    colnames(sumeR1)=c("Proposed","Proposed1","Bai","CSProposed","CSProposed1","CSBai")
    rownames(sumeR1)=c("mean","sd")
    
    

    CSeR2=eR2[CSind,]
    
    colnames(eR2)=c("Proposed","Proposed1","Bai")
    colnames(CSeR2)=c("Proposed","Proposed1","Bai")
    sumeR2=cbind(rbind(colMeans(eR2),apply(eR2,2,sd)),
                     rbind(colMeans(CSeR2),apply(CSeR2,2,sd)))
    colnames(sumeR2)=c("Proposed","Proposed1","Bai","CSProposed","CSProposed1","CSBai")
    rownames(sumeR2)=c("mean","sd")
    
    
    #######################
  }
  
  ####################### Price of Risk summary ######################
  
  sumProfrp=sumProfrp/indfrp1
  sumPro1frp=sumPro1frp/indfrp1
  sumBaifrp=sumBaifrp/indfrp1
  
  TProfrp=TProfrp/indfrp2
  TPro1frp=TPro1frp/indfrp2
  TBaifrp=TBaifrp/indfrp2
  
  ####################################################################
  
  
  ############################# Time summary #######################
  
  colnames(Pro_time)=c("User","System","Elapse")
  Pro_timems=rbind(colMeans(Pro_time),apply(Pro_time,2,sd))
  colnames(Pro_timems)=c("User","System","Elapse")
  rownames(Pro_timems)=c("Pro_time","sd")
  
  
  colnames(Pro1_time)=c("User","System","Elapse")
  Pro1_timems=rbind(colMeans(Pro1_time),apply(Pro1_time,2,sd))
  colnames(Pro1_timems)=c("User","System","Elapse")
  rownames(Pro1_timems)=c("Pro1_time","sd")
  
  
  colnames(Bai_time)=c("User","System","Elapse")
  Bai_timems=rbind(colMeans(Bai_time),apply(Bai_time,2,sd))
  colnames(Bai_timems)=c("User","System","Elapse")
  rownames(Bai_timems)=c("Bai_time","sd")
  
  
  
  ################## Weighting coefficients summary ##################
  
  colnames(saveweia)=1:nuA
  sumweia=rbind(colMeans(saveweia),apply(saveweia,2,sd))
  colnames(sumweia)=1:nuA
  rownames(sumweia)=c("mean","sd")
  

  
  ########################### Output #############################
  out=list()
  out$reMSE=reMSE
  out$reTPR=reTPR
  out$reTNR=reTNR
  
  out$sumMSE=sumMSE
  out$sumTPR=sumTPR
  out$sumTNR=sumTNR
  
  out$rePRO=rePRO
  out$sumPRO=sumPRO
  
  out$FcLcMSE=FcLcMSE
  out$FcLcUE=FcLcUE
  
  out$sumFcLcMSE=sumFcLcMSE
  out$sumFcLcUE=sumFcLcUE

  out$FjLjMSE=FjLjMSE
  out$sumFjLjMSE=sumFjLjMSE

  
  ## Save weighting coefficients estimation ##    
  out$saveprea=saveprea
  out$saveweia=saveweia
  out$sumweia=sumweia
  out$saveA=saveA
  
  
  ## Save R1 and R2 ##
  out$eR1=eR1
  out$eR2=eR2
  
  out$sumeR1=sumeR1
  out$sumeR2=sumeR2
  
  ## Save Price of Risk ##
  out$ProPoR=ProPoR
  out$Pro1PoR=Pro1PoR
  out$BaiPoR=BaiPoR
  
  out$sumProfrp=sumProfrp
  out$sumPro1frp=sumPro1frp
  out$sumBaifrp=sumBaifrp
  
  out$TProfrp=TProfrp
  out$TPro1frp=TPro1frp
  out$TBaifrp=TBaifrp
  
  out$indfrp1=indfrp1
  out$indfrp2=indfrp2
  
  #########################
  
  out$Pro_time=Pro_time
  out$Pro1_time=Pro1_time
  out$Bai_time=Bai_time
  
  out$Pro_timems=Pro_timems
  out$Pro1_timems=Pro1_timems
  out$Bai_timems=Bai_timems
  
  out$saveprog=saveprog
  out$saveprobeta=saveprobeta
  out$saveproFc=saveproFc
  out$saveproLamc=saveproLamc
  out$saveproFj=saveproFj
  out$saveproLamj=saveproLamj
  out$saveproobj=saveproobj
  out$saveprogmem=saveprogmem
  
  
  out$savepro1g=savepro1g
  out$savepro1beta=savepro1beta
  out$savepro1Fc=savepro1Fc
  out$savepro1Lamc=savepro1Lamc
  out$savepro1Fj=savepro1Fj
  out$savepro1Lamj=savepro1Lamj
  out$savepro1obj=savepro1obj
  out$savepro1gmem=savepro1gmem

  
  out$saveBaig=saveBaig
  out$saveBaibeta=saveBaibeta
  out$saveBaiFc=saveBaiFc
  out$saveBaiLamc=saveBaiLamc
  out$saveBaiFj=saveBaiFj
  out$saveBaiLamj=saveBaiLamj
  out$saveBaiobj=saveBaiobj
  out$saveBaigmem=saveBaigmem
  
  
  out$saveIngmem=saveIngmem
  out$saveInbeta=saveInbeta
  out$saveInFc=saveInFc
  out$saveInLamc=saveInLamc
  out$saveInFj=saveInFj
  out$saveInLamj=saveInLamj
  
  
  
  ## Save optimal tuning parameters ##

  out$saveopS=saveopS
  out$saveopits=saveopits
  out$saveopkapN=saveopkapN
  out$saveopk=saveopk
  out$saveopkS=saveopkS
  out$saveopalp=saveopalp
  out$saveopalp1=saveopalp1
  out$savefiPIC=savefiPIC
  out$savefiPIC1=savefiPIC1
  
  out$SUCO=SUCO
  out$rUCO=rUCO
  out$rjUCO=rjUCO
  out$torjUCO=torjUCO
  

  
  return(out)
  

}
























