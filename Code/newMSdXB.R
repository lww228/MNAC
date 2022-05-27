

Ini=function(X,Y,N,T,S,r,rS,kapN){
  ############################## Initial Parameter Values ##########################
  cen=matrix(0,S,T)
  for(j in 1:S){
    inr=(j-1)*N/S+N/(2*S)
    cen[j,]=t(Y)[inr,]
  }
  Ingmem=kmeans(t(Y),cen)$cluster
  
  
  Inbeta=vector("list",N)
  for(i in 1:N){
    Inbeta[[i]]=solve(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%Y[,i]
  }
  
  
  InZc=matrix(0,T,N)
  for(i in 1:N){
    InZc[,i]=Y[,i]-X[[i]]%*%Inbeta[[i]]
  }
  InFc=sqrt(T)*eigen(InZc%*%t(InZc))$vectors[,1:r]
  InLamc=t(InZc)%*%InFc/T
  
  
  InFj=vector("list",S)
  InLamj=vector("list",S)
  
  for(j in 1:S){
    IniIndex=which(Ingmem==j)
    InNj=length(IniIndex)
    InZj=matrix(0,T,InNj)
    ll=1
    for(l in IniIndex){
      InZj[,ll]=Y[,l]-X[[l]]%*%Inbeta[[l]]-InFc%*%as.matrix(InLamc[l,])
      ll=ll+1
    }
    InFj[[j]]=sqrt(T)*eigen(InZj%*%t(InZj))$vectors[,1:rS[j]]
    InLamj[[j]]=t(InZj)%*%InFj[[j]]/T
  }
  
  

  out=list()
  out$Ingmem=Ingmem
  out$Inbeta=Inbeta
  out$InFc=InFc
  out$InLamc=InLamc
  out$InFj=InFj
  out$InLamj=InLamj
  
  return(out)
  
  ################# End for Initial Parameter Values ##########################
}




###################### Model Selection ########################

PICfun=function(Y,X,beta,gmem,Fc,Lamc,Fj,T,N,C,sighat2,S,k,kS){
  PIC=0
  
  for(i in 1:N){
    gi=gmem[i]
    
    prei=Y[,i]-X[[i]]%*%beta[[i]]-Fc%*%as.matrix(Lamc[i,])
    Lgii=t(Fj[[gi]])%*%prei/T
    pre_eps1=prei-Fj[[gi]]%*%Lgii
    eps1=t(pre_eps1)%*%pre_eps1
    
    PIC=PIC+eps1
  }
  
  PIC=PIC/(N*T)
  
  com1=PIC
  com2=C*k*sighat2*(T+N)*log(T*N)/(T*N)
  com3=0

  
  PIC=PIC+com2
  
  for(j in 1:S){
    Nj=sum(gmem==j)
    
    if(Nj>0){
      preca=C*kS[j]*sighat2*(T+Nj)*log(T*Nj)/(T*Nj)
      com3=com3+preca
      
      PIC=PIC+preca
    }
    
  }
  
  out=list()
  out$com=c(com1,com2,com3,PIC)
  out$PIC=PIC
  return(out)


}


################# Calculate Vc2 #######################

Vc2=function(Smax,sar,sarj){
  prer=(sar-mean(sar))^2
  vc2=mean(prer)
  
  for(j in 1:Smax){
    prerj=(sarj[,j]-mean(sarj[,j]))^2
    vc2=vc2+mean(prerj)
  }
  
  return(vc2)
}


OCfun=function(x){
  C=canC[x]
  
  saS=rep(0,lena)
  sar=rep(0,lena)
  sarj=matrix(0,lena,Smax)
  sakapN=matrix(0,lena,N)
  saits=rep(0,lena)
  saPICcom=vector("list",lena)
  
  for(aa in 1:lena){
    NN=Na[aa]
    TT=Ta[aa]
    
    XX=saXX[[aa]]
    YY=saYY[[aa]]

    
    #start_time=proc.time()
    cS=length(canS)
    comPIC=rep(0,cS)
    difPIC=vector("list",cS)
    
    PICcom=matrix(0,cS,4)
    
    
    for(gs in 1:cS){
      ss=canS[gs]
      
      #start_time=proc.time()
      BaiSPIC=BaifixS(XX,YY,ss,TT,NN,cankapi,cank,cankj,
                      gamma,max_step,show_step,tol0,C,Baisighat2,tolMS,maxMS)
      
      #end_time=proc.time()-start_time
      #print(end_time)
      
      its=length(BaiSPIC$savePIC)
      comPIC[gs]=BaiSPIC$savePIC[its]
      difPIC[[gs]]=BaiSPIC
      
      PICcom[gs,]=tail(BaiSPIC$savePICcom,1)
      
    }
    
    optind=which.min(comPIC)
    opS=canS[optind]
    opPIC=difPIC[[optind]]
    
    opits=length(opPIC$savek)
    opkapN=opPIC$savekapN[opits,]
    opk=opPIC$savek[opits]
    opkS=opPIC$savekS[opits,]
    
    ###########################################################################
    
    saS[aa]=opS
    sar[aa]=opk
    sarj[aa,1:opS]=opkS
    sakapN[aa,1:NN]=opkapN
    saits[aa]=opits
    
    saPICcom[[aa]]=PICcom
  }
  
  rVc2=Vc2(Smax,sar,sarj)
  tunl=list()
  tunl$saS=saS
  tunl$sar=sar
  tunl$sarj=sarj
  tunl$sakapN=sakapN
  tunl$saits=saits
  tunl$saPICcom=saPICcom
  
  out=list()
  out$rVc2=rVc2
  out$tunl=tunl
  return(out)
}



BPfun=function(x){
  ss=canS[x]
  
  #start_time=proc.time()
  BaiSPIC=BaifixS(XX,YY,ss,TT,NN,cankapi,cank,cankj,
                  gamma,max_step,show_step,tol0,C,Baisighat2,tolMS,maxMS)
  
  #end_time=proc.time()-start_time
  #print(end_time)
  
  its=length(BaiSPIC$savePIC)
  cPIC=BaiSPIC$savePIC[its]
  dPIC=BaiSPIC
  tPIC=tail(BaiSPIC$savePICcom,1)
  
  out=list()
  out$cPIC=cPIC
  out$dPIC=dPIC
  out$tPIC=tPIC
  return(out)
}





##################################### Add Baisighat2 #######################################



Baiscasig2=function(X,Y,canS,T,N,cankapi,cank,cankj,gamma,max_step,show_step,tol0,NS){
  
  upS=max(canS)
  upkapN=rep(min(cankapi),N)
  upr=max(cank)
  uprS=rep(max(cankj),upS)
  
  
  ############################## Initial Parameter Values ##########################
  Init=Ini(X,Y,N,T,upS,upr,uprS,upkapN)
  Ingmem=Init$Ingmem
  Inbeta=Init$Inbeta
  InFc=Init$InFc
  InLamc=Init$InLamc
  InFj=Init$InFj
  InLamj=Init$InLamj
  
  ################# End for Initial Parameter Values ##########################
  
  
  sigres=BaiOutloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,upS,T,N,upkapN,upr,uprS,gamma,max_step,show_step,tol0)
  ## Proposed ##
  signewg=sigres$newg
  signewbeta=sigres$newbeta
  signewFc=sigres$newFc
  signewLamc=sigres$newLamc
  signewFj=sigres$newFj
  signewLamj=sigres$newLamj
  
  sighat2=0
  
  for(i in 1:N){
    gi=signewg[i]
    
    prei=Y[,i]-X[[i]]%*%signewbeta[[i]]-signewFc%*%as.matrix(signewLamc[i,])
    Lgii=t(signewFj[[gi]])%*%prei/T
    pre_eps1=prei-signewFj[[gi]]%*%Lgii
    eps1=t(pre_eps1)%*%pre_eps1
    
    sighat2=sighat2+eps1
  }
  sighat2=sighat2/(N*T)
  
  return(sighat2)
}



######################################## Update Bai kapi #######################################

Baiupdatekapi=function(X,Y,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0,
                       C,sighat2,i,cankapi){
  
  lkapi=length(cankapi)
  rekapi=rep(0,lkapi)
  
  reloop=vector("list",lkapi)
  for(ll in 1:lkapi){
    InkapN[i]=cankapi[ll]
    
    ############################## Initial Parameter Values ##########################
    Init=Ini(X,Y,N,T,InS,Ink,InkS,InkapN)
    Ingmem=Init$Ingmem
    Inbeta=Init$Inbeta
    InFc=Init$InFc
    InLamc=Init$InLamc
    InFj=Init$InFj
    InLamj=Init$InLamj
    
    ################# End for Initial Parameter Values ##########################
    
    
    upres=BaiOutloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0)
    reloop[[ll]]=upres
    
    upnewg=upres$newg
    upnewbeta=upres$newbeta
    upnewFc=upres$newFc
    upnewLamc=upres$newLamc
    upnewFj=upres$newFj
    upnewLamj=upres$newLamj
    upobj_fun=upres$obj_fun
    upiter_gmem=upres$iter_gmem
    
    rekapi[ll]=PICfun(Y,X,upnewbeta,upnewg,upnewFc,upnewLamc,upnewFj,T,N,C,sighat2,InS,Ink,InkS)
  }
  
  out=list()
  out$PIC=rekapi
  out$minPIC=min(rekapi)
  
  opmin=which.min(rekapi)
  out$opt=cankapi[opmin]
  out$loop=reloop[[opmin]]
  
  return(out)
  
}

######################################## End for Update Bai kapi #######################################


######################################## Bai Update k ##################################################


Baiupdatek=function(X,Y,InS,T,N,InkapN,InkS,gamma,max_step,show_step,tol0,
                    C,sighat2,cank){
  
  lk=length(cank)
  rek=rep(0,lk)
  
  reloop=vector("list",lk)
  
  rePICcom=matrix(0,lk,4)
  
  for(ll in 1:lk){
    Ink=cank[ll]
    
    ############################## Initial Parameter Values ##########################
    Init=Ini(X,Y,N,T,InS,Ink,InkS,InkapN)
    Ingmem=Init$Ingmem
    Inbeta=Init$Inbeta
    InFc=Init$InFc
    InLamc=Init$InLamc
    InFj=Init$InFj
    InLamj=Init$InLamj
    
    ################# End for Initial Parameter Values ##########################
    
    upres=BaiOutloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0)
    reloop[[ll]]=upres
    
    upnewg=upres$newg
    upnewbeta=upres$newbeta
    upnewFc=upres$newFc
    upnewLamc=upres$newLamc
    upnewFj=upres$newFj
    upnewLamj=upres$newLamj
    upobj_fun=upres$obj_fun
    upiter_gmem=upres$iter_gmem
    
    prePIC=PICfun(Y,X,upnewbeta,upnewg,upnewFc,upnewLamc,upnewFj,T,N,C,sighat2,InS,Ink,InkS)
    rek[ll]=prePIC$PIC
    rePICcom[ll,]=prePIC$com
  }
  
  out=list()
  out$PIC=rek
  out$minPIC=min(rek)
  
  opmin=which.min(rek)
  out$opt=cank[opmin]
  out$loop=reloop[[opmin]]
  
  out$rePICcom=rePICcom
  out$opPICcom=rePICcom[opmin,]
  
  return(out)
  
}

######################################## End for Bai Update k ################################################


############################################ Bai Update kj ##################################################

Baiupdatekj=function(X,Y,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0,
                     C,sighat2,j,cankj){
  
  lkj=length(cankj)
  rekj=rep(0,lkj)
  
  reloop=vector("list",lkj)
  
  rePICcom=matrix(0,lkj,4)
  
  for(ll in 1:lkj){
    InkS[j]=cankj[ll]
    
    ############################## Initial Parameter Values ##########################
    Init=Ini(X,Y,N,T,InS,Ink,InkS,InkapN)
    Ingmem=Init$Ingmem
    Inbeta=Init$Inbeta
    InFc=Init$InFc
    InLamc=Init$InLamc
    InFj=Init$InFj
    InLamj=Init$InLamj
    
    ################# End for Initial Parameter Values ##########################
    
    upres=BaiOutloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0)
    reloop[[ll]]=upres
    
    upnewg=upres$newg
    upnewbeta=upres$newbeta
    upnewFc=upres$newFc
    upnewLamc=upres$newLamc
    upnewFj=upres$newFj
    upnewLamj=upres$newLamj
    upobj_fun=upres$obj_fun
    upiter_gmem=upres$iter_gmem
    
    prePIC=PICfun(Y,X,upnewbeta,upnewg,upnewFc,upnewLamc,upnewFj,T,N,C,sighat2,InS,Ink,InkS)
    rekj[ll]=prePIC$PIC
    rePICcom[ll,]=prePIC$com
  }
  
  out=list()
  out$PIC=rekj
  out$minPIC=min(rekj)
  
  opmin=which.min(rekj)
  out$opt=cankj[opmin]
  out$loop=reloop[[opmin]]
  
  out$rePICcom=rePICcom
  out$opPICcom=rePICcom[opmin,]
  
  return(out)
  
}

############################################ End for Bai Update kj ##############################################

BaifixS=function(X,Y,s,T,N,cankapi,cank,cankj,
                 gamma,max_step,show_step,tol0,C,sighat2,tolMS,maxMS){
  InS=s
  
  savekapN=NULL
  savek=c()
  savekS=NULL
  savePIC=c()
  
  savePICcom=NULL
  
  InkapN=rep(cankapi[1],N)
  InkapN0=InkapN
  Ink=cank[2]
  InkS=rep(cankj[2],InS)
  InkS0=InkS
  
  k=1
  
  
  newkapN=rep(0,N)

  newkS=rep(0,InS)
  for(j in 1:InS){
    newkS[j]=Baiupdatekj(X,Y,InS,T,N,newkapN,Ink,InkS,
                         gamma,max_step,show_step,tol0,C,sighat2,j,cankj)$opt
    InkS[j]=newkS[j]
  }
  
  upk=Baiupdatek(X,Y,InS,T,N,newkapN,newkS,
                 gamma,max_step,show_step,tol0,C,sighat2,cank)
  
  newk=upk$opt
  newPIC=upk$minPIC
  
  newPICcom=upk$opPICcom
  savePICcom=rbind(savePICcom,newPICcom)
  
  savekapN=rbind(savekapN,newkapN)
  savek=c(savek,newk)
  savekS=rbind(savekS,newkS)
  savePIC=c(savePIC,newPIC)
  
  tol=sum(abs(newkapN-InkapN0))+abs(newk-Ink)+sum(abs(newkS-InkS0))
  
  while((tol>tolMS)&(k<maxMS)){
    k=k+1
    
    InkapN=newkapN
    InkapN0=InkapN
    Ink=newk
    InkS=newkS
    InkS0=InkS
    
    newkapN=rep(0,N)

    newkS=rep(0,InS)
    for(j in 1:InS){
      newkS[j]=Baiupdatekj(X,Y,InS,T,N,newkapN,Ink,InkS,
                           gamma,max_step,show_step,tol0,C,sighat2,j,cankj)$opt
      InkS[j]=newkS[j]
    }
    
    upk=Baiupdatek(X,Y,InS,T,N,newkapN,newkS,
                   gamma,max_step,show_step,tol0,C,sighat2,cank)
    
    newk=upk$opt
    newPIC=upk$minPIC
    
    newPICcom=upk$opPICcom
    savePICcom=rbind(savePICcom,newPICcom)
    
    savekapN=rbind(savekapN,newkapN)
    savek=c(savek,newk)
    savekS=rbind(savekS,newkS)
    savePIC=c(savePIC,newPIC)
    
    tol=sum(abs(newkapN-InkapN0))+abs(newk-Ink)+sum(abs(newkS-InkS0))
  }
  
  out=list()
  out$savekapN=savekapN
  out$savek=savek
  out$savekS=savekS
  out$savePIC=savePIC
  
  out$savePICcom=savePICcom
  
  return(out)
  
}






############################################## Update alpha ################################################

updatealp=function(X,Y,A,A2,Am12,InS,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0,
                   C,sighat2,canalp){
  
  ############################## Initial Parameter Values ##########################
  Init=Ini(X,Y,N,T,InS,Ink,InkS,InkapN)
  Ingmem=Init$Ingmem
  Inbeta=Init$Inbeta
  InFc=Init$InFc
  InLamc=Init$InLamc
  InFj=Init$InFj
  InLamj=Init$InLamj
  
  ################# End for Initial Parameter Values ##########################
  
  
  lalp=length(canalp)
  realp=rep(0,lalp)
  
  reloop=vector("list",lalp)
  for(ll in 1:lalp){
    Inalpha=canalp[ll]
    upres=Outloop(Ingmem,Inbeta,InFc,InLamc,InFj,InLamj,X,Y,A,A2,Am12,InS,Inalpha,T,N,InkapN,Ink,InkS,gamma,max_step,show_step,tol0)
    reloop[[ll]]=upres
    
    upnewg=upres$newg
    upnewbeta=upres$newbeta
    upnewFc=upres$newFc
    upnewLamc=upres$newLamc
    upnewFj=upres$newFj
    upnewLamj=upres$newLamj
    upobj_fun=upres$obj_fun
    upiter_gmem=upres$iter_gmem
    
    realp[ll]=PICfun(Y,X,upnewbeta,upnewg,upnewFc,upnewLamc,upnewFj,T,N,C,sighat2,InS,Ink,InkS)$PIC
  }
  
  out=list()
  out$PIC=realp
  out$minPIC=min(realp)
  
  opmin=which.min(realp)
  out$opt=canalp[opmin]
  out$loop=reloop[[opmin]]
  
  return(out)
  
}

########################################## End for Update alpha ##############################################








































