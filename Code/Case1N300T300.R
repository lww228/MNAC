library(MASS)
library(ncvreg)
library(gtools)


setwd('')



source("./Simulation.R")
source("./Iter.R")
source("./newupdate.R")
source("./newMSdXB.R")

#tmp = tempfile(fileext = ".out")
#Rprof(tmp) 


DGP=1
n=50
T=300
S=3

N1=100; N2=100; N3=100
N=N1+N2+N3
NS=c(N1,N2,N3)


G=c(rep(1,N1),rep(2,N2),rep(3,N3))
Gperm=permutations(S,S,v=1:S)
nG=nrow(Gperm)
Gcom=matrix(0,nG,N)
for(gg in 1:nG){
  a=NULL
  for(j in 1:S){
    ele=Gperm[gg,j]
    a=c(a,rep(ele,NS[ele]))
  }
  Gcom[gg,]=a
}



G1=which(G==1)
G2=which(G==2)
G3=which(G==3)
Group=list(G1,G2,G3)


r=2
r1=3; r2=3; r3=3
rS=c(r1,r2,r3)


p=rep(80,N)
q=rep(3,N)


pre_beta=vector("list",S)
pre_beta[[1]]=c(3,3.5,3)
pre_beta[[2]]=c(-2,-2,-2.5)
pre_beta[[3]]=c(1,0.5,1.5)



gamma=3.7
max_step=1000
show_step=5
tol0=1e-3

A0=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(G[i]==G[j]){A0[i,j]=A0[j,i]=1}
  }
}


## A080

GA=c(rep(1,N1-N*0.1),rep(2,N2+N*0.2),rep(3,N3-N*0.1))

A080=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A080[i,j]=A080[j,i]=1}
  }
}


## A060

GA=c(rep(1,N1-N*0.2),rep(2,N2+N*0.4),rep(3,N3-N*0.2))

A060=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A060[i,j]=A060[j,i]=1}
  }
}



## A052

GA=c(rep(1,N1-N*0.24),rep(2,N2+N*0.48),rep(3,N3-N*0.24))

A052=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A052[i,j]=A052[j,i]=1}
  }
}



## A040

GA=c(rep(1,N1-N*0.3),rep(2,N2+N*0.6),rep(3,N3-N*0.3))

A040=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A040[i,j]=A040[j,i]=1}
  }
}



## A0ran

GA=sample(c(1,2,3),N,replace = T)

A0ran=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A0ran[i,j]=A0ran[j,i]=1}
  }
}


## A0ran1

GA=sample(c(1,2,3),N,replace = T)

A0ran1=matrix(0,N,N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(GA[i]==GA[j]){A0ran1[i,j]=A0ran1[j,i]=1}
  }
}


################

zA=matrix(0,N,N)
zA2=zA^2
zAm12=(1-zA)^2


## observed networks ##

ep=1e-7
nuA=2
obsA=vector("list",nuA)
obsA[[1]]=A0
obsA[[2]]=A0ran

############################## Candidate Parameter Values ##########################

canalp=c(0.01,0.03,0.5,1,2,4)


C=2
tolMS=0
maxMS=10


######## Ando and Bai (2017) ###########

#cana=c(0,1:5)
#canNT=cbind(N-cana*10,T-cana*10)

#canC=0.1*c(0,1:20)

#kapk=c(0,1:12)
#cankapi=10^{-3+0.25*kapk}

#canS=2:4

#cank=c(0,1:8)
#cankj=c(0,1:8)

#canalp=2:5

############################## End for Candidate Parameter Values ##########################


start_time=proc.time()
com=comparison(DGP,n,T,N,NS,r,rS,p,q,pre_beta,A,A2,Am12,S,gamma,max_step,show_step,tol0,
               Gperm,nG,Gcom,G,cankapi,canS,cank,cankj,canalp,C,tolMS,maxMS,ep,nuA,obsA,zA,zA2,zAm12)
end_time=proc.time()-start_time
print(end_time)



#Rprof(NULL)
#summaryRprof(tmp)



save(com,file='newRdata/G3nuA2/Case1N300T300.RData')

























