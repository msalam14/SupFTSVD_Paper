rm(list=ls())
# Parallel computing
library(parallel)
library(snow)
library(tidyverse)
library(dplyr)
library(rsvd)
library(matrixcalc)
library(doParallel)
library(doSNOW)
library(foreach)
library(doParallel)
library(doRNG)
#library(doFuture)

# Required function
source(paste(getwd(),'/ftsvd-main/FTSVD.R',sep = ""))

#Simulation settings
SS<-c(20,50,100,200)
pdim<-c(500)
Tau2<-c(1,4) # tensor error
Eta2<-list(c(1),c(2),c(5)) # subject loading error variance

#Model components
lr<-1
r<-lr

# True parameters
lmd_val<-c(80) 
sig<-0 # weight for the reminder
gam<-list(c(1.5,3))#partial slopes for subject loading model


# Grid of Time points
nres<-51
Time<-seq(0,1,length.out=nres)



# Grid scenario

Mi<-c(3)

# Cluster formation
#cls<-makeCluster( detectCores()-1)
cls<-parallel::makeCluster(5)
registerDoParallel(cl=cls)
#cls<-makeCluster( (mpi.universe.size()-1) , type='MPI' )
#registerDoSNOW(cls)
#registerDoParallel(cls)



dirP<-"/share/astaicu/malam3/SupFTSVD/"


source("SupFTSVD_SimRCodes/SupFTSVD_SimRCodes_Feb15_24/supFTSVD_iter.R")
NSim<-10
sitn<-expand.grid(SS,pdim,Mi,Eta2,Tau2)
k<-1
sam_s<<-SS[k]
p<-pdim[1]
tau2<<-Tau2[1]
eta2<<-Eta2[[1]]
m<-Mi[1]          
kfold<-5
n<-sam_s

set.seed(sam_s,kind="L'Ecuyer-CMRG")
Vmat<<-cbind(round(runif(sam_s),2),round(rbeta(sam_s,1,1),2))
sbjMEM<-lapply(1:NSim,function(iter){
  sapply(1:r,function(k){as.matrix(rnorm(n,mean = 0,sd = sqrt(eta2[k])))})
})
# Generating b_j
Bval<<-sapply(1:r, function(b){runif(p)})
bval<<-Bval*outer(rep(1,p),1/apply(Bval,2,norm,type="2"))
# Singular Function
PhiFunc<-lapply(1:r, function(x){
  function(s){
    a<-t(sapply(s,function(j){
      c(1,sapply(2:10, function(i){sqrt(2)*cos((i-1)*(pi)*j)}))
    }))%*%matrix(sapply(1:10,function(j){runif(1,-1/j,1/j)}))
    (a/norm(a,type="2"))*sqrt(length(s))
  }
})
PhiF<-sapply(1:r,function(r){PhiFunc[[r]](Time)})
##
iterT<-lapply(1:NSim,function(iter){
  lapply(1:n, function(i){sort(sample(1:length(Time),sample(m:(m+5),1)))})
})
tenERR<-lapply(1:NSim,function(iter){
  lapply(1:n, function(i){
    matrix(rnorm(length(iterT[[iter]][[i]])*p,mean = 0,sd = sqrt(tau2)),nrow=p)
  })
})

iterKF<-lapply(1:NSim,function(iter){
  tM<-length(do.call(c,iterT[[iter]]))
  int<-tM%%kfold
  if(int==0){
    Kfold<-rep(1:kfold,each=tM/kfold)
  } else{
    Kfold<-c(rep(1:kfold,each=tM%/%kfold),1:(tM%%kfold))
  }
  KInd<-sample(Kfold,tM,replace = FALSE)
})

## Subject loading expected value
EAval<-as.matrix(sapply(1:r, function(k){
  (Vmat%*%matrix(gam[[k]],ncol=1))
}))

iter_aval<-lapply(1:NSim, function(iter){
  EAval+sbjMEM[[iter]]
})

#Tpos<-lapply(1:n, function(u){sort(sample(1:length(Time),m))})
#Tpos<<-lapply(1:sam_s, function(u){sort(sample(1:length(Time),sample(m:(m+5),1)))})
###
i<-1

# Removed from iteration function
#seedN<-n+iter+p+tau2+i+m+200
#set.seed(seedN)
n<-sam_s
#set.seed(45)
sim_obj<-lapply(1:NSim, function(iter){
  aval<-EAval+sbjMEM[[iter]]
  Tpos<-iterT[[iter]]
  Tij<-lapply(1:n, function(u){Time[Tpos[[u]]]})
  mi<-sapply(Tij,length)
  id<-lapply(1:n,function(u){rep(u,mi[u])})
  xiT<-lapply(Tpos,FUN=function(u){
    as.matrix(PhiF[u,])
  })
  # Generating the smooth tensor
  lowr_comp<-lapply(1:n, function(u){
    Reduce(`+`,lapply(1:r, function(v){
      if(v<=lr){
        lmd_val[v]*aval[u,v]*(outer(as.numeric(bval[,v]),xiT[[u]][,v]))
      } else{
        sig*lmd_val[v]*aval[u,v]*(outer(as.numeric(bval[,v]),xiT[[u]][,v]))
      }
    }))
  })
  # Functional tensor with noise
  obs_dat<-lapply(1:n, function(i){
    lowr_comp[[i]]+tenERR[[iter]][[i]]
  })
  
  # data structure for FTSVD
  dat_list<-lapply(1:n,FUN=function(u){
    rbind(Tij[[u]],obs_dat[[u]])
  })
  
  dat_list
})


############


#set.seed(23)
pt<-proc.time()
#plan(cluster,workers=5)
par_res<-foreach(j=1:NSim,
                 .combine = 'rbind',.packages = c("tidyverse",
                                                  "dplyr",
                                                  "rsvd",
                                                  "matrixcalc"),
                 .export = c("ftsvd",
                             "SUPftsvd",
                             "supFTSVD",
                             "bernoulli_kernel",
                             "freg_rkhs",
                             "cv_freg_rkhs",
                             "lr",
                             "r",
                             "lmd_val",
                             "sig",
                             "gam",
                             "PhiF",
                             "Time",
                             "nres",
                             "bval",
                             "i",
                             "tau2",
                             "eta2",
                             "sam_s",
                             "p",
                             "m",
                             "Vmat",
                             "sbjMEM",
                             "iterT",
                             "tenERR",
                             "iterKF",
                             "sim_obj",
                             "iter_aval",
                             "iter_sim_SupFTSVD")) %dopar%
  iter_sim_SupFTSVD(iter=j,n=sam_s)
pt1<-proc.time()-pt



par_res


#write.csv(Res2,file=paste(dirP,"SupFTSVD",pdim,"_DSCN1_3.csv",sep=""),row.names = FALSE)?
stopCluster(cls)
#mpi.exit()




