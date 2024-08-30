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
library(Rmpi)



# Required function
source(paste(getwd(),'/ftsvd-main/FTSVD.R',sep = ""))

# feature dimension 

pdim<-c(500)

# No of rank-1 components
lr<-1
r<-lr

# Basis dimension
lmd_val<-c(80) #lmd_min*(r:1)

# weight for reminder
sig<-0

# noise variance in tensor
Tau2<-c(1,4)

# Parameters for supervised component
Eta2<-list(c(1),c(2),c(5)) # error variance
gam<-list(c(0,0))#,c(2,3.4)) # partial slopes


# Grid of Time points
nres<-51
Time<-seq(0,1,length.out=nres)
#Pht<-PhiF[[1]](Tg)
set.seed(nres)
PhiFunc<-lapply(1:r, function(x){
  function(s){
    a<-t(sapply(s,function(j){
      c(1,sapply(2:10, function(i){sqrt(2)*cos((i-1)*(pi)*j)}))
    }))%*%matrix(sapply(1:10,function(j){runif(1,-1/j,1/j)}))
    a/norm(a,type="2")
  }
})
PhiF<-sapply(1:r,function(r){sqrt(length(Time))*PhiFunc[[r]](Time)})


# Grid scenario

Mi<-c(3,5,8,10)

# Cluster formation
#cls<-makeCluster( detectCores()-1)
#cls<-parallel::makeCluster(25)
#registerDoParallel(cl=cls)
cls<-makeCluster( (mpi.universe.size()-1) , type='MPI' )
registerDoSNOW(cls)
registerDoParallel(cls)



dirP<-"/share/astaicu/malam3/SupFTSVD/"


source("IterSimSupFTSVD1GDPE.R")
NSim<-100
SS<-c(100)

Res2<-do.call(rbind,lapply(seq_len(length(SS)),function(k){
  sam_s<<-SS[k]
  set.seed(sam_s)
  Vmat<<-cbind(round(runif(sam_s),2),round(rbeta(sam_s,1,1),2))
  do.call(rbind, lapply(pdim,function(p){
    # Generating b_j
    set.seed(p)
    Bval<<-sapply(1:r, function(b){runif(p)})
    bval<<-Bval*outer(rep(1,p),1/apply(Bval,2,norm,type="2"))
    
    do.call(rbind,lapply(Tau2,function(g){
      tau2<<-g
      do.call(rbind,lapply(seq_len(length(Eta2)), function(i){
        eta2<<-Eta2[[i]]
        do.call(rbind,lapply(Mi,function(m){
          
          #Tpos<-lapply(1:n, function(u){sort(sample(1:length(Time),m))})
          #Tpos<<-lapply(1:sam_s, function(u){sort(sample(1:length(Time),sample(m:(m+5),1)))})
          ###
          par_res<-foreach(j=1:NSim,
                           .combine = 'rbind',
                           .packages = c("tidyverse",
                                         "dplyr",
                                         "rsvd",
                                         "matrixcalc"),
                           .export = c("ftsvd",
                                       "SUPftsvd",
                                       "supFTSVD",
                                       "predict.supFTSVD",
                                       "predict.ftsvd",
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
                                       "iter_sim_SupFTSVD")) %dopar%
            iter_sim_SupFTSVD(iter=j,n=sam_s)
          
          print(paste("Simulation is done for", "(", sam_s,",",p,",",m,")", "with", eta2, "when tensor noise is", tau2))
          
          par_res
        }))
      }))
    }))
  }))
})
)

write.csv(Res2,file=paste(dirP,"SupFTSVD",pdim,"USPE_DGDSCN1_3.csv",sep=""),row.names = FALSE)

stopCluster(cls)
mpi.exit()




