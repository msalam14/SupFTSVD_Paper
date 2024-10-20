iter_sim_SupFTSVD<-function(iter,n){
  EAval<-sapply(1:length(gam), function(u){
    (Vmat%*%matrix(gam[[u]],ncol=1))
  })
  
  seedN<-n+iter+p+tau2+i+m+200
  set.seed(seedN)
  aval<-sapply(1:length(gam), function(u){
    (Vmat%*%matrix(gam[[u]],ncol=1))+rnorm(n,mean = 0,sd = sqrt(eta2[u]))
  })

  # Generating observation time points
  #Tpos<-lapply(1:n, function(u){sort(sample(1:length(Time),sample(m,1)))})
  Tpos<-lapply(1:n, function(u){sort(sample(1:length(Time),sample(m:(m+5),1)))})
  Tij<-lapply(1:n, function(u){Time[Tpos[[u]]]})
  mi<-sapply(Tij,length)
  id<-lapply(1:n,function(u){rep(u,mi[u])})
  ## index for cross-validation
  kfold<-5
  tM<-sum(mi)
  int<-tM%%kfold
  if(int==0){
    Kfold<-rep(1:kfold,each=tM/kfold)
  } else{
    Kfold<-c(rep(1:kfold,each=tM%/%kfold),1:(tM%%kfold))
  }
  KInd<-sample(Kfold,tM,replace = FALSE)
  
  
  ### Storing observation points
  TimeGen<-cbind(do.call(c,Tij),do.call(c,id),iter)
  MiP<-cbind(1:n,mi,iter)
  
  # Generating the xi functions
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
  obs_dat<-lapply(lowr_comp, function(u){
    u+matrix(rnorm(prod(dim(u)),mean = 0,sd=sqrt(tau2)),nrow = nrow(u),ncol = ncol(u))
  })
  
  # True componentwise R2
  yy<-do.call(c,lapply(obs_dat,as.numeric))
  tR2<-sapply(1:lr, function(l){
    xx<-do.call(c,lapply(1:n, function(i){
      as.numeric(lmd_val[l]*EAval[i,l]*outer(bval[,l],xiT[[i]][,l]))
    }))
    tfit<-lm(yy~xx-1)
    yy<-tfit$residuals
    summary(tfit)$r.squared
  })
  
  
  
  # data structure for FTSVD
  dat_list<-lapply(1:n,FUN=function(u){
    rbind(Tij[[u]],obs_dat[[u]])
  })
  
  # Fitting unsupervised method
  res_ftsvd <- ftsvd(datlist = dat_list, interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=kfold,smooth=round(exp(seq(-5,2,length.out=25)),3),maxiter = 200,epsilon = 1e-6,KInd=KInd)


  # Fitting supervised FTSVD
  sup_res_ftsvd <- supFTSVD(datlist = dat_list,response = Vmat,interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=kfold,smooth=round(exp(seq(-5,2,length.out=25)),3),maxiter = 200,epsilon = 1e-4,KInd=KInd)
  

  # all results in a data frame
  do.call(rbind,lapply(seq_len(lr),function(u){
    data.frame("FTSVD1"=min(sqrt(sum((((res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*aval[,u]))^2)))/n,
                            sqrt(sum((((-1*res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*aval[,u]))^2)))/n),
               "FTSVD1M"=min(sqrt(sum((((res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*EAval[,u]))^2)))/n,
                             sqrt(sum((((-1*res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*EAval[,u]))^2)))/n),
               "FTSVD1R"=min(sqrt(sum((((res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*aval[,u]))^2)/((lmd_val[u]*aval[,u])^2)))/n,
                             sqrt(sum((((-1*res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*aval[,u]))^2)/((lmd_val[u]*aval[,u])^2)))/n),
               "FTSVD1MR"=min(sqrt(sum((((res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n,
                              sqrt(sum((((-1*res_ftsvd$Lambda[u]*res_ftsvd$A.hat[,u])-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n),
               "SupFTSVD1"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*aval[,u]))^2)))/n,
                               sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*aval[,u]))^2)))/n),
               "SupFTSVD1M"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*EAval[,u]))^2)))/n,
                                sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*EAval[,u]))^2)))/n),
               "SupFTSVD1CM"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*(sup_res_ftsvd$Xb[,u])-(lmd_val[u]*EAval[,u]))^2)))/n,
                               sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*(sup_res_ftsvd$Xb[,u])-(lmd_val[u]*EAval[,u]))^2)))/n),
               "SupFTSVD1R"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*aval[,u]))^2)/((lmd_val[u]*aval[,u])^2)))/n,
                                sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*aval[,u]))^2)/((lmd_val[u]*aval[,u])^2)))/n),
               "SupFTSVD1MR"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n,
                                 sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*sup_res_ftsvd$A.hat[,u]-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n),
               "SupFTSVD1CMR"=min(sqrt(sum(((sup_res_ftsvd$Lambda[u]*(sup_res_ftsvd$Xb[,u])-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n,
                                sqrt(sum(((-1*sup_res_ftsvd$Lambda[u]*(sup_res_ftsvd$Xb[,u])-(lmd_val[u]*EAval[,u]))^2)/((lmd_val[u]*EAval[,u])^2)))/n),
               "bFTSVD1"=min(norm(matrix(res_ftsvd$B.hat[,u]-bval[,u],ncol=1),type="F"),norm(matrix(-res_ftsvd$B.hat[,u]-bval[,u],ncol=1),type="F")),
               "bSupFTSVD1"=min(norm(matrix(sup_res_ftsvd$B.hat[,u]-bval[,u],ncol=1),type="F"),norm(matrix(-sup_res_ftsvd$B.hat[,u]-bval[,u],ncol=1),type="F")),
               "xiFTSVD1"=min(mean((res_ftsvd$Phi.hat[,u]-PhiF[,u])^2),mean((-res_ftsvd$Phi.hat[,u]-PhiF[,u])^2)),
               "xiSupFTSVD1"=min(mean((sup_res_ftsvd$Phi.hat[,u]-PhiF[,u])^2),mean((-sup_res_ftsvd$Phi.hat[,u]-PhiF[,u])^2)),
               "lamFTSVD1"=res_ftsvd$Lambda[u],
               "lamSupFTSVD1"=sup_res_ftsvd$Lambda[u],
               "Gamma1"=sup_res_ftsvd$Gamma[1,u],
               "Gamma2"=sup_res_ftsvd$Gamma[2,u],
               "CWR2_FTSVD1"=res_ftsvd$r.square[u],
               "CWR2_SupFTSVD1"=sup_res_ftsvd$r.square[u],
               "CWR2_TRUE"=tR2[u],
               "AR2_FTSVD1"=res_ftsvd$accum.r.square[u],
               "AR2_SupFTSVD1"=sup_res_ftsvd$accum.r.square[u],
               "sigRcomp1"=(sup_res_ftsvd$Sigma2R[u]-(lmd_val[u]^2)*eta2[u])^2,
               "sigRcomp1R"=((sup_res_ftsvd$Sigma2R[u]-(lmd_val[u]^2)*eta2[u])^2)/((lmd_val[u]^2)*eta2[u])^2,
               "Sig2B"=ifelse(u==1,sup_res_ftsvd$Sigma2-tau2,NA),
               "Sig2M"=ifelse(u==1,(sup_res_ftsvd$Sigma2-tau2)^2,NA),
               "RSig2B"=ifelse(u==1,(sup_res_ftsvd$Sigma2-tau2)/tau2,NA),
               "RSig2M"=ifelse(u==1,((sup_res_ftsvd$Sigma2-tau2)^2)/tau2^2,NA),
               "Cov2R"= summary(lm(res_ftsvd$A.hat[,u]~Vmat-1))$adj.r.squared,
               "SupCov2R"= summary(lm(sup_res_ftsvd$A.hat[,u]~Vmat-1))$adj.r.squared,
               "FTSVD1MP"=min(sqrt(sum((((res_ftsvd$Lambda[u]*as.numeric(fitted(lm(res_ftsvd$A.hat[,u]~Vmat-1))))-(lmd_val[u]*EAval[,u]))^2)))/n,
                             sqrt(sum((((-1*res_ftsvd$Lambda[u]*as.numeric(fitted(lm(res_ftsvd$A.hat[,u]~Vmat-1))))-(lmd_val[u]*EAval[,u]))^2)))/n),
               "Iter"=iter,
               "Sig2E"=eta2[u],
               "comp"=u,
               "SeedNumber"=seedN,
               "Tau"=tau2,
               "n"=n,
               "Grid"=m,
               "nFeat"=p)
  }))
}
