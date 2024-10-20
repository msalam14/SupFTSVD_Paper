iter_sim_SupFTSVD<-function(iter,n){
  EAval<-sapply(1:length(gam), function(u){
    (Vmat%*%matrix(gam[[u]],ncol=1))
  })
  
  seedN<-n+iter+p+tau2+i+m+100
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
  
  # Test data generation
  tVmat<-cbind(round(runif(100),2),round(rbeta(100,1,1),2))
  tEAval<-sapply(1:length(gam), function(u){
    (tVmat%*%matrix(gam[[u]],ncol=1))
  })
  
  taval<-sapply(1:length(gam), function(u){
    (tVmat%*%matrix(gam[[u]],ncol=1))+rnorm(100,mean = 0,sd = sqrt(eta2[u]))
  })
  
  testY<-lapply(1:100, function(u){
    Reduce(`+`,lapply(1:r, function(v){
      if(v<=lr){
        lmd_val[v]*taval[u,v]*(outer(as.numeric(bval[,v]),PhiF[,v]))
      } else{
        sig*lmd_val[v]*taval[u,v]*(outer(as.numeric(bval[,v]),PhiF[,v]))
      }
    }))+matrix(rnorm((p*nrow(PhiF)),mean = 0,sd=sqrt(tau2)),nrow = p,ncol = nrow(PhiF))
  })
  
  true_testY<-lapply(1:100, function(u){
    Reduce(`+`,lapply(1:r, function(v){
      if(v<=lr){
        lmd_val[v]*taval[u,v]*(outer(as.numeric(bval[,v]),PhiF[,v]))
      } else{
        sig*lmd_val[v]*taval[u,v]*(outer(as.numeric(bval[,v]),PhiF[,v]))
      }
    }))
  })
  
  testTpos<-lapply(1:100, function(u){sort(sample(1:length(Time),sample(m:(m+5),1)))})
  testTij<-lapply(1:100, function(u){Time[testTpos[[u]]]})
  
  test_dat<-lapply(1:100, function(i){
    testY[[i]][,testTpos[[i]]]
  })
  
  
  # Fitting unsupervised method
  res_ftsvd <- ftsvd(datlist = dat_list, interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=kfold,smooth=round(exp(seq(-5,2,length.out=25)),3),maxiter = 200,epsilon = 1e-6,KInd=KInd)
  
  # Estimation and prediction error
  ftsvd_supF<-predict.ftsvd(res_ftsvd,new_dat = test_dat,newT = testTij)
  
  ftsvdPERR<-sqrt(colMeans(t(sapply(1:length(testY), function(i){
    c(mean(rowMeans((testY[[i]]-ftsvd_supF[[i]])^2)),
      mean(rowMeans((true_testY[[i]]-ftsvd_supF[[i]])^2)))
  }))))

  # Fitting supervised FTSVD
  sup_res_ftsvd <- supFTSVD(datlist = dat_list,response = Vmat,interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=kfold,smooth=round(exp(seq(-5,2,length.out=25)),3),maxiter = 200,epsilon = 1e-4,KInd=KInd)
  
  # Estimation and prediction error
  pred_supF<-predict.supFTSVD(sup_res_ftsvd,designM = tVmat,new_dat = test_dat,newT = testTij)

  perr<-sqrt(colMeans(t(sapply(1:length(testY), function(i){
    c(mean(rowMeans((testY[[i]]-pred_supF$subTRJ[[i]])^2)),
    mean(rowMeans((testY[[i]]-pred_supF$meanTRJ[[i]])^2)),
    mean(rowMeans((true_testY[[i]]-pred_supF$subTRJ[[i]])^2)),
    mean(rowMeans((true_testY[[i]]-pred_supF$meanTRJ[[i]])^2)))
  }))))
  data.frame("subTR"=perr[1],
             "meanTR"=perr[2],
             "TsubTR"=perr[3],
             "TmeanTR"=perr[4],
             "ftsvdTR"=ftsvdPERR[1],
             "TftsvdTR"=ftsvdPERR[2],
             "Iter"=iter,
             "Sig2E1"=eta2[1],
             "Sig2E2"=eta2[2],             
             "SeedNumber"=seedN,
             "Tau"=tau2,
             "n"=n,
             "Grid"=m,
             "nFeat"=p)
}
