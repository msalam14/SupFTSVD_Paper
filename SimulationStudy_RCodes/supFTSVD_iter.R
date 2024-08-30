iter_sim_SupFTSVD<-function(iter,n){
  # Fitting unsupervised method
  res_ftsvd <- ftsvd(datlist = sim_obj[[iter]], interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=5,smooth=round(exp(seq(-5,5,length.out=25)),3),maxiter = 200,epsilon = 1e-6,KInd=iterKF[[iter]])
  
  
  # Fitting supervised FTSVD
  sup_res_ftsvd <- supFTSVD(datlist = sim_obj[[iter]],response = Vmat,interval = c(Time[1],Time[length(Time)]), r = lr, resolution = nres,CVPhi = TRUE,K=5,smooth=round(exp(seq(-5,5,length.out=25)),3),maxiter = 200,epsilon = 1e-4,KInd=iterKF[[iter]])
  
  
  # all results in a data frame
  aval<-iter_aval[[iter]]
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
               "Tau"=tau2,
               "n"=n,
               "Grid"=m,
               "nFeat"=p)
  }))
}
