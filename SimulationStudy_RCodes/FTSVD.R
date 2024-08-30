library(fdapace)
library(rTensor)
library(KernSmooth)
library(locfit)
library(vegan)
library(pROC)


#' Temporal Functional SVD
#' @param datlist A list of n matrices data, each matrix represents a subject;
#' the first row represents the sampling timestamp; 
#' following rows represent the feature measurements.
#' @param r CP-rank/Number of principle components. Default: 3.
#' @param interval The range of time points. Default: range of all observed data. 
#' @param resolution Resolution for the output eigencurves. Default: 100.
#' @param smooth Smoothing parameter for RKHS norm. Larger means smoother functions. Default: 1e-8.
#' @param maxiter Maximum number of itereation. Default: 20.
#' @param epsilon Convergence criteria for a and b. Default: 1e-4.
#' @return The estimations of the principle components in suject/feature/time;
#'       Var.prop: Explained variances.
ftsvd <- function(datlist, interval = NULL, r = 3, resolution = 251, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                  maxiter=20, epsilon=1e-4){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  
  
  # Calculate range.
  timestamps.all = NULL
  for (i in 1:n){
    timestamps.all = c(timestamps.all, datlist[[i]][1,])
  }
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res = NULL
  Lambda = rep(0, r)
  X = NULL
  y0 <- NULL
  Rsq <- accumRsq <- rep(0, r)
  
  ti <- vector(mode='list', length=n)
  for (i in 1:n){
    temp = 1 + round((resolution-1) * (datlist[[i]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    ti[[i]] <- temp
  }
  
  tipos <- vector(mode='list', length=n)
  for (i in 1:n){
    keep <- ti[[i]]>0
    tipos[[i]] <- keep
    y0 <- c(y0, as.vector(t(datlist[[i]][2:(p+1),keep])))
  }
  
  Lt = list()
  ind_vec <- NULL
  for (i in 1:n){
    Lt = c(Lt, list(datlist[[i]][1,]))
    ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
  }
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  for (s in 1:r){ 
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    print(sprintf("Calculate the %dth Component", s))
    
    # intialization of b
    data.unfold = NULL
    y <- NULL
    for (i in 1:n){
      data.unfold = cbind(data.unfold, datlist[[i]][2:(p+1),])
      y <- c(y, as.vector(t(datlist[[i]][2:(p+1),tipos[[i]]])))
    }
    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    b.hat = b.initials[,1]
    # initialization of a
    a.hat <- rep(1,n)/sqrt(n)
    
    # iteratively update a,b,phi
    t <- 0
    dif <- 1
    while(t<=(maxiter+cvT) & dif>epsilon){
      # update phi:
      Ly = list()
      for (i in 1:n){
        Ly = c(Ly, list(a.hat[i]*as.numeric(b.hat%*%datlist[[i]][2:(p+1),])))
      }
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K)
        Smv<-cvfit[[1]]
        phi.hat = cvfit[[2]]
      } else{
        phi.hat = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
      }
      
      #phi.hat = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi.hat = phi.hat / sqrt(sum(phi.hat^2))
      
      # update a:
      a.tilde <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        a.tilde[i] <- b.hat %*% datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] 
        a.tilde[i] <- a.tilde[i] / sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      a.new <- a.tilde / sqrt(sum(a.tilde^2))
      if(t>cvT)
        dif <- sum((a.hat - a.new)^2)
      a.hat <- a.new
      
      # update b:
      temp.num <- matrix(0,p,n)
      temp.denom <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        temp.num[,i] <- datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]]
        temp.denom[i] <-sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      b.tilde <- as.numeric(temp.num%*%a.hat) / as.numeric(temp.denom%*%(a.hat^2))
      b.new <- b.tilde / sqrt(sum(b.tilde^2))
      if(t>cvT)
        dif <- max(dif, sum((b.hat - b.new)^2))
      b.hat <- b.new
      
      t <- t+1
    }
    
    # calculate lambda
    x = NULL
    for (i in 1:n){
      t.temp = ti[[i]]
      t.temp <- t.temp[t.temp>0]
      x <- c(x,as.vector(t(a.hat[i]*b.hat%o%phi.hat[t.temp])))
    }
    X = cbind(X, x)
    l.fit = lm(y~x-1)
    lambda = as.numeric(l.fit$coefficients)
    A[,s] = a.hat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    Phi[,s] = t(phi.hat)
    Lambda[s] = lambda
    Rsq[s] <- summary(l.fit)$r.squared
    accumRsq[s] <- summary(lm(y0~X-1))$r.squared
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] = datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
    }
    print(paste0("Convergence reached at dif=", dif, ', iter=', t-cvT))
  }
  l.fit = lm(y0~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  
  # revise the sign of Lambda
  for (r in length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] = -Lambda[r]
      B[,r] = -B[,r]
    }
  }
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results = list("A.hat" = A, "B.hat" = B, 
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq)
  return(results)
}

#' Supervised Functional SVD
#' @param datlist A list of n matrices data, each matrix represents a subject;
#' the first row represents the sampling timestamp; 
#' following rows represent the feature measurements.
#' @param r CP-rank/Number of principle components. Default: 3.
#' @param interval The range of time points. Default: range of all observed data. 
#' @param resolution Resolution for the output eigencurves. Default: 100.
#' @param smooth Smoothing parameter for RKHS norm. Larger means smoother functions. Default: 1e-8.
#' @param maxiter Maximum number of itereation. Default: 20.
#' @param epsilon Convergence criteria for a and b. Default: 1e-4.
#' @return The estimations of the principle components in suject/feature/time;
#'       Var.prop: Explained variances.
supFTSVD <- function(datlist, response, interval = NULL, r = 3,resolution=100, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                     maxiter=20, epsilon=1e-4){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  # response in matrix form
  if(!is.matrix(response)){
    response<-matrix(response,ncol=1)
  }
  
  # storing original response
  ires<-response
  
  
  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  ti = lapply(1:n,FUN=function(u){
    temp = 1 + round((resolution-1) * (datlist[[u]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    temp
  })
  
  y0=do.call(c,lapply(1:n,FUN=function(u){
    as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
  }))
  
  y<-y0
  
  tipos<-lapply(ti,function(u){u>0})
  
  Lt<-lapply(datlist,FUN = function(u){u[1,]})
  ind_vec<-rep(1:n,sapply(Lt,length))
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  # Initial values
  ## initialization of b
  data.unfold = do.call(cbind,lapply(datlist,FUN=function(u){u[-1,]}))
  
  b.initials <- rsvd(data.unfold, k=r)$u
  b.hat = as.matrix(b.initials[,1:r])
  
  ## initialization of a
  if(r==1){
    a.hat<-as.matrix(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  } else{
    a.hat<-t(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  }
  
  
  ## initialize U
  Ahat<-a.hat
  
  # Initialize xi function
  # compress the data and apply function PCA.
  phi.hat<-sapply(1:r, function(k){
    Ly<-lapply(1:n,FUN = function(i){
      Ahat[i,k]*as.numeric(b.hat[,k]%*%datlist[[i]][2:(p+1),])
    })
    
    if(CVPhi){
      phiH = cv_freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth, kfold = K)[[2]]
    } else{
      phiH = freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth)
    }
    #phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    phiH / sqrt(sum(phiH^2))
  })
  
  t_arg<-seq(interval[1],interval[2],length.out = resolution)
  
  allT<-do.call(c,Lt)
  pXI<-apply(phi.hat,2,function(u){
    approx(t_arg,u,xout = allT,rule=2)$y
  })
  
  
  # initial values for error variance
  mi<-as.numeric(sapply(Lt,length))
  Mi<-cumsum(mi)
  
  prdXI<-lapply(1:n,function(w){
    if(w ==1){
      as.matrix(pXI[1:Mi[1],])
    } else{
      as.matrix(pXI[(Mi[w-1]+1):(Mi[w]),])
    }
  })
  
  cwY2<-lapply(1:r, function(k){
    do.call(cbind,lapply(1:n, function(i){
      outer(b.hat[,k],prdXI[[i]][,k])*Ahat[i,k]
    }))
  })
  
  prdY2<-Reduce(`+`,cwY2)
  
  sig<-mean(as.numeric((data.unfold-prdY2)^2))
  
  # Initial values for sigma^2_k
  # lambda adjusted a_i
  LAhat<-sapply(1:r,function(k){
    (as.numeric(coef(lm(as.numeric(data.unfold)~as.numeric(cwY2[[k]])-1))))*Ahat[,k]
  })
  
  Ahat<-LAhat
  
  
  ## Initialize for gamma
  Resp<-response #cbind(1,response)
  gammaP<-sapply(1:r,function(k){
    as.numeric((solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(Ahat[,k],ncol = 1)))
  })
  
  
  # for all subject, rows correspond to subjects
  Xb<-Resp%*%gammaP
  
  
  ## initialize for sigma_r
  
  sigR<-colMeans((Ahat-Xb)^2)
  
  
  # Necessary matrices
  Hmat<-lapply(1:n, function(i){
    sapply(1:r, function(k){
      as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
    })
  })
  
  # for all subject, rows correspond to subjects
  HiB<-Resp%*%gammaP
  
  
  ## Subject wise super X matrix
  Hbeta<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
  })
  
  # Subject wise vectorized data 
  Yvec<-lapply(1:n, function(i){
    as.numeric(datlist[[i]][-1,])
  })
  
  
  # Starting of the EM based estimation
  t<-1
  iter_dif<-NULL
  while(t<=(maxiter+cvT)){
    dif<-NULL
    # E step for EM algorithm
    # Conditional variance
    sw_vcov<-lapply(1:n, function(i){
      solve(((t(Hmat[[i]])%*%Hmat[[i]])*(1/sig))+(diag(r)*(1/sigR)))
    })
    
    if(r==1){
      varU<-as.matrix(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    } else{
      varU<-t(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    }
    
    # conditional mean for all subjects
    if(r==1){
      conU<-as.matrix(sapply(1:n,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    } else{
      conU<-t(sapply(1:n,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    }
    
    
    # another necessary matrix
    HUmat<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
    })
    
    
    
    ## objective function value
    if(t==1){
      obj_val<-(-(sum(sapply(Yvec,length)))*log(sig))+
        (-sum(sapply(1:n,function(i){
          sum(diag((t(Hmat[[i]])%*%Hmat[[i]])%*%sw_vcov[[i]]))
        })))/(2*sig)+
        (-sum(sapply(1:n,FUN = function(i){
          sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2)
        }))/(2*sig))+
        (-sum(colSums(conU^2+varU)/sigR))+
        (-(n/2)*log(sig))
    }
    # M step
    # update of beta_k for all k
    eta_val<-NULL
    for(k in 1:r){
      #Xb<-Resp%*%gammaP
      Anew<-conU+Xb
      if(r>1){
        indx<-c(1:r)[-k]
      }
      
      # residual data
      if(r==1){
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]
        })
      } else{
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]-Reduce(`+`,lapply(indx, function(l){
            Anew[i,l]*outer(b.hat[,l],prdXI[[i]][,l])
          }))
        })
      }
      
      
      # design matrix for beta update
      Xs<-do.call(rbind,lapply(1:n, function(i){
        do.call(rbind,lapply(prdXI[[i]][,k], function(j){
          t(outer(Resp[i,],b.hat[,k])*j)
        }))
      }))
      
      ## Creating Z vector
      Zs<-matrix(do.call(c,lapply(1:n, function(i){
        as.numeric(rY[[i]]-(conU[i,k]*outer(b.hat[,k],prdXI[[i]][,k])))
      })))
      
      ngamma<-solve(t(Xs)%*%(Xs))%*%(t(Xs)%*%Zs)
      dif<-c(dif,sum((ngamma-gammaP[,k])^2)/sum(gammaP[,k]^2))
      gammaP[,k]<-ngamma
      
      # updated Xb
      Xb<-Resp%*%gammaP
      Anew<-conU+Xb
      
      # Scaling factor
      scD<-sqrt((Anew[,k]^2)+sapply(sw_vcov,function(u){u[k,k]}))
      
      if(r==1){
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]]))*(1/scD[i])
        })
      } else{
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]])-Reduce(`+`,lapply(indx,function(l){
            outer(b.hat[,l],prdXI[[i]][,l])*sw_vcov[[i]][k,l]
          })))*(1/scD[i])
        })
      }
      
      # update of xi_b k
      
      # scaled Residual data
      sY<-sapply(1:p,function(b){
        do.call(c,lapply(srY, function(u){u[b,]}))
      })
      
      
      nbhat<-NULL
      for(b in 1:p){
        bX<-do.call(c,lapply(1:n,function(i){scD[i]*prdXI[[i]][,k]}))
        
        nbhat<-c(nbhat,sum(sY[,b]*bX)/sum(bX^2))
      }
      n_bhat<-nbhat/sqrt(sum(nbhat^2))
      dif<-c(dif,sum((n_bhat-b.hat[,k])^2))
      b.hat[,k]<-n_bhat
      
      # Update of psi function
      Lty = lapply(1:n,function(i){
        list(scD[i]*as.numeric(b.hat[,k]%*%srY[[i]]))
        #(matrixcalc::hadamard.prod(scF[[i]],rY[[i]]))))
      })
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K)
        Smv<-cvfit[[1]]
        phiH = cvfit[[2]]
        eta_val<-c(eta_val,Smv)
      } else{
        phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
        eta_val<-c(eta_val,ifelse(CVPhi,Smv,smooth))
      }
      phi.new = phiH / sqrt(sum(phiH^2))
      dif<-c(dif,sum((phi.hat[,k] - phi.new)^2))
      phi.hat[,k] <- phi.new
      
      kpXI<-approx(t_arg,phi.new,xout = allT,rule=2)$y
      
      for(i in 1:n){
        if(i ==1){
          prdXI[[i]][,k]<-kpXI[1:Mi[1]]
        } else{
          prdXI[[i]][,k]<-kpXI[(Mi[i-1]+1):(Mi[i])]
        }
      }
    }
    
    # update of sigma2_k
    sk2_new<-colMeans(conU^2)
    
    dif<-c(dif,((sk2_new-sigR)^2)/sigR^2)
    
    sigR<-sk2_new
    
    # to avoid sigma2k=0
    if(any(sk2_new<1e-10)){
      sk2_new[sk2_new<1e-10]<-1e-10
    }
    
    
    # update of sigma2
    # Necessary matrices (updated)
    Hmat<-lapply(1:n, function(i){
      sapply(1:r, function(k){
        as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
      })
    })
    
    # for all subject, rows correspond to subjects
    HiB<-Resp%*%gammaP
    
    
    ## Subject wise super X matrix
    Hbeta<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
    })
    
    HUmat<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
    })
    
    nsig2<-sum(sapply(1:n, function(i){
      sum(diag((t(Hmat[[i]])%*%Hmat[[i]])%*%sw_vcov[[i]]))+
        sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2)
    }))/sum(sapply(Yvec,length))
    
    dif<-c(dif,((nsig2-sig)^2)/sig^2)
    sig<-nsig2
    
    # objective function
    nobj_val<-(-(sum(sapply(Yvec,length)))*log(sig))+
      (-sum(sapply(1:n,function(i){
        sum(diag((t(Hmat[[i]])%*%Hmat[[i]])%*%sw_vcov[[i]]))
      })))/(2*sig)+
      (-sum(sapply(1:n,FUN = function(i){
        sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2)
      }))/(2*sig))+
      (-sum(colSums(conU^2+varU)/sigR))+
      (-(n/2)*log(sig))+sum(eta_val)
    
    dif<-c(dif,nobj_val)
    iter_dif<-rbind(iter_dif,dif)
    
    if(t<=(cvT+1)){
      relC<-epsilon+0.05
    } else{
      relC<-(nobj_val-obj_val)/abs(obj_val)
    }
    
    obj_val<-nobj_val
    
    t <- t+1
    
    if(t>cvT & relC<=epsilon) 
      break
  }
  
  # calculate lambda
  tenY<-do.call(c,Yvec)
  
  # using Ahat
  cX<-do.call(rbind,lapply(1:n,function(i){
    sapply(1:r, function(k){
      as.numeric(Anew[i,k]*outer(b.hat[,k],prdXI[[i]][,k]))
    })
  }))
  
  cl.fit = lm(tenY~cX-1)
  cLambda = as.numeric(cl.fit$coefficients)
  cR2<-sapply(1:r, function(k){
    summary(lm(tenY~cX[,k]-1))$r.squared
  })
  
  cRsq<-sapply(1:r, function(k){
    summary(lm(tenY~cX[,1:k]-1))$r.squared
  })
  
  # using Xb
  X<-do.call(rbind,lapply(1:n,function(i){
    sapply(1:r, function(k){
      as.numeric(Xb[i,k]*outer(b.hat[,k],prdXI[[i]][,k]))
    })
  }))
  
  
  
  l.fit = lm(tenY~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  mRsq=sapply(1:r, function(k){
    summary(lm(tenY~X[,k]-1))$r.squared
  })
  
  Rsq<-sapply(1:r, function(k){
    summary(lm(tenY~X[,1:k]-1))$r.squared
  })
  
  
  # elements of returning object
  PCname <- paste('Component', 1:r)
  colnames(Anew) = PCname
  colnames(Xb) = PCname
  colnames(b.hat) = PCname
  colnames(phi.hat) = PCname
  rownames(Anew) = names(datlist)
  rownames(Xb) = names(datlist)
  rownames(b.hat) = rownames(datlist[[1]])[-1]
  
  ## Time point where phi functions are estiamted
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  
  # return object
  results = list("A.hat" = Anew, "B.hat" = b.hat, 
                 "Phi.hat" = phi.hat, "time" = time.return,
                 "Lambda" = cLambda, "r.square" = cR2, "accum.r.square" = cRsq, 
                 "r.square.Xb" = mRsq, "accum.r.square.Xb" = Rsq,"supR2"=cRsq[r],"supR2.Xb"=Rsq[r],
                 "Gamma"=gammaP,"Xb"=Xb,"Sigma2R"=sigR,"Sigma2"=sig,cnvrgnt=iter_dif)
  return(results)
}


#' Prediction using the fitted objects from supFTSVD
#' 
#' @param obj fitted object by supFTSVD
#' @param designM matrix of covariates for all subjects (rows represent subject)
#' @param new_dat observed data for subjects; only needed when subject specific 
#' prediction is the interest. The default is NULL
#' @param newT a list with time points where data for new subjects are observed
#' @param TimeG Grid in the time domain where prediction will be made.
#' @returns a list with predicted values. Elements correspond to the subjects
#' @export
predict.supFTSVD<-function(obj,designM,new_dat=NULL,newT=NULL,TimeG=NULL){
  # number of components in the fitted model
  r<-ncol(obj$A.hat)
  # number of subjects for which prediction will be made
  ns<-nrow(designM)
  # Prediction of subject loading
  Xb<-designM%*%obj$Gamma
  # Predicted singular functions
  if(!is.null(new_dat)){
    if(is.null(newT)){
      stop("Provide value for newT")
    }
    newSF<-lapply(newT, function(nT){
      sapply(1:r, function(k){
        approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
      })
    })
    
    Hmat<-lapply(1:ns, function(i){
      as.matrix(sapply(1:r,function(k){
        as.numeric(outer(obj$B.hat[,k],newSF[[i]][,k]))
      }))
    })
    
    Hbeta<-lapply(1:ns, function(i){
      as.numeric(Hmat[[i]]%*%matrix(Xb[i,],ncol=1))
    })
    
    sw_vcov<-lapply(1:ns, function(i){
      solve(((t(Hmat[[i]])%*%Hmat[[i]])*(1/obj$Sigma2))+(diag(r)*(1/obj$Sigma2R)))
    })
    
    ## Vectorization of data
    Yvec<-lapply(1:ns, function(i){
      as.numeric(new_dat[[i]])
    })
    
    # conditional mean for all subjects
    if(r==1){
      conU<-as.matrix(sapply(1:ns,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    } else{
      conU<-t(sapply(1:ns,FUN = function(i){
        sw_vcov[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))
    }
    Ahat<-Xb+conU
  } else{
    Ahat<-Xb
  }
  
  list("meanTRJ"=lapply(1:ns,function(i){
    Reduce(`+`,lapply(1:r, function(k){
      outer(obj$B.hat[,k],obj$Phi.hat[,k])*Xb[i,k]*obj$Lambda[k]
    }))
  }),
  "subTRJ"=lapply(1:ns,function(i){
    Reduce(`+`,lapply(1:r, function(k){
      outer(obj$B.hat[,k],obj$Phi.hat[,k])*Ahat[i,k]*obj$Lambda[k]
    }))
  }))
}


#' Prediction using the fitted objects from supFTSVD
#' 
#' @param obj fitted object by supFTSVD
#' @param new_dat observed data for subjects for which we want to make prediction
#' @param newT a list with time points where data for new subjects are observed
#' @param TimeG Grid in the time domain where prediction will be made.
#' @returns a list with predicted values. Elements correspond to the subjects
#' @export
predict.ftsvd<-function(obj,new_dat,newT,TimeG){
  # number of components in the fitted model
  r<-ncol(obj$A.hat)
  # number of subjects for which prediction will be made
  ns<-length(new_dat)
  # Predicted singular functions
  newSF<-lapply(newT, function(nT){
    sapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
    })
  })
  
  # matrix for estimating the subject loading
  Hmat<-lapply(1:ns, function(i){
    as.matrix(sapply(1:r,function(k){
      as.numeric(outer(obj$B.hat[,k],newSF[[i]][,k]))
    }))
  })
  
  ## Vectorization of data
  Yvec<-lapply(1:ns, function(i){
    as.numeric(new_dat[[i]])
  })
  
  ## Subject loading
  if(r==1){
    Ahat<-as.matrix((sapply(1:ns, function(i){
      as.numeric((solve(t(Hmat[[i]])%*%Hmat[[i]]))%*%(t(Hmat[[i]])%*%Yvec[[i]]))
    })))
  } else{
    Ahat<-t(sapply(1:ns, function(i){
      as.numeric((solve(t(Hmat[[i]])%*%Hmat[[i]]))%*%(t(Hmat[[i]])%*%Yvec[[i]]))
    }))
  }
  
  # Predicted data
  lapply(1:ns,function(i){
    Reduce(`+`,lapply(1:r, function(k){
      outer(obj$B.hat[,k],obj$Phi.hat[,k])*Ahat[i,k]
    }))
  })
}



#' Supervised Functional SVD
#' @param datlist A list of n matrices data, each matrix represents a subject;
#' the first row represents the sampling timestamp; 
#' following rows represent the feature measurements.
#' @param r CP-rank/Number of principle components. Default: 3.
#' @param interval The range of time points. Default: range of all observed data. 
#' @param resolution Resolution for the output eigencurves. Default: 100.
#' @param smooth Smoothing parameter for RKHS norm. Larger means smoother functions. Default: 1e-8.
#' @param maxiter Maximum number of itereation. Default: 20.
#' @param epsilon Convergence criteria for a and b. Default: 1e-4.
#' @return The estimations of the principle components in suject/feature/time;
#'       Var.prop: Explained variances.
SUPftsvd <- function(datlist, response, interval = NULL, r = 3,resolution=100, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                     maxiter=20, epsilon=1e-4){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  # response in matrix form
  if(!is.matrix(response)){
    response<-matrix(response,ncol=1)
  }
  
  
  
  Gamma_par=matrix(0,r,ncol(response))
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  SAhat = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(SAhat) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  
  # storing original response
  ires<-response
  
  
  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res = NULL
  Lambda = rep(0, r)
  X = NULL
  cX = NULL
  Rsq <- accumRsq <- rep(0, r)
  cRsq <- caccumRsq <- rep(0, r)
  Sig2<-NULL
  Sig2R<-NULL
  ASigma2<-NULL
  supR2<-NULL
  
  ti = lapply(1:n,FUN=function(u){
    temp = 1 + round((resolution-1) * (datlist[[u]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    temp
  })
  
  y0=do.call(c,lapply(1:n,FUN=function(u){
    as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
  }))
  
  y<-y0
  
  tipos<-lapply(ti,function(u){u>0})
  
  Lt<-lapply(datlist,FUN = function(u){u[1,]})
  ind_vec<-rep(1:n,sapply(Lt,length))
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  
  
  for (s in 1:r){ 
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    #print(sprintf("Calculate the %dth Component", s))
    print(paste("Calculating the Component", s))
    
    
    ## initialize U
    a.hat<-sapply(Lt,length) #rep(1,length(Lt))/sqrt(length(Lt))
    a.hat<-a.hat/sqrt(sum(a.hat^2))
    
    ## 
    data.unfold = do.call(cbind,lapply(datlist,FUN=function(u){u[-1,]}))
    
    #    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    #    b.hat = b.initials[,1]
    b.initials <- rsvd(data.unfold, k=r)$u
    b.hat = b.initials[,1]
    
    # Initialize xi function
    # compress the data and apply function PCA.
    Ly<-lapply(1:n,FUN = function(i){
      a.hat[i]*as.numeric(b.hat%*%datlist[[i]][2:(p+1),])
    })
    
    if(CVPhi){
      phiH = cv_freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth, kfold = K)[[2]]
    } else{
      phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    }
    #phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    phi.hat = phiH / sqrt(sum(phiH^2))
    
    t_arg<-seq(interval[1],interval[2],length.out = resolution)
    
    allT<-do.call(c,Lt)
    pXI<-approx(t_arg,phi.hat,xout = allT,rule=2)$y
    mi<-as.numeric(sapply(Lt,length))
    Mi<-cumsum(mi)
    
    prdXI<-lapply(1:n,function(w){
      if(w ==1){
        pXI[1:Mi[1]]
      } else{
        pXI[(Mi[w-1]+1):(Mi[w])]
      }
    })
    
    
    # initial values for error variance
    AB<-outer(b.hat,a.hat)
    prdY2<-do.call(cbind,lapply(1:n,function(w){
      outer(AB[,w],prdXI[[w]])
    }))
    
    sig<-mean(as.numeric((data.unfold-prdY2)^2))
    # lambda adjusted a_i
    a.hat<-(as.numeric(coef(lm(as.numeric(data.unfold)~as.numeric(prdY2)-1))))*a.hat
    
    ## Initialize for gamma
    Resp<-response#cbind(1,response)
    gammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(a.hat,ncol = 1))
    
    
    ## initialize for sigma_r
    
    sigR<-mean((a.hat-(Resp%*%gammaP))^2)
    
    
    # iteratively update a,b,phi
    t <- 1
    dif <- 1
    while(t<=(maxiter+cvT) & dif>epsilon){
      #sdif<-NULL
      # E-step
      # conditional mean
      dFAC<-lapply(1:n,function(w){
        sigR*outer(b.hat,prdXI[[w]])
      })
      
      
      conV<-sapply(1:n,function(w){
        qA2<-sum(as.numeric(((1/sigR)*hadamard.prod(dFAC[[w]],dFAC[[w]]))))+sig
        (sig*sigR)/as.numeric(qA2)
      })
      
      sconF<-conV/(sig*sigR)
      
      conM<-sapply(1:n,function(w){
        qA<-dFAC[[w]]
        qB<-sum(as.numeric(hadamard.prod(datlist[[w]][-1,],qA)))+as.numeric((sig*(t(Resp[w,])%*%gammaP)))
        as.numeric(qB*sconF[w])
      })
      
      #conM<-conM/sqrt(sum(conM^2))
      
      #dif <- sum((a.hat - conM)^2)/sum(a.hat^2)
      #a.hat<-conM
      #sdif<-c(sdif,dif)
      
      
      
      # M-Step
      
      # update of b:
      b.tilde<-apply(data.unfold,1,function(w){
        sum(w*rep(conM,mi)*pXI)/sum(rep((conV+conM^2),mi)*pXI^2)
      })
      
      b.new<-b.tilde/sqrt(sum(b.tilde^2))
      if(t>cvT)
        dif <- max(dif, sum((b.hat - b.new)^2))
      #dif <- sum((b.hat - b.new)^2)
      b.hat <- b.new
      #sdif<-c(sdif,dif)
      
      # update phi:
      scD<-sqrt(conV+conM^2)
      scF<-conM/scD
      Lty = list()
      for (i in 1:n){
        Lty = c(Lty, list(scD[i]*as.numeric(b.hat%*%(scF[i]*datlist[[i]][2:(p+1),]))))
      }
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K)
        Smv<-cvfit[[1]]
        phiH = cvfit[[2]]
      } else{
        phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
      }
      #phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi.new = phiH / sqrt(sum(phiH^2))
      if(t>cvT)
        dif <- max(dif, sum((phi.hat - phi.new)^2))
      phi.hat <- phi.new
      #sdif<-c(sdif,dif)
      
      # update gamma
      #ngammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(conM,ncol = 1))
      subAF<-lm(conM~Resp-1)
      ngammaP<-as.numeric(coef(subAF))
      
      if(t>cvT)
        dif <- max(dif, sum((gammaP - ngammaP)^2))
      gammaP <- ngammaP
      #sdif<-c(sdif,dif)
      
      ## update of sigma2_r
      #Ahat<-Resp%*%gammaP
      Ahat<-matrix(as.numeric(fitted(subAF)),ncol=1)
      if(t>cvT)
        dif <- sum((a.hat - Ahat)^2)/sum(a.hat^2)
      a.hat<-Ahat 
      ## Residual in sub_sing_vector
      Ehat<-matrix(as.numeric(resid(subAF)),ncol=1)#conM-Ahat
      sR2<-summary(subAF)$r.squared
      
      nsigR<-mean(conV)+mean((conM-Ahat)^2)
      if(t>cvT)
        dif <- max(dif, ((sigR - nsigR)^2)/(sigR^2))
      sigR <- nsigR # variance of sub_singV
      #sdif<-c(sdif,dif)
      
      ## update of sigma2
      pXI<-approx(t_arg,phi.hat,xout = allT,rule=2)$y
      prdXI<-lapply(1:n,function(w){
        if(w ==1){
          pXI[1:Mi[1]]
        } else{
          pXI[(Mi[w-1]+1):(Mi[w])]
        }
      })
      
      AB<-outer(b.hat,conM)
      prdY2<-do.call(cbind,lapply(1:n,function(w){
        outer(AB[,w],prdXI[[w]])
      }))
      
      adF<-sum(sapply(1:length(prdXI), function(w){
        sum(conV[w]*(prdXI[[w]]^2))
      }))/(p*sum(sapply(Lt,length)))
      
      nsig<-mean(as.numeric((data.unfold-prdY2)^2))+adF
      
      if(t>cvT)
        dif <- max(dif, ((sig - nsig)^2)/(sig^2))
      sig <- nsig # component wise sigma^2
      #sdif<-c(sdif,dif)
      #print(sdif)
      
      t <- t+1
    }
    
    # calculate lambda
    x = do.call(c,lapply(1:n,FUN=function(u){
      t.ind<-which(ti[[u]]>0)
      outer(a.hat[u],as.numeric(t(outer(b.hat,prdXI[[u]]))))
    }))

    cx = do.call(c,lapply(1:n,FUN=function(u){
      t.ind<-which(ti[[u]]>0)
      outer(conM[u],as.numeric(t(outer(b.hat,prdXI[[u]]))))
    }))
        
    #    x = NULL
    #    for (i in 1:n){
    #      t.temp = ti[[i]]
    #      t.temp <- t.temp[t.temp>0]
    #      x <- c(x,as.vector(t(a.hat[i]*b.hat%o%phi.hat[t.temp])))
    #    }
    X = cbind(X, x)
    cX = cbind(cX, cx)
    l.fit = lm(y~x-1)
    lambda = as.numeric(l.fit$coefficients)
    A[,s] = a.hat
    SAhat[,s] = Ehat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    Phi[,s] = t(phi.hat)
    Lambda[s] = lambda
    Rsq[s] <- summary(l.fit)$r.squared
    acF<-summary(lm(y0~X-1))
    accumRsq[s] <- acF$r.squared
    supR2[s]<-sR2
    cRsq[s]<-summary(lm(y~cx-1))$r.squared
    caccumRsq[s] <-summary(lm(y0~cX-1))$r.squared

    Gamma_par[s,] <- gammaP
    Sig2R[s] <- sigR
    Sig2[s] <- sig #
    ASigma2[s]<-acF$sigma^2
    ## update of theta
    #tfit<-lm(y~X-1)
    #vfit<-lm(ires~A-1)
    #a1q<-summary(tfit)$r.square
    #a2q<-summary(vfit)$r.square
    #theta<-(1-a2q)/(2-a1q-a2q)
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] = datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * conM[i] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
      #response[i]<-response[i]-(a.hat[i]*gam_par)  
    }
    
    y=do.call(c,lapply(1:n,FUN=function(u){
      as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
    }))
    
    print(paste0("Convergence reached at dif=", dif, ', iter=', t-cvT))
  }
  l.fit = lm(y0~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  #v.fit<-lm(ires~A-1)
  
  # revise the sign of Lambda
  for (r in length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] = -Lambda[r]
      B[,r] = -B[,r]
    }
  }
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results = list("A.hat" = A, "B.hat" = B, 
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq, "r.square1" = cRsq, "accum.r.square1" = caccumRsq,"supR2"=supR2,
                 "Gamma"=Gamma_par,"CSigma2"=Sig2,"ASigma2"=ASigma2,"Sigma2"=summary(l.fit)$sigma^2,"Sigma2R"=Sig2R,"SAhat"=SAhat)
  return(results)
}



#' Supervised Functional SVD
#' @param datlist A list of n matrices data, each matrix represents a subject;
#' the first row represents the sampling timestamp; 
#' following rows represent the feature measurements.
#' @param r CP-rank/Number of principle components. Default: 3.
#' @param interval The range of time points. Default: range of all observed data. 
#' @param resolution Resolution for the output eigencurves. Default: 100.
#' @param smooth Smoothing parameter for RKHS norm. Larger means smoother functions. Default: 1e-8.
#' @param maxiter Maximum number of itereation. Default: 20.
#' @param epsilon Convergence criteria for a and b. Default: 1e-4.
#' @return The estimations of the principle components in suject/feature/time;
#'       Var.prop: Explained variances.
SUPftsvdOAT <- function(datlist, response, interval = NULL, r = 3,resolution=100, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                     maxiter=20, epsilon=1e-4){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  # response in matrix form
  if(!is.matrix(response)){
    response<-matrix(response,ncol=1)
  }
  
  
  
  Gamma_par=matrix(0,r,ncol(response))
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  SAhat = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(SAhat) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  
  # storing original response
  ires<-response
  
  
  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res = NULL
  Lambda = rep(0, r)
  X = NULL
  Rsq <- accumRsq <- rep(0, r)
  Sig2<-NULL
  Sig2R<-NULL
  ASigma2<-NULL
  supR2<-NULL
  
  ti = lapply(1:n,FUN=function(u){
    temp = 1 + round((resolution-1) * (datlist[[u]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    temp
  })
  
  y0=do.call(c,lapply(1:n,FUN=function(u){
    as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
  }))
  
  y<-y0
  
  tipos<-lapply(ti,function(u){u>0})
  
  Lt<-lapply(datlist,FUN = function(u){u[1,]})
  ind_vec<-rep(1:n,sapply(Lt,length))
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  
  
  for (s in 1:r){ 
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    #print(sprintf("Calculate the %dth Component", s))
    print(paste("Calculating the Component", s))
    
    
    ## initialize U
    a.hat<-sapply(Lt,length) #rep(1,length(Lt))/sqrt(length(Lt))
    a.hat<-a.hat/sqrt(sum(a.hat^2))
    
    ## 
    data.unfold = do.call(cbind,lapply(datlist,FUN=function(u){u[-1,]}))
    
    #    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    #    b.hat = b.initials[,1]
    b.initials <- rsvd(data.unfold, k=r)$u
    b.hat = b.initials[,1]
    
    # Initialize xi function
    # compress the data and apply function PCA.
    Ly<-lapply(1:n,FUN = function(i){
      a.hat[i]*as.numeric(b.hat%*%datlist[[i]][2:(p+1),])
    })
    
    if(CVPhi){
      phiH = cv_freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth, kfold = K)[[2]]
    } else{
      phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    }
    #phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    phi.hat = phiH / sqrt(sum(phiH^2))
    
    t_arg<-seq(interval[1],interval[2],length.out = resolution)
    
    allT<-do.call(c,Lt)
    pXI<-approx(t_arg,phi.hat,xout = allT,rule=2)$y
    mi<-as.numeric(sapply(Lt,length))
    Mi<-cumsum(mi)
    
    prdXI<-lapply(1:n,function(w){
      if(w ==1){
        pXI[1:Mi[1]]
      } else{
        pXI[(Mi[w-1]+1):(Mi[w])]
      }
    })
    
    
    # initial values for error variance
    AB<-outer(b.hat,a.hat)
    prdY2<-do.call(cbind,lapply(1:n,function(w){
      outer(AB[,w],prdXI[[w]])
    }))
    
    sig<-mean(as.numeric((data.unfold-prdY2)^2))
    # lambda adjusted a_i
    a.hat<-(as.numeric(coef(lm(as.numeric(data.unfold)~as.numeric(prdY2)-1))))*a.hat
    
    ## Initialize for gamma
    if(is.null(response)|s>ncol(response)){
      Resp<-matrix(1,nrow=n,ncol=1)
      gammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(a.hat,ncol = 1))
    } else{
      Resp<-response[,s]#cbind(1,response)
      gammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(a.hat,ncol = 1))
    }
    #Resp<-response#cbind(1,response)
    #gammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(a.hat,ncol = 1))
    
    ## initialize for sigma_r
    
    sigR<-mean((a.hat-(Resp%*%gammaP))^2)
    
    print(paste("initialization is done for component",s))
    # iteratively update a,b,phi
    t <- 1
    dif <- 1
    while(t<=(maxiter+cvT) & dif>epsilon){
      #sdif<-NULL
      # E-step
      # conditional mean
      dFAC<-lapply(1:n,function(w){
        sigR*outer(b.hat,prdXI[[w]])
      })
      
      
      conV<-sapply(1:n,function(w){
        qA2<-sum(as.numeric(((1/sigR)*hadamard.prod(dFAC[[w]],dFAC[[w]]))))+sig
        (sig*sigR)/as.numeric(qA2)
      })
      
      sconF<-conV/(sig*sigR)
      
      conM<-sapply(1:n,function(w){
        qA<-dFAC[[w]]
        qB<-sum(as.numeric(hadamard.prod(datlist[[w]][-1,],qA)))+as.numeric((sig*(Resp[w]%*%gammaP)))
        as.numeric(qB*sconF[w])
      })
      
      #conM<-conM/sqrt(sum(conM^2))
      
      #dif <- sum((a.hat - conM)^2)/sum(a.hat^2)
      #a.hat<-conM
      #sdif<-c(sdif,dif)
      
      
      
      # M-Step
      
      # update of b:
      b.tilde<-apply(data.unfold,1,function(w){
        sum(w*rep(conM,mi)*pXI)/sum(rep((conV+conM^2),mi)*pXI^2)
      })
      
      b.new<-b.tilde/sqrt(sum(b.tilde^2))
      if(t>cvT)
        dif <- max(dif, sum((b.hat - b.new)^2))
      #dif <- sum((b.hat - b.new)^2)
      b.hat <- b.new
      #sdif<-c(sdif,dif)
      
      # update phi:
      scD<-sqrt(conV+conM^2)
      scF<-conM/scD
      Lty = list()
      for (i in 1:n){
        Lty = c(Lty, list(scD[i]*as.numeric(b.hat%*%(scF[i]*datlist[[i]][2:(p+1),]))))
      }
      
      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K)
        Smv<-cvfit[[1]]
        phiH = cvfit[[2]]
      } else{
        phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
      }
      #phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi.new = phiH / sqrt(sum(phiH^2))
      if(t>cvT)
        dif <- max(dif, sum((phi.hat - phi.new)^2))
      phi.hat <- phi.new
      #sdif<-c(sdif,dif)
      
      # update gamma
      #ngammaP<-(solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(conM,ncol = 1))
      subAF<-lm(conM~Resp-1)
      ngammaP<-as.numeric(coef(subAF))
      
      if(t>cvT)
        dif <- max(dif, sum((gammaP - ngammaP)^2))
      gammaP <- ngammaP
      #sdif<-c(sdif,dif)
      
      ## update of sigma2_r
      #Ahat<-Resp%*%gammaP
      Ahat<-matrix(as.numeric(fitted(subAF)),ncol=1)
      if(t>cvT)
        dif <- sum((a.hat - Ahat)^2)/sum(a.hat^2)
      a.hat<-Ahat
      ## Residual in sub_sing_vector
      Ehat<-matrix(as.numeric(resid(subAF)),ncol=1)#conM-Ahat
      sR2<-summary(subAF)$r.squared
      
      nsigR<-mean(conV)+mean((conM-Ahat)^2)
      if(t>cvT)
        dif <- max(dif, ((sigR - nsigR)^2)/(sigR^2))
      sigR <- nsigR # variance of sub_singV
      #sdif<-c(sdif,dif)
      
      ## update of sigma2
      pXI<-approx(t_arg,phi.hat,xout = allT,rule=2)$y
      prdXI<-lapply(1:n,function(w){
        if(w ==1){
          pXI[1:Mi[1]]
        } else{
          pXI[(Mi[w-1]+1):(Mi[w])]
        }
      })
      
      AB<-outer(b.hat,conM)
      prdY2<-do.call(cbind,lapply(1:n,function(w){
        outer(AB[,w],prdXI[[w]])
      }))
      
      adF<-sum(sapply(1:length(prdXI), function(w){
        sum(conV[w]*(prdXI[[w]]^2))
      }))/(p*sum(sapply(Lt,length)))
      
      nsig<-mean(as.numeric((data.unfold-prdY2)^2))+adF
      
      if(t>cvT)
        dif <- max(dif, ((sig - nsig)^2)/(sig^2))
      sig <- nsig # component wise sigma^2
      #sdif<-c(sdif,dif)
      #print(sdif)
      
      t <- t+1
    }
    
    # calculate lambda
    x = do.call(c,lapply(1:n,FUN=function(u){
      t.ind<-which(ti[[u]]>0)
      outer(a.hat[u],as.numeric(t(outer(b.hat,prdXI[[u]]))))
    }))
    
    #    x = NULL
    #    for (i in 1:n){
    #      t.temp = ti[[i]]
    #      t.temp <- t.temp[t.temp>0]
    #      x <- c(x,as.vector(t(a.hat[i]*b.hat%o%phi.hat[t.temp])))
    #    }
    X = cbind(X, x)
    l.fit = lm(y~x-1)
    lambda = as.numeric(l.fit$coefficients)
    A[,s] = a.hat
    SAhat[,s] = Ehat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    Phi[,s] = t(phi.hat)
    Lambda[s] = lambda
    Rsq[s] <- summary(l.fit)$r.squared
    acF<-summary(lm(y0~X-1))
    accumRsq[s] <- acF$r.squared
    supR2[s]<-sR2
    
    Gamma_par[s,] <- gammaP
    Sig2R[s] <- sigR
    Sig2[s] <- sig #
    ASigma2[s]<-acF$sigma^2
    ## update of theta
    #tfit<-lm(y~X-1)
    #vfit<-lm(ires~A-1)
    #a1q<-summary(tfit)$r.square
    #a2q<-summary(vfit)$r.square
    #theta<-(1-a2q)/(2-a1q-a2q)
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] = datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * conM[i] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
      #response[i]<-response[i]-(a.hat[i]*gam_par)  
    }
    
    y=do.call(c,lapply(1:n,FUN=function(u){
      as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
    }))
    
    print(paste0("Convergence reached at dif=", dif, ', iter=', t-cvT))
  }
  l.fit = lm(y0~X-1)
  Lambda = as.numeric(l.fit$coefficients)
  #v.fit<-lm(ires~A-1)
  
  # revise the sign of Lambda
  for (r in length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] = -Lambda[r]
      B[,r] = -B[,r]
    }
  }
  
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results = list("A.hat" = A, "B.hat" = B, 
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq,"supR2"=supR2,
                 "Gamma"=Gamma_par,"CSigma2"=Sig2,"ASigma2"=ASigma2,"Sigma2"=summary(l.fit)$sigma^2,"Sigma2R"=Sig2R,"SAhat"=SAhat)
  return(results)
}



log_comp_centralize <- function(datlist, r = 1){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  log.comp = matrix(0,n,p)
  for (i in 1:length(datlist)){
    log.comp[i,] = apply(datlist[[i]][-1,], 1, mean)
  }
  log.comp.svd = svd(log.comp, nu=r, nv=r)
  log.comp.mean = log.comp.svd$u %*% t(log.comp.svd$v * log.comp.svd$d[1:r])
  mf.new = datlist
  for (i in 1:length(datlist)){
    mf.new[[i]][-1,] = datlist[[i]][-1,] - log.comp.mean[i,]
  }
  results = list("datlist" = mf.new, "A.tilde" = log.comp.svd$u, 
                 "B.tilde" = log.comp.svd$v, "lambda.tilde" = log.comp.svd$d[1:r])
  return(results)
}


#' Format data table into input of ftsvd
#' @param taxon_table A table of read counts, with n rows for samples and p columns for taxa.
#' @param time_point The time stamp of each sample, relative to the start of the study. 
#' A length n vector.
#' @param subjectID The subject ID of each sample. A length n vector.
#' @param threshold A threshold for taxon filtering. 
#' Taxa with zero counts percentage >= threshold will be excluded.
#' @param pseudo_count A small number to add to all the counts before 
#' normalizing into proportions and log transformation.
#' @return Input for ftsvd. A list of matrices, each representing the 
format_ftsvd <- function(taxon_table, time_point, subjectID, threshold=0.95, 
                         feature_names=NULL, 
                         pseudo_count=0.5, transform="log_comp"){
  # format data table into a list as input for ftsvd(), 
  # read counts all have 1/2 added, before being normalized into proportions and log transformation
  # check length of subID and time_point
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:', 
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(taxon_table)) 
    stop('length of subjectID does not match taxon_table!')
  if (length(time_point)!=nrow(taxon_table)) 
    stop('length of time_point does not match taxon_table!')
  # keep taxon that has non-zeros in >=1-threshold samples
  if (is.null(feature_names)){
    taxon_table <- taxon_table[,colMeans(taxon_table==0)<threshold]
  }else{
    taxon_table <- taxon_table[,feature_names]
  }
  if(transform=='log_comp'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }else if(transform=='comp'){
    taxon_table <- taxon_table
    taxon_table <- t(taxon_table/rowSums(taxon_table))
  }else if(transform=='ast'){
    taxon_table <- taxon_table
    taxon_table <- t(asin(sqrt(taxon_table/rowSums(taxon_table))))
  }else if(transform=='clr'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- log(taxon_table/rowSums(taxon_table))
    taxon_table <- t(taxon_table-rowMeans(taxon_table))
  }else if(transform=='logit'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(taxon_table/rowSums(taxon_table))
    taxon_table <- log(taxon_table/(1-taxon_table))
  }else if(transform=='none'){
    taxon_table <- t(taxon_table)
  }else{
    print('Input transformation method is wrong! log_comp is applied instead')
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }
  taxon_table <- rbind(time_point, taxon_table)
  rownames(taxon_table)[1] <- 'time_point'
  subID <- unique(subjectID)
  nsub <- length(subID)
  
  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID
  
  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    # print(i)
    datlist[[i]] <- taxon_table[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
  }
  return(datlist)
}




##########################
# RKHS functional regression
#########################
bernoulli_kernel <- function(x, y){
  k1.x = x-0.5
  k1.y = y-0.5
  k2.x = 0.5*(k1.x^2-1/12)
  k2.y = 0.5*(k1.y^2-1/12)
  xy = abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy = 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy = k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

freg_rkhs <- function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=1e-8){
  A <- Kmat
  for (i in 1:length(Ly)){
    #A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]
    A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]^2
  }
  #cvec <- Kmat%*%unlist(Ly)
  cvec <- unlist(Ly)
  
  #A.temp = crossprod(A) + smooth * Kmat
  #A.temp.eig <- eigen(A.temp, symmetric = TRUE)
  #A.d <- A.temp.eig$value
  #A.d[A.d<1e-10] <- 1e-10
  #beta <- ( (A.temp.eig$vector)%*%(t(A.temp.eig$vector)/A.d) ) %*% cvec
  
  A.temp <- A + smooth*diag(ncol(A))
  beta <- solve(A.temp)%*%cvec
  
  phi.est <- Kmat_output %*% beta
  return(phi.est)
}

cv_freg_rkhs<-function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth,kfold=5){
  A <- Kmat
  for (i in 1:length(Ly)){
    A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]^2
  }
  cvec <- unlist(Ly)
  
  int<-length(cvec)%%kfold
  if(int==0){
    Kfold<-rep(1:kfold,each=length(cvec)/kfold)
  } else{
    Kfold<-c(rep(1:kfold,each=length(cvec)%/%kfold),1:(length(cvec)%%kfold))
  }
  
  KInd<-sample(Kfold,length(cvec),replace = FALSE)
  
  KFres<-sapply(1:length(smooth),function(j){
    sapply(1:kfold,function(i){
      aind<-which(KInd!=i)
      A.temp <- A[aind,aind] + smooth[j]*diag(length(aind))
      beta <- solve(A.temp)%*%cvec[aind]
      cor(cvec[KInd==i],as.numeric((Kmat[KInd==i,aind]%*%beta)))^2
      #mean((cvec[KInd==i]-as.numeric((Kmat[KInd==i,aind]%*%beta)))^2)
    })
  })
  Sm<-smooth[which.max(colMeans(KFres))]
  print(Sm)
  #CV fit
  A.temp <- A + Sm*diag(ncol(A))
  beta <- solve(A.temp)%*%cvec
  
  phi.est <- Kmat_output %*% beta
  list(Sm,phi.est)
}



#' Multiply loadings from res_ftsvd into denoised tensor
#' @param res_ftsvd Output of ftsvd
#' @param mean_svd Outpuf of log_comp_centralize
#' @return The denoised functional tensor
tdenoise <- function(res_ftsvd, mean_svd=NULL){
  n <- nrow(res_ftsvd$A.hat)
  p <- nrow(res_ftsvd$B.hat)
  resol <- nrow(res_ftsvd$Phi.hat)
  tensor.est <- array(0,dim=c(n,p,resol))
  if (!is.null(mean_svd))
    tensor.est <- (mean_svd$A.tilde %*% t(mean_svd$B.tilde * mean_svd$lambda.tilde)) %o%
    rep(1, resol)
  for (i in 1:ncol(res_ftsvd$A.hat)){
    tensor.est <- tensor.est+res_ftsvd$A.hat[,i]%o%res_ftsvd$B.hat[,i]%o%res_ftsvd$Phi.hat[,i]*res_ftsvd$Lambda[i]
  }
  dimnames(tensor.est)[[3]] <- res_ftsvd$time
  return(tensor.est)
}




plot_time_loading <- function(res){
  Phi.data <- res$Phi.hat
  npc <- ncol(Phi.data)
  ntime <- nrow(Phi.data)
  Phi.data <- data.frame(time=res$time, value=as.vector(Phi.data), 
                         component=as.factor(as.vector(t(matrix(rep(1:npc,ntime),npc,)))))
  ptime <- ggplot(data=Phi.data, aes(x=time, y=value, color=component)) + geom_line()
  return(ptime)
}

aggregate_feature <- function(res_ftsvd, res_svd=NULL, datlist){
  # estimated
  B.data <- as.data.frame(res_ftsvd$B.hat)
  tensor.est <- tdenoise(res_ftsvd, res_svd)
  tensor.est.agg <- apply(tensor.est, c(1,3), function(x){(t(B.data)%*%x)})
  dim(tensor.est.agg)
  npc <- ncol(B.data)
  metafeature.est <- NULL
  for (i in 1:npc){ 
    tmp <- data.frame(value=as.vector(tensor.est.agg[i,,]),
                      subID=rep(dimnames(tensor.est.agg)[[2]], dim(tensor.est.agg)[3]),
                      timepoint=as.vector(t(matrix(res_ftsvd$time,length(res_ftsvd$time),dim(tensor.est.agg)[2]))),
                      PC=paste0('Component ',i))
    metafeature.est <- rbind(metafeature.est,tmp)
  }
  metafeature.est$type <- 'estimated'
  # observed
  datlist.agg <- sapply(datlist, function(x){t(B.data)%*%x[-1,]}, simplify=F)
  metafeature.obs <- NULL
  for (i in 1:length(datlist.agg)){
    tmp <- data.frame(value=as.vector(datlist.agg[[i]]), 
                      subID=names(datlist.agg)[i],
                      timepoint=as.vector(t(matrix(datlist[[i]][1,], ncol(datlist[[i]]), npc))),
                      PC=rep(rownames(datlist.agg[[i]]), ncol(datlist.agg[[i]])))
    
    metafeature.obs <- rbind(metafeature.obs, tmp)
  }
  metafeature.obs$type <- 'observed'
  return(list(metafeature.obs=metafeature.obs, 
              metafeature.est=metafeature.est))
}


aggregate_subject <- function(res_ftsvd, res_svd=NULL, datlist){
  # estimated
  A.data <- as.data.frame(res_ftsvd$A.hat)
  tensor.est <- tdenoise(res_ftsvd, res_svd)
  tensor.est.agg <- apply(tensor.est, c(2,3), function(x){(t(A.data)%*%x)})
  dim(tensor.est.agg)
  npc <- ncol(A.data)
  metasubject.est <- NULL
  for (i in 1:npc){ 
    tmp <- data.frame(value=as.vector(tensor.est.agg[i,,]),
                      feature=rep(dimnames(tensor.est.agg)[[2]], dim(tensor.est.agg)[3]),
                      timepoint=as.vector(t(matrix(res_ftsvd$time,length(res_ftsvd$time),dim(tensor.est.agg)[2]))),
                      PC=paste0('Component ',i))
    metasubject.est <- rbind(metasubject.est,tmp)
  }
  
  return(metasubject.est)
}


tabular_feature_est <- function(tensor.est, feature_sel){
  tensor.est <- tensor.est[,feature_sel,]
  tab_est <- NULL
  tm <- as.numeric(dimnames(tensor.est)[[3]])
  for (i in 1:dim(tensor.est)[1]){
    tmp <- data.frame(subID=rownames(tensor.est)[i],
                      time=as.vector(t(matrix(tm, length(tm), ncol(tensor.est)))),
                      feature=rep(colnames(tensor.est), length(time)),
                      value=as.vector(tensor.est[i,,]))
    tab_est <- rbind(tab_est, tmp)
  }
  tab_est$type <- 'estimated'
  return(tab_est)
}


tabular_feature_obs <- function(datlist, feature_sel){
  tab_obs <- NULL
  for (i in 1:length(feature_sel)){
    value <- unlist(sapply(datlist, function(x){x[feature_sel[i],]}, simplify=F))
    time_point <- unlist(sapply(datlist, function(x){x['time_point',]}, simplify=F))
    nobs <- sapply(datlist, function(x){ncol(x)}, simplify=F)
    subID <- unlist(mapply(function(i){rep(names(datlist)[i], nobs[i])}, 
                           1:length(nobs), SIMPLIFY=F))
    tmp <- data.frame(subID=subID, time=time_point, feature=feature_sel[i], value=value)
    rownames(tmp) <- NULL
    tab_obs <- rbind(tab_obs, tmp)
  }
  tab_obs$type <- 'observed'
  return(tab_obs)
}



#' Calculate ROC using logistic regression
#' @param clust the observed labels
#' @param Xdata the predictors
#' @return output of function roc 
roc_logistic <- function(xtrain, ytrain, xtest, ytest){
  dftrain <- data.frame(y=ytrain, x=xtrain)
  dftest <- data.frame(y=ytest, x=xtest)
  fit <- glm(y ~ ., data = dftrain, family = "binomial")
  prob <- predict(fit, newdata=dftest, type = c("response"))
  g <- roc(dftest$y~prob)
  return(g$auc)
}

#' Estimate A of testing data based on B and Phi from training data
#' @param datlist testing data
#' @param res_ftsvd ftsvd result from training data
#' @return estimated A of testing data
est_A <- function(datlist, res_ftsvd){
  B <- res_ftsvd$B.hat
  Phi <- res_ftsvd$Phi.hat
  Lambda <- res_ftsvd$Lambda
  time.return <- res_ftsvd$time
  n <- length(datlist)
  p <- nrow(B)
  r <- ncol(B)
  resolution <- length(time.return)
  A_test <- matrix(0,n,r)
  y <- NULL
  ti <- vector(mode = "list", length = n)
  # get the coordinate of observed time points in the returned time grid
  for (i in 1:n){
    ti[[i]] <- sapply(datlist[[i]][1,], function(x){which.min(abs(x-time.return))})
    y <- c(y, as.numeric(t(datlist[[i]][-1,ti[[i]]>0])))
  }
  
  for (s in 1:r){
    for (i in 1:n){
      t.temp <- ti[[i]]>0
      A_test[i,s] <- B[,s] %*% datlist[[i]][2:(p+1),t.temp] %*% Phi[ti[[i]][t.temp],s] 
      A_test[i,s] <- A_test[i,s] / sum((Phi[ti[[i]][t.temp],s])^2) / Lambda[s]
      datlist[[i]][2:(p+1),t.temp] <- datlist[[i]][2:(p+1),t.temp] - 
        Lambda[s] * A_test[i,s] * (B[,s] %*% t(Phi[ti[[i]][t.temp],s])) 
    }
  }
  rownames(A_test) <- names(datlist)
  colnames(A_test) <- paste('Component', 1:r)
  return(A_test)
}


#' Conditional mean for new subjects using results of supervised FTSVD
#' 
#' @param fitOBJ a fitted object from SUPftsvd
#' @export
sub_pred<-function(fitOBJ,test_datlist){
  n<-length(test_datlist)
  
  # conditional mean
  dFAC<-lapply(1:n,function(w){
    sigR*outer(b.hat,prdXI[[w]])
  })
  
  
  conV<-sapply(1:n,function(w){
    qA2<-sum(as.numeric(((1/sigR)*hadamard.prod(dFAC[[w]],dFAC[[w]]))))+sig
    (sig*sigR)/as.numeric(qA2)
  })
  
  sconF<-conV/(sig*sigR)
  
  conM<-sapply(1:n,function(w){
    qA<-dFAC[[w]]
    qB<-sum(as.numeric(hadamard.prod(datlist[[w]][-1,],qA)))+as.numeric((sig*(t(Resp[w,])%*%gammaP)))
    as.numeric(qB*sconF[w])
  })
}
