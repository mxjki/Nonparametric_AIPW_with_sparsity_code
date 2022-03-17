#' This is based on revision.r, but the code are docomposed into parts.
#' This is based on the discussion on 1-22-2021.

library(quantreg)
library(MASS)
library('kernlab')
library('foreach')
library('parallel')##################################################### DC ###################################################################################
#------------------------------------------------------------------------#
##################################################### DC ###################################################################################
#------------------------------------------------------------------------#


min.fun<-function(min.x,min.c){
  
  min.temp<-rep(NA,length(min.x))
  for(i in 1:length(min.x)){
    min.temp[i]<-min(min.c,min.x[i])
  }
  
  
  return(min.temp)
}

##################################################### General Variable selection Method ##############################
############################################################################################################################


update.large<-function(up.x,up.y,up.K,up.sigma,up.lambda){
  up.n<-nrow(up.x)
  up.p<-ncol(up.x)
  ############################ Calculate Lipshitz Constant ############################################
  up.LipTemp<-matrix(rep(0),up.n,up.n)
  result<-rep(NA,up.p)
  # for(k in 1:up.n){
  #   up.LipTemp<- up.LipTemp+ 1/up.n*up.K[,k]%*%t(up.K[,k])
  # }
  #solve(up.LipTemp+up.lambda*up.K)%*%  rowSums(t(up.y)%x%rep(1,up.n)*up.K/up.n)
  ############################ Algorithm ###########################################################################
  up.alpha<-ginv(up.K+up.n*up.lambda*diag(up.n), 0)%*%up.y
 # for(j in 1:(up.p/up.paral)){
  #  result[((j-1)*up.paral+1):(j*up.paral)]<-sqrt(colSums(matrix(t(up.alpha)%*%t(rep(1,up.paral)%x%up.K * -matrix(as.vector((up.x[,((j-1)*up.paral+1):(j*up.paral)]%x%rep(1,up.n)-rep(1,up.n)%x%up.x[,((j-1)*up.paral+1):(j*up.paral)])),up.n*up.paral,up.n,byrow = T)/up.sigma),up.n,up.paral,byrow=F)^2)/up.n)
    result<-sqrt(colSums(matrix(t(up.alpha)%*%t(rep(1,up.p)%x%up.K * -matrix(as.vector((up.x%x%rep(1,up.n)-rep(1,up.n)%x%up.x)),up.n*up.p,up.n,byrow = T)/up.sigma),up.n,up.p,byrow=F)^2)/up.n)
 # }
  return(result)
}


#############################################################################################################################
############################################################################################################################
count<-function(cou.x,cou.thresh){
  cou.n<-length(cou.x)
  final<-rep(1,cou.n)
  for(i in 1:cou.n){
    if(cou.x[i]<=cou.thresh){
      final[i]<-0
    }
    
  }
  return(final)
}


max.fun<-function(max.x,max.c){
  
  max.temp<-rep(NA,length(max.x))
  for(i in 1:length(max.x)){
    max.temp[i]<-max(max.c,max.x[i])
  }
  
  
  return(max.temp)
}

############################################################################################################################
###############################################################################################################################################

calculate<-function(result1,result2){
  cou.p<-length(result1)
  n11<-0
  n12<-0
  n21<-0
  n22<-0
  Kapp<-0
  for(i in 1:cou.p){
    if(result1[i]==1 & result2[i]==1){
      n11<-n11+1
    }else if(result1[i]==1 & result2[i]==0){
      n12<-n12+1
    }else if(result1[i]==0&result2[i]==1){
      n21<-n21+1
    }
    else if(result1[i]==0&result2[i]==0){
      n22<-n22+1
    }
  }
  if((n11==cou.p)|(n22==cou.p) ){
    Kapp<--1
  }else{
    
    pra<-(n11+n22)/cou.p
    pre<-(n11+n12)*(n11+n21)/(cou.p^2)+(n12+n22)*(n21+n22)/(cou.p^2)
    
    Kapp<-(pra-pre)/(1-pre)
  }
  return(Kapp)
}
############################################################################################################################

################################################################  Kappa############################################################

##########################################

Kapp.cal<-function(kap.re1, kap.re2, kap.lam){
  
  
  
  kap.value<-0
  
  kap.re1[(kap.re1-kap.lam)<=0]<-0
  kap.re1[(kap.re1-kap.lam)>0]<-1
  kap.re2[(kap.re2-kap.lam)<=0]<-0
  kap.re2[(kap.re2-kap.lam)>0]<-1
  
  for(k1 in 1:nrow(kap.re1)){
    kap.value<-kap.value+calculate(kap.re1[k1,], kap.re2[k1,])
  }
  return(kap.value)
}






Kappa.fun<-function(kapp.x,kapp.y,kapp.lambda, kappa.tune){
  n.kapp<-length(kapp.y)
  p.kapp<-ncol(kapp.x)
  
  repeattime<-1
  Kappa<-rep(0,length(kapp.lambda))
  kappa.time<-5
  result1<-matrix(rep(0), kappa.time, p.kapp)
  result2<-matrix(rep(0), kappa.time, p.kapp)
  while(repeattime<=kappa.time){
    
    set.seed(repeattime)
    index1<-sample(1:n.kapp,floor(n.kapp/2),replace = F)
    #print(sort(index1,decreasing = T))
    x.1<-kapp.x[index1,]
    y.1<-kapp.y[index1]
    x.2<-kapp.x[-index1,]
    y.2<-kapp.y[-index1]
    ################# Group 1 ###############################
    
    sigma1<-1/sigest(x.1,scaled=F,frac=1)[2]
    K1<-kernelMatrix(rbfdot(sigma=0.5/sigma1),x.1)
   
    
    
    
    result1[repeattime,]<- update.large(x.1,y.1,K1,sigma1,kappa.tune)
    
    ################# Group 2 ###############################
    #  set.seed(kapp.seed)
    sigma2<-1/sigest(x.2,scaled=F,frac=1)[2]
    K2<-kernelMatrix(rbfdot(sigma=0.5/sigma2),x.2)
    
    
    
    
    result2[repeattime,]<-update.large(x.2,y.2,K2,sigma2, kappa.tune)
    
    repeattime<-repeattime+1
  }
  ################### Comparing ############################
  
  
  Kappa<- sapply(kapp.lambda, Kapp.cal, kap.re1=result1, kap.re2=result2)
  # print(calculate( count( result1-kapp.lambda,0), count(result2-kapp.lambda,0)))
  remove(K1,K2,result2,result1,x.1,y.1,x.2,y.2)
  
  #print(Kappa)
  return(Kappa)
}



Kernel.test<-function(Kx.tr,Kx.te){
  
  

  
  
  kn.tr<-nrow(Kx.tr)
  kn.te<-nrow(Kx.te)
  kp.tr<-ncol(Kx.tr)
  Kern.te<-kern.temp<-matrix(rep(0),kn.tr,kn.te)
  for(i in 1:kn.te){
    kern.temp[,i]<-rowSums((Kx.tr-rep(1,kn.tr)%x%t(Kx.te[i,]))^2)
  }
  kern.sigma<-median(kern.temp)
  Kern.te<-exp(kern.temp/(-2*kern.sigma))
  return(Kern.te)
}


main.fun<-function(main.x,main.y,main.x.test){
  

#### 默认设置 ##########
  main.lambda<-10^seq(-3,3,0.1)
  main.tune<-1e-3
  main.n<-nrow(main.x)
  ################################



  kappa.store<-Kappa.fun(main.x,main.y,main.lambda, main.tune)
  lam.sel<-main.lambda[  min(      which(kappa.store>=( 0.95*max(kappa.store))  )           )  ]
  
  #set.seed(1)
  main.sigma<-1/sigest(main.x,scaled=F,frac=1)[2]
  main.K<-kernelMatrix(rbfdot(sigma=0.5/main.sigma),main.x)
  result.store<-update.large(main.x,main.y,main.K,main.sigma, main.tune)
  final.store<-count(result.store-lam.sel,0)
  #print(result.store)
 ################ compute the fitted value for training sample ####################
  final.x<-main.x[,which(final.store!=0)]
  final.sigma<-1/sigest(final.x,scaled=F,frac=1)[2]
  final.K<-kernelMatrix(rbfdot(sigma=0.5/final.sigma),final.x)
  main.paral<-ncol(final.x)
  
  final.alpha<-ginv(final.K+main.n*1e-3*diag(main.n), 0)%*%main.y
  final.fit<-  final.K%*%final.alpha
  #final.fit.error<-  main.y-final.K%*%final.alpha
  ############### compute the fitted value for missing covarites ################
  final.alpha.test<-solve(final.K+main.n*1e-3*diag(main.n))%*%main.y
  final.x.test<-main.x.test[,which(final.store!=0)]
  final.K.test<-Kernel.test(final.x,final.x.test)
 final.fit.test<-  t(final.K.test)%*%final.alpha.test

  ##################################################
   return(list(final.store,final.fit,final.fit.test,lam_sel = lam.sel))

}
    
    

