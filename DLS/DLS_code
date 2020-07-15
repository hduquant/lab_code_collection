library(lavaan)
library(expm)
library(doParallel)



# There are two options for DLS: DLS_M which uses gammaN.model (based on model-implied covariance)
# and DLS_S which uses gammaN.sample (based on sample covariance)
methodlist <- c("gammaN.model","gammaN.sample")

# Choose a method
id <- 1
method <- methodlist[id]

# Call the Grant-White data from lavaan
data <- HolzingerSwineford1939
data <- data[data$school=="Grant-White",7:15]
xb <- data
colnames(xb) <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")

# specify the number of factors and number of variables
m <- 3
p <- 9
ps <- p*(p+1)/2
ms <- m*(m-1)/2
q <- 2*p+ms+1
df <- ps-q
n <- 145

# Get sample covariance
ps <- p*(p+1)/2
Scov <-cov(xb)*(n-1)/n
vscov <- vech(Scov)

# Defind an important function in DLS which provides the model-implied covariance
# and first order derivative
sigdsig <- function(p,m,theta0){
  
  ps <- p*(p+1)/2
  p_h <- p/m
  p2 <- p*2
  
  theta.important<-theta0[4]
  theta0<-  theta0[-4]
  
  lamb <- matrix(0,p,m)
  for (j in 1:m){
    p_hj <- (j-1)*p_h+1
    p_hj1 <- j*p_h
    lamb[p_hj:p_hj1,j] <- theta0[p_hj:p_hj1]
  }
  
  lamb[9,1]<- theta.important
  
  psi_vec <- theta0[(p+1):p2] 
  phi <- matrix(1,m,m)
  k <- 0
  if (m>1){
    for (i in 1:(m-1)){
      for (j in (i+1):m){
        k <- k+1
        phi[j,i] <- theta0[p2+k] 
        phi[i,j] <- theta0[p2+k]
      }
    }
  }
  # Calculate model-implied covariance
  Sigma0 <- lamb%*%phi%*%t(lamb) + diag(psi_vec) #model-implied covariance
  
  # Derivative with Lambda
  DL <- matrix(0,ps,(p+1)) 
  lambphi <- lamb%*%phi
  lambphi_t <- t(lambphi)
  p_hj1=0
  for (j in 1:m){
    p_hj <-  p_hj1+1
    p_hj1 <-  p_hj+2
    if (j==1){ p_hj1=4}
    for (k in p_hj:p_hj1){
      
      if (j==1){
      tt <- matrix(0,p,m)
      tt[k,j] <- 1
      ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
      DL[,k] <- vech(ttt)} else {
        tt <- matrix(0,p,m)
        tt[k-1,j] <- 1
        ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
        DL[,k] <- vech(ttt)
      }
      
      if (k==4){
        tt <- matrix(0,p,m)
        tt[9,1] <- 1
        ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
        DL[,k] <- vech(ttt)
      }
     # print(list(j,k,tt))
    }
  }
  
  #Derivative with Psi
  Dps <- matrix(0,ps,p)
  for (j in 1:p){
    tt <- matrix(0,p,p)
    tt[j,j] <- 1
    Dps[,j] <- vech(tt)
  }
  
  #Derivative with Phi
  if (m>1){
    ms <- m*(m-1)/2 
    Dphi <- matrix(0,ps,ms)  
    k <- 0 
    for (i in 1:(m-1)){
      for (j in (i+1):m){
        k <- k+1 
        tt <- matrix(0,m,m) 
        tt[j,i] <- 1 ; tt[i,j] <- 1 
        ttt <- lamb%*%tt%*%t(lamb)
        Dphi[,k] <- vech(ttt)
      }
    } 
    vdsig <- cbind(DL,Dps,Dphi)
  } else {  vdsig <- cbind(DL,Dps)}
  
  out <- list('Sigma0'=Sigma0, 'vdsig'=vdsig)
  return(out)
}

# calculate the sample covariance matrix
  cov.s <-cov(xb)*(n-1)/n
vscov <- vech(cov.s)

# Calculate W matrix
# Calculate gamma_adf (the asymptotically distribution free gamma)

  mean_xt <- apply(xb,2,mean)
  x_c <- xb-matrix(1,n,1)%*%mean_xt
  k=1
  sigmaele<-matrix(NA,nrow=ps,ncol=n)
  for(i in 1:p){
    for(j in i:p){
      sigmaele[k,]<- x_c[,i]*x_c[,j]
      k=k+1
    }
  }
  sigmaij=sigmakl=rowSums( sigmaele)/n 
  
  sigmaijkl=c()
  gammaadf=matrix(NA,nrow=ps,ncol=ps)
  k=1
  for(i in 1:ps){
    for(j in i:ps){
      sigmaijkl[k]<-sum(sigmaele[i,]*sigmaele[j,])/n #n
      gammaadf[i,j]<-sigmaijkl[k]-sum(sigmaele[i,])*sum(sigmaele[j,])/(n)^2
      gammaadf[j,i]<- gammaadf[i,j]
      k=k+1
    }
  }
  
  kk=1
  index<-matrix(NA,nrow=p,ncol=p)
  for(i in 1:p){
    for(j in i:p){
      index[i,j]=index[j,i]=kk
      kk=kk+1
    }}
  
# Calculate gamma_N.s (the normal theory based gamma with sample covariance)
  kk=1
  index<-matrix(NA,nrow=p,ncol=p)
  for(i in 1:p){
    for(j in i:p){
      index[i,j]=index[j,i]=kk
      kk=kk+1
    }}
  
  gammaNs<-matrix(NA,nrow=ps,ncol=ps)
  for(i in 1:p){
    for(j in i:p){
      for (k in  1:p){
        for (l in 1:p){
          if ( index[k,l]>=index[i,j] & l>=k){
            #    print(c(i,j,k,l))
            gammaNs[index[k,l],index[i,j]]<-  Scov[i,k]*Scov[j,l]+Scov[i,l]*Scov[j,k]
            gammaNs[index[i,j],index[k,l]]<-gammaNs[index[k,l],index[i,j]]
          }
        }
      }
    }
  }

# Select a tunning parameter value  
    a<-0.75

    if (method=="gammaN.sample"){
      weight<-try(solve( (1-a)* gammaadf+ a*gammaNs))
      if( class( weight)=="try-error"){
        diconverg <- 1   # not converge
      }
    }
# Starting value (MLE estiamtes)
theta0=c(0.82, 0.54, 0.69, 0.46, 0.97, 0.96, 0.93, 0.71, 0.90, 0.45, 0.65, 0.93, 0.60,
        0.48, 0.31, 0.42, 0.41, 0.56, 0.29,0.55, 0.39, 0.24)
   # theta0=theta.s

    sig <- sigdsig(p,m, theta0)  
    vdsig<-sig$vdsig  
    Sigma0<-sig$Sigma0  
    vsig0 <- vech(Sigma0)
    
    diconverg=0
    
# Stop when the update is smaller than ep 
    ep <- 0.0001 
    
#Begin iterations
    
    for(t in 1:300){
# Calculate gamma_N.model (the normal theory based gamma with model-implied covariance); 
# updated theta0 (point estiamte) and sigma0 (model-implied covariance)
      if (method=="gammaN.model"){
        gammaNm<-matrix(NA,nrow=ps,ncol=ps)
        for(i in 1:p){
          for(j in i:p){
            for (k in  1:p){
              for (l in 1:p){
                if ( index[k,l]>=index[i,j] & l>=k){
                  #  print(c(i,j,k,l))
                  gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                  gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                }
              }
            }
          }
        }
        
        weight<-try(solve((1-a)* gammaadf+ a*gammaNm))
        if( class( weight)=="try-error"){
          diconverg <- 1   # not converge
          break
        }
      }
      
      stdi <- try(solve(t(vdsig) %*% weight%*%vdsig))
      if( class(  stdi)=="try-error"){
 # Give 1 if it does not converge
        diconverg <- 1   # not converge
        break
      }
      eresdu <- vscov-vsig0 
      dtheta <- t( eresdu) %*% weight %*%  vdsig%*%  stdi
      theta0 <- theta0 + dtheta
      delta <- max(abs(dtheta))
      
      sig <- sigdsig(p,m, theta0)  
      vdsig<-sig$vdsig 
      Sigma0<-sig$Sigma0
      vsig0 <- vech(Sigma0)
      
      if(delta<=ep) {
# updated theta0 and sigma0 for the last time
        if (method=="gammaN.model"){
          gammaNm<-matrix(NA,nrow=ps,ncol=ps)
          for(i in 1:p){
            for(j in i:p){
              for (k in  1:p){
                for (l in 1:p){
                  if ( index[k,l]>=index[i,j] & l>=k){
                    #  print(c(i,j,k,l))
                    gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                    gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                  }
                }
              }
            }
          }
          
          weight<-try(solve((1-a)* gammaadf+a*gammaNm))
          if( class( weight)=="try-error"){
            diconverg <- 1   # not converge
            break
          }
        }
        break};
    }
    
t<300
diconverg==0
        
        dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
        dwe<-t(vdsig)%*%weight
        U<- weight-t(dwe)%*%dtd%*%dwe
        ugamma <-  U %*% gammaadf
      
        r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
        # Robust SE
        r.SE <- sqrt(diag(r.Var)/(n-1))
        r.SE
        
        Fstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))
        
        Tstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))*(n-1)
 
        # Calculate the Jiang-Yuan rank adjusted statistic 
        c<-sum(diag(ugamma ))/qr(ugamma)$rank
        rTstats <- Tstats/c
        
        # DLS estiamtes of parameters
        theta0
        
   
