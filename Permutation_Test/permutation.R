library(doParallel)
library("psych")
library(Matrix)

permutation.p<-function(B=10^3, N, t, p, r, target=c(1:r), x,z,y, ncores=4) {
  # y reads in the outcome, which is a t by N matrix. 
  # t is the level-1 sample size, such as the number of measurements
  # N is the level-2 sample size, such as the number of individuals.
  # p reads in the number of fixed effects.
  # r reads in the number of random effects.
  # x reads in the design matrix for the fixed effects, which is a t*N by p matrix.
  #    And in x, the first individual's t*p matrix is presented first, then the second individual's,
  #    and so on. The level-2 covariates should be after the level-1 covariates.
  # z reads in the design matrix for the random effects, which is a t*N by r matrix.
  #    And in z, the first individual's t*r matrix is presented first, then the second individual's,
  #    and so on.
  # target specifies the tested variance components, such as the first one (target=1 
  #    or all of them target=c(1:r).
  # ncores reads the number of parallel processes that run simultaneously to increase efficiency .
##################
y.list<-split(t(y), seq(ncol(y)))
grp_vect<- rep(1:N,each=t*p) #p dimension of x
x<-split(t(x), as.factor( grp_vect))
x<-lapply (x,function(x) matrix(x,byrow=T,nrow=t,ncol=p))
grp_vect<- rep(1:N,each=t*r) #r dimension of z
z<-split(t(z), as.factor( grp_vect))
z<-lapply (z,function(z) matrix(z,byrow=T,nrow=t,ncol=r))
#####calculate observed statistic

fun1<-function(xi){
  t(xi)%*%xi
}
w<-solve(Reduce('+', lapply(x,fun1)))

fun2<-function(xi,yi){
  t(xi)%*%yi
}
beta_ols<- w%*%apply(mapply(fun2,x,y.list),1,sum)

fun3<-function(zi,xi){
  t(zi)%*%zi-t(zi)%*%xi%*%w%*%t(xi)%*%zi
}
cv<- matrix(Reduce('+', mapply(fun3,z,x, SIMPLIFY=FALSE)))

fun4<-function(zi,xi){
  kronecker(t(zi)%*%zi,t(zi)%*%zi) - kronecker(t(zi)%*%zi,t(zi)%*%xi%*%w%*%t(xi)%*%zi)-
    kronecker(t(zi)%*%xi%*%w%*%t(xi)%*%zi,t(zi)%*%zi)
}
fun5<-function(zi, xi){
  kronecker(t(zi)%*%xi%*%w,t(zi)%*%xi%*%w) 
}
fun6<-function(zi,xi){
  kronecker(t(xi)%*%zi,t(xi)%*%zi) 
}
h<- Reduce('+', mapply(fun4,z,x,SIMPLIFY=FALSE))+Reduce('+', mapply(fun5,z,x,SIMPLIFY=FALSE))%*%
  Reduce('+', mapply(fun6,z,x,SIMPLIFY=FALSE))

q<-c(N*t-p-t(cv)%*%solve(h)%*%cv)

e.list<-mapply('-',y.list,lapply(x,function(xi) xi%*%beta_ols),SIMPLIFY=FALSE)

A<-matrix(NA,nrow=r*r+1,ncol=r*r+1)
A[1]<-N*t-p
A[2:(r*r+1),1]<-cv
A[1,2:(r*r+1)]<-t(cv)
A[2:(r*r+1),2:(r*r+1)]<-h

fun7<-function(zi,ei){
  kronecker(t(zi)%*%ei,t(zi)%*%ei) 
}

B.m<- matrix(c(sum(sapply(e.list,function(ei) t(ei)%*%ei )),
               apply(mapply(fun7,z,e.list),1,sum)))
vls<-solve(A)%*%B.m
U_vls<-vls[2:(r*r+1),1]
D_vls<- matrix(U_vls,nrow=r,ncol=r,byrow=F)

sigma_h<-vls[1,1]

fun8<-function(xi,zi,yi){
  D_vls%*%t(zi)%*%solve(diag(sigma_h, nrow(zi)) +
          zi%*%D_vls%*%t(zi))%*%(yi-xi%*%beta_ols)
}
b_h<-mapply(fun8,x,z,y.list)

##### y under null
b_i<-b_h
b_i[target,]<-0
b_i<-split(b_i,col(b_i))

fun9<-function(xi,yi,zi,bi){
  yi-xi%*%beta_ols-zi%*%bi
}

y1<-mapply(fun9,x,y.list,z,b_i)

y1.list<-split(t(y1), seq(ncol(y1)))

fun2<-function(xi,yi){
  t(xi)%*%yi
}
beta_ols_1<- w%*%apply(mapply(fun2,x,y1.list),1,sum)

e.list_1<-mapply('-',y1.list,lapply(x,function(xi) xi%*%beta_ols_1),SIMPLIFY=FALSE)

B.m_1<- matrix(c(sum(sapply(e.list_1,function(ei) t(ei)%*%ei )),apply(mapply(fun7,z,e.list_1),1,sum)))
vls_1<-solve(A)%*%B.m_1
U_vls_1<-vls_1[2:(r*r+1),1]
D_vls_1<- matrix(U_vls_1,nrow=r,ncol=r,byrow=F)
o_D_1<- D_vls_1

sigma_h_1<-vls_1[1,1]

z_star<-as.matrix(.bdiag(lapply(z,function(a) matrix(a[,target],nrow=t)))) 
funtb<-function(i){
  t_s<- 1/N*tr(z_star%*% kronecker(diag(1,nrow=N,ncol=N),o_D_1[target,target]) %*%t(z_star))
  return(t_s)}
ot_s<-mclapply(1:1,funtb, mc.cores=ncores)
ot_s<- ot_s[[1]]

########permutation
big.fun<-function (z,x,w,y.list,beta_ols) {
  fun3<-function(zi,xi){
    t(zi)%*%zi-t(zi)%*%xi%*%w%*%t(xi)%*%zi
  }
  cv<- matrix(Reduce('+', mapply(fun3,z,x, SIMPLIFY=FALSE)))
  
  fun4<-function(zi,xi){
    kronecker(t(zi)%*%zi,t(zi)%*%zi) - kronecker(t(zi)%*%zi,t(zi)%*%xi%*%w%*%t(xi)%*%zi)-
      kronecker(t(zi)%*%xi%*%w%*%t(xi)%*%zi,t(zi)%*%zi)
  }
  fun5<-function(zi, xi){
    kronecker(t(zi)%*%xi%*%w,t(zi)%*%xi%*%w) 
  }
  fun6<-function(zi,xi){
    kronecker(t(xi)%*%zi,t(xi)%*%zi) 
  }
  h<- Reduce('+', mapply(fun4,z,x,SIMPLIFY=FALSE))+Reduce('+', mapply(fun5,z,x,SIMPLIFY=FALSE))%*%
    Reduce('+', mapply(fun6,z,x,SIMPLIFY=FALSE))
  
  q<-c(N*t-p-t(cv)%*%solve(h)%*%cv)
  
  e.list<-mapply('-',y.list,lapply(x,function(xi) xi%*%beta_ols),SIMPLIFY=FALSE)
  
  return (list(cv,h,e.list))
}

estimate<-function(j) {
  xa<-x.m
  grp_vect<- rep(rep(1:t,each=p),N)
  xa<-split(t(xa), as.factor( grp_vect))
  xa<-lapply (xa,function(xa) matrix(xa,byrow=TRUE,nrow=N,ncol=p))
  
  za<-z.m
  grp_vect<- rep(rep(1:t,each=r),N)
  za<-split(t(za), as.factor( grp_vect))
  za<-lapply (za,function(za) matrix(za,byrow=TRUE,nrow=N,ncol=r))
  
  fun5<-function(i){
    index_i<-sample(c(1:N), N, replace = FALSE)
    yb<-y1[i,index_i]
    xb<-xa[[i]][index_i,]
    zb<-za[[i]][index_i,]
    return (list(xb,yb,zb))
  }
  set.seed(201888)
  xyb<-lapply(c(1:t),FUN=fun5)
  xb<-c()
  yb<-c()
  zb<-c()
  for (q in 1:t){
    xb<-rbind(xb,xyb[[q]][[1]])
    yb<-rbind(yb,xyb[[q]][[2]])
    zb<-rbind(zb,xyb[[q]][[3]])
  }
  
  yb.list<-split(yb, col(yb))
  
  grp_vect<- rep(rep(1:N,each=p),t)
  xb<-split(t(xb), as.factor( grp_vect))
  xb<-lapply (xb,function(xb) matrix(xb,byrow=TRUE,nrow=t,ncol=p))
  
  grp_vect<- rep(rep(1:N,each=r),t)
  zb<-split(t(zb), as.factor( grp_vect))
  zb<-lapply (zb,function(zb) matrix(zb,byrow=TRUE,nrow=t,ncol=r))
  
  fun1<-function(xi){
    t(xi)%*%xi
  }
  wb<-solve(Reduce('+', lapply(xb,fun1)))
  
  fun2<-function(xi,yi){
    t(xi)%*%yi
  }
  beta_ols_b<- wb%*%apply(mapply(fun2,xb,yb.list),1,sum)
  
  
  cvb<-big.fun (zb,xb,wb,yb.list,beta_ols_b)[[1]]
  hb<-big.fun (zb,xb,wb,yb.list,beta_ols_b)[[2]]
  eb.list<-big.fun (zb,xb,wb,yb.list,beta_ols_b)[[3]]
  
  Ab<-matrix(NA,nrow=r*r+1,ncol=r*r+1)
  Ab[1]<-N*t-p
  Ab[2:(r*r+1),1]<-cvb
  Ab[1,2:(r*r+1)]<-t(cvb)
  Ab[2:(r*r+1),2:(r*r+1)]<-hb
  
  fun7<-function(zi,ei){
    kronecker(t(zi)%*%ei,t(zi)%*%ei) 
  }
  
  Bb.m<- matrix(c(sum(sapply(eb.list,function(ei) t(ei)%*%ei )),
                  apply(mapply(fun7,zb,eb.list),1,sum)))
  vlsb<-solve(Ab)%*%Bb.m
  U_vlsb<-vlsb[2:(r*r+1),1]
  D_vlsb<- matrix(U_vlsb,nrow=r,ncol=r,byrow=F)
  
  z_star<-as.matrix(.bdiag(lapply(zb,function(a) matrix(a[,target],nrow=t)))) 
  funtb3<-function(i){
    t_s<- 1/N*tr(z_star%*% kronecker(diag(1,nrow=N,ncol=N),D_vlsb[target,target]) %*%t(z_star))
    return(t_s)}
  t_sb<-mclapply(1:1,funtb3, mc.cores=ncores)
  t_s_b<- t_sb[[1]]
  
  return(c(t_s_b))
  
}
set.seed(201888)
niter=B
results <-mclapply(1:niter,estimate, mc.cores=ncores)

results2<-do.call(rbind, results)

t_s_b<-results2

p_s<-sum(t_s_b >= ot_s )/niter 

return(paste("The p-value of testing the variance components (",paste(target,collapse=", "), ") is", p_s))
}

#######EXAMPLE########
## A growth curve model with two variance components
########################
library("mvtnorm")
set.seed(201888)
N=100;t=10
beta0=1;beta1=2;
beta<-matrix(c(beta0,beta1),nrow=1)
D11<-1; D_off<-0 ;D_d<-0
psi<- matrix(D_off,2,2)
psi[1,1]<-D11
psi[2,2]=D_d
x0<-rep(rep(1,t),N)
x1<-rep(c(1:t),N)
x<-cbind(x0,x1)
x.m<-x
p=2;r=2
z=x
z.m<-z
grp_vect<- rep(1:N,each=t*p) #p dimension of x
x<-split(t(x), as.factor( grp_vect))
x<-lapply (x,function(x) matrix(x,byrow=T,nrow=t,ncol=p))
grp_vect<- rep(1:N,each=t*r) #r dimension of z
z<-split(t(z), as.factor( grp_vect))
z<-lapply (z,function(z) matrix(z,byrow=T,nrow=t,ncol=r))
b<- rmvnorm(N, mean = rep(0, nrow(psi)),sigma = psi)
b<-split(b, seq(nrow(b)))
gen.y<- function(bi,xi,zi){
  xi%*%t(beta) +zi%*%as.matrix(bi) + rnorm(t,0,1)
}
y<-mapply(gen.y,b,x,z)

###testing both the intercept and slope variances
permutation.p(B=1000, N, t, p=2, r=2, target=c(1:r), x=x.m,z=z.m,y)
###testing only the slope variance
permutation.p(B=1000, N, t, p=2, r=2, target=c(2:2), x=x.m,z=z.m,y)


