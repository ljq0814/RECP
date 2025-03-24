datageneration <- function(n,d,cploc,delta,existcp = TRUE, distri = "Case1"){
    muvec = rep(0,d)
    covmat = diag(1,d)
    #have no finite first moments#
    #Case I, null Cauchy(0,1) to Cauchy(delta, 1)#
    if(distri == 'Case1'){
        if(existcp){
            x = rbind(LaplacesDemon::rmvc(cploc,muvec,covmat),LaplacesDemon::rmvc(n-cploc,muvec+delta,covmat))
        }
        else{
            x = LaplacesDemon::rmvc(n,muvec,covmat)
        }
    }
    #case 2 mean change of log-normal variable#
    if(distri == 'Case2'){
      if(existcp){
        v <- LaplacesDemon::rmvn(cploc,muvec,3*covmat) #with variance 3
        w <- LaplacesDemon::rmvn(n-cploc,muvec+delta,3*covmat) #delta = 0.5
        x1 <- sapply(1:d, function(i){exp(v[,i])})
        x2 <- sapply(1:d, function(i){exp(w[,i])})
        x <- rbind(x1,x2)
      }
      else{
        v <- LaplacesDemon::rmvn(n,muvec,3*covmat) #with variance 3
        x <- sapply(1:d, function(i){exp(v[,i])})
      }
    }
    #case 3 mean change of Laplace variable#
    if(distri == "Case3"){
      if(existcp){
        x1 <- sapply(1:d, function(i){rlaplace(cploc,0,1)})
        x2 <- sapply(1:d, function(i){rlaplace(n-cploc,delta,1)}) #there is no change point
        x <- rbind(x1,x2)
      }
      else{
        x <- sapply(1:d, function(i){rlaplace(n,0,1)})
      }
    }
    #Case 4, multivariate normal mean change#
    if(distri == 'Case4'){
      if(existcp){
        x = rbind(LaplacesDemon::rmvn(cploc,muvec,3*covmat),LaplacesDemon::rmvn(n-cploc,muvec+delta,3*covmat))
      }
      else{
        x = LaplacesDemon::rmvn(n,muvec,3*covmat)
      }
    }
    #Case 5, multivariate normal mixture change#
    if(distri == "Case5"){
      if(existcp){
        x1 <- LaplacesDemon::rmvn(cploc,muvec,2*arcov(d,0.5))
        x2 <- rmixnorm(n-cploc,muvec,2*arcov(d,0.5),muvec+delta,2*arcov(d,0.5),0.5)
        x <- rbind(x1,x2)
      }
      else{
        x <- LaplacesDemon::rmvn(n,muvec,2*arcov(d,0.5))
      }
    }
    #Case 6, covariance change under normal distribution, cov(0.1) to cov(0.6)# AR
    if(distri == 'Case6'){
      if(existcp){
        cov1 <- toeplitz(0.1^(0:(d-1)))
        cov2 <- toeplitz((0.1+delta)^(0:(d-1)))
        x1 <- LaplacesDemon::rmvn(cploc,muvec,cov1)
        x2 <- LaplacesDemon::rmvn(n-cploc,muvec,cov2)
        x <- rbind(x1,x2)
      }
      else{
        cov1 <- toeplitz(0.1^(0:(d-1)))
        x <- LaplacesDemon::rmvn(n,muvec,cov1)
      }
    }
    #log-normal, cov change# AR
    #Case 7, delta = 0.5#
    if(distri == 'Case7'){
        if(existcp){
            cov1 <- toeplitz(0.1^(0:(d-1)))
            cov2 <- toeplitz((0.1 + delta)^(0:(d-1)))
            v <- LaplacesDemon::rmvn(cploc,muvec,cov1)
            w <- LaplacesDemon::rmvn(n-cploc,muvec,cov2)
            x1 <- exp(v)
            x2 <- exp(w)
            x <- rbind(x1,x2)
        }
        else{
            cov1 <- toeplitz(0.1^(0:(d-1)))
            v <- LaplacesDemon::rmvn(n,muvec,cov1)
            x <- exp(v)
        }
    }
    #Case 8, covariance change under normal distribution, rho from 0.1 to 0.6,delta=0.5# dense
    if(distri == 'Case8'){
      if(existcp){
        cov1 <- matrix(0.1,d,d)
        diag(cov1) = 1
        cov2 <- matrix((0.1 + delta),d,d)
        diag(cov2) = 1
        x1 <- LaplacesDemon::rmvn(cploc,muvec,cov1)
        x2 <- LaplacesDemon::rmvn(n-cploc,muvec,cov2)
        x <- rbind(x1,x2)
      }
      else{
        cov1 <- matrix(0.1,d,d)
        diag(cov1) = 1
        x <- LaplacesDemon::rmvn(n,muvec,cov1)
      }
    }
    #log-normal, cov change# dense
    #Case 9
    if(distri == 'Case9'){
        if(existcp){
            cov1 <- matrix(0.1,d,d)
            diag(cov1) = 1
            cov2 <- matrix((0.1 + delta),d,d)
            diag(cov2) = 1
            v <- LaplacesDemon::rmvn(cploc,muvec,cov1)
            w <- LaplacesDemon::rmvn(n-cploc,muvec,cov2)
            x1 <- exp(v)
            x2 <- exp(w)
            x <- rbind(x1,x2)
        }
        else{
            cov1 <- matrix(0.1,d,d)
            diag(cov1) = 1
            v <- LaplacesDemon::rmvn(n,muvec,cov1)
            x <- exp(v)
        }
    }
    x = as.matrix(x)
    return(x)
}


rspher <- function(n,d){
  tmp1 <- matrix(rnorm(n*d,0,1),n,d)
  nor <- apply(tmp1,1,norm,type = "2")
  return(tmp1/(nor%*%t(rep(1,d))))
}
recp <- function(n,d,mean,Sigma){
  u <- rspher(n,d)
  f <- rf(n,d,1)
  return(mean+u%*%halfsigma(Sigma)*(f%*%t(rep(1,d))))
}
halfsigma <- function(Sigma){
  p <- ncol(Sigma)
  tmpeig <- eigen(Sigma)
  return(tmpeig$vectors%*%diag(tmpeig$values^(0.5))%*%t(tmpeig$vectors))
}
rmixnorm <- function(n,mean1,Sigma1,mean2,Sigma2,prob){
  x1 <- LaplacesDemon::rmvn(n,mean1,Sigma1)
  x2 <- LaplacesDemon::rmvn(n,mean2,Sigma2)
  p <- rbinom(n,1,prob)
  return(x1*(1-p)+x2*p)
}
arcov <- function(p,rho = 0.5,var = 1){
  res <- matrix(1,p,p)
  for (i in 1:p){
    for (j in 1:p){
      res[i,j] <- rho^(abs(i-j))
    }
  }
  return(res*var)
}
cscov <- function(p,rho = 0.5,var = 1){
  res <- matrix(rho,p,p)
  for (i in 1:p){
    res[i,i] <- 1
  }
  return(res*var)
}
