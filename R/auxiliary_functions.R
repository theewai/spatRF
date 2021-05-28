#######################
# Auxiliary Functions #
#######################
.fastSolve <- function(mat){
  out <- try(chol2inv(chol(mat)),silent=TRUE)
  if(class(out)=="try-error") {
    out <- try(solve(mat),silent=TRUE)
    if(class(out)=="try-error"){
      out <- svd(mat)
      out <- out$u %*% diag(1/out$d) %*% t(out$v)
    }
  }
  out
}

.makeCmat <- function(mod,X,domain){
  num.branches <- length(mod$covar)
  n <- dim(X)[1]
  d <- length(unique(domain))
  Ck <- matrix(0,n,d*(num.branches+1))
  Ck[cbind(1:n,domain)] <- 1
  leaf.mem <- rep(1,n)
  for (a in 1:num.branches){
    Ca <-  as.numeric(as.numeric(X[,mod$covar[a]]) <= mod$split[a] & leaf.mem==mod$leaf[a])
    Ck[cbind(which(Ca==1),a*d+domain[which(Ca==1)])] <- 1 
    leaf.mem[which(Ca==1)] <- a+1
  }
  return(Ck)
}

.makeSpatBas <- function(pars,coords.train,coords.test=NULL, coords.vldt=NULL,cov.type, 
                         cov.opts=list(NULL),out.Z = FALSE){
  n <- dim(coords.train)[1]
  n.test <- dim(coords.test)[1]
  if(pars$psill!=0){
    if(cov.type == "exp"){
      ## Full Rank Kriging with Specified Range
      train.dist.mat <- as.matrix(crossDist(coords.train))
      Sigma <-  exp(-train.dist.mat/cov.opts$range) 
      i.sig <- .fastSolve(pars$psill * Sigma + pars$nugget * diag(n))
      if(!is.null(coords.test)){
        test.dist.mat <- as.matrix(crossDist(coords.test,coords.train))
        Sig.test <- pars$psill * exp(-test.dist.mat/cov.opts$range)
      }
      if(!is.null(coords.vldt)){
        vldt.dist.mat <- as.matrix(crossDist(coords.vldt,coords.train))
        Sig.vldt <- pars$psill * exp(-vldt.dist.mat/cov.opts$range)
      }
      if(out.Z){
        Zb <- .sqrtMat(Sigma)
      } 
    } else if ( cov.type == "LRK" ){
      ## Low Rank Kriging with Specified Range
      train.dist.mat <- as.matrix(crossDist(coords.train,cov.opts$knots))
      Z1 <- exp(-train.dist.mat/cov.opts$range)
      Om <- exp(-crossDist(cov.opts$knots)/cov.opts$range)
      i.sqrt.Om <- .iSqrtMat(Om,by.chol=TRUE)
      Zb <- Z1 %*% i.sqrt.Om
      i.sig <- .blockSolve(pars,Zb)
      if(!is.null(coords.test)){
        test.dist.mat <- as.matrix(crossDist(coords.test,cov.opts$knots))
        Z2 <- exp(-test.dist.mat/cov.opts$range)
        Sig.test <- pars$psill*Z2 %*% i.sqrt.Om %*% t(Z1 %*% i.sqrt.Om)
      }
    } else if(cov.type == "TPRS"){
      
      if(cov.opts$k==0){
        cov.opts$k <- n
      }
      K <- cov.opts$k
      if(is.null(cov.opts$m)) cov.opts$m <- 2
      m <- cov.opts$m
      ## Thin Plate Regression Splines
      Xbas <- cbind(1, coords.train)
      if(m > 2){
        for(jj in 2:(m-1)){
          temp.T <- coords.train^jj
          Xbas <- cbind(Xbas, temp.T)
          rm(temp.T)
        }
      }
      Xbas <- as.matrix(Xbas)
      M <- choose(m+2-1, 2) #d=2
      E <- .etaFunc(r=crossDist(coords.train), d=2, m=m)
      Ek <- eigen(E)
      Uk <- Ek$vectors[, 1:K]
      Dk <- diag(Ek$values[1:K])
      tUkT <- t(Uk)%*%Xbas
      QR <- qr(tUkT)
      Zk <- qr.Q(QR, complete=TRUE)[, (M+1):K]
      O <- t(Zk)%*%Dk%*%Zk
      O <- svd(O)
      O$d <- sqrt(O$d)
      O <- solve(O$u%*%diag(O$d)%*%t(O$v))
      Zb <- Uk%*%Dk%*%Zk%*%O
      i.sig <- .blockSolve(pars,Zb)
      if(!is.null(coords.test)){
        E1 <- .etaFunc(r=crossDist(coords.test, coords.train), d=2, m=m)
        E2 <- .etaFunc(r=crossDist(coords.train, coords.train), d=2, m=m)
        Sig.test <- pars$psill * E1%*%Uk%*%Zk%*%O%*%O%*%t(Zk)%*%t(Uk)%*%t(E2)
      }
    }
  } else{
    i.sig <- diag(rep(1/pars$nugget,n))
    if(!is.null(coords.test)){
      Sig.test <- matrix(0,n.test,n)
    }
  }
  out <- list(i.sig = i.sig)
  if(!is.null(coords.test)){
    out$sig.test <- Sig.test
  }
  if(!is.null(coords.vldt)){
    out$sig.vldt <- Sig.vldt
  }
  if(out.Z){
    out$Zb <- Zb
  }
  
  return(out)
}

.etaFunc <- function(m, d, r){
  if(d%%2 == 0){
    out <- (-1)^(m+1+d/2)/
      (2^(2*m-1)*pi^(d/2)*factorial(m-1)*factorial(m-d/2))*r^(2*m-d)*log(r)
  } else {
    out <- gamma(d/2-m)/(2^(2*m)*pi^(d/2)*factorial(m-1))*r^(2*m-d)
  }
  
  out[is.na(out)] <- 0
  return(out)
}

.iSqrtMat <- function(mat,by.chol=TRUE){
  if(by.chol){
    out <- try(chol(mat),silent=TRUE)
    if(class(out)=="try-error") {
      by.chol=FALSE
    } else{
      out <- (chol(chol2inv(out)))
    }
  }
  if(!by.chol){
    out <- svd(mat)
    out <- out$u %*% diag(1/sqrt(out$d)) %*% t(out$v)
  } 
  out
}

.sqrtMat <- function(mat,by.chol=TRUE){
  if(by.chol){
    out <- try(t(chol(mat)))
    if(class(out)=="try-error"){
      by.chol = FALSE
    }
  }
  if(!by.chol){
    out <- svd(mat)
    out <- out$u %*% diag(sqrt(out$d)) %*% t(out$v)
  }
  out
}

.blockSolve <- function(pars,Z){
  n <- dim(Z)[1]
  k <- dim(Z)[2]
  out <- try(diag(n) - Z %*% 
               .fastSolve (t(Z) %*% Z  + pars$nugget/pars$psill* diag(k)) %*% 
               t(Z)) / pars$nugget
  if(class(out) == "try-error"){
    out <- .fastSolve(pars$psill* Z %*% t(Z) + pars$nugget * diag(n))
  }
  out
}

.sampleBlock <- function(domain,group,replace=FALSE){
  ds <- sort(unique(domain))
  out <- NULL
  for(i in 1:length(ds)){
    groups <- unique(group[which(domain==ds[i])])
    if(replace){
      out <- c(out,sample(groups,replace=replace))
    } else{
      out <- c(out,sample(groups,replace=replace,size=0.632 * length(groups)))
    }
  }
  return(out)
}

.vecToBlock <- function(X,domain){
  ds <- unique(domain)
  d <- length(ds)
  out <- matrix(0,dim(X)[1],dim(X)[2]*length(ds))
  for( i in 1:dim(X)[2]){
    out[cbind(1:dim(X)[1],d*(i-1) + match(domain,ds))] <- X[,i]
  }
  return(out)
}

##' Computed the Euclidian distance matrix between to sets of points.
##'
##' @title Computed the Euclidian Distance Matrix
##' @param coord1,coord2 Matrices with the coordinates of locations, between
##'   which distances are to be computed.
##' @return A \code{dim(coord1)[1]}-by-\code{dim(coord2)[1]} distance matrix.
##' 
##' 
##' @author Johan Lindström
##' @family covariance functions
##' @family basic linear algebra
##' @export
##' @useDynLib spatRF C_dist

crossDist <- function(coord1, coord2=coord1){
  if( missing(coord2) ){
    symmetric <- TRUE
  }else{
    symmetric <- FALSE
  }
  ##cast to double
  coord1 <- as.matrix(coord1)
  storage.mode(coord1) <- "double"
  coord2 <- as.matrix(coord2)
  storage.mode(coord2) <- "double"
  
  ##call C-code, internal error-checking
  return( .Call(C_dist, coord1, coord2, as.integer(symmetric)) )
}