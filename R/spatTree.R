#' spatTree
#'
#' @param Y Observed values
#' @param X Covariates used to construct the spatially adjusted regression tree
#' @param domain Used to construct multiple regression trees simultaneously. For
#'               details, see Hee Wai et al., 2020
#' @param X.fix Covariates to fit as linear fixed effects in addition to the
#'              spatially regression tree (Optional)
#' @param coords x and y coordinates for locations of the points
#' @param cov.type Type of covariance function to use.  Currently supported are 
#'                 'exp' (Exponential), 'LRK' (low rank kriging), and 'TPRS',
#'                 thin plate regression spline basis functions
#' @param cov.opts Parameters for the covariance function.  See additional
#'                 documentation for more details
#' @param i.Sig Alternative to specifying coords, cov.type, and cov.opts allows
#'              specification of an inverse covariance matrix.  If all are not 
#'              specified assumes independence
#' @param m Number of covariates to randomly sample for each split of
#'          the spatially adjusted regression tree. Applicable only for use in 
#'          Random Spatial Forests
#' @param r Minimum Leaf Size
#' @param max.depth Maximum number of splits to construct (if applicable)
#' @param lambda Penalty parameter for L. Applicable only is multiple domains 
#'               are used.
#' @param pen.fix logical indicating whether or not to penalize the fixed 
#'                effects in addition to the tree estimates.  Applicable only if
#'                multiple domains are used.
#' @param L Penalty matrix encouraging similarity among estimates in different 
#'          domains.  Applicable only is multiple domains are used.
#'
#' @return
#' @export
#' @useDynLib spatRF C_splits
#' @useDynLib spatRF C_splits1

spatTree <- function(Y, X, domain=rep(1,length(Y)), X.fix = NULL, coords=NULL, 
                     cov.type=NULL, cov.opts=NULL, i.Sig = NULL, 
                     m = ceiling(dim(X)[2]/3), r = 5, max.depth = NULL, 
                     lambda = 0, pen.fix=TRUE, 
                     L=matrix(0, length(unique(domain)),
                              length(unique(domain))) ){
  ## Error Checking
  if(length(Y) != dim(X)[1]) {
    stop("Number of Observations in X and Y must match")
  }
  if(length(Y) != length(domain)) {
    stop("Arguments Y and domain have different lengths")
  }
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  
  ## Sort observations into their respective domains
  ds <- sort(unique(domain))
  d <- length(ds)
  num_ds <- as.numeric(table(domain)[which(names(table(domain))==ds)])
  sorted.ind <- order(domain,1:N)
  
  ## Re-order X and Y to maintain block structure
  X <- X[sorted.ind,,drop=FALSE]
  Y <- Y[sorted.ind]
  coords <- coords[sorted.ind,,drop=FALSE]
  domain <- domain[sorted.ind]
  
  ## Check if there are any fixed covariates for initial adjustment
  if(!is.null(X.fix)){
    num.fix <- dim(X.fix)[2]
    X.fix.old <- X.fix[sorted.ind,,drop=FALSE]
    X.fix <- matrix(0,N,d*num.fix)
    strt.pen <- vector(mode = "list", length = num.fix+1)
    for(i in 1:num.fix){
      if(pen.fix){
        strt.pen[[i]] <- L
      }else{
        strt.pen[[i]] <- matrix(0,d,d)
      }
      X.fix[cbind(1:N,d*(i-1)+match(domain,ds))] <- X.fix.old[,i]
    }
    strt.pen[[i+1]] <- matrix(0,d,d)
  }
  
  
  strt.ind <- cumsum(c(1,num_ds))
  stop.ind <- cumsum(num_ds)
  
  if(dim(X)[2] != 1){
    sortind <- t(t(matrix(order(rep(1:p,each=N),c(X)),N,p)) - c(0:(p-1)*N))
  } else{
    sortind <- t(t(matrix(order(rep(1:p,each=N),c(X)),N,p)))
  }
  
  if(is.null(max.depth)) max.depth <- floor(max(num_ds)/r)
  
  
  if(is.null(i.Sig)){
    block_sig <- vector(mode = "list", length = d)
    
    for(i in 1:d){
      if(!is.null(coords)){
        domain.coords <- coords[sorted.ind[strt.ind[i]:stop.ind[i]],,drop=FALSE]
        obj <- .makeSpatBas(pars=cov.opts$pars,
                            coords.train = domain.coords,
                            cov.type=cov.type,
                            cov.opts=cov.opts)
        block_sig[[i]] <- obj$i.sig
      } else{
        if(i==1) warning("No coordinates input - assuming independence")
        block_sig[[i]] <- diag(num_ds[i])
      }
    }
    i.Sig <- as.matrix(Matrix::bdiag(block_sig))
  }
  
  Ck <- matrix(0,N,d*(max.depth+1))
  for(i in 1:d){
    Ck[strt.ind[i]:stop.ind[i],i] <- 1
  }
  
  if(is.null(X.fix)){
    b <- NULL
  } else{
    p.fix <- dim(X.fix)[2]
    b <- rep(NA,p.fix)
  }
  
  
  
  ## First, Start with the intercept estimate
  ## Calculate starting kernel estimate
  if(!is.null(X.fix)) {
    X.start <- cbind(X.fix,Ck[,1:d])
    Kern.C <- try(i.Sig - i.Sig %*% X.start %*% 
                    solve( as.matrix(Matrix::bdiag(strt.pen)) + t(X.start) %*% i.Sig %*% X.start,
                           t(X.start) %*% i.Sig),silent=TRUE)
    if(class(Kern.C)=="try-error"){
      Kern.C <- i.Sig - i.Sig %*% X.start %*% 
        .fastSolve( as.matrix(Matrix::bdiag(strt.pen)) + t(X.start) %*% i.Sig %*% X.start) %*%
        t(X.start) %*% i.Sig
    }
  } else {
    X.start <- Ck[,1:d]
    Kern.C <- i.Sig - i.Sig %*% X.start %*% 
      solve( t(X.start) %*% i.Sig %*% X.start,
             t(X.start) %*% i.Sig)
  }
  
  loss.full <- t(Y) %*% Kern.C %*% Y
  omegaY <- as.numeric(Kern.C %*% Y)
  
  L.full <- vector(mode = "list", length = max.depth)
  covar.num <- split.val <- branch.loss <- leaf.split <- rep(NA,max.depth)
  num.per.leaf <- rep(NA,max.depth)
  num.per.leaf[1] <- N
  num.per.group <- matrix(NA,d,(max.depth+1))
  num.per.group[,1] <- num_ds
  leaf.mem <- leaf.mem.use<-rep(1,N)
  split.leaf <- a <- 1
  
  while( max(num.per.group[,c(1:a)]) > (2 * r - 1)  & a <= max.depth ) {
    
    smpl <- sample(1:p,m)
    
     if(d==1){
       split_C<- .Call("C_splits1",as.numeric(X),
                                   as.integer(num.per.leaf[split.leaf]),
                                   as.integer(1:length(split.leaf)),
                                   as.integer(leaf.mem.use),
                                   as.integer(smpl),
                                   as.integer(r),
                                   as.integer(sortind),
                                   as.numeric(Kern.C),
                                   as.numeric(omegaY)
       ) 
     } else{
      split_C<- .Call("C_splits",as.numeric(X),
                                  as.integer(num.per.leaf[split.leaf]),
                                  as.integer(1:length(split.leaf)),
                                  as.integer(leaf.mem.use),
                                  as.integer(smpl),
                                  as.integer(r),
                                  as.integer(sortind),
                                  as.numeric(Kern.C),
                                  as.numeric(omegaY),
                                  as.numeric(L),
                                  as.integer(num.per.group[,split.leaf]),
                                  as.integer(domain)
      ) 
    }
    
    tot.splits <- split_C[[4]]
    if(tot.splits!=0){
      
      if(max(split_C[[5]][1:tot.splits],na.rm=T) <= 0){
        break
      } else{
        ind <- which.max( split_C[[5]][1:tot.splits] )
      }
      
      covar.num[a] <- split_C[[2]][ind]
      split.val[a] <- split_C[[1]][ind]
      leaf.split[a] <- split.leaf[split_C[[3]][ind]]
      branch.loss[a] <- split_C[[5]][ind]
      
      Ca <-  as.numeric(as.numeric(X[,covar.num[a]]) <= split.val[a] & leaf.mem==leaf.split[a])
      num.per.group[,a+1] <- as.matrix(table(domain,Ca))[,2]
      num.per.group[,leaf.split[a]] <- num.per.group[,leaf.split[a]]- num.per.group[,a+1]
      
      if(any(apply(num.per.group[,1:a+1,drop=FALSE],2,max) < 5)){
        num.per.group
      }
      
      num.per.leaf[a+1] <- sum(Ca)
      num.per.leaf[leaf.split[a]] <- num.per.leaf[leaf.split[a]]-num.per.leaf[a+1]
      leaf.mem[which(Ca==1)] <- a+1
      if(d>1){
        split.leaf <- unique(which(num.per.group[,1:(a+1)] > (2*r-1),arr.ind=T)[,2])
      } else{
        split.leaf <- unique(which(num.per.group[,1:(a+1)] > (2*r-1)))
      }
      leaf.mem.use <- match(leaf.mem,split.leaf)
      leaf.mem.use <- ifelse(is.na(leaf.mem.use),0,leaf.mem.use)
      Ck[cbind(which(Ca==1),a*d+domain[which(Ca==1)])] <- 1  
      
      ## Update Kern.C and OmegaY
      if (d > 1){
        
        rmv <- which(colSums(Ck[,(a*d+1):((a+1)*d),drop=FALSE])==0)
        d.rmv <- length(rmv)
        if(length(rmv)>0){
          C.add <- Ck[,(a*d+1):((a+1)*d),drop=FALSE][,-rmv,drop=FALSE]
          L.add <- -L[1,2]*((d-d.rmv) *diag((d-d.rmv)) - rep(1,d-d.rmv) %*% t(rep(1,d-d.rmv))) #try(L[-rmv,,drop=FALSE][,-rmv,drop=FALSE])
        } else{
          C.add <- Ck[,(a*d+1):((a+1)*d),drop=FALSE]
          L.add <- L
        }
      } else{
        C.add <- Ck[,a+1]
        L.add <- 0
      }
      v <- crossprod(Kern.C,C.add)
      Kern.C <- Kern.C - v %*% solve(t(v) %*% C.add +L.add, t(v))
      omegaY <- crossprod(Kern.C,Y)
      
      L.full[[a]] <- L.add
      
      a <- a + 1
      
    } else{
      break
    }
  }
  
  Ck <- Ck[,1:(a*d),drop=FALSE]
  
  if(is.null(X.fix)){
    if(d ==1){
      muk <- solve(t(Ck) %*% i.Sig %*% Ck ,t(Ck) %*% i.Sig %*% Y)
    } else{
      muk <- try(solve(t(Ck) %*% i.Sig %*% Ck +Matrix::bdiag(append(list(matrix(0,d,d)),rep(list(L),a-1))),t(Ck) %*% i.Sig %*% Y))
      
      if(class(muk)=="try-error"){
        rmv <- which(colSums(Ck)==0)
        muk <- rep(NA,a*d)
        Ck.use <- Ck[,-rmv]
        muk[-rmv] <- (solve(t(Ck.use) %*% i.Sig %*% Ck.use + Matrix::bdiag(append(list(matrix(0,d,d)),L.full[1:a-1])),t(Ck.use) %*% i.Sig %*% Y))
        for(i in 1:rmv){
          muk[i] <- mean(muk[(floor(i/d)*d+1):(ceiling(i/d)*d)])
        }
      }
    }
  } else{
    X.out <- cbind(X.fix,Ck)
    muk <- try(solve(t(X.out) %*% i.Sig %*% X.out + 
                       Matrix::bdiag(append(strt.pen,rep(list(L),a-1))),
                     t(X.out) %*% i.Sig %*% Y),silent=TRUE)
    if(class(muk)=="try-error"){
      muk <- try(.fastSolve(t(X.out) %*% i.Sig %*% X.out + 
                              Matrix::bdiag(append(strt.pen,rep(list(L),a-1)))) %*%
                   t(X.out) %*% i.Sig %*% Y)
    }
    if(class(muk)=="try-error"){
      if(any(colSums(X.out)==0)){
        rmv <- which(colSums(X.out)==0)
        Ck.use <- cbind(X.fix,Ck)[,-rmv]
        muk <- rep(NA,dim(X.out)[2])
        muk.nz <- try(solve(t(Ck.use) %*% i.Sig %*% Ck.use + as.matrix(Matrix::bdiag(append(strt.pen,L.full[1:a-1])))),
                      t(Ck.use) %*% i.Sig %*% Y)
        if(class(muk.nz)=="try-error"){
          muk.nz <- .fastSolve(t(Ck.use) %*% i.Sig %*% Ck.use + 
                                 as.matrix(Matrix::bdiag(append(strt.pen,L.full[1:a-1])))) %*%
            t(Ck.use) %*% i.Sig %*% Y
        }
        muk[-rmv] <- muk.nz
        
        for(i in 1:rmv){
          muk[i] <- mean(muk[(floor(i/d)*d+1):(ceiling(i/d)*d)])
        }
      } else{
        stop("Problem Calculating muk for current design matrix")
      }
    } 
    
    b <- muk[1:(num.fix*d)]
    muk <- muk[(num.fix*d+1):length(muk)]
  }
  
  
  ## Return the tree 
  return(list(covar=covar.num[1:(a-1)],split=split.val[1:(a-1)],
              leaf=leaf.split[1:(a-1)],loss=branch.loss[1:(a-1)],
              Ck=Ck,muk=muk,X=X.fix,b=b,max.loss = loss.full))
}