utils::globalVariables(c("num_ds"))

#' Title
#'
#' @param pars 
#' @param Y Observed values
#' @param X Covariates used to construct the spatially adjusted regression tree
#' @param domain Used to construct multiple regression trees simultaneously. For
#'               details, see Hee Wai et al., 2020
#' @param coords x and y coordinates for locations of the points
#' @param group 
#' @param X.fix 
#' @param cov.type 
#' @param cov.opts 
#' @param X.test 
#' @param X.fix.test 
#' @param test.coords 
#' @param test.domain 
#' @param t 
#' @param m 
#' @param r 
#' @param k 
#' @param replace 
#' @param s 
#' @param var.imp 
#' @param imp.msr 
#' @param pen.fix 
#' @param L 
#' @param scale 
#' @param ... 
#'
#' @return A list of result. With...
#' @export

spatRF <- function( pars, Y, X, domain, coords, group = 1:length(Y),
                    X.fix = NULL, cov.type = "exp", cov.opts = NULL, 
                    X.test = NULL, X.fix.test = NULL, test.coords = NULL, 
                    test.domain = NULL, t = 500, m = floor(dim(X)[2]/3),  r = 5, 
                    k = floor(length(Y)/r), replace = TRUE,
                    s = ifelse (replace,length(Y), 
                                ceiling(.632*length(Y))), 
                    var.imp=FALSE,imp.msr = c("perm","imp"), pen.fix=TRUE,
                    L=matrix(0,length(unique(domain)),length(unique(domain))), 
                    scale=FALSE,... ){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  ds <- unique(domain)
  d <- length(ds)
  
  use.pars <- vector(mode = "list", length = d)
  if(length(pars)==1){
    if(!is.list(pars)){
      for(i in 1:d){
        use.pars[[i]]<- list(psill=pars,nugget=1-pars)
      }
    } else{
      for(i in 1:d){
        use.pars[[i]]<- pars
      }
    }
  } else if (!is.list(pars)){
    for(i in 1:d){
      use.pars[[i]]<- list(psill=pars[i],nugget=1-pars[i])
    }
  }
  
  use.cov.opts <- vector(mode = "list", length = d)
  if(length(cov.opts)==1){
    for(i in 1:d){
      use.cov.opts[[i]]<- cov.opts[i]
    }
  } else{
    use.cov.opts <- cov.opts
  }
  
  ### Sort by domain
  domain <- match(domain,ds)
  
  ordr <- order(domain)
  Y <- Y[ordr]
  X <- X[ordr,,drop=FALSE]
  coords <- coords[ordr,,drop=FALSE]
  domain <- domain[ordr]
  group <- group[ordr]
  X.fix.train <- X.fix.vldt <- NULL
  
  if(!is.null(X.fix)){
    X.fix <- X.fix[ordr,,drop=FALSE]
  }
  
  times.selected <- times.trained <- rep(0,n)
  Yoob <- foob <- zoob <- rep(0,n)
  
  if(cov.type=="TPRS" ){
    X.fix <- cbind(X.fix,as.matrix(coords))
    X.fix.test <- cbind(X.fix.test,as.matrix(test.coords))
  } 
  
  if(!is.null(X.test)) {
    n.test <- dim(X.test)[1]
    Ytest <- ftest <- ztest <- rep(0,n.test)
    test.domain <- match(test.domain,ds)
    ordr.test <- order(test.domain)
    X.test <- X.test[ordr.test,,drop=FALSE]
    if(!is.null(X.fix.test)){
      X.fix.test <- X.fix.test[ordr.test,,drop=FALSE]
    }
    test.coords <- test.coords[ordr.test,,drop=FALSE]
    test.domain <- test.domain[ordr.test]
  }
  
  if(var.imp){
    var.imp.msr <- rep(0,p)
    imp.msr <- imp.msr[1]
    if(imp.msr=="perm"){
      perm.ind <- replicate(p,sample(1:n))
      perm.mat <- matrix(0,nrow=n,ncol=p)
    }
  }
  
  ## Build t Spatial trees
  for(i in 1:t){
    
    ## For each tree, subsample the data into training and test
    vldt.sample <- replace
    smpl <- .sampleBlock(domain,group,replace=replace)
    train <- NULL
    for(i in 1:length(smpl)){
      train <- c(train,which(group==smpl[i]))
    }
    #train <- which(group %in% smpl)
    vldt <- (1:n)[-which(1:n %in% train)]
    Y.train <- Y[train]
    X.train <- X[train,,drop=FALSE]
    domain.train <- domain[train]
    domain.vldt <- domain[vldt]
    coords.train <- coords[train,,drop=FALSE]
    coords.vldt <- coords[vldt,,drop=FALSE]
    
    if(!is.null(X.fix)){
      X.fix.train <- X.fix[train,,drop=FALSE]
      X.fix.vldt <- X.fix[vldt,,drop=FALSE]
    }
    
    ## If TPRS, generate the fixed basis functions
    if(cov.type=="TPRS") {
      X.fix.train <- cbind(X.fix.train,as.matrix(coords[train,,drop=FALSE]))
      X.fix.vldt <- cbind(X.fix.vldt,as.matrix(coords[vldt,,drop=FALSE]))
    } 
    
    ## Generate the training and test covariance matrices
    block_sig <- vector(mode = "list", length = d)
    block_vldt <- vector(mode="list",length=d)
    if(!is.null(X.test)) block_test <- vector(mode="list",length=d)
    
    for(i in 1:d){
      if(!is.null(coords)){
        domain.coords <- coords.train[which(domain.train==i),,drop=FALSE]
        domain.vldt.coords <- coords.vldt[which(domain.vldt==i),,drop=FALSE]
        obj <- .makeSpatBas(pars=use.pars[[i]],
                            coords.train = domain.coords, coords.vldt=domain.vldt.coords,
                            coords.test = test.coords[which(test.domain==i),,drop=FALSE],cov.type=cov.type,
                            cov.opts=use.cov.opts[[i]])
        if(scale){
          block_sig[[i]] <- obj$i.sig / dim(domain.coords)[1]
          if(!is.null(X.test)) block_test[[i]] <- obj$sig.test * dim(domain.coords)[1]
          block_vldt[[i]] <- obj$sig.vldt * dim(domain.coords)[1]
        } else{
          block_sig[[i]] <- obj$i.sig 
          if(!is.null(X.test)) block_test[[i]] <- obj$sig.test
          block_vldt[[i]] <- obj$sig.vldt
        }
      } else{
        if(i==1) warning("No coordinates input - assuming independence")
        block_sig[[i]] <- diag(num_ds[i])
      }
    }
    i.Sig <- as.matrix(Matrix::bdiag(block_sig))
    sig.vldt <- as.matrix(Matrix::bdiag(block_vldt))
    if(!is.null(X.test)) sig.test <- as.matrix(Matrix::bdiag(block_test))
    
    ## Run the tree
    mod <- spatTree(Y=Y.train,X=X.train,domain=domain.train,X.fix = X.fix.train,
                                i.Sig =i.Sig,r=r,max.depth=k,pen.fix=pen.fix,L=L)
    
    ## Compile predictions on held out samples
    C.vldt <- .makeCmat(mod,X[vldt,,drop=FALSE],domain=domain.vldt)
    times.selected[vldt] <- times.selected[vldt] + 1
    
    if(dim(C.vldt)[2] != length(mod$muk)){
      C.vldt <- C.vldt[,1:length(mod$muk)]
    }
    
    if(!is.null(X.fix)){
      foob[vldt] <- foob[vldt] + C.vldt%*% mod$muk + .vecToBlock(X.fix.vldt,domain.vldt) %*% mod$b 
      zoob[vldt] <- zoob[vldt] +
        sig.vldt %*% i.Sig%*%
        (Y.train - cbind(mod$Ck,mod$X) %*% c(mod$muk,mod$b))
    } else{
      foob[vldt] <- foob[vldt] + C.vldt %*% mod$muk
      zoob[vldt] <- zoob[vldt] + sig.vldt %*% i.Sig %*% 
        (Y.train - mod$Ck %*% mod$muk )
    }
    
    ## Generate predictions at held out sites using the current tree
    if( !is.null(X.test) ) {
      C.test <- .makeCmat(mod,X.test,domain=test.domain)
      if(dim(C.test)[2] != length(mod$muk)){
        C.test <- C.test[,1:length(mod$muk)]
      }
      
      if(!is.null(X.fix.test)){
        
        num.fix <- dim(X.fix.test)[2]
        
        ftest <- ftest + C.test %*% mod$muk + .vecToBlock(X.fix.test,test.domain) %*% mod$b 
        ztest <- ztest + sig.test %*% i.Sig %*% 
          (Y.train - cbind(mod$Ck,mod$X) %*% c(mod$muk,mod$b))
      } else{
        ftest <- ftest + C.test %*% mod$muk
        ztest <- ztest + sig.test %*% i.Sig %*% 
          (Y.train - mod$Ck %*% mod$muk) 
      }
    }
    
    ## Generate variable importance measures if requested
    if(var.imp){
      if(imp.msr=="imp"){
        loss.explain <- rowsum(mod$loss/c(mod$max.loss),mod$covar)
        var.imp.msr[as.numeric(rownames(loss.explain))] <- 
          var.imp.msr[as.numeric(rownames(loss.explain))] + c(loss.explain)
      } else if(imp.msr=="perm"){
        if(!all(is.na(mod$covar))){
          for (w in unique(mod$covar)){
            X.perm <- X
            X.perm[,w] <- X[perm.ind[,w],w]
            C.perm <- .makeCmat(mod,X.perm[vldt,,drop=FALSE],domain=domain.vldt)
            if(!is.null(X.fix)){
              perm.mat[vldt,w] <- perm.mat[vldt,w] + C.perm%*% mod$muk + .vecToBlock(X.fix.vldt,domain.vldt) %*% mod$b 
            } else{
              perm.mat[vldt,w] <- perm.mat[vldt,w] + C.perm%*% mod$muk 
            }
          }
          notused <- (1:p)[-unique(mod$covar)]
          if(!is.null(X.fix)){
            perm.mat[vldt,notused] <- perm.mat[vldt,notused] + c(C.vldt%*% mod$muk + .vecToBlock(X.fix.vldt,domain.vldt) %*% mod$b )
          } else{
            perm.mat[vldt,notused] <- perm.mat[vldt,notused] + c(C.vldt%*% mod$muk)
          }
        } else{
          if(!is.null(X.fix)){
            perm.mat[vldt,] <- perm.mat[vldt,] + c(C.vldt%*% mod$muk + .vecToBlock(X.fix.vldt,domain.vldt) %*% mod$b )
          } else{
            perm.mat[vldt,] <- perm.mat[vldt,] + c(C.vldt%*% mod$muk)
          }
        }
      }
      
    }
  }
  
  
  out <-list(fpredicted=foob/times.selected[match(1:n,ordr)],zpredicted = zoob/times.selected[match(1:n,ordr)])
  
  if(!is.null(X.test)) {
    out$ftest <- (ftest/t)[match(1:n.test,ordr.test)]
    out$ztest <- (ztest/t)[match(1:n.test,ordr.test)]
  }
  
  if(var.imp){
    if(imp.msr=="perm"){
      for(w in 1:p){
        var.imp.msr[w] <- (sum((Y-(perm.mat[,w]+zoob)/times.selected)^2)/n - 
                             sum((Y-(foob + zoob)/times.selected)^2)/n) /
          sum((Y-(foob + zoob)/times.selected)^2)/n
      }
    }
    out$var.imp <- var.imp.msr
  }
  
  return(out)
}