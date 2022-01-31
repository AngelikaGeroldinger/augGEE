# Supplementary file 2: R functions for single-step and iterated augmented GEE
# Author: Angelika Geroldinger
# Date: 28 January 2022


library(logistf)
library(mmmgee) # this package uses "scale" weights as in proc gee and proc genmod in SAS;
library(Matrix) # to create block diagonal matrix R; 
library(expm) # to get the square root of a (positive definite) matrix; 


# The function getSW computes the Sandwich variance covariance matrix according to the formula in Morel JG, Bokossa MC and Neerchal NK. Small sample correction for the variance of GEE estimators. Biometrical J 2003;
# note that already the uncorrected sandwich matrix is different from the one defined in Molenberghs G and Verbeke G. Models for Discrete Longitudinal Data. New York: Springer, 2006. ; 
# getSW is used in the augGEE functions to compute the sandwich estimate using only the original observations, not the pseudo observations;
getSW <- function (YY, XX, eta, R.alpha.inv, clusterid){
  mu <- (1+exp(-eta))^-1
  w <- as.vector(mu*(1-mu))
  hess <- Matrix::crossprod(diag(sqrt(w)) %*% XX, R.alpha.inv %*% diag(sqrt(w)) %*% XX)
  hessinv <- solve(as.matrix(hess))
  n <- length(unique(clusterid)) # number of clusters
  d <- sapply(unique(clusterid), function(x) t(YY[clusterid==x]- mu[clusterid==x]) %*% Diagonal(x=1/sqrt(w[clusterid==x])) %*% R.alpha.inv[clusterid==x,clusterid==x, drop=FALSE] %*% Diagonal(x=sqrt(w[clusterid==x])) %*% XX[clusterid==x,, drop=FALSE])
  d_bar <- Reduce("+", d)/n  
  I1 <- sapply(d, function(x) t(x-d_bar) %*% (x-d_bar))
  I1 <- Reduce('+', I1)
  sw <- I1 %*% hessinv
  sw <- hessinv %*% sw 
  # apply the small sample correction by Morel;
  k <- ncol(XX) # number of parameters (including intercept)
  nstar <- nrow(XX) # number of observations
  mult <- ((nstar-1)/(nstar-k)) * (n/(n-1))
  phi <- max(1, sum(diag(mult * as.matrix(hessinv %*% I1))) /k)
  delta <- min(0.5, k/(n-k))
  swcor <- as.matrix( mult* sw + delta * phi * hessinv)
  return(list(sw=sw, swcor=swcor, hess=hess))
}


# the cluster-id is expected to be stored in a column named "id" in data;
augGEE <- function(form, data, corstr, initGEE=NULL, tol=0.001){
  nobs <- nrow(data)
  sample_size <- length(unique(data$id))
  resfl <- logistf(formula=form, data=data, firth=TRUE, pl=FALSE)  # set pl=FALSE for faster computation; 
  y <- data[, all.vars(form[[2]])]
  data_aug <- data.frame(rbind(data, data, data), pseudo=c(rep(0,nobs), rep(1, 2*nobs)))
  data_aug[,all.vars(form[[2]])] <- c(y, y, 1-y)
  data_aug$neww <- c(rep(1,nobs), 0.5* resfl$hat.diag, 0.5* resfl$hat.diag)
  data_aug$id <- c(data$id, data$id + max(data$id), data$id + 2*max(data$id)) 
  data_aug <- data_aug[order(data_aug$id),] #geeglm requires data to be ordered by clusters! 
  if (is.null(initGEE)){
    initGEE <- resfl$coefficients
  }
  res <- geem2(formula=form, family = binomial("logit"), id = id, data = data_aug, weights=neww, init.beta=initGEE, 
               scale.fix=TRUE, init.phi=1, corstr=corstr, useP=FALSE, maxit=30, tol=tol, 
               conv.criterion="difference") 
  ind <- 1:nrow(data)
  varorig_list <- getSW(YY=res$sandw.args$Y[ind], XX=res$sandw.args$X[ind,], eta=res$sandw.args$eta[ind], 
                        R.alpha.inv=res$sandw.args$R.alpha.inv[ind, ind], clusterid=res$sandw.args$id[ind]) # function getSW defined in augGEEit.R, small-sample corrected sandwich estimator;  
  varorig <- varorig_list$sw
  varorigcor <- varorig_list$swcor
  return(list(res=res, varorig=varorig, varorigcor=varorigcor, data_aug=data_aug, form=form))
}


# currently only exchangeable correlation structure implemented in the iterated augmented GEE algorithm; 
# the cluster-id is expected to be stored in a column named "id" in data;
augGEEit <- function(form,  data, corstr, tol=0.001, maxit=20){
  data <- data[order(data$id),]
  nobs <- nrow(data)
  nclust <- table(data$id)
  sample_size <- length(unique(data$id))
  y <- data[, all.vars(form[[2]])]
  # use coef from independent fit as starting values; 
  resfl <- logistf(formula=form, data=data, firth=TRUE, pl=FALSE)  # set pl=FALSE for faster computation; 
  if ("(Intercept)" %in% names(resfl$coefficients)){
    X <- cbind(1,as.matrix(data[,all.vars(form[[3]])]))
  } else  X <- as.matrix(data[,all.vars(form[[3]])])
  pi<- resfl$predict
  correst <- 0
  betachange <- 1
  beta <- resfl$coefficients
  nit <- 0
  
  while(max(abs(betachange)) > tol & nit <= maxit){
    Wroot <- sqrt(pi * (1-pi))
    Wroot <- split(x=Wroot, f=data$id)
    Rinv <- lapply(nclust, function(x) solve( toeplitz(c(1, rep(correst, x-1))) ) )
    omega <- mapply(function(x,y) Diagonal(x=x) %*% y %*% Diagonal(x=x), Wroot, Rinv)
    omegaroot <- lapply(omega, function(x) sqrtm(x))
    I <- (t(X) %*% bdiag(omega) %*% X)
    Iinv <- solve(I)
    hat <- bdiag(omegaroot) %*% X %*% Iinv %*% t(X) %*% bdiag(omegaroot)
    data_aug <- data.frame(rbind(data, data, data), pseudo=c(rep(0,nobs), rep(1, 2*nobs)))
    data_aug[,all.vars(form[[2]])] <- c(y, y, 1-y)
    data_aug$neww <- c(rep(1,nobs), 0.5* diag(hat), 0.5* diag(hat))
    data_aug$id <- c(data$id, data$id + max(data$id), data$id + 2*max(data$id)) 
    #data_aug <- data_aug[order(data_aug$id),] 
    res <- geem2(formula=form, family = binomial("logit"), id = id, data = data_aug, weights=neww, scale.fix=TRUE, 
                 corstr=corstr, useP=FALSE, maxit=30, tol=1e-03, init.beta=beta, conv.criterion="difference")  
    pi <- (1+ exp(-res$eta[1:nobs]))^(-1)
    correst <- res$alpha
    betachange <- res$beta - beta
    beta <- res$beta
    nit <- nit +1 
    if(res$converge==FALSE){
      nit <- maxit + 1
    }
  }
  ind <- 1:nrow(data)
  varorig_list <- getSW(YY=res$sandw.args$Y[ind], XX=res$sandw.args$X[ind,], eta=res$sandw.args$eta[ind], 
                        R.alpha.inv=res$sandw.args$R.alpha.inv[ind, ind], clusterid=res$sandw.args$id[ind])
  varorig <- varorig_list$sw
  varorigcor <- varorig_list$swcor
  if(max(abs(betachange)) > tol){
    warning("Did not converge.")
  }
  return(list(res=res, varorig=varorig, varorigcor=varorigcor, data_aug=data_aug, form=form, nit=nit))
}


#############################################################################################
### Example
#############################################################################################

data(datasim) # from package mmmgee;

fit1 <- augGEE(form=Y.bin~gr.lang+x3, data=datasim, corstr="exchangeable", initGEE=NULL)
fit1$res

fit2 <- augGEEit(form=Y.bin~gr.lang+x3, data=datasim, corstr="exchangeable")
fit2$res

