# R script for ICA of Pacific and Indian Ocean deep-sea sediments.
# Originally written by K.Yasukawa on August 22, 2014. (Modified on August 12, 2016)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

numIC <- 9   # Number of ICs

# Data loading
path <- getwd()

datapath <- list.files(path = path, pattern = "csv")
ALLDATA <-
  read.csv("762C_only_NEW.csv",
           stringsAsFactors = F,
           fileEncoding = "utf-8")
REY <- ALLDATA[, colnames(ALLDATA) == "SREY"]
ALLDATA <- ALLDATA[, colnames(ALLDATA) != "SREY"]
X <- ALLDATA

index_elements_start <- 7 # index of the first element in the dataset (Mg) of X Data (without Na)
index_elements_end <- 46 # index of the last element in the dataset (CaCO3) of X

Y <-
  X[, index_elements_start:index_elements_end]		# For calculation, numerical data matrix Y is defined here

# Standardization of the data

mean <- apply(Y, 2, mean)
stdev <- apply(Y, 2, sd)
YA <- Y
Y0 <- t(Y)
Ys <- t((Y0 - mean) / stdev)
Y <- Ys

# Estimation of the number of ICs with PCA
#pca and visualization
pca <- function(df) {
  df <- na.omit(df)
  ans <-
    prcomp(df, scale = TRUE)	# prcomp is a function to implement PCA, which is originally incorporated in R
  ans$loadings <- t(t(ans$rotation) * ans$sdev)
  ans$eigenvalues <- ans$sdev ^ 2
  invisible(ans)
}

par(mfrow = c(1, 1))
Z <- pca(Y)
summary(Z)	# Display the summarized PCA result
screeplot(Z)	# Draw the bar chart of eigenvalues (variances)
summary(Y)	# Display the statistical summary of the original data
pcaResult <- summary(Z)
write.csv(pcaResult$importance,paste(path, "/Result_762C/ICnum=",numIC,"/Result/pcaResult.txt",sep=""))
# Variable Parameters

numDim <-
  index_elements_end - index_elements_start + 1        # Number of sample variables
numData <- unlist(length(X$Site))
#numEM <-  unlist(length(em$Endmember))  # Number of data in Dataset of reference materials


# Independent Component Analysis with fastICA algorithm (Hyvarinen et al., 2001)


# === FastICA function modified from R Package "fastICA" ver 1.1-16 by Marchini & Heaton (2012) ===
fastICA2 <- function(X, 
                     n.comp, 
                     alg.typ = c("parallel", "deflation"), 
                     fun = c("kurtosis", "skewness", "logcosh", "exp"),
                     a1 = 1, 
                     a2 = 1, 
                     maxit = 200,
                     tol = 1e-04, 
                     verbose = TRUE, 
                     w.init = NULL){
  X <- as.matrix(X)
  XT <- t(X)
  
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  if(n.comp>min(n,p)){
    message("'n.comp' is too large: reset to ", min(n,p))
    n.comp <- min(n,p)
  }
  
  if(is.null(w.init))
    w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  else{
    if(!is.matrix(w.init) || length(w.init) != (n.comp^2))
      stop("w.init is not a matrix or is the wrong size.")
  }
  
  if(verbose)
    message("Centering")
  mean <- apply(X, 2, mean)
  X0 <- t(X)
  Xc <- t(X0-mean)
  XcT <- t(Xc)                       # Q-mode-like calculation: transpose of input data matrix X
  
  if(verbose)
    message("Whitening")
  V <- XcT %*% Xc/n
  s <- La.svd(V)
  D <- diag(c(1/sqrt(s$d)))
  KT <- D %*% t(s$u)
  KT <- matrix(KT[1:n.comp,], n.comp, p)
  ZT <- KT %*% XcT
  
  if(alg.typ == "deflation")
    WT <- as.matrix(
      ica2.def(
        ZT, 
        n.comp, 
        tol = tol, 
        fun = fun, 
        a1 = a1, 
        a2 = a2, 
        maxit = maxit, 
        verbose = verbose, 
        w.init = w.init))
  
  else if(alg.typ == "parallel")
    WT <- as.matrix(
      ica2.par(
        ZT, 
        n.comp, 
        tol = tol, 
        fun = fun, 
        a1 = a1, 
        a2 = a2, 
        maxit = maxit, 
        verbose = verbose, 
        w.init = w.init))
  
  w <- WT %*% KT
  ST <- w %*% XcT
  AT <- t(w) %*% solve(w %*% t(w))
  
  return(list(X = t(XT), K = t(KT), W = t(WT), A = t(AT), S = t(ST)))
}

#  fastICA algorithm using a deflation scheme for orthogonalization
ica2.def <- function(X, n.comp, tol, fun, a1, a2, maxit, verbose, w.init){
  
  if(verbose && fun ==  "kurtosis")
    message("Deflation fastICA using kurtosis approximation to non-gaussianity")
  if(verbose && fun == "skewness")
    message("Deflation fastICA using skewness approximation to non-gaussianity")
  if(verbose && fun ==  "logcosh")
    message("Deflation fastICA using logcosh approximation to neg-entropy")
  if(verbose && fun == "exp")
    message("Deflation fastICA using exponential approximation to neg-entropy")
  
  n <- nrow(X)
  p <- ncol(X)
  W <- matrix(0, n.comp, n.comp)
  
  for(i in 1:n.comp){
    if(verbose)
      message("Component ", i)
    w <- matrix(w.init[i,], n.comp, 1)    # w is i-th row vector of matrix W
    if(i>1){
      tmp <- w
      tmp[1:length(tmp)] <- 0
      for(u in 1:(i-1)){
        k <- sum(w*W[u,])
        tmp <- tmp + k*W[u,]
      }
      w <- w - tmp
    }
    w <- w/sqrt(sum(w^2))
    
    lim <- rep(1000, maxit)
    it <- 1
    
    if(fun == "kurtosis"){
      while(lim[it] > tol && it < maxit){
        wx <- t(w) %*% X
        gwx <- wx^3
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X*gwx
        V1 <- apply(xgwx, 1, FUN = mean)
        
        gwx2 <- 3
        gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
        V2 <- mean(gwx2)*w
        
        wp <- V1 - V2
        wp <- matrix(wp, n.comp, 1)
        it <- it + 1
        
        if(i>1){
          tmp <- wp
          tmp[1:length(tmp)] <- 0
          for(u in 1:(i-1)){
            k <- sum(wp*W[u,])
            tmp <- tmp + k*W[u,]
          }
          wp <- wp - tmp
        }
        wp <- wp/sqrt(sum(wp^2))
        
        lim[it] <- Mod(Mod(sum((wp*w)))-1)
        if(verbose)
          message("Iteration ", it-1, "  tol = ", format(lim[it]))
        
        w <- matrix(wp, n.comp, 1)
      }
    }
    
    if(fun == "skewness"){
      while(lim[it] > tol && it < maxit){
        wx <- t(w) %*% X
        gwx <- wx^2
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X*gwx
        V1 <- apply(xgwx, 1, FUN = mean)
        
        gwx2 <- 0
        gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
        V2 <- mean(gwx2)*w
        
        wp <- V1 - V2
        wp <- matrix(wp, n.comp, 1)
        it <- it + 1
        
        if(i>1){
          tmp <- wp
          tmp[1:length(tmp)] <- 0
          for(u in 1:(i-1)){
            k <- sum(wp*W[u,])
            tmp <- tmp + k*W[u,]
          }
          wp <- wp - tmp
        }
        wp <- wp/sqrt(sum(wp^2))
        
        lim[it] <- Mod(Mod(sum((wp*w)))-1)
        if(verbose)
          message("Iteration ", it-1, "  tol = ", format(lim[it]))
        
        w <- matrix(wp, n.comp, 1)      
      }
    }
    
    if(fun == "logcosh"){
      while(lim[it] > tol && it < maxit){
        wx <- t(w) %*% X
        gwx <- tanh(a1*wx)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X*gwx
        V1 <- apply(xgwx, 1, FUN = mean)
        
        gwx2 <- a1*(1 - (tanh(a1*wx))^2)
        gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
        V2 <- mean(gwx2)*w
        
        wp <- V1 - V2
        wp <- matrix(wp, n.comp, 1)
        it <- it + 1
        
        if(i>1){
          tmp <- wp
          tmp[1:length(tmp)] <- 0
          for(u in 1:(i-1)){
            k <- sum(wp*W[u,])
            tmp <- tmp + k*W[u,]
          }
          wp <- wp - tmp
        }
        wp <- wp/sqrt(sum(wp^2))
        
        lim[it] <- Mod(Mod(sum((wp*w)))-1)
        if(verbose)
          message("Iteration ", it-1, "  tol = ", format(lim[it]))
        
        w <- matrix(wp, n.comp, 1)
      }
    }
    
    if(fun == "exp"){
      while(lim[it] > tol && it < maxit){
        wx <- t(w) %*% X
        gwx <- wx*exp(-a2*(wx^2)/2)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X*gwx
        V1 <- apply(xgwx, 1, FUN = mean)
        
        gwx2 <- (1 - a2*wx^2)*exp(-a2*(wx^2)/2)
        gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
        V2 <- mean(gwx2)*w
        
        wp <- V1 - V2
        wp <- matrix(wp, n.comp, 1)
        it <- it + 1
        
        if(i>1){
          tmp <- wp
          tmp[1:length(tmp)] <- 0
          for(u in 1:(i-1)){
            k <- sum(wp*W[u,])
            tmp <- tmp + k*W[u,]
          }
          wp <- wp - tmp
        }
        wp <- wp/sqrt(sum(wp^2))
        
        lim[it] <- Mod(Mod(sum((wp*w)))-1)
        if(verbose)
          message("Iteration ", it-1, "  tol = ", format(lim[it]))
        
        w <- matrix(wp, n.comp, 1)
      }
    } 
    W[i,] <- w 
  }  
  W
}

#  fastICA algorithm using a symmetric scheme for orthogonalization
ica2.par <- function(X, n.comp, tol, fun, a1, a2, maxit, verbose, w.init){
  
  Diag <- function(d) if (length(d) > 1L)
    diag(d)
  else as.matrix(d)
  
  n <- nrow(X)
  p <- ncol(X)
  
  W <- w.init
  sW <- La.svd(W)
  
  W <- sW$u %*% Diag(1/sW$d) %*% W
  Wp <- W
  
  lim <- rep(1000, maxit)
  it <- 1
  
  if(fun == "kurtosis"){
    message("Symmetric fastICA using kurtosis approximation to non-gaussianity")
    while(lim[it] > tol && it < maxit){
      wx <- W %*% X
      gwx <- wx^3
      V1 <- gwx %*% t(X)/p
      gwx2 <- matrix(3, n.comp, p)
      V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
      Wp <- V1 - V2
      sWp <- La.svd(Wp)
      Wp <- sWp$u %*% Diag(1/sWp$d) %*% t(sWp$u) %*% Wp
      lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(W))) - 1))
      
      W <- Wp
      
      if(verbose)
        message("Iteration ", it, "  tol = ", format(lim[it + 1]))
      
      it <- it + 1
    }
  }
  
  if(fun == "skewness"){
    message("Symmetric fastICA using skewness approximation to non-gaussianity")
    while(lim[it] > tol && it < maxit){
      wx <- W %*% X
      gwx <- wx^2
      V1 <- gwx %*% t(X)/p
      gwx2 <- matrix(0, n.comp, p)
      V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
      Wp <- V1 - V2
      sWp <- La.svd(Wp)
      Wp <- sWp$u %*% Diag(1/sWp$d) %*% t(sWp$u) %*% Wp
      lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(W))) - 1))
      
      W <- Wp
      
      if(verbose)
        message("Iteration ", it, "  tol = ", format(lim[it + 1]))
      
      it <- it + 1
    }
  }
  
  if(fun == "logcosh"){
    message("Symmetric fastICA using logcosh approximation to neg-entropy")
    V1 <- matrix(0, nrow=n.comp, ncol=n.comp, byrow=TRUE)
    while(lim[it] > tol && it < maxit){
      wx <- W %*% X
      gwx <- tanh(a1*wx)
      V1 <- gwx %*% t(X)/p
      gwx2 <- a1*(1 - (gwx)^2)
      V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
      Wp <- V1 - V2
      sWp <- La.svd(Wp)
      Wp <- sWp$u %*% Diag(1/sWp$d) %*% t(sWp$u) %*% Wp
      lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(W))) - 1))
      
      W <- Wp
      
      if(verbose)
        message("Iteration ", it, "  tol = ", format(lim[it + 1]))
      
      it <- it + 1
    }
  }
  
  if(fun == "exp"){
    message("Symmetric fastICA using exponential approximation to neg-entropy")
    while(lim[it] > tol && it < maxit){
      wx <- W %*% X
      gwx <- wx*exp(-a2*(wx^2)/2)
      V1 <- gwx %*% t(X)/p
      gwx2 <- (1 - a2*wx^2)*exp(-a2*(wx^2)/2)
      V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
      Wp <- V1 - V2
      sWp <- La.svd(Wp)
      Wp <- sWp$u %*% Diag(1/sWp$d) %*% t(sWp$u) %*% Wp
      lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(W))) - 1))
      
      W <- Wp
      
      if(verbose)
        message("Iteration ", it, "  tol = ", format(lim[it + 1]))
      
      it <- it + 1
    }
  }
  W
}

#

# ICA loop for robustness check

loop <- 100

Sp <- array(0, dim=c(numIC, numDim, loop))
Sn <- array(0, dim=c(numIC, numDim, loop))

numPair <- choose(numDim,2)
slope <- array(0, dim=c(numIC, numPair, loop))

ICinRealSpaceall <- array(dim=c(1,numDim,1))


for(i in 1:loop) {
  
  message("i-loop: ", i)
  
  
  Y1 <- Y		# Repetitive
  
  #A <- fastICA2(Y1, numIC, "parallel", "kurtosis", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  #A <- fastICA2(Y1, numIC, "parallel", "skewness", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  #A <- fastICA2(Y1, numIC, "parallel", "logcosh", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  A <- fastICA2(Y1, numIC, "parallel", "exp", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  
  #A <- fastICA2(Y1, numIC, "deflation", "kurtosis", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  #A <- fastICA2(Y1, numIC, "deflation", "skewness", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  #A <- fastICA2(Y1, numIC, "deflation", "logcosh", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
  #A <- fastICA2(Y1, numIC, "deflation", "exp", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
 
  center <- apply(Y, 2, mean)
  

  
  # Projection of ICs to compositional space as vectors
  
  sp <- diag(100, numIC, numIC)
  sn <- diag(-100, numIC, numIC)
  
  Sp_temp <- sp%*%A$A
  Sn_temp <- sn%*%A$A
  
  
  mean0 <- t(c(apply(Y, 2, mean)))
  mean1 <- mean0
  
  for(j in 1:(numIC-1)){ 
    mean1 <- rbind(mean1, mean0) 
    }
  
  Sp_temp <- Sp_temp + mean1
  Sn_temp <- Sn_temp + mean1
  

  
  stdICvecR_dmy <- rbind(Sp_temp, center, Sn_temp)
  
  
  ICinRealSpace <-  stdICvecR_dmy

  
  for (i in 1:(2*numIC +1)){
    ICinRealSpace[i, ] <- (stdICvecR_dmy[i, ]*stdev+mean)
  }

  

  ICinRealSpaceall <- rbind(ICinRealSpaceall, ICinRealSpace)
  
  write.csv(
    ICinRealSpaceall,
    paste(path,"/Result_762C/ICnum=",numIC,"/Result_Robustness/ICvecReal_Repetitive.csv", sep = ""),
    row.names = F,
    quote = F
  ) 
 

    
  Sp_dmy <- Sp_temp%o%rep(1,2)
  Sn_dmy <- Sn_temp%o%rep(1,2)
  
  Sp[,,i] <- Sp_dmy[,,1]
  Sn[,,i] <- Sn_dmy[,,1]
  

  
  
  i <- i + 1
  
}


# ----- End of the script -----