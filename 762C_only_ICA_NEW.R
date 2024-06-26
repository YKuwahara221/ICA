# R script for ICA of ODP Site 762C in the Indian ocean
# Originally written by K.Yasukawa on August 22, 2014. (Modified on August 12, 2016)
# Refactored by M.Hirako on September 21, 2020, and modified by Y. Kuwahara on Jan 24, 2022
# For Site 762C paper, modified by Y. Kuwahara on Jan 6, 2023

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

numIC <- 3   # Number of ICs

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

#numEM <-  unlist(length(em$Endmember))  # Number of data in Dataset of reference materials



# Independent Component Analysis with FastICA algorithm (Hyvarinen et al., 2001)

#

# === FastICA function modified from R Package "fastICA" ver 1.1-16 by Marchini & Heaton (2012) ===
fastICA2 <-
  function(X,
           n.comp,
           alg.typ = c("parallel", "deflation"),
           fun = c("kurtosis", "skewness", "logcosh", "exp"),
           a1 = 1,
           a2 = 1,
           maxit = 200,
           tol = 1e-04,
           verbose = TRUE,
           w.init = NULL) {
    X <- as.matrix(X)
    XT <- t(X)
    
    alg.typ <- match.arg(alg.typ)
    fun <- match.arg(fun)
    n <- nrow(X)
    p <- ncol(X)
    if (n.comp > min(n, p)) {
      message("'n.comp' is too large: reset to ", min(n, p))
      n.comp <- min(n, p)
    }
    
    if (is.null(w.init))
      w.init <- matrix(rnorm(n.comp ^ 2), n.comp, n.comp)
    else{
      if (!is.matrix(w.init) || length(w.init) != (n.comp ^ 2))
        stop("w.init is not a matrix or is the wrong size.")
    }
    
    if (verbose)
      message("Centering")
    mean <- apply(X, 2, mean)
    X0 <- t(X)
    Xc <- t(X0 - mean)
    XcT <-t(Xc)                       # Q-mode-like calculation: transpose of input data matrix X
    
    if (verbose)
      message("Whitening")
    V <- XcT %*% Xc / n
    s <- La.svd(V)
    D <- diag(c(1 / sqrt(s$d)))
    KT <- D %*% t(s$u)
    KT <- matrix(KT[1:n.comp, ], n.comp, p)
    ZT <- KT %*% XcT
    
    if (alg.typ == "deflation")
      WT <-
      as.matrix(
        ica2.def(
          ZT,
          n.comp,
          tol = tol,
          fun = fun,
          a1 = a1,
          a2 = a2,
          maxit = maxit,
          verbose = verbose,
          w.init = w.init
        )
      )
    
    else if (alg.typ == "parallel")
      WT <-
      as.matrix(
        ica2.par(
          ZT,
          n.comp,
          tol = tol,
          fun = fun,
          a1 = a1,
          a2 = a2,
          maxit = maxit,
          verbose = verbose,
          w.init = w.init
        )
      )
    
    w <- WT %*% KT
    ST <- w %*% XcT
    AT <- t(w) %*% solve(w %*% t(w))
    
    return(list(
      X = t(XT),
      K = t(KT),
      W = t(WT),
      A = t(AT),
      S = t(ST)
    ))
  }

#  fastICA algorithm using a deflation scheme for orthogonalization
ica2.def <-
  function(X,
           n.comp,
           tol,
           fun,
           a1,
           a2,
           maxit,
           verbose,
           w.init) {
    if (verbose && fun ==  "kurtosis")
      message("Deflation fastICA using kurtosis approximation to non-gaussianity")
    if (verbose && fun == "skewness")
      message("Deflation fastICA using skewness approximation to non-gaussianity")
    if (verbose && fun ==  "logcosh")
      message("Deflation fastICA using logcosh approximation to neg-entropy")
    if (verbose && fun == "exp")
      message("Deflation fastICA using exponential approximation to neg-entropy")
    
    n <- nrow(X)
    p <- ncol(X)
    W <- matrix(0, n.comp, n.comp)
    
    for (i in 1:n.comp) {
      if (verbose)
        message("Component ", i)
      w <-
        matrix(w.init[i, ], n.comp, 1)    # w is i-th row vector of matrix W
      if (i > 1) {
        tmp <- w
        tmp[1:length(tmp)] <- 0
        for (u in 1:(i - 1)) {
          k <- sum(w * W[u, ])
          tmp <- tmp + k * W[u, ]
        }
        w <- w - tmp
      }
      w <- w / sqrt(sum(w ^ 2))
      
      lim <- rep(1000, maxit)
      it <- 1
      
      if (fun == "kurtosis") {
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- wx ^ 3
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          V1 <- apply(xgwx, 1, FUN = mean)
          
          gwx2 <- 3
          gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
          V2 <- mean(gwx2) * w
          
          wp <- V1 - V2
          wp <- matrix(wp, n.comp, 1)
          it <- it + 1
          
          if (i > 1) {
            tmp <- wp
            tmp[1:length(tmp)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(wp * W[u, ])
              tmp <- tmp + k * W[u, ]
            }
            wp <- wp - tmp
          }
          wp <- wp / sqrt(sum(wp ^ 2))
          
          lim[it] <- Mod(Mod(sum((wp * w))) - 1)
          if (verbose)
            message("Iteration ", it - 1, "  tol = ", format(lim[it]))
          
          w <- matrix(wp, n.comp, 1)
        }
      }
      
      if (fun == "skewness") {
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- wx ^ 2
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          V1 <- apply(xgwx, 1, FUN = mean)
          
          gwx2 <- 0
          gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
          V2 <- mean(gwx2) * w
          
          wp <- V1 - V2
          wp <- matrix(wp, n.comp, 1)
          it <- it + 1
          
          if (i > 1) {
            tmp <- wp
            tmp[1:length(tmp)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(wp * W[u, ])
              tmp <- tmp + k * W[u, ]
            }
            wp <- wp - tmp
          }
          wp <- wp / sqrt(sum(wp ^ 2))
          
          lim[it] <- Mod(Mod(sum((wp * w))) - 1)
          if (verbose)
            message("Iteration ", it - 1, "  tol = ", format(lim[it]))
          
          w <- matrix(wp, n.comp, 1)
        }
      }
      
      if (fun == "logcosh") {
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- tanh(a1 * wx)
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          V1 <- apply(xgwx, 1, FUN = mean)
          
          gwx2 <- a1 * (1 - (tanh(a1 * wx)) ^ 2)
          gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
          V2 <- mean(gwx2) * w
          
          wp <- V1 - V2
          wp <- matrix(wp, n.comp, 1)
          it <- it + 1
          
          if (i > 1) {
            tmp <- wp
            tmp[1:length(tmp)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(wp * W[u, ])
              tmp <- tmp + k * W[u, ]
            }
            wp <- wp - tmp
          }
          wp <- wp / sqrt(sum(wp ^ 2))
          
          lim[it] <- Mod(Mod(sum((wp * w))) - 1)
          if (verbose)
            message("Iteration ", it - 1, "  tol = ", format(lim[it]))
          
          w <- matrix(wp, n.comp, 1)
        }
      }
      
      if (fun == "exp") {
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- wx * exp(-a2 * (wx ^ 2) / 2)
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          V1 <- apply(xgwx, 1, FUN = mean)
          
          gwx2 <- (1 - a2 * wx ^ 2) * exp(-a2 * (wx ^ 2) / 2)
          gwx2 <- matrix(gwx2, n.comp, p, byrow = TRUE)
          V2 <- mean(gwx2) * w
          
          wp <- V1 - V2
          wp <- matrix(wp, n.comp, 1)
          it <- it + 1
          
          if (i > 1) {
            tmp <- wp
            tmp[1:length(tmp)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(wp * W[u, ])
              tmp <- tmp + k * W[u, ]
            }
            wp <- wp - tmp
          }
          wp <- wp / sqrt(sum(wp ^ 2))
          
          lim[it] <- Mod(Mod(sum((wp * w))) - 1)
          if (verbose)
            message("Iteration ", it - 1, "  tol = ", format(lim[it]))
          
          w <- matrix(wp, n.comp, 1)
        }
      }
      W[i, ] <- w
    }
    W
  }

#  fastICA algorithm using a symmetric scheme for orthogonalization
ica2.par <-
  function(X,
           n.comp,
           tol,
           fun,
           a1,
           a2,
           maxit,
           verbose,
           w.init) {
    Diag <- function(d)
      if (length(d) > 1L)
        diag(d)
    else
      as.matrix(d)
    
    n <- nrow(X)
    p <- ncol(X)
    
    W <- w.init
    sW <- La.svd(W)
    
    W <- sW$u %*% Diag(1 / sW$d) %*% W
    Wp <- W
    
    lim <- rep(1000, maxit)
    it <- 1
    
    if (fun == "kurtosis") {
      message("Symmetric fastICA using kurtosis approximation to non-gaussianity")
      while (lim[it] > tol && it < maxit) {
        wx <- W %*% X
        gwx <- wx ^ 3
        V1 <- gwx %*% t(X) / p
        gwx2 <- matrix(3, n.comp, p)
        V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
        Wp <- V1 - V2
        sWp <- La.svd(Wp)
        Wp <- sWp$u %*% Diag(1 / sWp$d) %*% t(sWp$u) %*% Wp
        lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(
          W
        ))) - 1))
        
        W <- Wp
        
        if (verbose)
          message("Iteration ", it, "  tol = ", format(lim[it + 1]))
        
        it <- it + 1
      }
    }
    
    if (fun == "skewness") {
      message("Symmetric fastICA using skewness approximation to non-gaussianity")
      while (lim[it] > tol && it < maxit) {
        wx <- W %*% X
        gwx <- wx ^ 2
        V1 <- gwx %*% t(X) / p
        gwx2 <- matrix(0, n.comp, p)
        V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
        Wp <- V1 - V2
        sWp <- La.svd(Wp)
        Wp <- sWp$u %*% Diag(1 / sWp$d) %*% t(sWp$u) %*% Wp
        lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(
          W
        ))) - 1))
        
        W <- Wp
        
        if (verbose)
          message("Iteration ", it, "  tol = ", format(lim[it + 1]))
        
        it <- it + 1
      }
    }
    
    if (fun == "logcosh") {
      message("Symmetric fastICA using logcosh approximation to neg-entropy")
      V1 <- matrix(0,
                   nrow = n.comp,
                   ncol = n.comp,
                   byrow = TRUE)
      while (lim[it] > tol && it < maxit) {
        wx <- W %*% X
        gwx <- tanh(a1 * wx)
        V1 <- gwx %*% t(X) / p
        gwx2 <- a1 * (1 - (gwx) ^ 2)
        V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
        Wp <- V1 - V2
        sWp <- La.svd(Wp)
        Wp <- sWp$u %*% Diag(1 / sWp$d) %*% t(sWp$u) %*% Wp
        lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(
          W
        ))) - 1))
        
        W <- Wp
        
        if (verbose)
          message("Iteration ", it, "  tol = ", format(lim[it + 1]))
        
        it <- it + 1
      }
    }
    
    if (fun == "exp") {
      message("Symmetric fastICA using exponential approximation to neg-entropy")
      while (lim[it] > tol && it < maxit) {
        wx <- W %*% X
        gwx <- wx * exp(-a2 * (wx ^ 2) / 2)
        V1 <- gwx %*% t(X) / p
        gwx2 <- (1 - a2 * wx ^ 2) * exp(-a2 * (wx ^ 2) / 2)
        V2 <- Diag(apply(gwx2, 1, FUN = mean)) %*% W
        Wp <- V1 - V2
        sWp <- La.svd(Wp)
        Wp <- sWp$u %*% Diag(1 / sWp$d) %*% t(sWp$u) %*% Wp
        lim[it + 1] <- max(Mod(Mod(diag(Wp %*% t(
          W
        ))) - 1))
        
        W <- Wp
        
        if (verbose)
          message("Iteration ", it, "  tol = ", format(lim[it + 1]))
        
        it <- it + 1
      }
    }
    W
  }

# Performing fastICA2 function  (Remove "#" to select an arbitral combination of orthogonalization and G-function)

#A <- fastICA2(Y, numIC, "parallel", "kurtosis", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
#A <- fastICA2(Y, numIC, "parallel", "skewness", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
#A <- fastICA2(Y, numIC, "parallel", "logcosh", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
A <- fastICA2(Y, numIC, "parallel", "exp", 1.0, 1.0, 200, 1e-04, TRUE, NULL)

#A <- fastICA2(Y, numIC, "deflation", "kurtosis", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
#A <- fastICA2(Y, numIC, "deflation", "skewness", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
#A <- fastICA2(Y, numIC, "deflation", "logcosh", 1.0, 1.0, 200, 1e-04, TRUE, NULL)
#A <- fastICA2(Y, numIC, "deflation", "exp", 1.0, 1.0, 200, 1e-04, TRUE, NULL)

center <- apply(Y, 2, mean)

# Projection of ICs to Real Space

sp <-
  diag(100, numIC, numIC)       # diagonal matrix of numIC x numIC (base vectors in the positive side of IC space)
sn <-
  diag(-100, numIC, numIC)      # diagonal matrix of numIC x numIC (base vectors in the negative side of IC space)

Sp <-
  sp %*% A$A      # Sp is a projected composition corresponding to IC score of 100 (positive side)
Sn <-
  sn %*% A$A      # Sn is a projected composition corresponding to IC score of -100 (negative side)

mean0 <- t(c(apply(Y, 2, mean)))
mean1 <- mean0
for (i in 1:(numIC - 1)) {
  mean1 <- rbind(mean1, mean0)
}

Sp <-
  Sp + mean1        # Parallel translation of Sp to center in the compositional data
Sn <-
  Sn + mean1        # Parallel translation of Sn to center in the compositional data

stdICvec <- rbind(Sp, center, Sn)

# Exporting results as text files
write.csv(
  stdICvec,
  paste(path, "/Result_762C/ICnum=",numIC,"/Result/StdICvec.csv", sep = ""),
  row.names = F,
  quote = F
)

ICinRealSpace <- stdICvec
for (i in 1:(2 * numIC + 1)) {
  ICinRealSpace[i, ] <- (stdICvec[i, ] * stdev + mean)
}

write.csv(
  ICinRealSpace,
  paste(path, "/Result_762C/ICnum=",numIC,"/Result/RealICvec.csv", sep = ""),
  row.names = F,
  quote = F
)


# Data in IC space
ICscore <- rbind(A$S)
sampleID_all <- X$Sample.ID
ICscore <- cbind(sampleID_all, ICscore)
nameICscore <- rep("numIC")
for (i in 1:numIC) {
  nameICscore[i] <-  paste("IC", i, sep = "")
}
colnames(ICscore) <- nameICscore
write.csv(
  ICscore,
  paste(path, "/Result_762C/ICnum=",numIC,"/Result/ICscores.csv", sep = ""),
  row.names = F,
  quote = F
)


# Data in PC space
Xpc <- as.matrix(X[, index_elements_start:index_elements_end])
XpcT <- t(Xpc)
Xpc_centered <- t(XpcT - center)
Xpc_whitened <- Xpc_centered %*% A$K
PCscore <- cbind(sampleID_all, Xpc_whitened)

# Projection of endmember data to IC space

centerT <- as.vector(center)
center <- t(centerT)


# Projection of PCs to IC space
zp <- diag(10, numIC, numIC)
zn <- diag(-10, numIC, numIC)

Zp <- zp %*% A$W
Zn <- zn %*% A$W

center0 <- diag(0, 1, numIC)


# Relative loadings

Loadings <- A$A
Xsd <- apply(Y, 2, sd)
Xsd <- as.vector(Xsd)
colors <- diag(0, numIC, numDim)

tmp <- array(0, dim = numIC)

for (i in 1:numIC) {
  for (j in 1:numDim) {
    Loadings[i, j] <- Loadings[i, j] / Xsd[j]
    tmp[i] <- tmp[i] + abs(Loadings[i, j])
    colors[i, j] <- ifelse(Loadings[i, j] >= 0, "blue", "red")
  }
}


ICnames <- rep("", numIC)
for (i in 1:numIC) {
  ICnames[i] <- paste("IC", i, sep = "")
}
elements <- names(as.data.frame(Y))
par(mfrow = c(3, 1))
for (i in 1:numIC) {
  barplot(Loadings[i, ],
          col = colors[i, ],
          names = elements,
          main = paste("IC", i))
}

Loadings <- as.data.frame(Loadings)
Loadings <- rbind(elements, Loadings)
Loadings <- t(as.matrix(Loadings))
Loadings <- as.data.frame(Loadings)
names(Loadings) <- append(ICnames, "Elements", after = 0)

write.csv(
  Loadings,
  file = paste(path, "/Result_762C/ICnum=",numIC,"/Result/Loadings.csv", sep = ""),
  row.names = F,
  quote = F
)


#histgram of IC scores

ICscores <- as.data.frame(A$S)

names(ICscores) <- ICnames
ResultAll <- cbind(X,REY, ICscores)
write.csv(
  ResultAll,
  file = paste(path, "/Result_762C/ICnum=",numIC,"/Result/ResultAll.csv", sep = ""),
  row.names = F,
  quote = F
)
#StdResult <- cbind(X[, 1:43], Y,REY, ICscores) #modified by Y. Kuwahara
StdResult <- cbind(Y,REY, ICscores)
write.csv(
  StdResult,
  paste(path, "/Result_762C/ICnum=",numIC,"/Result/StdResult.csv", sep = ""),
  row.names = F,
  quote = F
)


par(mfrow = c(3, 3))
for (i in 1:numIC) {
  name <- paste("IC", i, sep = "")
  hist(ICscores[, i], 100, main = name, xlab = "")
}


# ----- End of the script -----
