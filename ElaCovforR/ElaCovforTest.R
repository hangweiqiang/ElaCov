ElaCovforTest <- function(X, Cr, Cv){
  # Compute the estimated correlation coefficient matrix 
  # and covariance matrix according to observation X
  # 
  # Args:
  #   X : nxp observation matrix
  #   Cr: penalty parameter for the estimation of correlation coefficient
  #   Cv: penalty parameter for the estimation of variances
  # 
  # Returns:
  #   List$R: The estimated correlation coefficient matrix of pxp
  #   List$S: The estimated covariance matrix of pxp
  # 
  # Error handling
  if (TRUE %in% is.na(X)){
  	stop("Arguments X must not have missing values.")
  }
  if (length(Cv)>1){
  	stop("Arguments Cv must not contain more than one constant")
  }

  size <- dim(X);
  n <- size[1];
  p <- size[2];
  x <- X-matrix(colMeans(X), n, p, byrow=TRUE);

  nCr <- length(Cr);
  SS0 <- myMatrix2tmpR(x, Cr + 1.0e-10)  
  R   <- SS0$SS;

  V <- matrix(diag(var(X)), 1, p);
  Vk <- matrix(0, p, 1);
  W <- sqrt(diag(var(X^2))) / sqrt(n);
  for (i in 1:p){
  	Ii <- (V <= V[i] + Cv * W[i]) & (V >= V[i] - Cv * W[i]);
  	Vk[i] <- mean(V[Ii])
  }
  V <- sqrt(Vk);
  S <- R;
  for (i in 1:nCr){
    S[ , , i] <- V %*% t(V) * R[ , , i];
  }
  ans <- list(R = R, S = S);
  return(ans);
}

myMatrix2tmpR=function(X,C){
	dimX <- dim(X);
	n <- dimX[1];p=dimX[2];

	Xm <- matrix(rep(colMeans(X), n), n, p, byrow=TRUE);
	X <- X-Xm;
	Xm <- matrix(rep(apply(X,2,sd), n), n, p, byrow=TRUE);
	X <- X/Xm;
  X <- X*sqrt(n/(n-1));

	sd <- diag(p);
  S0 <- diag(p)-100;
	for(i in 1:(p-1)){
	  for(j in (i+1):p){
	    e <- X[ ,i]*X[ ,j];
	    mij <- mean(e);
	    S0[i,j] <- mij;
	    sd[i,j] <- sqrt(sum((e-mij)^2)/(n-1));
	  }
	}
	sd <- sd/sqrt(n);

	stmp <- S0[S0>-10];
	SS <- rep(1, p*p*length(C));
	dim(SS) <- c(p, p, length(C));

	for(ic in 1:length(C)){
	  Csd <- C[ic]*sd;
	  for(i in 1:(p-1)){
		for(j in (i+1):p){
		  Iij <- (abs(stmp-S0[i,j])<=Csd[i,j]);
		  SS[i,j,ic] <- sum(stmp*Iij)/sum(Iij);
		  SS[j,i,ic] <- SS[i,j,ic];
		}
	  }
	}
	ans <- list(S0=S0,SS=SS);
	return(ans)
}