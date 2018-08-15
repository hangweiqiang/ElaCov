ElaCov <- function(X, Cr, Cv){
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
  # Error handling end
  size <- dim(X);
  n <- size[1];
  p <- size[2];
  x <- X-matrix(colMeans(X), n, p, byrow=TRUE);

  nCr <- length(Cr);
  SS0 <- myMatrix2tmpRC(X, Cr)
  R   <- SS0$SS;

  V <- matrix(diag(var(X)), 1, p);
  Vs <- sort(V);
  csum <- cumsum(Vs);                  # cumulative sum of Vs
  rsum <- cumsum(Vs[length(Vs):1]);   
  rsum <- rsum[length(rsum):1];        # inverse cummulative sum of Vs
  
  Vk <- matrix(0, p, 1);
  W <- sqrt(diag(var(X^2))) / sqrt(n);
  for (i in 1:p){
    # 
    # The code in the for-loop is a faster way, using binary search method, 
    # to realize the codes below 
    # Ii <- (V <= V[i] + Cv * W[i]) & (V >= V[i] - Cv * W[i]);
    # Vk[i] <- mean(V[Ii])
    # 
    First <- V[i]-Cv*W[i];
    Second<- V[i]+Cv*W[i];
    pl <- 0;
    pr <- length(Vs);
    while (pl<pr){
      mid <- floor((pl+pr+1)/2);
      if (Vs[mid]<=Second){
        pl <- mid;
      }else{
        pr <- mid-1;
      }
    }
    SecondPos <- pl;
    pl <- 1;
    pr <- length(Vs)+1;
    while (pl<pr){
      mid <- floor((pl+pr)/2);
      if (Vs[mid]>=First){
        pr <- mid;
      }else{
        pl <- mid+1;
      }
    }
    pl <- 1;
    pr <- length(Vs)+1;
    while (pl<pr){
      mid <- floor((pl+pr)/2);
      if (Vs[mid]>First){
        pr <- mid;
      }else{
        pl <- mid+1;
      }
    }
    FirstPos <- pl;
    if (FirstPos != length(Vs) + 1 && SecondPos != 0){
      SumIi <- SecondPos - FirstPos + 1;
      SumVIi <- rsum[FirstPos] + csum[SecondPos] - rsum[1];
      Vk[i] <- SumVIi / SumIi;
    }else{
      Vk[i] <- 0;
    }
  }
  V <- sqrt(Vk);
  S <- R;
  for (i in 1:nCr){
    S[ , , i] <- V %*% t(V) * R[, , i];
  }
  ans <- list(R = R, S = S);
  return(ans);  
}

myMatrix2tmpRC <- function(X,Cr)
{
  # Compute the estimated correlation coefficient matrix of pxp according to X and 
  # it is a faster way to realize the function myMatrix2tmpR in ElaCovforTest.R by 
  # using binary search method and the core part is written by C++ language.
  # 
  # The working directory should be set to the directory of myMatrix2tmpCppCoreForR.dll
  #
  # Args:
  #   X : nxp observation matrix
  #   Cr: penalty parameter for the estimation of correlation coefficient
  # 
  # Returns:
  #   List$SS: The estimated correlation coefficient matrix of pxp
  #   List$S0: 
  # 
	dyn.load('myMatrix2tmpCppCoreForR.dll');  # load myMatrix2tmpCppCoreForR.dll
	X <- data.matrix(X);
	Cr <- data.matrix(Cr);
	ans <- .Call("myMatrix2tmpC",X,Cr+1.0e-10);
	dimX <- dim(X);
	dimCr <- dim(Cr);
	n <- dimX[1];
	p <- dimX[2];
  c <- dimCr[1];
	dim(ans$S0) <- c(p,p);
	dim(ans$SS) <- c(p,p,c);
	dyn.unload("myMatrix2tmpCppCoreForR.dll")  # unload myMatrix2tmpCppCoreForR.dll
	return(ans);
}

