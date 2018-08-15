# 
# The working directory should be set to "~\\ElaCovforR"
#
source("ElaCOV.R")
source("ElaCOVforTest.R")

n <- 100;
p <- 100;
X <- matrix(rnorm(n*p, 0, 1), n, p)
Cv <- .3
Cr <- c(.3,.5,1)

ans1 <- ElaCov(X,Cr,Cv)
ans2 <- ElaCovforTest(X,Cr,Cv)
max(abs(ans1$R - ans2$R))
max(abs(ans1$S - ans2$S))