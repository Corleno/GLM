#!/Library/Frameworks/R.framework/Versions/3.2/Resources/bin/Rscript
rm(list = ls())
data = read.csv("fabric_1.csv", header = FALSE, sep = "")
y = data[,2]
x = data[,1]
glm.poisson = glm(y~x, family = poisson)
summary(glm.poisson)
glm.quasipoisson = glm(y~x, family = quasipoisson)
summary(glm.quasipoisson)

###Compute the expected Fisher matrix
n = length(x)
X = cbind(rep(1,n),x)

beta = glm.poisson$coefficients
mu = exp(X %*% beta)

Jf <- function(X,mu,k,j,n){
  Jf = 0
  for (i in 1:n){
    Jf = Jf + X[i,k]*X[i,j]*mu[i]
  }
  return(Jf)
}

Jf_M <- function(X,mu,n){
  Jf_M = matrix(c(Jf(X,mu,1,1,n),Jf(X,mu,2,1,n),Jf(X,mu,1,2,n),Jf(X,mu,2,2,n)),nrow = 2, ncol = 2)
  return(Jf_M)
}

FM = Jf_M(X,mu,n)
FM_Q = FM/2.238357
