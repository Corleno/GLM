rm(list = ls())
time = seq(1,20)
x = log(time)
y = c(1,6,16,23,27,39,31,30,43,51,63,70,88,97,91,104,110,113,149,159)
X = rbind(rep(1,20),x)

Lik_func <- function(beta, y, x, X){
  val = sum(y*log(t(beta)%*%X)-t(beta)%*%X-log(factorial(y)))
  return(val)
}

Sco_func <- function(beta, y, x, X){
  val = c(sum(y-exp(t(beta)%*%X)),sum(y*x-exp(t(beta)%*%X)*x))
  return(val)
}

f1 <- function(s, x){
  val = matrix(0, nrow = 2, ncol = 2)
  for (i in 1:length(s)){
    val = val + s[i] * matrix(c(1,x[i],x[i],x[i]^2), nrow = 2, ncol = 2)
  }
  return(val)
}

Fish_func <- function(beta, x, X){
  val = f1(exp(t(beta)%*%X),x)
  return(val)
}

###Initialization
beta_0 = c(0,0)
beta_1 = c(1,-1)
k = 1
###Iteration
while (norm(as.matrix(beta_1-beta_0),"f") > 0.01){
  k = k+1
  beta_0 = beta_1
  beta_1 = beta_0 + solve(Fish_func(beta_0, x, X))%*%Sco_func(beta_0, y, x, X)
}


###Beta VS Likelihood
beta = c(2,2)
beta = c(0.5,0.5)
Lik_func(beta, y, x, X)


  