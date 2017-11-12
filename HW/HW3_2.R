########a
rm(list = ls())
x = c(1.6907, 1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
m = c(59,60,62,56,63,59,62,60)
y = c(6,13,18,28,52,53,61,60)

n = length(x)
beta = c(0,0)
u = rep(0,n)

times = 0
tot = 20
beta1 = rep(0,tot)
beta2 = rep(0,tot)


# ###Method 1 (Slice gibb sampling)
# ###Normal prior with N(0,I)
# while (times < tot){
#   for (i in 1:n){
#     u[i] = runif(1,0,choose(m[i],y[i])*(1-exp(-exp(beta[1]+beta[2]*x[i])))^y[i]*(exp(-exp(beta[1]+beta[2]*x[i])))^(m[i]-y[i]))
#   }
#   beta = rnorm(2)
#   flag = TRUE
#   beta_list = list()
#   while (flag) {
#     flag = FALSE
#     for (i in 1:n){
#       if (u[i] >= choose(m[i],y[i])*(1-exp(-exp(beta[1]+beta[2]*x[i])))^y[i]*(exp(-exp(beta[1]+beta[2]*x[i])))^(m[i]-y[i])){
#         flag = TRUE
#         beta = rnorm(2)
#         break
#       }
#     }
#   }
#   times = times + 1
#   beta1[times] = beta[1]
#   beta2[times] = beta[2]
# }
# 
# ###flat prior
# while (times < tot){
#   for (i in 1:n){
#     u[i] = runif(1,0,choose(m[i],y[i])*(1-exp(-exp(beta[1]+beta[2]*x[i])))^y[i]*(exp(-exp(beta[1]+beta[2]*x[i])))^(m[i]-y[i]))
#   }
#   beta = runif(2,-1000,1000)
#   flag = TRUE
#   beta_list = list()
#   while (flag) {
#     flag = FALSE
#     for (i in 1:n){
#       if (u[i] >= choose(m[i],y[i])*(1-exp(-exp(beta[1]+beta[2]*x[i])))^y[i]*(exp(-exp(beta[1]+beta[2]*x[i])))^(m[i]-y[i])){
#         flag = TRUE
#         beta = runif(2,-1000,1000)
#         break
#       }
#     }
#   }
#   times = times + 1
#   beta1[times] = beta[1]
#   beta2[times] = beta[2]
# }

###Method 2 (Metroplolis-Hasting algorithm)
rm(list = ls())
x = c(1.6907, 1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
m = c(59,60,62,56,63,59,62,60)
y = c(6,13,18,28,52,53,61,60)

n = length(x)
beta = c(0,0)
u = rep(0,n)

times = 0
tot = 1000
beta1 = rep(0,tot)
beta2 = rep(0,tot)

###flat prior
alpha_func <- function(x,m,y,beta,new_beta){
  n = length(y)
  X = rbind(rep(1,n),x)
  Eta_0 = t(beta) %*% X
  Eta_1 = t(new_beta) %*% X
  out = exp(sum(y*log(1-exp(-exp(Eta_1)))+(m-y)*log(exp(-exp(Eta_1)))) - sum(y*log(1-exp(-exp(Eta_0)))+(m-y)*log(exp(-exp(Eta_0)))))
  if (is.nan(out)){
    out = 0
  }
  return(out)
}

###normal prior
alpha_func_n <- function(x,m,y,beta,new_beta){
  n = length(y)
  X = rbind(rep(1,n),x)
  Eta_0 = t(beta) %*% X
  Eta_1 = t(new_beta) %*% X
  out = pnorm(new_beta[1])*pnorm(new_beta[2])/(pnorm(beta[1])*pnorm(beta[2]))*exp(sum(y*log(1-exp(-exp(Eta_1)))+(m-y)*log(exp(-exp(Eta_1)))) - sum(y*log(1-exp(-exp(Eta_0)))+(m-y)*log(exp(-exp(Eta_0)))))
  if (is.nan(out)){
    out = 0
  }
  return(out)
}

###posterior dist for flat prior
while (times<tot){
  flag = FALSE
  while (!flag){
    new_beta = beta + rnorm(2)
    u = runif(1)
    p = min(1,alpha_func(x,m,y,beta,new_beta))
    if (u< p){
      flag = TRUE
      beta = new_beta
      times = times + 1
      beta1[times] = beta[1]
      beta2[times] = beta[2]
    }
  }
}
###posterior dist for the median
med = (log(log(2))-beta1)/beta2

pdf(file = "pic/HW3_2/a1.pdf")
par(mfrow=c(1,3))
hist(beta1[100:tot], main = "Histogram of beta_1", xlab = "Beta_1")
hist(beta2[100:tot], main = "Histogram of beta_2", xlab = "Beta_2")
hist(med[100:tot], main = "Histogram of median lethal dose", xlab = "x")
dev.off()

###point and interval estimates for the dose-reponse curver
grid = seq(1.65,1.9,length.out = 50)
n = length(grid)
G = rbind(rep(1,n),grid)
Beta = cbind(beta1,beta2)
Eta = Beta %*% G
Pi = 1-exp(-exp(Eta))
mean_vec = colMeans(Pi)
Q_0 = apply(Pi, 2, function(x) quantile(x,0.025))
Q_1 = apply(Pi, 2, function(x) quantile(x,0.975))
pdf(file = "pic/HW3_2/a3.pdf")
plot(grid, mean_vec, type = 'l', xlab = "x", ylab = expression(pi), main = "Estimated dose-response curve")
lines(grid, Q_0, type = 'l', lty = 2, col = 2)
lines(grid, Q_1, type = 'l', lty = 2, col = 2)
dev.off()

RG = rbind(rep(1,length(x)), x)
REta = Beta %*% RG
RPi = colMeans(1-exp(-exp(REta)))
res1 = y/m - RPi
round(res1, 2)
var_f <- function(X){
  v = rep(0,dim(X)[2])
  for (i in 1:dim(X)[2]){
    v[i] = var(X[,i])
  }
  return(v)
}
var1 = sum(var_f(1-exp(-exp(REta))))


###posterior dist for normal prior
n = length(x)
beta = c(0,0)
u = rep(0,n)

times = 0
tot = 1000
beta1 = rep(0,tot)
beta2 = rep(0,tot)

while (times<tot){
  flag = FALSE
  while (!flag){
    new_beta = beta + rnorm(2)
    u = runif(1)
    p = min(1,alpha_func_n(x,m,y,beta,new_beta))
    if (u< p){
      flag = TRUE
      beta = new_beta
      times = times + 1
      beta1[times] = beta[1]
      beta2[times] = beta[2]
    }
  }
}

pdf(file = "pic/HW3_2/a2.pdf")
par(mfrow=c(1,2))
hist(beta1[100:tot], main = "Histogram of beta_1", xlab = "Beta_1")
hist(beta2[100:tot], main = "Histogram of beta_2", xlab = "Beta_2")
dev.off()


########b
rm(list = ls())
x = c(1.6907, 1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
m = c(59,60,62,56,63,59,62,60)
y = c(6,13,18,28,52,53,61,60)

n = length(x)
beta = c(0,0)
u = rep(0,n)

times = 0
tot = 1000
beta1 = rep(0,tot)
beta2 = rep(0,tot)

###flat prior
alpha_func <- function(x,m,y,beta,new_beta){
  n = length(y)
  X = rbind(rep(1,n),x)
  Eta_0 = t(beta) %*% X
  Eta_1 = t(new_beta) %*% X
  out = exp(sum(y*log(exp(Eta_1)/(1+exp(Eta_1)))+(m-y)*log(1/(1+exp(Eta_1)))) - sum(y*log(exp(Eta_0)/(1+exp(Eta_0)))+(m-y)*log(1/(1+exp(Eta_0)))))
  if (is.nan(out)){
    out = 0
  }
  return(out)
}

###posterior dist for flat prior
while (times<tot){
  flag = FALSE
  while (!flag){
    new_beta = beta + rnorm(2)
    u = runif(1)
    p = min(1,alpha_func(x,m,y,beta,new_beta))
    if (u< p){
      flag = TRUE
      beta = new_beta
      times = times + 1
      beta1[times] = beta[1]
      beta2[times] = beta[2]
    }
  }
}

###posterior dist for the median
med = (-beta1)/beta2

pdf(file = "pic/HW3_2/b1.pdf")
par(mfrow=c(1,3))
hist(beta1[100:tot], main = "Histogram of beta_1", xlab = "Beta_1")
hist(beta2[100:tot], main = "Histogram of beta_2", xlab = "Beta_2")
hist(med[100:tot], main = "Histogram of median lethal dose", xlab = "x")
dev.off()

###point and interval estimates for the dose-reponse curver
grid = seq(1.65,1.9,length.out = 50)
n = length(grid)
G = rbind(rep(1,n),grid)
Beta = cbind(beta1,beta2)
Eta = Beta %*% G
Pi = exp(Eta)/(1+exp(Eta))
mean_vec = colMeans(Pi)
Q_0 = apply(Pi, 2, function(x) quantile(x,0.025))
Q_1 = apply(Pi, 2, function(x) quantile(x,0.975))
pdf(file = "pic/HW3_2/b2.pdf")
plot(grid, mean_vec, type = 'l', xlab = "x", ylab = expression(pi), main = "Estimated dose-response curve")
lines(grid, Q_0, type = 'l', lty = 2, col = 2)
lines(grid, Q_1, type = 'l', lty = 2, col = 2)
dev.off()

RG = rbind(rep(1,length(x)), x)
REta = Beta %*% RG
RPi = colMeans(exp(REta)/(1+exp(REta)))
res2 = y/m - RPi
round(res2,2)
var_f <- function(X){
  v = rep(0,dim(X)[2])
  for (i in 1:dim(X)[2]){
    v[i] = var(X[,i])
  }
  return(v)
}
var2 = sum(var_f(exp(REta)/(1+exp(REta))))

########c
rm(list = ls())
x = c(1.6907, 1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
m = c(59,60,62,56,63,59,62,60)
y = c(6,13,18,28,52,53,61,60)

n = length(x)
alpha = 1
beta = c(0,0)
u = rep(0,n)

times = 0
tot = 1000
alpha1 = rep(0,tot)
beta1 = rep(0,tot)
beta2 = rep(0,tot)

###flat prior
alpha_func <- function(x,m,y,alpha, new_alpha,beta,new_beta){
  n = length(y)
  X = rbind(rep(1,n),x)
  Eta_0 = t(beta) %*% X
  Eta_1 = t(new_beta) %*% X
  out = exp(sum(y*log(exp(new_alpha*Eta_1)/((1+exp(Eta_1))^new_alpha))+(m-y)*log(1-exp(new_alpha*Eta_1)/(1+exp(Eta_1))^new_alpha)) - sum(y*log(exp(alpha*Eta_0)/((1+exp(Eta_0))^alpha))+(m-y)*log(1-exp(alpha*Eta_0)/(1+exp(Eta_0))^alpha)))
  if (is.nan(out)){
    out = 0
  }
  return(out)
}

###posterior dist for flat prior
while (times<tot){
  flag = FALSE
  while (!flag){
    new_alpha = alpha + rnorm(1)
    new_beta = beta + rnorm(2)
    u = runif(1)
    p = min(1,alpha_func(x,m,y,alpha,new_alpha,beta,new_beta))
    if (u< p){
      flag = TRUE
      beta = new_beta
      alpha = new_alpha
      times = times + 1
      beta1[times] = beta[1]
      beta2[times] = beta[2]
      alpha1[times] = alpha
    }
  }
}

###posterior dist for the median
med = (log(1/(2^(1/alpha1)-1)) - beta1)/beta2

pdf(file = "pic/HW3_2/c1.pdf")
par(mfrow=c(1,4))
hist(alpha1[100:tot], main = "Histogram of alpha", xlab = "Alpha")
hist(beta1[100:tot], main = "Histogram of beta_1", xlab = "Beta_1")
hist(beta2[100:tot], main = "Histogram of beta_2", xlab = "Beta_2")
hist(med[100:tot], main = "Histogram of median lethal dose", xlab = "x")
dev.off()

###point and interval estimates for the dose-reponse curver
grid = seq(1.65,1.9,length.out = 50)
n = length(grid)
G = rbind(rep(1,n),grid)
Beta = cbind(beta1,beta2)
Eta = Beta %*% G
Pi = exp(alpha1*Eta)/(1+exp(Eta))^alpha1
mean_vec = colMeans(Pi)
Q_0 = apply(Pi, 2, function(x) quantile(x,0.025))
Q_1 = apply(Pi, 2, function(x) quantile(x,0.975))
pdf(file = "pic/HW3_2/c2.pdf")
plot(grid, mean_vec, type = 'l', xlab = "x", ylab = expression(pi), main = "Estimated dose-response curve")
lines(grid, Q_0, type = 'l', lty = 2, col = 2)
lines(grid, Q_1, type = 'l', lty = 2, col = 2)
dev.off()

RG = rbind(rep(1,length(x)), x)
REta = Beta %*% RG
RPi = colMeans(exp(alpha1*REta)/(1+exp(REta))^alpha1)
res3 = y/m - RPi
var_f <- function(X){
  v = rep(0,dim(X)[2])
  for (i in 1:dim(X)[2]){
    v[i] = var(X[,i])
  }
  return(v)
}
var3 = sum(var_f(exp(alpha1*REta)/(1+exp(REta))^alpha1))