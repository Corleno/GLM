# #!/Library/Frameworks/R.framework/Versions/3.2/Resources/bin/Rscript
# rm(list = ls())
# data = read.csv("fabric_1.csv", header = FALSE, sep = "")
# y = data[,2]
# x = data[,1]
# y = y[order(x)]
# x = x[order(x)]
# 
# 
# ###MCMC
# tot = 5000
# beta1 = beta2 = rep(0,tot)
# times = 0
# beta1_c = beta2_c = 0
# beta1_n = beta2_n = 0
# 
# ratio_func <- function(beta1_n, beta2_n, beta1_c, beta2_c, x, y){
#   n = length(x)
#   ita_c = t(c(beta1_c, beta2_c)) %*% rbind(rep(1,n),x)
#   lambda_c = exp(ita_c)
#   cur = sum(y*ita_c - lambda_c)
#   ita_n = t(c(beta1_n, beta2_n)) %*% rbind(rep(1,n),x)
#   lambda_n = exp(ita_n)
#   new = sum(y*ita_n - lambda_n)
#   if (new == -Inf){
#     ratio = 0
#   }else{
#     r = exp(new - cur)
#     ratio = min(1,r)
#   }
# }
# 
# while (times < tot){
#   times = times + 1
#   flag = FALSE
#   nflag = 0
#   while (!flag){
#     nflag = nflag + 1
#     beta1_n = beta1_c + rnorm(1, sd = 1/nflag)
#     beta2_n = beta2_c + rnorm(1, sd = 1/nflag)
#     r = runif(1)
#     acc = ratio_func(beta1_n, beta2_n, beta1_c, beta2_c, x, y)
#     if (r <= acc){
#       flag=TRUE
#       beta1[times] = beta1_c = beta1_n
#       beta2[times] = beta2_c = beta2_n
#     }
#   }
# }
# 
# ###Plot the parameter curve
# #plot(beta1, type = "l")
# #plot(beta2, type = "l")
# pdf(file = "pic/HW5_1/pos.pdf")
# par(mfrow = c(1,2))
# post_beta1 = beta1[1001:tot]
# post_beta2 = beta2[1001:tot]
# h1 = hist(post_beta1, freq = FALSE, xlab = expression(beta[1]), ylab = "Density")
# lines(h1$mids, h1$density, type = "l", lty = 2, col="red")
# #abline(h = 1, col = "blue")
# h2 = hist(post_beta2, freq = FALSE, xlab = expression(beta[1]), ylab = "Density")
# lines(h2$mids, h2$density, type = "l", lty = 2, col="red")
# #abline(h = 1, col = "blue")
# dev.off()
# 
# ###Response mean
# inv = seq(min(x), max(x), length.out = 1000)
# p_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,1000),inv)
# p_mu = exp(p_ita)
# p_mu_inv = apply(p_mu, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))
# pdf(file = "pic/HW5_1/RMC.pdf")
# plot(inv, p_mu_inv[1,], col = 1, type = "l", xlab = "Length", ylab = "Number of fault")
# lines(inv, p_mu_inv[2,], col = 1, lty = 2)
# lines(inv, p_mu_inv[3,], col = 1, lty = 2)
# dev.off()
# 
# ###Posterior predictive residuals
# n = length(x)
# p0_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,n),x)
# p0_mu = exp(p0_ita)
# p0_y = apply(p0_mu, 1:2, function(x) rpois(1,x))
# p0_e = p0_y - matrix(rep(y,4000), ncol = 32, byrow = TRUE)
# quan = apply(p0_e, 2, function(x) c(quantile(x,c(0.2,0.4,0.6,0.8))))
# pdf(file = "pic/HW5_1/box.pdf")
# boxplot(p0_e)
# abline(h = 0, col="red")
# dev.off()






# ###Hierarchical Model
# rm(list = ls())
# data = read.csv("fabric_1.csv", header = FALSE, sep = "")
# y = data[,2]
# x = data[,1]
# y = y[order(x)]
# x = x[order(x)]
# 
# acc_func1 <- function(mu_c, gamma_c, log_lambda_c, log_lambda_n, y){
#   n = length(y)
#   cur = new = 0
#   lambda_c = exp(log_lambda_c)
#   lambda_n = exp(log_lambda_n)
#   for (i in 1:n){
#     cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_c[i]/lambda_c, log = TRUE)
#     new = new + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_n, scale = gamma_c[i]/lambda_n, log = TRUE)
#   }
#   if (new == -Inf){
#     acc = 0
#   }else{
#     acc = min(1, exp(new - cur))
#   }
#   return(acc)
# }
# 
# acc_func2 <- function(mu_c, beta1_c, beta2_c, beta1_n, lambda_c, beta2_n, x, y){
#   n = length(x)
#   gamma_c = exp(t(c(beta1_c, beta2_c)) %*% rbind(rep(1,n), x))
#   gamma_n = exp(t(c(beta1_n, beta2_n)) %*% rbind(rep(1,n), x))
#   cur = new = 0
#   for (i in 1:n){
#     cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_c[i]/lambda_c, log = TRUE)
#     new = new + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_n[i]/lambda_c, log = TRUE)
#   }
#   if (new == -Inf){
#     acc = 0
#   }else{
#     acc = min(1, exp(new - cur))
#   }
# }
# 
# tot = 10000
# ###Test: 
# ###tot = 10
# beta1 = beta2 = lambda = rep(0,tot)
# n = length(x)
# 
# times = 0
# beta1_c = beta2_c = lambda_c = 0.01
# beta1_n = beta2_n = lambda_n = 0.01
# mu_c = rep(0,n)
# 
# while (times < tot){
#   if (times %% 1000 == 0){
#     print(times)
#   }
#   times = times + 1
#   ###Step 1: Sample mu 
#   gamma_c = exp(t(c(beta1_c, beta2_c)) %*% rbind(rep(1,n),x))
#   for (i in 1:n){
#      mu_c[i] = rgamma(1,shape = (lambda_c+y[i]) ,scale = ((lambda_c/gamma_c[i]+1)^(-1)))
#   }
# #   cur = -Inf
# #   while (cur == -Inf){
# #      for (i in 1:n){
# #        mu_c[i] = rgamma(1,shape = (lambda_c+y[i]) ,scale = (lambda_c/gamma_c[i]+1))
# #      }
# #      cur = 0
# #      for (i in 1:n){
# #        cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = lambda_c/gamma_c[i], log = TRUE)
# #      }
# #   }
#   ###Step 2: Sample lambda
#   flag_lambda = F
#   nflag_lambda = 0
#   log_lambda_c = log(lambda_c)
#   while (!flag_lambda){
#     nflag_lambda = nflag_lambda + 1
#     log_lambda_n = log_lambda_c + rnorm(1, sd = 0.1/nflag_lambda)
#     acc_labmda = acc_func1(mu_c, gamma_c, log_lambda_c, log_lambda_n, y)
#     r = runif(1)
#     if (r < acc_labmda){
#       flag_lambda = T
#       lambda[times] = lambda_c = exp(log_lambda_n)
#     }
#   }
#   ###Step 3: Sample beta
#   flag_beta = F
#   nflag_beta = 0
#   while (!flag_beta){
#     nflag_beta = nflag_beta + 1
#     beta1_n = beta1_c + rnorm(1, sd = 1/nflag_beta)
#     beta2_n = beta2_c + rnorm(1, sd = 0.1/nflag_beta)
#     acc_beta = acc_func2(mu_c, beta1_c, beta2_c, beta1_n, lambda_c, beta2_n, x, y)
#     r = runif(1)
#     if (r < acc_beta){
#       flag_beta = T
#       beta1[times] = beta1_c = beta1_n
#       beta2[times] = beta2_c = beta2_n
#     }
#   }
# }
# 
# ##Plot the parameter curve
# par(mfrow = c(1,3))
# plot(beta1, type = "l")
# plot(beta2, type = "l")
# plot(lambda, type = "l")
# pdf(file = "pic/HW5_1/pos0.pdf")
# par(mfrow = c(1,3))
# post_beta1 = beta1[2001:tot]
# post_beta2 = beta2[2001:tot]
# post_lambda = lambda[2001:tot]
# h1 = hist(post_beta1, freq = FALSE, xlab = expression(beta[1]), ylab = "Density")
# lines(h1$mids, h1$density, type = "l", lty = 2, col="red")
# #abline(h = 1, col = "blue")
# h2 = hist(post_beta2, freq = FALSE, xlab = expression(beta[2]), ylab = "Density")
# lines(h2$mids, h2$density, type = "l", lty = 2, col="red")
# #abline(h = 1, col = "blue")
# h3 = hist(post_lambda, freq = FALSE, breaks = 100, xlab = expression(lambda), ylab = "Density")
# lines(h3$mids, h3$density, type = "l", lty = 2, col="red")
# #abline(h = 1, col = "blue")
# dev.off()
# 
# ###Response mean curve
# inv = seq(min(x), max(x), length.out = 1000)
# p_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,1000),inv)
# p_mu = exp(p_ita)
# p_mu_inv = apply(p_mu, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))
# pdf(file = "pic/HW5_1/RMC0.pdf")
# plot(inv, p_mu_inv[1,], col = 1, type = "l", xlab = "Length", ylab = "Number of fault")
# lines(inv, p_mu_inv[2,], col = 1, lty = 2)
# lines(inv, p_mu_inv[3,], col = 1, lty = 2)
# dev.off()
# 
# ###Posterior predictive residuals
# n = length(x)
# p0_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,n),x)
# p0_gamma = exp(p0_ita)
# p0_lambda = post_lambda
# p0_mu = matrix(rep(0,n*(tot-2000)), ncol = 32)
# for (i in 1:8000){
#   for (j in 1:n){
#     p0_mu[i,j] = rgamma(1, shape = p0_lambda[i], scale = p0_gamma[i,j]/p0_lambda[i])
#   }
# }
# p0_y = apply(p0_mu, 1:2, function(x) rpois(1,x))
# p0_e = p0_y - matrix(rep(y,8000), ncol = 32, byrow = TRUE)
# pdf(file = "pic/HW5_1/box0.pdf")
# boxplot(p0_e)
# abline(h = 0, col="red")
# dev.off()
# 
# ###quadratic loss L measure
# p = 1
# var0 = apply(p0_y, 2, function(x) var(x))
# sqd = apply(p0_y, 2, function(x) (y-mean(x))^2)
# L = sum(var0) + p*sum(sqd)








###Hierarchical Model with different priors
rm(list = ls())
data = read.csv("fabric_1.csv", header = FALSE, sep = "")
y = data[,2]
x = data[,1]
y = y[order(x)]
x = x[order(x)]

acc_func1 <- function(mu_c, gamma_c, log_lambda_c, log_lambda_n, y){
  n = length(y)
  cur = new = 0
  lambda_c = exp(log_lambda_c)
  lambda_n = exp(log_lambda_n)
  for (i in 1:n){
    cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_c[i]/lambda_c, log = TRUE)
    new = new + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_n, scale = gamma_c[i]/lambda_n, log = TRUE)
  }
  cur = cur - 2*log(lambda_c+1)
  new = new - 2*log(lambda_n+1)
  if (new == -Inf){
    acc = 0
  }else{
    acc = min(1, exp(new - cur))
  }
  return(acc)
}

acc_func2 <- function(mu_c, beta1_c, beta2_c, beta1_n, lambda_c, beta2_n, x, y){
  n = length(x)
  gamma_c = exp(t(c(beta1_c, beta2_c)) %*% rbind(rep(1,n), x))
  gamma_n = exp(t(c(beta1_n, beta2_n)) %*% rbind(rep(1,n), x))
  cur = new = 0
  for (i in 1:n){
    cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_c[i]/lambda_c, log = TRUE)
    new = new + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = gamma_n[i]/lambda_c, log = TRUE)
  }
  if (new == -Inf){
    acc = 0
  }else{
    acc = min(1, exp(new - cur))
  }
}

tot = 10000
###Test: 
###tot = 10
beta1 = beta2 = lambda = rep(0,tot)
n = length(x)

times = 0
beta1_c = beta2_c = lambda_c = 0.01
beta1_n = beta2_n = lambda_n = 0.01
mu_c = rep(0,n)

while (times < tot){
  if (times %% 1000 == 0){
    print(times)
  }
  times = times + 1
  ###Step 1: Sample mu 
  gamma_c = exp(t(c(beta1_c, beta2_c)) %*% rbind(rep(1,n),x))
  for (i in 1:n){
    mu_c[i] = rgamma(1,shape = (lambda_c+y[i]) ,scale = ((lambda_c/gamma_c[i]+1)^(-1)))
  }
  #   cur = -Inf
  #   while (cur == -Inf){
  #      for (i in 1:n){
  #        mu_c[i] = rgamma(1,shape = (lambda_c+y[i]) ,scale = (lambda_c/gamma_c[i]+1))
  #      }
  #      cur = 0
  #      for (i in 1:n){
  #        cur = cur + dpois(y[i], lambda = mu_c[i], log = TRUE) + dgamma(mu_c[i], shape = lambda_c, scale = lambda_c/gamma_c[i], log = TRUE)
  #      }
  #   }
  ###Step 2: Sample lambda
  flag_lambda = F
  nflag_lambda = 0
  log_lambda_c = log(lambda_c)
  while (!flag_lambda){
    nflag_lambda = nflag_lambda + 1
    log_lambda_n = log_lambda_c + rnorm(1, sd = 0.1/nflag_lambda)
    acc_labmda = acc_func1(mu_c, gamma_c, log_lambda_c, log_lambda_n, y)
    r = runif(1)
    if (r < acc_labmda){
      flag_lambda = T
      lambda[times] = lambda_c = exp(log_lambda_n)
    }
  }
  ###Step 3: Sample beta
  flag_beta = F
  nflag_beta = 0
  while (!flag_beta){
    nflag_beta = nflag_beta + 1
    beta1_n = beta1_c + rnorm(1, sd = 1/nflag_beta)
    beta2_n = beta2_c + rnorm(1, sd = 0.1/nflag_beta)
    acc_beta = acc_func2(mu_c, beta1_c, beta2_c, beta1_n, lambda_c, beta2_n, x, y)
    r = runif(1)
    if (r < acc_beta){
      flag_beta = T
      beta1[times] = beta1_c = beta1_n
      beta2[times] = beta2_c = beta2_n
    }
  }
}

##Plot the parameter curve
par(mfrow = c(1,3))
plot(beta1, type = "l")
plot(beta2, type = "l")
plot(lambda, type = "l")
pdf(file = "pic/HW5_1/pos1.pdf")
par(mfrow = c(1,3))
post_beta1 = beta1[2001:tot]
post_beta2 = beta2[2001:tot]
post_lambda = lambda[2001:tot]
h1 = hist(post_beta1, freq = FALSE, xlab = expression(beta[1]), ylab = "Density")
lines(h1$mids, h1$density, type = "l", lty = 2, col="red")
#abline(h = 1, col = "blue")
h2 = hist(post_beta2, freq = FALSE, xlab = expression(beta[2]), ylab = "Density")
lines(h2$mids, h2$density, type = "l", lty = 2, col="red")
#abline(h = 1, col = "blue")
h3 = hist(post_lambda, freq = FALSE, xlab = expression(lambda), ylab = "Density")
lines(h3$mids, h3$density, type = "l", lty = 2, col="red")
#abline(h = 1, col = "blue")
dev.off()

###Response mean curve
inv = seq(min(x), max(x), length.out = 1000)
p_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,1000),inv)
p_mu = exp(p_ita)
p_mu_inv = apply(p_mu, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))
pdf(file = "pic/HW5_1/RMC1.pdf")
plot(inv, p_mu_inv[1,], col = 1, type = "l", xlab = "Length", ylab = "Number of fault")
lines(inv, p_mu_inv[2,], col = 1, lty = 2)
lines(inv, p_mu_inv[3,], col = 1, lty = 2)
dev.off()

###Posterior predictive residuals
n = length(x)
p0_ita = cbind(post_beta1, post_beta2) %*% rbind(rep(1,n),x)
p0_gamma = exp(p0_ita)
p0_lambda = post_lambda
p0_mu = matrix(rep(0,n*(tot-2000)), ncol = 32)
for (i in 1:8000){
  for (j in 1:n){
    p0_mu[i,j] = rgamma(1, shape = p0_lambda[i], scale = p0_gamma[i,j]/p0_lambda[i])
  }
}
p0_y = apply(p0_mu, 1:2, function(x) rpois(1,x))
p0_e = p0_y - matrix(rep(y,8000), ncol = 32, byrow = TRUE)
pdf(file = "pic/HW5_1/box1.pdf")
boxplot(p0_e)
abline(h = 0, col="red")
dev.off()

###quadratic loss L measure
p = 1
var0 = apply(p0_y, 2, function(x) var(x))
sqd = apply(p0_y, 2, function(x) (y-mean(x))^2)
L = sum(var0) + p*sum(sqd)