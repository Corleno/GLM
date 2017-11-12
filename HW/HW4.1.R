rm(list = ls())
x = c(0, 62.5, 125, 250, 500)
y1 = c(15, 17, 22, 38, 144)
y2 = c(1, 0, 7, 59, 132)
y3 = c(281, 225, 283, 202, 9)
m = c(297, 242, 312, 299, 285)
# 
# g1 = glm(cbind(y1,m-y1) ~ x, family = binomial)
# g2 = glm(cbind(y2,y3) ~ x, family = binomial)
# hat_coe1 = g1$coefficients
# hat_coe2 = g2$coefficients
# 
# 
# 
# inv = seq(0,500,length.out = 1000)
# ita1 = t(hat_coe1) %*% rbind(rep(1,1000),inv)
# ita2 = t(hat_coe2) %*% rbind(rep(1,1000),inv)
# hat_pi1 = 1/(1+exp(-ita1))
# hat_pi2 = 1/(1+exp(-ita2))*(1-hat_pi1)
# hat_pi3 = 1- hat_pi1 - hat_pi2
# 
# ###Plot the response curves
# pdf(file = "pic/HW4_1/ERC.pdf")
# plot(inv, hat_pi1, col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated response curves")
# lines(inv, hat_pi2, col = 2)
# lines(inv, hat_pi3, col = 3)
# legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))
# dev.off()

###Ratio
logistic_func <- function(x){
  return(1/(1+exp(-x)))
}

ratio_func <- function(alpha_c, beta_c, alpha_n, beta_n, x, y1, y2, y3, m){
  n = length(x)
  ita1 = t(c(alpha_c[1], beta_c[1])) %*% rbind(rep(1,n),x)
  ita2 = t(c(alpha_c[2], beta_c[2])) %*% rbind(rep(1,n),x)
  p1 = logistic_func(ita1)
  p2 = logistic_func(ita2)
  cur = (sum(y1*log(p1)+(m-y1)*log(1-p1)+y2*log(p2)+y3*log(1-p2)))
  ita1_n = t(c(alpha_n[1], beta_n[1])) %*% rbind(rep(1,n),x)
  ita2_n = t(c(alpha_n[2], beta_n[2])) %*% rbind(rep(1,n),x)
  p1_n = logistic_func(ita1_n)
  p2_n = logistic_func(ita2_n)
  new = (sum(y1*log(p1_n)+(m-y1)*log(1-p1_n)+y2*log(p2_n)+y3*log(1-p2_n)))
  if (new == -Inf){
    ratio = 0
  }else{
    r = exp(new - cur)
    ratio = min(1,r)
  }
  return(ratio)
}

###MCMC
tot = 5000
alpha1 = alpha2 = beta1 = beta2 = rep(0,tot)

times = 0
alpha_c1 = alpha_c2 = beta_c1 = beta_c2 = 0
alpha_n1 = alpha_n2 = beta_n1 = beta_n2 = 0
while (times<tot){
  times = times+1
  flag = FALSE
  nflag = 0
  while (!flag){
    nflag = nflag+1
    variance = 1/flag
    alpha_n1 = alpha_c1 + rnorm(1,sd = 1/nflag)
    alpha_n2 = alpha_c2 + rnorm(1,sd = 1/nflag)
    beta_n1 = beta_c1 + rnorm(1,sd = 1/nflag)
    beta_n2 = beta_c2 + rnorm(1,sd = 1/nflag)
    r = runif(1)
    acc = ratio_func(c(alpha_c1, alpha_c2), c(beta_c1,beta_c2), c(alpha_n1,alpha_n2), c(beta_n1, beta_n2), x, y1, y2, y3, m)
    if (r <= acc){
      flag=TRUE
      alpha1[times] = alpha_c1 = alpha_n1
      alpha2[times] = alpha_c2 = alpha_n2
      beta1[times] = beta_c1 = beta_n1
      beta2[times] = beta_c2 = beta_n2
    }
  }
}

###
# hat_alpha1 = mean(alpha1[101:1000])
# hat_alpha2 = mean(alpha2[101:1000])
# hat_beta1 = mean(beta1[101:1000])
# hat_beta2 = mean(beta2[101:1000])
inv = seq(0,500,length.out = 1000)

# ita1 = t(c(hat_alpha1, hat_beta1)) %*% rbind(rep(1,1000),inv)
# ita2 = t(c(hat_beta1, hat_beta2)) %*% rbind(rep(1,1000),inv)
# hat_pi1 = 1/(1+exp(-ita1))
# hat_pi2 = 1/(1+exp(-ita2))*(1-hat_pi1)
# hat_pi3 = 1- hat_pi1 - hat_pi2

ita1_m = cbind(alpha1[101:tot], beta1[101:tot])%*%rbind(rep(1,1000),inv)
ita2_m = cbind(alpha2[101:tot], beta2[101:tot])%*%rbind(rep(1,1000),inv)
hat_pi1_v = 1/(1+exp(-ita1_m))
hat_pi2_v = 1/(1+exp(-ita2_m))*(1-hat_pi1_v)
hat_pi3_v = 1- hat_pi1_v - hat_pi2_v
hat_pi1_inv = apply(hat_pi1_v, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))
hat_pi2_inv = apply(hat_pi2_v, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))
hat_pi3_inv = apply(hat_pi3_v, 2, function(x) c(mean(x),quantile(x,0.025),quantile(x,0.975)))


###Plot the response curves
pdf(file = "pic/HW4_1/BERC.pdf")
plot(inv, hat_pi1_inv[1,], col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated response curves in Bayesion")
lines(inv, hat_pi1_inv[2,], col = 1, lty = 2)
lines(inv, hat_pi1_inv[3,], col = 1, lty = 2)
lines(inv, hat_pi2_inv[1,], col = 2)
lines(inv, hat_pi2_inv[2,], col = 2, lty = 2)
lines(inv, hat_pi2_inv[3,], col = 2, lty = 2)
lines(inv, hat_pi3_inv[1,], col = 3)
lines(inv, hat_pi3_inv[2,], col = 3, lty = 2)
lines(inv, hat_pi3_inv[3,], col = 3, lty = 2)
legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))
dev.off()