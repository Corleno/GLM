rm(list = ls())
M_len = c(1.30,1.32,1.32,1.40,1.42,1.42,1.47,1.47,1.50,1.52,1.63,1.65,1.65,1.65,1.65,1.68,1.70,1.73,1.78,1.78,1.80,1.85,1.93,1.93,1.98,2.03,2.03,2.31,2.36,2.46,3.25,3.28,3.33,3.56,3.58,3.66,3.68,3.71,3.89)
M_cat = c("I","F","F","F","I","F","I","F","I","I","I","O","O","I","F","F","I","O","F","O","F","F","I","F","I","F","F","F","F","F","O","O","F","F","F","F","O","F","F")
F_len = c(1.24,1.30,1.45,1.45,1.55,1.60,1.60,1.65,1.78,1.78,1.80,1.88,2.16,2.26,2.31,2.36,2.39,2.41,2.44,2.56,2.67,2.72,2.79,2.84)
F_cat = c("I","I","I","O","I","I","I","F","I","O","I","I","F","F","F","F","F","F","F","O","F","I","F","F")
x = len = sort(unique(c(M_len, F_len)))
n = length(len)
In = Fi = Ot = rep(0,n)
for (i in 1:n){
  In[i] = sum(M_cat[M_len == len[i]] == "I") + sum(F_cat[F_len == len[i]] == "I")
  Fi[i] = sum(M_cat[M_len == len[i]] == "F") + sum(F_cat[F_len == len[i]] == "F")
  Ot[i] = sum(M_cat[M_len == len[i]] == "O") + sum(F_cat[F_len == len[i]] == "O")
}
M = In + Fi + Ot

tot = 5000
alpha1 = alpha2 = beta1 = beta2 = rep(0,tot)
alpha_c1 = alpha_c2 = beta_c1 = beta_c2 = 0
alpha_n1 = alpha_n2 = beta_n1 = beta_n2 = 0
times = 0

ratio_func <- function(alpha_c1, alpha_c2, beta_c1, beta_c2, alpha_n1, alpha_n2, beta_n1, beta_n2, x, In, Fi, Ot, M){
  n = length(x)
  ita1 = t(c(alpha_c1, beta_c1)) %*% rbind(rep(1, n),x)
  ita2 = t(c(alpha_c2, beta_c2)) %*% rbind(rep(1, n),x)
  cur = (sum(In*ita1 +  Fi*ita2 - M*log(1 + exp(ita1) + exp(ita2))))
  ita1_n = t(c(alpha_n1, beta_n1)) %*% rbind(rep(1, n),x)
  ita2_n = t(c(alpha_n2, beta_n2)) %*% rbind(rep(1, n),x)
  new = (sum(In*ita1_n +  Fi*ita2_n - M*log(1 + exp(ita1_n) + exp(ita2_n))))
  if (new == -Inf){
    ratio = 0
  }else{
    r = exp(new - cur)
    #r = exp(new-cur)*dnorm(alpha2)*dnorm(beta2)/dnorm(alpha1)/dnorm(beta1)
    ratio = min(1, r)
  }
  return(ratio)
}

while (times < tot){
  times = times + 1
  flag = FALSE
  nflag = 0
  while (!flag){
    nflag = nflag + 1
    alpha_n1 = alpha_c1 + rnorm(1, sd = 1/nflag)
    alpha_n2 = alpha_c2 + rnorm(1, sd = 1/nflag)
    beta_n1 = beta_c1 + rnorm(1, sd = 1/nflag)
    beta_n2 = beta_c2 + rnorm(1, sd = 1/nflag)
    r = runif(1)
    acc = ratio_func(alpha_c1, alpha_c2, beta_c1, beta_c2, alpha_n1, alpha_n2, beta_n1, beta_n2, x, In, Fi, Ot, M)
    if (r < acc){
      flag = TRUE
      alpha1[times] = alpha_c1 = alpha_n1
      alpha2[times] = alpha_c2 = alpha_n2
      beta1[times] = beta_c1 = beta_n1
      beta2[times] = beta_c2 = beta_n2
    }
  }
}

###Point and interval estimates of response probability
inv = seq(1,3,length.out = 100)
ita1 = cbind(alpha1[101:tot], beta1[101:tot]) %*% rbind(rep(1, 100),inv)
ita2 = cbind(alpha2[101:tot], beta2[101:tot]) %*% rbind(rep(1, 100),inv)

pi_3 = (1+exp(ita1)+exp(ita2))^(-1)
pi_1 = exp(ita1)*pi_3
pi_2 = exp(ita2)*pi_3

curve1 = apply(pi_1, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve2 = apply(pi_2, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve3 = apply(pi_3, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))

# pdf(file = "pic/HW4_2/BERC.pdf")
# plot(inv, curve1[1,], col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated response curves in Bayesion")
# lines(inv, curve1[2,], col = 1, lty = 2)
# lines(inv, curve1[3,], col = 1, lty = 2)
# lines(inv, curve2[1,], col = 2)
# lines(inv, curve2[2,], col = 2, lty = 2)
# lines(inv, curve2[3,], col = 2, lty = 2)
# lines(inv, curve3[1,], col = 3)
# lines(inv, curve3[2,], col = 3, lty = 2)
# lines(inv, curve3[3,], col = 3, lty = 2)
# legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))
# dev.off()
# pdf(file = "pic/HW4_2/BERC_N.pdf")
# plot(inv, curve1[1,], col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated response curves in Bayesion")
# lines(inv, curve1[2,], col = 1, lty = 2)
# lines(inv, curve1[3,], col = 1, lty = 2)
# lines(inv, curve2[1,], col = 2)
# lines(inv, curve2[2,], col = 2, lty = 2)
# lines(inv, curve2[3,], col = 2, lty = 2)
# lines(inv, curve3[1,], col = 3)
# lines(inv, curve3[2,], col = 3, lty = 2)
# lines(inv, curve3[3,], col = 3, lty = 2)
# legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))
# dev.off()

###Curves for parameters
pdf(file = "pic/HW4_2/pc.pdf")
par(mfrow = c(2,2))
plot(alpha1, type = 'l')
plot(alpha2, type = 'l')
plot(beta1, type = 'l')
plot(beta2,type = 'l')
dev.off()
# pdf(file = "pic/HW4_2/pc_N.pdf")
# par(mfrow = c(2,2))
# plot(alpha1, type = 'l')
# plot(alpha2, type = 'l')
# plot(beta1, type = 'l')
# plot(beta2,type = 'l')
# dev.off()



























###Consider both length and Gender
rm(list = ls())
M_len = c(1.30,1.32,1.32,1.40,1.42,1.42,1.47,1.47,1.50,1.52,1.63,1.65,1.65,1.65,1.65,1.68,1.70,1.73,1.78,1.78,1.80,1.85,1.93,1.93,1.98,2.03,2.03,2.31,2.36,2.46,3.25,3.28,3.33,3.56,3.58,3.66,3.68,3.71,3.89)
M_cat = c("I","F","F","F","I","F","I","F","I","I","I","O","O","I","F","F","I","O","F","O","F","F","I","F","I","F","F","F","F","F","O","O","F","F","F","F","O","F","F")
F_len = c(1.24,1.30,1.45,1.45,1.55,1.60,1.60,1.65,1.78,1.78,1.80,1.88,2.16,2.26,2.31,2.36,2.39,2.41,2.44,2.56,2.67,2.72,2.79,2.84)
F_cat = c("I","I","I","O","I","I","I","F","I","O","I","I","F","F","F","F","F","F","F","O","F","I","F","F")
x = len = sort(unique(c(M_len, F_len)))
n = length(len)
In_M = In_F = Fi_M = Fi_F = Ot_M = Ot_F = rep(0,n)
for (i in 1:n){
  In_M[i] = sum(M_cat[M_len == len[i]] == "I")
  In_F[i] = sum(F_cat[F_len == len[i]] == "I")
  Fi_M[i] = sum(M_cat[M_len == len[i]] == "F")
  Fi_F[i] = sum(F_cat[F_len == len[i]] == "F")
  Ot_M[i] = sum(M_cat[M_len == len[i]] == "O")
  Ot_F[i] = sum(F_cat[F_len == len[i]] == "O")
}
M_M = In_M + Fi_M + Ot_M
M_F = In_F + Fi_F + Ot_F

tot = 5000
alpha1 = alpha2 = beta1 = beta2 = gamma1 = gamma2 = rep(0,tot)
alpha_c1 = alpha_c2 = beta_c1 = beta_c2 = gamma_c1 = gamma_c2 = 0
alpha_n1 = alpha_n2 = beta_n1 = beta_n2 = gamma_n1 = gamma_n2 = 0
times = 0

ratio_func <- function(alpha_c1, alpha_c2, beta_c1, beta_c2, gamma_c1, gamma_c2, alpha_n1, alpha_n2, beta_n1, beta_n2, gamma_n1, gamma_n2, x, In_M, In_F, Fi_M, Fi_F, Ot_M, Ot_F, M_M, M_F){
  n = length(x)
  ita10 = t(c(alpha_c1, beta_c1)) %*% rbind(rep(1, n),x)
  ita11 = t(c(alpha_c1, beta_c1, gamma_c1)) %*% rbind(rep(1, n),x,rep(1, n))
  ita20 = t(c(alpha_c2, beta_c2)) %*% rbind(rep(1, n),x)
  ita21 = t(c(alpha_c2, beta_c2, gamma_c2)) %*% rbind(rep(1, n),x,rep(1, n))
  cur = (sum(In_M*ita10+In_F*ita11+Fi_M*ita20+Fi_F*ita21-M_M*log(1+exp(ita10)+exp(ita20))-M_F*log(1+exp(ita11)+exp(ita21))))
  ita10_n = t(c(alpha_n1, beta_n1)) %*% rbind(rep(1, n),x)
  ita11_n = t(c(alpha_n1, beta_n1, gamma_n1)) %*% rbind(rep(1, n),x,rep(1, n))
  ita20_n = t(c(alpha_n2, beta_n2)) %*% rbind(rep(1, n),x)
  ita21_n = t(c(alpha_n2, beta_n2, gamma_n2)) %*% rbind(rep(1, n),x,rep(1, n))
  new = (sum(In_M*ita10_n+In_F*ita11_n+Fi_M*ita20_n+Fi_F*ita21_n-M_M*log(1+exp(ita10_n)+exp(ita20_n))-M_F*log(1+exp(ita11_n)+exp(ita21_n))))
  if (new == -Inf){
    ratio = 0
  }else{
    r = exp(new - cur)
    ratio = min(1, r)
  }
  return(ratio)
}

while (times < tot){
  times = times + 1
  flag = FALSE
  nflag = 0
  while (!flag){
    nflag = nflag + 1
    alpha_n1 = alpha_c1 + rnorm(1, sd = 1/nflag)
    alpha_n2 = alpha_c2 + rnorm(1, sd = 1/nflag)
    beta_n1 = beta_c1 + rnorm(1, sd = 1/nflag)
    beta_n2 = beta_c2 + rnorm(1, sd = 1/nflag)
    gamma_n1 = gamma_c1 + rnorm(1, sd = 1/nflag)
    gamma_n2 = gamma_c2 + rnorm(1, sd = 1/nflag)
    r = runif(1)
    acc = ratio_func(alpha_c1, alpha_c2, beta_c1, beta_c2, gamma_c1, gamma_c2, alpha_n1, alpha_n2, beta_n1, beta_n2, gamma_n1, gamma_n2, x, In_M, In_F, Fi_M, Fi_F, Ot_M, Ot_F, M_M, M_F)
    if (r < acc){
      flag = TRUE
      alpha1[times] = alpha_c1 = alpha_n1
      alpha2[times] = alpha_c2 = alpha_n2
      beta1[times] = beta_c1 = beta_n1
      beta2[times] = beta_c2 = beta_n2
      gamma1[times] = gamma_c1 = gamma_n1
      gamma2[times] = gamma_c2 = gamma_n2
    }
  }
}

###Point and interval estimates of response probability
inv = seq(1,3,length.out = 100)
ita10 = cbind(alpha1[101:tot], beta1[101:tot]) %*% rbind(rep(1, 100),inv)
ita11 = cbind(alpha1[101:tot], beta1[101:tot], gamma1[101:tot]) %*% rbind(rep(1, 100),inv,rep(1,100))
ita20 = cbind(alpha2[101:tot], beta2[101:tot]) %*% rbind(rep(1, 100),inv)
ita21 = cbind(alpha2[101:tot], beta2[101:tot], gamma2[101:tot]) %*% rbind(rep(1, 100),inv,rep(1, 100))


pi_30 = (1+exp(ita10)+exp(ita20))^(-1)
pi_31 = (1+exp(ita11)+exp(ita21))^(-1)
pi_10 = exp(ita10)*pi_30
pi_20 = exp(ita20)*pi_30
pi_11 = exp(ita11)*pi_31
pi_21 = exp(ita21)*pi_31

curve10 = apply(pi_10, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve20 = apply(pi_20, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve30 = apply(pi_30, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve11 = apply(pi_11, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve21 = apply(pi_21, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))
curve31 = apply(pi_31, 2, function(x) c(mean(x), quantile(x,0.025), quantile(x,0.975)))



pdf(file = "pic/HW4_2/BERC_1.pdf")
par(mfrow = c(1,2))
plot(inv, curve10[1,], col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated male response curves", cex.main = 1)
lines(inv, curve10[2,], col = 1, lty = 2)
lines(inv, curve10[3,], col = 1, lty = 2)
lines(inv, curve20[1,], col = 2)
lines(inv, curve20[2,], col = 2, lty = 2)
lines(inv, curve20[3,], col = 2, lty = 2)
lines(inv, curve30[1,], col = 3)
lines(inv, curve30[2,], col = 3, lty = 2)
lines(inv, curve30[3,], col = 3, lty = 2)
legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))


plot(inv, curve11[1,], col = 1, type = "l", ylim = c(0,1), xlab = "Concentration", ylab = "Response probability", main = "Estimated female response curves", cex.main = 1)
lines(inv, curve11[2,], col = 1, lty = 2)
lines(inv, curve11[3,], col = 1, lty = 2)
lines(inv, curve21[1,], col = 2)
lines(inv, curve21[2,], col = 2, lty = 2)
lines(inv, curve21[3,], col = 2, lty = 2)
lines(inv, curve31[1,], col = 3)
lines(inv, curve31[2,], col = 3, lty = 2)
lines(inv, curve31[3,], col = 3, lty = 2)
legend("topright", legend = c(expression(hat(pi)[1]), expression(hat(pi)[2]), expression(hat(pi)[3])), lty = c(1,1,1), col = c(1,2,3))
dev.off()