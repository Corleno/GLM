rm(list = ls())
x = c(1.6907, 1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
m = c(59,60,62,56,63,59,62,60)
y = c(6,13,18,28,52,53,61,60)

###(a)
round(y/m,2)

glm_1 = glm(cbind(y,m-y)~x, family = binomial(link = logit))
summary(glm_1)
glm_1$deviance
round(resid(glm_1),2)
round(glm_1$fitted.values,2)
beta = glm_1$coefficients
link_1 <- function(x, beta){
  n = length(x)
  eta = t(beta) %*% rbind(rep(1,n), x)
  out = exp(eta)/(1+exp(eta))
  return(out)
}
inv = seq(1.65,1.9,length.out = 50)
pred_y = link_1(inv, beta)
pdf(file = "pic/HW3_1/a1.pdf")
plot(inv, pred_y, type = 'l', main = "Dose response curve with logit link")
points(x, y/m, col = 2)
points(x, glm_1$fitted.values, col = 3)
dev.off()

glm_2 = glm(cbind(y,m-y)~x, family = binomial(link = probit))
summary(glm_2)
glm_2$deviance
round(resid(glm_2),2)
round(glm_2$fitted.values,2)
beta = glm_2$coefficients
link_2 <- function(x, beta){
  n = length(x)
  eta = t(beta) %*% rbind(rep(1,n), x)
  out = pnorm(eta)
  return(out)
}
inv = seq(1.65,1.9,length.out = 50)
pred_y = link_2(inv, beta)
pdf(file = "pic/HW3_1/a2.pdf")
plot(inv,pred_y,type = 'l', main = "Dose response curve with probit link")
points(x, y/m, col = '2')
points(x, glm_2$fitted.values, col = 3)
dev.off()

glm_3 = glm(cbind(y,m-y)~x, family = binomial(link = cloglog))
summary(glm_3)
glm_3$deviance
round(resid(glm_3),2)
round(glm_3$fitted.values,2)
beta = glm_3$coefficients
link_3 <- function(x, beta){
  n = length(x)
  eta = t(beta) %*% rbind(rep(1,n), x)
  out = 1 - exp(-exp(eta))
  return(out)
}
inv = seq(1.65,1.9,length.out = 50)
pred_y = link_3(inv, beta)
pdf("pic/HW3_1/a3.pdf")
plot(inv,pred_y,type = 'l', main = "Dose response curve with complementary log-log link")
points(x, y/m, col = '2')
points(x, glm_3$fitted.values, col = 3)
dev.off()

###(b)
rm(list = ls())
f_alpha <- function(x, alpha){
  y = exp(alpha * x)/(1 + exp(x))^alpha
}
inv = seq(-10,10,length.out = 100) 
alpha = 0.5
pdf(file = "pic/HW3_1/b.pdf")
y1 = f_alpha(inv, alpha)
plot(inv, y1, type = 'l', xlab = expression(eta), ylab = expression(pi), main = paste0("Comparison with different ", expression(alpha)))
alpha = 1
y2 = f_alpha(inv,alpha)
lines(inv, y2, type = 'l', col = 2)
alpha = 2
y3 = f_alpha(inv,alpha)
lines(inv, y3, type = 'l', col = 3)
legend("topleft",legend = c(expression(alpha == 0.5), expression(alpha == 1), expression(alpha == 2)), lty = c(1,1,1), col = c(1,2,3), cex = 1)
dev.off()

###(c)

beta = c(-113.625, 62.5)
alpha = 0.279

f_ba <- function(x, beta ,alpha){
  n = length(x)
  X = rbind(rep(1,n), x)
  eta = beta %*% X
  pi = exp(alpha * eta)/((1+exp(eta))^alpha)
  return(pi)
}
hat_pi = f_ba(x, beta, alpha)
round(hat_pi,2)
inv = seq(1.65,1.9,length.out = 50) 
pred_y = f_ba(inv,beta, alpha)
pdf(file = "pic/HW3_1/c.pdf")
plot(inv, pred_y, type = 'l', xlab = "x", ylab = "Estimated dose-response")
points(x, y/m)
dev.off()

Dev_f <- function(y ,mu, m){
  if (m == y){
    out = 2*(y*log(y/mu))
  }else{
    out = 2*(y*log(y/mu) + (m-y)*log((m-y)/(m-mu)))
  }
  return(out)
}

n = length(x)
mu = rep(0,n)
hat_mu = m * hat_pi
dev_res = rep(0,n)
for (i in 1:n){
  dev_res[i] = sqrt(Dev_f(y[i], hat_mu[i], m[i])) * sign(y[i] - hat_mu[i])
}
round(dev_res,2)
sum(dev_res^2)

log_lik <- function(y, m, pi){
  out = y*log(pi/(1-pi)) + m*log(1-pi) + log(choose(m,y))
  return(sum(out))
}
AIC = -2*log_lik(y,m,hat_pi) + 2*3
BIC = -2*log_lik(y,m,hat_pi) + 3*log(8)
