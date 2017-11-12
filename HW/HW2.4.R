#!/Library/Frameworks/R.framework/Versions/3.2/Resources/bin/Rscript
rm(list = ls())
data = read.table('ceriodaphnia.txt', header = FALSE, sep="")
y = data[,1]
x1 = data[,2]
x2 = data[,3]

library("plot3D")
###Plot the relation between the number of organisms and the concentration of jet fuel
col_index = x2 + 1
pdf(file = "pic/HW2_4/Identity_link.pdf")
plot(x1,y, col = col_index, xlab = "concentration of jet fuel", ylab = "g(the number of organisms)", main = "Identity link")
legend("topright", legend=c("strain = 0","strain = 1"), pch = c(1,1), col = c(1,2),  cex = 1)
dev.off()

pdf(file = "pic/HW2_4/Log_link.pdf")
plot(x1,log(y), col = col_index, xlab = "concentration of jet fuel", ylab = "g(the number of organisms)", main = "Log link")
legend("topright", legend=c("strain = 0","strain = 1"), pch = c(1,1), col = c(1,2),  cex = 1)
dev.off()

pdf(file = "pic/HW2_4/Inverse_link.pdf")
plot(x1,1/y, col = col_index, xlab = "concentration of jet fuel", ylab = "g(the number of organisms)", main = "Inverse link")
legend("topright", legend=c("strain = 0","strain = 1"), pch = c(1,1), col = c(1,2),  cex = 1)
dev.off()

pdf(file = "pic/HW2_4/Square_root_link.pdf")
plot(x1,sqrt(y), col = col_index, xlab = "concentration of jet fuel", ylab = "g(the number of organisms)", main = "Square root link")
legend("topright", legend=c("strain = 0","strain = 1"), pch = c(1,1), col = c(1,2),  cex = 1)
dev.off()

glm.poi.ide = glm(y ~ x1 + x2, family = poisson(link = "identity"))
glm.poi.ide$deviance
glm.poi.ide$aic
glm.poi.ide$aic - 2*3 + 2*log(32)
glm.poi.log = glm(y ~ x1 + x2, family = poisson(link = "log"))
glm.poi.log$deviance
glm.poi.log$aic
glm.poi.log$aic - 2*3 + 2*log(32)
glm.poi.inv = glm(y ~ x1 + x2, family = poisson(link = "inverse"))
glm.poi.inv$deviance
glm.poi.inv$aic
glm.poi.inv$aic - 2*3 + 2*log(32)
glm.poi.squ = glm(y ~ x1 + x2, family = poisson(link = "sqrt"))
glm.poi.squ$deviance
glm.poi.squ$aic
glm.poi.squ$aic - 2*3 + 2*log(32)

###Residual Computation
#Pearson residuals
pdf(file = "pic/HW2_4/p_res_ide.pdf")
plot(glm.poi.ide$fitted.values, residuals(glm.poi.ide, type = "pearson"), xlab = "Fitted value", ylab= "Pearson residual", main = "Identity link")
dev.off()
pdf(file = "pic/HW2_4/p_res_log.pdf")
plot(glm.poi.log$fitted.values, residuals(glm.poi.log, type = "pearson"), xlab = "Fitted value", ylab= "Pearson residual", main = "Log link")
dev.off()
pdf(file = "pic/HW2_4/p_res_inv.pdf")
plot(glm.poi.inv$fitted.values, residuals(glm.poi.inv, type = "pearson"), xlab = "Fitted value", ylab= "Pearson residual", main = "Inverse link")
dev.off()
pdf(file = "pic/HW2_4/p_res_squ.pdf")
plot(glm.poi.squ$fitted.values, residuals(glm.poi.squ, type = "pearson"), xlab = "Fitted value", ylab= "Pearson residual", main = "Square root link")
dev.off()


#Deviance residuals
pdf(file = "pic/HW2_4/d_res_ide.pdf")
plot(glm.poi.ide$fitted.values, residuals(glm.poi.ide, type = "deviance"), xlab = "Fitted value", ylab= "Deviance residual", main = "Identity link")
dev.off()
pdf(file = "pic/HW2_4/d_res_log.pdf")
plot(glm.poi.log$fitted.values, residuals(glm.poi.log, type = "deviance"), xlab = "Fitted value", ylab= "Deviance residual", main = "Log link")
dev.off()
pdf(file = "pic/HW2_4/d_res_inv.pdf")
plot(glm.poi.inv$fitted.values, residuals(glm.poi.inv, type = "deviance"), xlab = "Fitted value", ylab= "Deviance residual", main = "Inverse link")
dev.off()
pdf(file = "pic/HW2_4/d_res_squ.pdf")
plot(glm.poi.squ$fitted.values, residuals(glm.poi.squ, type = "deviance"), xlab = "Fitted value", ylab= "Deviance residual", main = "Square root link")
dev.off()

###regression functions 
est_f <- function(beta, X){
  mu <- exp(beta^T%*%X)
  return(mu)
}

x_inv = seq(0,2,length.out = 20) 
X0 = rbind(rep(1,20), x_inv, rep(0,20))
X1 = rbind(rep(1,20), x_inv, rep(1,20))
beta = glm.poi.log$coefficients
Y0 = est_f(beta,X0)
Y1 = est_f(beta,X1)

pdf(file = "pic/HW2_4/est_func.pdf")
plot(x_inv, Y0, col = 1, type = "l", xlab = "concentation of jet fuel", ylab = "the number of organisms", main = "Estimated regression functions")
lines(x_inv, Y1, col = 2, type = "l")
legend("topright", lty = 1, col = c(1,2), legend = c("strain = 0","strain = 1"), cex = 1)
dev.off()