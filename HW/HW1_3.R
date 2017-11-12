###3(b)
rm(list = ls())
n=9
y = c(-0.774,0.597,7.575, 0.397,-0.865,-0.318,-0.125,0.961,1.039)

Lik_func <- function(theta,y){
  val = -n*log(pi)-sum(log(1+(y-theta)^2))
  return(val)
}

Sco_func <- function(theta,y){
  val = sum(-2*(y-theta)/(1+(y-theta)^2))
  return(val)
}

Ob_fish_func <- function(theta,y){
  val = -sum((2-2*(y-theta)^2)/(1+(y-theta)^2)^2)
  return(val)
}

Ex_fish_func <- function(theta,y){
  val = -2/length(y)
}

###Initialize theta
theta_0 = 0
theta_1 = 1
k = 1
#theta_arr = lik_arr = c()
###Newton Raphson
while (abs(theta_1-theta_0)>0.01){
#  theta_arr[k] = theta_1
#  lik_arr[k] = Lik_func(theta_1,y)
  theta_0 = theta_1
  k = k+1
  theta_1 = theta_0 + Sco_func(theta_0,y)/Ob_fish_func(theta_0,y)
}
print(paste0("k = ",k," beta = ",theta_1))
###"k = 4 theta = 0.17918627636346" suggests using the Newton Raphson method with threshold 0.01 and starting point 1, the theta can be estimated at the 4th iteration and equal to 0.179.

###Scoring Method
theta_0 = 0
theta_1 = 1
k = 1
#theta_arr = lik_arr = c()
###Scoring Method
while (abs(theta_1-theta_0)>0.01){
#  theta_arr[k] = theta_1
#  lik_arr[k] = Lik_func(theta_1,y)
  theta_0 = theta_1
  k = k+1
  theta_1 = theta_0 + Sco_func(theta_0,y)/Ex_fish_func(theta_0,y)
}
print(paste0("k = ",k," beta = ",theta_1))
###"k = 51978 theta = 0.187783745581135" suggests using the Scoring method with threshold 0.01 and starting point 1, the theta can be estimated at the 51978th iteration and equal to 0.188.

###Theta VS Likelihood
theta_arr = seq(from = -1, to = 6, by = 0.01)
lik_arr = c()
for (i in 1:length(theta_arr)){
  lik_arr[i] = Lik_func(theta_arr[i],y)
}
plot(theta_arr, lik_arr, type = "l", main = "Theta VS Likelihood", xlab = "Theta", ylab = "Likelihood")



###3(c)
rm(list = ls())
n = 3
y = c(0,5,9)

Lik_func <- function(theta,y){
  val = -n*log(pi)-sum(log(1+(y-theta)^2))
  return(val)
}

Sco_func <- function(theta,y){
  val = sum(-2*(y-theta)/(1+(y-theta)^2))
  return(val)
}

Ob_fish_func <- function(theta,y){
  val = -sum((2-2*(y-theta)^2)/(1+(y-theta)^2)^2)
  return(val)
}

Ex_fish_func <- function(theta,y){
  val = -2/length(y)
}

###Initialize theta
theta_0 = 0
theta_1 = 10  ###-1, 4.67, 10
k = 1
#theta_arr = lik_arr = c()
###Newton Raphson
while (abs(theta_1-theta_0)>0.01){
  #  theta_arr[k] = theta_1
  #  lik_arr[k] = Lik_func(theta_1,y)
  theta_0 = theta_1
  k = k+1
  theta_1 = theta_0 + Sco_func(theta_0,y)/Ob_fish_func(theta_0,y)
}
print(paste0("k = ",k," beta = ",theta_1))
###theta_0 = -1: it doesn't converge.
###theta_0 = 4.67: it converges to 5.047 at the 4th iteration.
###theta_0 = 10: it doesn't converge.

###Scoring Method
theta_0 = 0
theta_1 = 10 ###-1, 4.67, 10
k = 1
#theta_arr = lik_arr = c()
###Scoring Method
while (abs(theta_1-theta_0)>0.01){
  #  theta_arr[k] = theta_1
  #  lik_arr[k] = Lik_func(theta_1,y)
  print(theta_1)
  theta_0 = theta_1
  k = k+1
  theta_1 = theta_0 + Sco_func(theta_0,y)/Ex_fish_func(theta_0,y)
}
print(paste0("k = ",k," beta = ",theta_1))
###theta_0 = -1: it converges to 0.358 at the 27th iteration
###theta_0 = 4.67: it doesn't converge.
###theta_0 = 11: it converges to 8.548 at the 11th iteration.

###Theta VS Likelihood
theta_arr = seq(from = -1, to = 10, by = 0.01)
lik_arr = c()
for (i in 1:length(theta_arr)){
  lik_arr[i] = Lik_func(theta_arr[i],y)
}
plot(theta_arr, lik_arr, type = "l", main = "Theta VS Likelihood", xlab = "Theta", ylab = "Likelihood")

###Comments: The starting point selection is very important for both algorithm. The bad staring point cause the local maximization even convergence failure.
