library(tidyverse)

caloriesdata = read.csv("calories.txt", header = T)

n = nrow(caloriesdata)
y = caloriesdata$Calories
x = caloriesdata$Time
numerator = sum(y*x) - n*mean(x)*mean(y)
denominator = sum(x^2) - n*mean(x)^2
beta1hat = numerator/denominator

beta0hat = mean(y) - beta1hat*mean(x)

sigma2RSE = sum((y-(beta0hat+beta1hat*x))^2)/(n-2)

varhatbeta1hat = sigma2RSE/(sum(x^2)-n*mean(x)^2)

n = nrow(caloriesdata)

intercept = rep(1, n)
Xmatrix = cbind(intercept, Time = caloriesdata$Time)
Yvector = caloriesdata$Calories

XtX=t(Xmatrix)%*%Xmatrix

invXtX = solve(XtX)

XtY = t(Xmatrix)%*%Yvector

mles = invXtX %*% XtY

sigma2RSE = t(Yvector - Xmatrix %*% mles) %*% (Yvector - Xmatrix %*% mles) / (n-2)

se_beta_time <- sqrt(.7221977)
se_beta_intercept <- sqrt(862.725)

t_mult <- qt(.975, df = n - 2)

hatmatrix = Xmatrix %*% invXtX %*% t(Xmatrix)

hii = diag(hatmatrix)

xnew = c(1, 20)
yhat = t(xnew) %*% mles

#prediction interval for the AVERAGE calorie amount at time = 20
lowerlimitavg = yhat - qt(.975, df = n - 2)*sqrt(sigma2RSE*t(xnew) %*% invXtX %*% xnew)
upperlimitavg = yhat + qt(.975, df = n - 2)*sqrt(sigma2RSE*t(xnew) %*% invXtX %*% xnew)

#prediction interval for an INDIVIDUAL'S calorie amount at time = 20
lowerlimitind = yhat - qt(.975, df = n - 2)*sqrt(sigma2RSE*(1+t(xnew) %*% invXtX %*% xnew))
upperlimitind = yhat + qt(.975, df = n - 2)*sqrt(sigma2RSE*(1+t(xnew) %*% invXtX %*% xnew))

yhat
c(lowerlimitavg, upperlimitavg)
c(lowerlimitind, upperlimitind)
