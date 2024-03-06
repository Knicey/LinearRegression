## script from regression review lecture

#a useful package.  type "install.packages("tidyverse") if you have not installed it yet.
library(tidyverse)

#read in the calories data-- data on calories consumed at lunch by 20 toddlers
#goal is to see if there is a relationship between time at table and calories consumed.

caloriesdata = read.csv("calories.txt", header = T)

#this command allows you to see the first 5 rows of the dataset
head(caloriesdata)

#this command provides a quick summary of the distributions of the variables
summary(caloriesdata)

#plot the data (commands without tidyverse)
plot(x = caloriesdata$Time, y = caloriesdata$Calories, xlab = "Time", ylab = "Calories", main = "Calories versus time")

#same plot using ggplot (note: ggplot is more useful for multivariate datasets)

ggplot(data = caloriesdata) + geom_point(mapping = aes(x = Time, y = Calories))

#get correlation
cor(caloriesdata$Time, caloriesdata$Calories)

#fit a linear regression model
regConT = lm(Calories ~ Time, data = caloriesdata, x=T)

#view the results of the regression
summary(regConT)

#get 95% confidence intervals for the coefficients
confint(regConT)

### checking the fit of the model with residual plots

#residuals versus predictors (old plotting function)
plot(y = regConT$residual, x=caloriesdata$Time, ylab = "Residuals", xlab = "Time", main = "Plot of residuals versus Time")
abline(0,0)  #adds a horizontal line

#note: you also can make residual plots using ggplot.
#I find the plot command a little easier.

#normal quantile plot or histograms of residuals can help assess normal distribution assumption
qqnorm(regConT$residuals, ylab = "Residuals", main = "Normal quantile plot of residuals")

###  confirming that the MLE formula matches output of lm command

#beta_1
n = nrow(caloriesdata)
y = caloriesdata$Calories
x = caloriesdata$Time
numerator = sum(y*x) - n*mean(x)*mean(y)
denominator = sum(x^2) - n*mean(x)^2
beta1hat = numerator/denominator
beta1hat

#beta_0
beta0hat = mean(y) - beta1hat*mean(x)
beta0hat

#sigma^2 (RSE^2 estimator)
sigma2RSE = sum((y-(beta0hat+beta1hat*x))^2)/(n-2)
sigma2RSE  ## this is estimate of sigma^2.  R reports estimate of sigma
sqrt(sigma2RSE)

#var(beta1hat)
varhatbeta1hat = sigma2RSE/(sum(x^2)-n*mean(x)^2)
varhatbeta1hat  ## this is estimate of the variance.  R reports estimate of the standard error.
sqrt(varhatbeta1hat)


### let's show that the results of the matrix algebra match the MLEs
## Class 6
#################
n = nrow(caloriesdata)

# make the n by 2 matrix of predictors X, with a column of ones and a column of times.
intercept = rep(1, n)
Xmatrix = cbind(intercept, Time = caloriesdata$Time)
Xmatrix
Yvector = caloriesdata$Calories

##compute X^tX
#use t(X) as the transpose and %*% for matrix multiplication
XtX=t(Xmatrix)%*%Xmatrix
XtX

##take the inverse of XtX
invXtX = solve(XtX)
invXtX

#compute X^tY
XtY = t(Xmatrix)%*%Yvector

#compute MLE
mles = invXtX %*% XtY
mles

##regression variance
sigma2RSE = t(Yvector - Xmatrix %*% mles) %*% (Yvector - Xmatrix %*% mles) / (n-2)
sigma2RSE
sqrt(sigma2RSE)

##regression covariance matrix
as.numeric(sigma2RSE)*invXtX

#the elements on the diagonal of this matrix are the variances of the mles
#take their square roots to get the standard errors
se_beta_time <- sqrt(.7221977)
se_beta_intercept <- sqrt(862.725)

#the quantity on the off-diagonal is the covariance between beta_0.hat and beta_1.hat
#we don't use the covariance often, but we will see an example in the future.

t_mult <- qt(.975, df = n - 2)

mles[1] + c(-1,1)*t_mult*se_beta_intercept
mles[2] + c(-1,1)*t_mult*se_beta_time

confint(regConT)

########################
# Class 7
#######################

## Checks for Non-linearity:

#residuals versus fitted values (old plotting function)
plot(x = regConT$fitted.values, y=regConT$residual, xlab = "Fitted Values",
     ylab = "Residuals", main = "Residuals vs Fitted Values")
#Note that this line will always be horizontal.
abline(reg=lm(regConT$residual ~ regConT$fitted.values))

#### leverage values
hatmatrix = Xmatrix %*% invXtX %*% t(Xmatrix)

#pull off the diagonal elements
hii = diag(hatmatrix)

plot(x=caloriesdata$Time, y = hii, xlab = "Time", ylab = "Leverage")

#as you can see, h_ii is large for times that are far from the average time
#and h_ii is small for times that are near the average time.

#incidently, R has a function to get the leverage values
hatvalues(regConT)

## Scale the variances by their residuals

#residuals versus predictors old plot
plot(y = regConT$residual, x=caloriesdata$Time, ylab = "Residuals", xlab = "Time", main = "Plot of residuals versus Time")
abline(0,0)  #adds a horizontal line at zero to help with visualization

# Plot residuals versus predictors after scaling the residuals
# by their variances.
# The scaling makes all the residuals have a unit variance.
plot( y = rstudent(regConT),
       x=caloriesdata$Time,
       ylab = "Scaled Residuals",
       xlab = "Time", main = "Plot of residuals versus Time")
abline(0,0)  #adds a horizontal line at zero to help with visualization

###prediction of calories at new times, both average calories and for individual's calories

#predictions of calorie amount at time = 30
xnew = c(1, 20)
yhat = t(xnew) %*% mles

#prediction interval for the AVERAGE calorie amount at time = 30
lowerlimitavg = yhat - qt(.975, df = n - 2)*sqrt(sigma2RSE*t(xnew) %*% invXtX %*% xnew)
upperlimitavg = yhat + qt(.975, df = n - 2)*sqrt(sigma2RSE*t(xnew) %*% invXtX %*% xnew)

#prediction interval for an INDIVIDUAL'S calorie amount at time = 30
lowerlimitind = yhat - qt(.975, df = n - 2)*sqrt(sigma2RSE*(1+t(xnew) %*% invXtX %*% xnew))
upperlimitind = yhat + qt(.975, df = n - 2)*sqrt(sigma2RSE*(1+t(xnew) %*% invXtX %*% xnew))

yhat

c(lowerlimitavg, upperlimitavg)

c(lowerlimitind, upperlimitind)

##how about using R directly?  let's add a few more times, because, why not!
#first make the dataset of new toddlers and their times
#make sure the predictor in the dataset has the same name as the predictor used in the regression model

newtimes = c(30, 35, 40)
newdata = data.frame(Time = newtimes)

#if you want to make predictions for the AVERAGE calorie amount at these times (not an individual), use

#predict.lm(regConT, newdata, interval = "confidence")

#here are the predictions for INDIVIDUAL toddlers

#predict.lm(regConT, newdata, interval = "prediction")


