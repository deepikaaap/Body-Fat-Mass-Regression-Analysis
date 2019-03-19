bodyfatwomen <- read.csv(file="bodyfatwomen.csv")
View(bodyfatwomen)
plot(bodyfatwomen)

n <- nrow(bodyfatwomen)
train <- sample(1:n, 4 * n / 5) # 2/3 for training
test <- (-train)
bodyfatwomen.test <- bodyfatwomen[test, ]
bodyfatwomen.train <- bodyfatwomen[train, ]

# Finding best model
model <- lm((DEXfat) ~ age + waistcirc + (hipcirc) + (elbowbreadth) + (kneebreadth) + anthro3a + anthro3b + anthro3c + anthro4, data = bodyfatwomen)
summary(model1)
par(mfrow = c(2,2))
plot(model1)
# Multiple R-squared:  0.9454,	Adjusted R-squared:  0.9347  
mse_train_1 <- mean(((bodyfatwomen.train$DEXfat) - predict.lm(model1, bodyfatwomen)) ^ 2)
# 0.008572449
#mse_test_1 <- mean((log(bodyfatwomen.test$DEXfat) - predict.lm(model1, bodyfatwomen)) ^ 2)
# 0.00949866
mse_cv_1 <- cross_vaildation_error(4,model1,bodyfatwomen)


model1 <- lm(log(DEXfat) ~ age + waistcirc + (hipcirc) + (elbowbreadth) + (kneebreadth) + anthro3c + anthro4, data = bodyfatwomen)
summary(model1)
# Multiple R-squared:  0.9454,	Adjusted R-squared:  0.9347  
mse_train_1 <- mean((log(bodyfatwomen.train$DEXfat) - predict.lm(model1, bodyfatwomen)) ^ 2)
# 0.008572449
#mse_test_1 <- mean((log(bodyfatwomen.test$DEXfat) - predict.lm(model1, bodyfatwomen)) ^ 2)
# 0.00949866
mse_cv_1 <- cross_vaildation_error(4,model1,bodyfatwomen)
# 0.01070745
AIC(model1)
# -85.59419
BIC(model1)
# -63.31532

############### Removing  elbowbreadth with highest p-value ###########################

model2 <- lm(log(DEXfat) ~ age + waistcirc + (hipcirc) + (kneebreadth) + anthro3a + anthro3b + anthro3c + anthro4, data = bodyfatwomen )
summary(model2)
# Multiple R-squared:  0.9447,	Adjusted R-squared:  0.9353
mse_train_2 <- mean((log(bodyfatwomen.train$DEXfat) - predict.lm(model2, bodyfatwomen.train)) ^ 2)
# 0.008680507
mse_test_2 <- mean((log(bodyfatwomen.test$DEXfat) - predict.lm(model2, bodyfatwomen.test)) ^ 2)
#  0.00910053
mse_cv_2 <- cross_vaildation_error(4,model2,bodyfatwomen)
# 0.01055181
AIC(model2)
# -86.8927
BIC(model2)
# -66.63918

############### Removing  anthro3b with highest p-value ###########################

model3 <- lm(log(DEXfat) ~ age + waistcirc + (hipcirc) + (kneebreadth) + anthro3a + anthro3c + anthro4, data = bodyfatwomen )
summary(model3)
# Multiple R-squared:  0.9447,	Adjusted R-squared:  0.9353
mse_train_3 <- mean((log(bodyfatwomen.train$DEXfat) - predict.lm(model3, bodyfatwomen.train)) ^ 2)
# 0.008752149
mse_test_3 <- mean((log(bodyfatwomen.test$DEXfat) - predict.lm(model3, bodyfatwomen.test)) ^ 2)
# 0.009549161
mse_cv_3 <- cross_vaildation_error(4,model3,bodyfatwomen)
# 0.01050037
AIC(model3)
# -88.43242
BIC(model3)
# -70.20426

################ Removing Kneebreadth ############################################

model4 <- lm(log(DEXfat) ~ age + waistcirc + (hipcirc) + anthro3a + anthro3c + anthro4, data = bodyfatwomen )
summary(model4)
# Multiple R-squared:  0.9447,	Adjusted R-squared:  0.9353
mse_train_4 <- mean((log(bodyfatwomen.train$DEXfat) - predict.lm(model4, bodyfatwomen.train)) ^ 2)
# 0.008895439
mse_test_4 <- mean((log(bodyfatwomen.test$DEXfat) - predict.lm(model4, bodyfatwomen.test)) ^ 2)
# 0.01153477
mse_cv_4 <- cross_vaildation_error(4,model4,bodyfatwomen)
# 0.0106356
AIC(model4)
# -89.52302
BIC(model4)
# -73.3202

###########################################################################################################
#Removing anthro3b due to high correlation 
model5<- lm(log(DEXfat) ~ age + waistcirc + (hipcirc) + (elbowbreadth) + (kneebreadth) + anthro3a + anthro3b + anthro3c + anthro4, data = bodyfatwomen)
plot(model5)
ols_mallows_cp(model5, model1)


# Function for cross vaildation error
# AIC and Mallows Cp calculation
library(leaps)
cross_vaildation_subset<-function(folds,dataset){
  n_folds <- 4
  folds_i <- sample(rep(1:n_folds, length.out = 50))
  
  cv.errors <- matrix(NA, nrow = n_folds, ncol = 7)
  for (k in 1:n_folds) {
    best.fit=regsubsets(log(DEXfat)~age + waistcirc + hipcirc + elbowbreadth + kneebreadth + anthro3c + anthro4,data=bodyfatwomen[folds_i!=3,],  method = "backward")
    rs<-summary(best.fit)
    rs$which
    # Printing Cp
    print(rs$cp)
    n<-nrow(bodyfatwomen)
    p<-7
    print(n)
    print(p)
    ( AIC<-n*log(rs$rss/n)+(2*p) )
    print(rs$bic)
    print(rs$rsq)
    print(rs$adjr2)
    
    (BIC <- n*log(rs$rss/n)+(log(n)*p))
    print(AIC)
    par(mfrow = c(1,1))
    plot(AIC~I(1:(p)),ylab="AIC",xlab="Number of Predictors")
    par(mfrow = c(1,1))
    plot(1:(p),rs$cp,xlab="Number of predictors",ylab="Cp Mallows")
    abline( 0, 1)
    
    test.mat = model.matrix(log(DEXfat)~age + waistcirc + hipcirc + elbowbreadth + kneebreadth + anthro3c + anthro4, data = bodyfatwomen[folds_i==k,]) 
    print(bodyfatwomen.train[folds_i==k,]$DEXfat)
    for(i in 1:7){
      coefi= coef(best.fit, id = i)
      print(names(coefi))
      pred=test.mat[,names(coefi)]%*%coefi
      print(pred)
      # Predict on the hold out part of the fold for that subset
      pred=predict(best.fit, bodyfatwomen.train[folds_i==k,],id=i)
      # Get the mean squared error for the model trained on the fold with the subset
      cv.errors[k,i]=mean((log(bodyfatwomen[folds_i==k,]$DEXfat)-pred)^2)
    }
  }
  cv <- colMeans(cv.errors)
  
  return(cv)
}


best_model <- lm(log(DEXfat) ~ waistcirc + (hipcirc) + anthro3c + anthro4, data = bodyfatwomen)
library(devtools)
library(gistr)
summary(best_model)
lev = hatvalues(best_model)
highlev <- print(lev[lev >  0.1126761]) #To find out which points have a high hat diagonal
library(car)
outlierTest(best_model)
print(highlev[highlev > 0.1126761])
influence <- influence.measures(best_model)
covratios = covratio(best_model)
print(covratios)
hats = hatvalues(best_model)
cook = cooks.distance(best_model)
print(cook)
summary(cook)

ols_plot_cooksd_bar(best_model)
ols_plot_cooksd_chart(best_model)



# Detecting outliers
p = 4
n = 71
leverage.cutoff <- 2*p/n # Montgomery p. 213
cooks.cutoff <- qf(0.5, p, n - p, lower.tail = FALSE) # Montgomery p. 215
studres.cutoff <- qt(0.05/2, n - p, lower.tail = FALSE) # Montgomery p. 135
leverage.cutoff
cooks.cutoff
studres.cutoff
influence <- influence.measures(best_model)
bodyfatwomen[27,] <- NA


# Bootstrapping and finding confidence intervals on coefficients
ms <- c(500)
par(mfrow = c(2, 3))
for (m in ms) {
  coefs <- c()
  for (i in seq(m)) {
    n <- nrow(bodyfatwomen)
    indices <- sample(n, n, replace = TRUE)
    residual.boot <- sum$residuals[indices]
    y.boot <- exp(best_model$fitted.values) + residual.boot # New bootstrap samples, see (eq. 5.19)
    best_model.boot <- lm(log(y.boot) ~ waistcirc + (hipcirc) + anthro3c + anthro4,data=bodyfatwomen)
    coefs <- rbind(coefs, coef(best_model.boot))
  }
  
  param.sd.boot <- apply(coefs, 2, sd)
  print(param.sd.boot)
  # Also create confidence intervals accoring to the percentile method
  # presented in Section 15.4.2, Montgomery
  conf.ints <- c()
  for (k in seq( length(best_model$coefficients))) {
    hist(coefs[, k],main=paste('Histogram of ',names( coef(best_model))[k]))
    quants <- quantile(coefs[, k], probs = c(0.025, 0.975))
    beta.est <- coef(best_model)[k]
    D1 <- beta.est - quants[1]
    D2 <- quants[2] - beta.est
    conf.ints <- rbind(conf.ints, c(beta.est - D2, beta.est + D1, beta.est))
  }
  
  colnames(conf.ints) <- c( names(quants), "beta est")
  rownames(conf.ints) <- names( coef(best_model))
  conf.ints
  confint(best_model)
}




