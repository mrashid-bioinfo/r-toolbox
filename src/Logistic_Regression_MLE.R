## MLE, PS 206 Class 1

## Linear Regression Example

regressdata <- read.table("ps206data1a.txt", header=T, sep="\t")
regressdata <- na.omit(regressdata)
attach(regressdata)
head(regressdata)

startv <- c(0,0,1)
y <- meanhap
X <- cbind(1, lgdp)

ols.lf <- function(beta) {
  e <- (y - X%*%beta[1:2])/beta[3]
  logl <- sum(log(dnorm(e)/beta[3]))
  return(logl)
} 

install.packages("maxLik")
library(maxLik)

regress.model <- maxNR(ols.lf, grad=NULL, hess=NULL, startv, print.level=4)
summary(regress.model)

summary(lm(meanhap ~ lgdp))


## Logit Example


logitdata <- read.table("ps206data1b.txt", header=T, sep="\t")

logitdata <- as.data.frame(logitdata)
attach(logitdata)

# The built-in routine to estimate logits and probits (and other models) in R -- estimation method is IRLS

logitmodel <- glm(itn ~  altitude, family=binomial(link="logit"), data=logitdata)
summary(logitmodel)

## Example of analytic solution -- constant only model

logitmodel2 <- glm(itn ~  1, family=binomial(link="logit"), data=logitdata)
summary(logitmodel2)

log(mean(itn)) - log(1 - mean(itn))


## writing our own likelihood function


X <- cbind(1, altitude)
y <- itn

## Counts number of columns in X so we know how many independent variables we have

K <- as.numeric(ncol(X))

## here's a good trick: use OLS for better start values ##

startvalue.model <- lm(itn ~ altitude)
startv <- startvalue.model$coefficients

## These functions are what we maximize to estimate our logit ##

## logit log-likelihood function ##

logit.lf <- function(beta) {

exb <- exp(X%*%beta) 
prob1 <- exb/(1+exb) 

logexb <- log(prob1)

y0 <- 1 - y
logexb0 <- log(1 - prob1)

yt <- t(y)
y0t <- t(y0)

logl <- -sum(yt%*%logexb + y0t%*%logexb0)

return(logl)

}


## logit gradient function ##

logit.gr <- function(beta) {

grad <- beta*0


exb <- exp(X%*%beta) 
prob1 <- exb/(1+exb) 

for (k in 1:K) { 
  grad[k] <- sum(X[,k] * (y - prob1))
  }

return(-grad)

}

logitmodel <- optim(startv, logit.lf, gr=logit.gr, method="BFGS", control=list(trace=TRUE, REPORT=1), hessian=TRUE)
coeffs <- logitmodel$par
covmat <- solve(logitmodel$hessian)
stderr <- sqrt(diag(covmat))
zscore <- coeffs/stderr
pvalue <- 2*(1 - pnorm(abs(zscore)))
results <- cbind(coeffs,stderr,zscore,pvalue)
colnames(results) <- c("Coeff.", "Std. Err.", "z", "p value")
print(results)  


## Graph likelihood for age coefficient, holding all else at maximum

altcoeff <- seq(-2,0,0.01)
coeff2 <- coeffs
LL <- NULL

# LL
for (m in altcoeff) {
  coeff2[2] <- m
  LLt <- sum(y*log(exp(X%*%coeff2)/(1 + exp(X%*%coeff2))) + (1-y)*log(1/(1 + exp(X%*%coeff2))))
  LL <- rbind(LL,LLt)
}

plot(altcoeff,LL, type="l", xlab="Alt. Coeff.", ylab="LL")

abline(v=coeffs[2], col="red")
abline(h=-(logitmodel$value), col="red")



## Trying different optimization routines

logit.lf <- function(beta) y*log(exp(X%*%beta)/(1 + exp(X%*%beta))) + (1-y)*log(1/(1 + exp(X%*%beta)))


## BFGS estimation

logitmodel.1 <- maxBFGS(logit.lf, start=startv2, print.level=2)
summary(logitmodel.1)

## NR estimation

logitmodel.2 <- maxNR(logit.lf, start=startv2, print.level=2)
summary(logitmodel.2)

## BHHH estimation

logitmodel.3 <- maxBHHH(logit.lf, start=startv2, print.level=2)
summary(logitmodel.3)

## Compare the results

coeffs <- logitmodel.1$estimate
covmat <- solve(-(logitmodel.1$hessian))   ### note we're using the negative of the hessian here (not working with negative as in optim)
stderr <- sqrt(diag(covmat))
zscore <- coeffs/stderr
pvalue <- 2*(1 - pnorm(abs(zscore)))
results.1 <- cbind(coeffs,stderr,zscore,pvalue)
colnames(results.1) <- c("Coeff.", "Std. Err.", "z", "p value")

coeffs <- logitmodel.2$estimate
covmat <- solve(-(logitmodel.2$hessian))   ### note we're using the negative of the hessian here (not working with negative as in optim)
stderr <- sqrt(diag(covmat))
zscore <- coeffs/stderr
pvalue <- 2*(1 - pnorm(abs(zscore)))
results.2 <- cbind(coeffs,stderr,zscore,pvalue)
colnames(results.2) <- c("Coeff.", "Std. Err.", "z", "p value")

coeffs <- logitmodel.3$estimate
covmat <- solve(-(logitmodel.3$hessian))   ### note we're using the negative of the hessian here (not working with negative as in optim)
stderr <- sqrt(diag(covmat))
zscore <- coeffs/stderr
pvalue <- 2*(1 - pnorm(abs(zscore)))
results.3 <- cbind(coeffs,stderr,zscore,pvalue)
colnames(results.3) <- c("Coeff.", "Std. Err.", "z", "p value")

print(results.1) # BFGS #
print(results.2) # NR #
print(results.3) # BHHH #

