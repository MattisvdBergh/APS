rm(list = ls())

library(MASS)
library(MBESS)
library(tibble)

n = 1000
mu = c(100, 100)
varNames = c("Performance", "Predictor")
testValidity = .4
sdPerf = sdPred = 15

corMat = matrix(nrow = length(varNames), ncol = length(varNames))
corMat[1, -1] = corMat[-1, 1] = testValidity
diag(corMat) = 1

covMat = cor2cov(corMat, c(sdPerf, sdPred))
row.names(corMat) = colnames(corMat) = row.names(covMat) = colnames(covMat) = varNames

simulatePopulation = function(n, mu, Sigma, emperical = TRUE, var_names = NULL, New = TRUE){
  # simulate data matrix using covariance matrix
  d <- mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = emperical); 
  d <- as.tibble(d); # transform to dataframe
  if(!is.null(var_names)) names(d) <- var_names; # add names if provided
  d$New <- New; # add column to identify new employees
  return(d);
}

pop = simulatePopulation(n = n, mu = mu, Sigma = covMat, var_names = varNames) 

corTest = cor(pop$Performance, pop$Predictor)

# Got this from:
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
# returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x
getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {
  if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- mar.fun(n)
  X <- cbind(X1,X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  #all.equal(X1,X[,1])
  #cor(X)
  
  return(df)
}


predNew = getBiCop(1000, 0.4, x = pop$Predictor - mean(pop$Predictor))

cor(pop$Performance, predNew[,2] +  100)

