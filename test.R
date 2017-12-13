rm(list = ls())

library(MASS)
library(MBESS)
library(tibble)
library(effsize)

simulatePopulation = function(n, mu, Sigma, nTest, emperical = TRUE, var_names = NULL, New = TRUE){
  # simulate data matrix using covariance matrix
  d <- mvrnorm(n = n, mu = rep(mu, nTest), Sigma = Sigma, empirical = emperical); 
  d <- as.tibble(d); # transform to dataframe
  if(!is.null(var_names)) names(d) <- var_names; # add names if provided
  d$New <- New; # add column to identify new employees
  return(d);
}

updatePopulation = function(pop, nNew, nRemove, mu, Sigma, nTest, emperical = FALSE, New = TRUE){
  
  # remove random set from the population
  toRemove <- sample(1:nrow(pop), nRemove)
  newPopTemp = pop[-toRemove,]
  
  # simulate data matrix using covariance matrix
  d <- mvrnorm(n = nNew, mu = rep(mu, nTest), Sigma = Sigma, empirical = emperical); 
  d <- as.tibble(d); # transform to dataframe
  d$New <- New; # add column to identify new employees
  newPop <- rbind(newPopTemp, d)
  return(newPop);
}

# Got this from:
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
# returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x
getBiCop <- function(dat, rho, nTest) {
  n = nrow(dat)
  X1 = dat$Predictor
  X3 = dat$HC
  mean1 = mean(dat$Performance)
  sd1 = sd(dat$Performance)
  
  C = rho
  
  
  C <- chol(C)
  
  X2 <- rnorm(n, mean = mean1, sd = sd1)
  X <- cbind(X1,X2,X3)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  # all.equal(X1,X[,1])
  # cor(X)
  cor(df)
  
  return(df)
}


getTriCop = function(dat, desiredCorrelations){
  n = nrow(dat)
  k = ncol(desiredCorrelations)
  mean1 = mean(dat$Performance)
  sd1 = sd(dat$Performance)
  x = matrix( c(rnorm(n, mean1, sd),
                rnorm(n),
                rnorm(n)), nc=k)
  x[,1] = dat$Predictor
  x[,2] = dat$HC
  y <- x %*% solve(chol(var(x))) %*% chol(desiredCorrelations)
  # cor(y)      # Desired correlation matrix
  return(y)
}

apply(dat, 2, summary)
apply(y, 2, summary)





simCorData = function(pop, rho, x1, mu, sd){

n <- nrow(pop)                    # length of vector
theta <- acos(rho)             # corresponding angle
x2    <- rnorm(n, mu, sd)      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
if(cor(x1, x) == rho){
  return(x)
} else {
  print("Something went wrong with simulating the correlated data")
}                                    # check correlation = rho
}

















# function parameters
nPop = 10000
nNew = nRemove = 100 
nOrg = 1000
mu = 0
varNames = c("Performance", "Predictor", "HC")
testValidity = .4
sdPerf = sdPred = sdExtra = 1
startYear = 2017
years = 50
attritionRate = 0.10
applicantsPerPosition = 10



####
corMat = matrix(c(1, testValidity, 0, 
                  testValidity, 1, .3,
                  0, .3, 1), nrow = 3)


# Make a correlation matrix
# corMat = matrix(nrow = length(varNames), ncol = length(varNames))
# corMat[1, -1] = corMat[-1, 1] = testValidity
# diag(corMat) = 1

# Compute the covariance matrix based on the correlation matrix
covMat = cor2cov(corMat, c(sdPerf, sdPred, sdExtra))
row.names(corMat) = colnames(corMat) = row.names(covMat) = colnames(covMat) = varNames

# Simulate a population
pop = simulatePopulation(n = 10000, mu = mu, Sigma = covMat, nTest = 3, var_names = varNames)

# Draw employees of an organisation from the population
orgRows = sample(1:nrow(pop), nOrg)
org = cbind.data.frame(Year = startYear,
                        HiredIn = startYear, 
                        pop[orgRows,])

# Remove selected applicants from population
pop = pop[-orgRows,]

# Test whether there is a mean difference between the population and the organization
test = rbind(cbind(pop, org = 0),
             cbind(org[,c(varNames, "New")], org = 1))
resultsPerf = glm(Performance~org, data = test)
resultsHC = glm(HC~org, data = test)

lPerfResults = lCohdResults = lHCResults = lSM = list()

spreadMeasuresPop = rbind(apply(pop[,varNames], 2, mean),
                       apply(pop[,varNames], 2, sd))
spreadMeasuresOrg = rbind(apply(org[,varNames], 2, mean),
                          apply(org[,varNames], 2, sd))
rownames(spreadMeasuresPop) = rownames(spreadMeasuresOrg) = c("Mean", "SD")

lSM[[1]] = cbind(spreadMeasuresPop, spreadMeasuresOrg)

lHCResults[[1]] = resultsHC
lPerfResults[[1]] = resultsPerf
lCohdResults[[1]] = cohen.d(Performance~as.factor(org), data = test)$estimate



for(idxYear in seq(years)){

  # new employees are now old employees
  org$New = FALSE
  
  # change year
  currentYear = idxYear + startYear
  
  # Random employees leave the organisation
  leaversRows = sample(1:nrow(org), nOrg*attritionRate)
  org = org[-leaversRows,]
  
  # Update the population
  pop = updatePopulation(pop = pop, nNew = nNew, nRemove = nRemove, 
                         nTest = 3, mu = mu, Sigma = covMat)
  
  # Compute new correlated performance for the population
  predNewPop = getTriCop(dat = pop, desiredCorrelations = corMat)
  
  # Compute new correlated performance for the organization
  predNewOrg = getTriCop(org, desiredCorrelations = corMat)

  pop$Performance = predNewPop[,2]
  org$Performance = predNewOrg[,2]
  
  # retrieve applicant pool from population
  applicantRows = sample(1:nrow(pop), nOrg * attritionRate * applicantsPerPosition)
  applicants = pop[applicantRows,]
  
  # coefficients within the organization
  predictedPerformance = coef(glm(Performance~Predictor, data = org))[1] +
    coef(glm(Performance~Predictor, data = org))[2] * applicants$Predictor
  
  # determine top estimated performers
  topApplicants = order(predictedPerformance, decreasing = T)[1 : (nOrg * attritionRate)]
  
  # add top estimated performers to organization
  org = rbind.data.frame(org, 
                         cbind.data.frame(Year = NA,
                                          HiredIn = currentYear, 
                                          applicants[topApplicants,]))
  org$Year = currentYear
  
  # remove selected applicants from the population
  pop = pop[-applicantRows[topApplicants], ] # remove selected applicants from population
  
  test = rbind(cbind(pop, org = 0),
               cbind(org[,3:5], org = 1))
  resultsTemp = glm(Performance~org, data = test)
  lRegResults[[length(lRegResults) + 1]] = resultsTemp
  lCohdResults[[length(lCohdResults) + 1]] = cohen.d(Performance~as.factor(org), 
                                                     data = test)$estimate
  
  
  print(idxYear)
    
  
}

plot(unlist(lCohdResults))


