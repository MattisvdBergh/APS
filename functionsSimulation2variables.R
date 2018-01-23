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

updatePerformance = function(data, rho) {
  n = nrow(data)
  X1 = data$Predictor
  mean1 = mean(data$Performance)
  sd1 = sd(data$Performance)
  
  C =  matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- rnorm(n, mean = mean1, sd = sd1)
  X <- cbind(X1,X2)
  
  # induce correlation (does not change X1)
  newData <- X %*% C
  
  data[,1:2] = newData
  
  return(newData)
}


runSimulation = function(nPop, nNew, nOrg, testValidity = .5, 
                         years, startYear = 2018, 
                         corUnrelVar = .3, attritionRate, applicantsPerPosition, 
                         changedTestValidity = NULL, changePoint = NULL,
                         varNames = c("Performance", "Predictor")){
  
  
  paste(nPop, nNew, nOrg, testValidity, years, attritionRate, applicantsPerPosition, changedTestValidity, changePoint, sep = ";")
  
  # covariance
  corMat = matrix(c(1, testValidity, 
                    testValidity, 1), nrow = 2)
  
  # Compute the covariance matrix based on the correlation matrix
  covMat = cor2cov(corMat, c(sdPerf, sdPred))
  row.names(corMat) = colnames(corMat) = row.names(covMat) = colnames(covMat) = varNames

  # Compute a new covariance matrix if a change is required  
  if(!is.null(changePoint)){
  corMat2 = matrix(c(1, changedTestValidity,
                     changedTestValidity, 1), nrow = 2)
  covMat2 = cor2cov(corMat2, c(sdPerf, sdPred))
  row.names(corMat2) = colnames(corMat2) = row.names(covMat2) = colnames(covMat2) = varNames
  }  
  
  # Simulate a population
  pop = simulatePopulation(n = 10000, mu = mu, Sigma = covMat, nTest = 2, var_names = varNames)
  
  # Draw employees of an organisation from the population
  orgRows = sample(1:nrow(pop), nOrg)
  org = cbind.data.frame(Year = startYear,
                         HiredIn = startYear, 
                         pop[orgRows,])
  
  # coefficients within the organization
  resLMperfPredOrg = glm(Performance ~ Predictor, data = org)
  
  # Remove selected applicants from population
  pop = pop[-orgRows,]
  
  # Test whether there is a mean difference between the population and the organization
  testData = rbind(cbind(pop, org = 0),
                   cbind(org[,c(varNames, "New")], org = 1))
  resultsPerf = glm(Performance~org, data = testData)
  
  lResLmPerfPerdOrg = lRegResults = lPerfResults = lCohdResults = lSM = vector("list", years)
  
  
  ### collect info on pop and org
  spreadMeasuresPop = rbind(apply(pop[,varNames], 2, mean),
                            apply(pop[,varNames], 2, sd))
  spreadMeasuresOrg = rbind(apply(org[,varNames], 2, mean),
                            apply(org[,varNames], 2, sd))
  rownames(spreadMeasuresPop) = rownames(spreadMeasuresOrg) = c("Mean", "SD")
  
  lSM[[1]] = cbind(spreadMeasuresPop, spreadMeasuresOrg)
  lPerfResults[[1]] = resultsPerf
  lCohdResults[[1]] = cohen.d(Performance~as.factor(org), data = testData)$estimate
  lResLmPerfPerdOrg[[1]] = resLMperfPredOrg
  
  for(idxYear in seq(years)){
    
    if(!is.null(changePoint)){
      if(idxYear == changePoint){
        covMat = covMat2
        corMat = corMat2
      }
    }
    
    # new employees are now old employees
    org$New = FALSE
    
    # change year
    currentYear = idxYear + startYear
    
    # Random employees leave the organisation
    leaversRows = sample(1:nrow(org), nOrg*attritionRate)
    org = org[-leaversRows,]
    
    # Update the population
    pop = updatePopulation(pop = pop, nNew = nNew, nRemove = nRemove, 
                           nTest = 2, mu = mu, Sigma = covMat)
    
    # Compute new correlated performance for the population
    predNewPop = updatePerformance(dat = pop, rho = corMat[1:2,1:2])
    
    # Compute new correlated performance for the organization
    predNewOrg = updatePerformance(org, rho = corMat[1:2,1:2])
    
    pop$Performance = predNewPop[,2]
    org$Performance = predNewOrg[,2]
    
    # retrieve applicant pool from population
    applicantRows = sample(1:nrow(pop), nOrg * attritionRate * applicantsPerPosition)
    applicants = pop[applicantRows,]
    
    # coefficients within the organization
    resLMperfPredOrg = glm(Performance ~ Predictor, data = org)
    
    predictedPerformance = 
      colSums(coef(resLMperfPredOrg) *
                rbind(1, applicants$Predictor))
    
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
    
    # Test whether there is a difference between the pop and org
    # and collect the results
    testData = rbind(cbind(pop, org = 0),
                     cbind(org[,3:5], org = 1))
    
    
    
    
    spreadMeasuresPop = rbind(apply(pop[,varNames], 2, mean),
                              apply(pop[,varNames], 2, sd))
    spreadMeasuresOrg = rbind(apply(org[,varNames], 2, mean),
                              apply(org[,varNames], 2, sd))
    rownames(spreadMeasuresPop) = rownames(spreadMeasuresOrg) = c("Mean", "SD")
    
    lSM[[idxYear]] = cbind(spreadMeasuresPop, spreadMeasuresOrg)
    lPerfResults[[idxYear]] = glm(Performance~org, data = testData)
    lCohdResults[[idxYear]] = cohen.d(Performance~as.factor(org), 
                                                       data = testData)$estimate
    lResLmPerfPerdOrg[[idxYear]] = resLMperfPredOrg
    
    print(idxYear)
  }
  
  return(list(spreadMeasuresPopOrg = lSM,
              regressionPerformance = lPerfResults,
              CohenD = lCohdResults,
              PerfOrg = lResLmPerfPerdOrg))
}
