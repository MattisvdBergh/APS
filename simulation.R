rm(list = ls())

# directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source functions
source("functionsSimulation.R")


# function parameters
nPop = 10000
nNew = nRemove = 100 
nOrg = 1000
mu = 0
varNames = c("Performance", "Predictor")
testValidity = .4
changeValidity = 0
sdPerf = sdPred = 1
startYear = 2017
years = 50
attritionRate = 0.10
applicantsPerPosition = 10
changePoint = years/2

test = runSimulation(nPop = nPop, nNew = nNew, nOrg = nOrg,
                     testValidity = 0.5, startYear = startYear, 
                     years = years, attritionRate = attritionRate, 
                     applicantsPerPosition = applicantsPerPosition)
                     changedTestValidity = 0.1,
                     changePoint = 10)


runManySimulations = function(nRep, setSeeds){
  for(idxnRep in 1:nRep){
  set.seed(setSeeds[idxnRep])
  
    runSimulation(nPop = nPop, nNew = nNew, nOrg = nOrg,
                  testValidity = 0.5, startYear = startYear, 
                  years = years, attritionRate = attritionRate, 
                  applicantsPerPosition = applicantsPerPosition,
                  changedTestValidity = 0.1, changePoint = 10)
  
  
}
  
}




plot(unlist(sapply(sapply(test$PerfOrg, coef), function(x){x[[2]]})))

plot(unlist(test$CohenD))


param = expand.grid(
  nPop = c(250000, 500000),
  nOrg = c(1000),
  nNew = c(50),
  attritionRate = c(0.10),
  applicantsPerPosition = c(10),
  testValidity = c(0.5),
  corUnrelVar = c(0.3),
  changedTestValidity = c(0),
  changePoint = c(0, 25),
  years = c(50)
)

nrow(param)

test = mapply(runSimulation, 
       nPop = param$nPop, 
       nNew = param$nNew, 
       nOrg = param$nOrg,
       testValidity = param$testValidity,
       corUnrelVar = param$corUnrelVar,
       attritionRate = param$attritionRate,
       applicantsPerPosition = param$applicantsPerPosition,
       changedTestValidity = param$changedTestValidity,
       changePoint = param$changePoint,
       year = param$years)
       
