# numEx from Tuo 2015
# replacement F
# AC
# 

library(ggplot2)
library(reshape2)
library(splines)

# # set number of total data and generate data
# numSeq = 10^3
# xSeq = seq(0, 2*pi, length.out = numSeq)
# ySeq = exp(xSeq/10)*sin(xSeq) +rnorm(numSeq, mean = 0, sd = sqrt(0.1))
# # ySeq = exp(xSeq/10)*sin(xSeq)
# # for (i in c(1:length(ySeq))){
# #   ySeq[i] = ySeq[i]+rnorm(numSeq, mean = 0, sd = sqrt(0.1*xSeq[i]))
# # }
# data = data.frame(x=xSeq, y = ySeq)

startPercent = 0.05
spSeq = seq(5,10,length.out = 5)
repetition = 1000

testLoad = "1000TestData.csv"
trainLoad = "1000TrainData.csv"
imperfectSave = "im1000TrainASGDEx1.csv"
perfectSave = "p1000TrainASGDEx1.csv"

# testLoad = "1000TestDataVariedError2.csv"
# trainLoad = "1000TrainDataVariedError2.csv"
# imperfectSave = "im1000TrainASGDEx2Ver2.csv"
# perfectSave = "p1000TrainASGDEx2Ver2.csv"


# set number of total data and generate data
setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/data')
data = read.csv(file = trainLoad)
numSeq = dim(data)[1]
xSeq = data$x
ySeq = data$y

#####SGD with Adaptive sampling####

# the following function calculates the sum of error squre from imperfect model
sumPolyLossLocal1SGDnonAdap<- function(theta, xSeq.,ySeq.){ 
  powerHat = exp(xSeq./10)*sin(xSeq.)-sqrt(theta^2-theta+1)*(sin(theta*xSeq.)+ cos(theta*xSeq.))
  errorFarm = (powerHat - ySeq.)^2
  return(errorFarm)
}

# the following function calculates the sum of error squre from perfect model
sumPolyLossLocal2SGDnonAdap<- function(theta, xSeq.,ySeq.){ 
  powerHat = exp(xSeq./10)*sin(xSeq.)-abs(theta+1)*(sin(theta*xSeq.)+ cos(theta*xSeq.))
  errorFarm = (powerHat - ySeq.)^2
  return(errorFarm)
}

# the following function takes data as input and calibrates theta accordingly
directFindLocal1SGDnonAdap <- function(xSeq, ySeq, start){
  # set the limit of sample in SGD 
  sampleLimit = length(xSeq)
  
  Lk = 1*10^-3
  ita = 1.5
  testTheta = 0.9
  
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = 100
  # have a matrix to record the number of sample size at each iteration
  sampleSize = array(0,dim = iter)
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  # set starting sample number at each iteration in SGD
  samNum = ceiling(startPercent*length(xSeq))
  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter = j

    # loop for each parameter
    for (k in c(1:length(start))){
      #sample for the current step
      samInx = sample(c(1:length(xSeq)), samNum, replace = F)
      
      samxSeq = xSeq[samInx]
      samySeq = ySeq[samInx]
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = sumPolyLossLocal1SGDnonAdap(theta1,  xSeq. = samxSeq, ySeq.=samySeq)
      errorAvg1 = mean(errors1)
      errors2 = sumPolyLossLocal1SGDnonAdap(theta2,  xSeq. = samxSeq, ySeq.=samySeq)
      errorAvg2 = mean(errors2)
      errorIndi = (errors1-errors2)/(2*h[k])
      sumError = (errorAvg1-errorAvg2)/(2*h[k])
      
      innerLeft = var(errorIndi*sumError)/samNum
      innerRight = (testTheta^2)*(sumError^4)
      
      # conduct inner product test
      if (innerLeft>innerRight){
        # add additional samples
        newSamNum = min(innerLeft*samNum/innerRight, sampleLimit)
        addSamNum = newSamNum - samNum
        addSamInx = sample(c(1:length(xSeq))[!(c(1:length(xSeq)) %in% samInx)], addSamNum, replace = F)
        # update
        errors1 = c(errors1, sumPolyLossLocal1SGDnonAdap(theta1,  xSeq. = xSeq[addSamInx], ySeq.=ySeq[addSamInx]))
        errorAvg1 = mean(errors1)
        errors2 = c(errors2, sumPolyLossLocal1SGDnonAdap(theta2,  xSeq. = xSeq[addSamInx], ySeq.=ySeq[addSamInx]))
        errorAvg2 = mean(errors2)
        errorIndi = (errors1-errors2)/(2*h[k])
        sumError = (errorAvg1-errorAvg2)/(2*h[k])
        samNum = newSamNum
      }
      
      dire = -sumError
      # conduct backtracking line search
      ak = 1 + (var(errors1-errors2)/(samNum*sumError^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      newErrors = sumPolyLossLocal1SGDnonAdap(currPara[k]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal1SGDnonAdap(currPara[k],  xSeq. = samxSeq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sumError^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal1SGDnonAdap(currPara[k]+dire/currLk,  xSeq. = samxSeq, ySeq.=samySeq)
        newLoss = sum(newErrors)/length(newErrors)
      }
      
      # move to next step
      step = 1/currLk
      h = c(0,0)
      moving = dire*step
      h[k] = moving
      movingIndex = 0
      currk = currPara[k] + moving
      # give theta a limit range
      while (currk<=-100 || currk>=100){
        moving = moving/2
        h[1] = moving
        currk = currPara[k]+h[k]
        movingIndex  =  movingIndex+1
      }
      currPara[k] = currPara[k] + moving
    }
    paras1[j+1,] = currPara
    # record current sample size
    sampleSize[j] = samNum
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]
  newList <- list("thetaChoice" = thetaChoice, "sampleSize" = sampleSize, "paras" = paras1, "iters" = actualIter)
  return (newList)
}

directFindLocal2SGDnonAdap <- function(xSeq, ySeq, start){
  # set the limit of sample in SGD 
  sampleLimit = length(xSeq)
  
  Lk = 1*10^-3
  ita = 1.5
  testTheta = 0.9
  
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = 100
  # have a matrix to record the number of sample size at each iteration
  sampleSize = array(0,dim = iter)
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  # set starting sample number at each iteration in SGD
  samNum = ceiling(startPercent*length(xSeq))
  

  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter = j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #sample for the current step
      samInx = sample(c(1:length(xSeq)), samNum, replace = F)
      samxSeq = xSeq[samInx]
      samySeq = ySeq[samInx]
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = sumPolyLossLocal2SGDnonAdap(theta1,  xSeq. = samxSeq, ySeq.=samySeq)
      errorAvg1 = mean(errors1)
      errors2 = sumPolyLossLocal2SGDnonAdap(theta2,  xSeq. = samxSeq, ySeq.=samySeq)
      errorAvg2 = mean(errors2)
      errorIndi = (errors1-errors2)/(2*h[k])
      sumError = (errorAvg1-errorAvg2)/(2*h[k])
      
      innerLeft = var(errorIndi*sumError)/samNum
      innerRight = (testTheta^2)*(sumError^4)
      
      # conduct inner product test
      if (innerLeft>innerRight){
        # add additional samples
        newSamNum = min(innerLeft*samNum/innerRight, sampleLimit)
        addSamNum = newSamNum - samNum
        addSamInx = sample(c(1:length(xSeq))[!(c(1:length(xSeq)) %in% samInx)], addSamNum, replace = F)
        # update
        errors1 = c(errors1, sumPolyLossLocal2SGDnonAdap(theta1,  xSeq. = xSeq[addSamInx], ySeq.=ySeq[addSamInx]))
        errorAvg1 = mean(errors1)
        errors2 = c(errors2, sumPolyLossLocal2SGDnonAdap(theta2,  xSeq. = xSeq[addSamInx], ySeq.=ySeq[addSamInx]))
        errorAvg2 = mean(errors2)
        errorIndi = (errors1-errors2)/(2*h[k])
        sumError = (errorAvg1-errorAvg2)/(2*h[k])
        samNum = newSamNum
        samInx = c(samInx, addSamInx)
      }
      
      dire = -sumError
      # conduct backtracking line search
      ak = 1 + (var(errors1-errors2)/(samNum*sumError^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      newErrors = sumPolyLossLocal2SGDnonAdap(currPara[k]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal2SGDnonAdap(currPara[k],  xSeq. = samxSeq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sumError^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal2SGDnonAdap(currPara[k]+dire/currLk,  xSeq. = samxSeq, ySeq.=samySeq)
        newLoss = sum(newErrors)/length(newErrors)
      }
      
      # move to next step
      step = 1/currLk
      h = c(0,0)
      moving = dire*step
      h[k] = moving
      movingIndex = 0
      currk = currPara[k] + moving
      # give theta a limit range
      while (currk<=-100 || currk>=100){
        moving = moving/2
        h[1] = moving
        currk = currPara[k]+h[k]
        movingIndex  =  movingIndex+1
      }
      currPara[k] = currk
    }
    paras1[j+1,] = currPara
    # record current sample size
    sampleSize[j] = samNum
    
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]
  newList <- list("thetaChoice" = thetaChoice, "sampleSize" = sampleSize, "paras" = paras1, "iters" = actualIter)
  return (newList)
}

# # imperfect model
# # time the function
# ptm <- proc.time()
# SGDnonkChoiceWSAll = directFindLocal1SGDnonAdap(xSeq, ySeq, 10)
# ntmAdaptive = proc.time() - ptm
# 
# # perfect model
# # time the function
# ptm <- proc.time()
# SGDnonk1ChoiceWSAll = directFindLocal2SGDnonAdap(xSeq, ySeq, 10)
# ntmAdaptive1 = proc.time() - ptm
# 
# # calculate RMSE based on the result above
# theta = SGDnonkChoiceWSAll$thetaChoice
# powerHat = exp(xSeq/10)*sin(xSeq)-sqrt(theta^2-theta+1)*(sin(theta*xSeq)+ cos(theta*xSeq))
# errorFarm = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
# 
# theta = SGDnonk1ChoiceWSAll$thetaChoice
# powerHat1 = exp(xSeq/10)*sin(xSeq)-abs(theta+1)*(sin(theta*xSeq)+ cos(theta*xSeq))
# errorFarm1 = sqrt(sum((powerHat1 - ySeq)^2)/length(powerHat))
# 
# errorFarm
# errorFarm1
# SGDnonkChoiceWSAll$thetaChoice
# SGDnonk1ChoiceWSAll$thetaChoice
# SGDnonkChoiceWSAll$sampleSize[1:30]
# SGDnonk1ChoiceWSAll$sampleSize[1:30]
# ntmAdaptive
# ntmAdaptive1
# 



# #Randomly shuffle the data
# data<-data[sample(nrow(data)),]
# 
# #Create 10 equally size folds
# folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
# 
# spSeq = seq(5,10,length.out = 5)
# imIterM <- pIterM <- imThetaM <- pThetaM <- imRMSEM <- pRMSEM <- matrix(0, nrow = 10, ncol = 5)
# imIter <- pIter <- imTheta <- pTheta <- imRMSE <- pRMSE <- array(0, dim = 10)
# 
# ptm <- proc.time()
# #Perform 10 fold cross validation
# for(i in c(1:10)){
#   #Segement your data by fold using the which() function 
#   testIndexes <- which(folds==i,arr.ind=TRUE)
#   testData <- data[testIndexes, ]
#   trainData <- data[-testIndexes, ]
#   trainX = trainData$x
#   trainY = trainData$y
#   for (j in c(1:5)){
#     thetaA <- iterA <-  array(0, dim = 100)
#     for (k in c(1:100)){
#       SGDnonk1ChoiceWSAll = directFindLocal1SGDnonAdap(trainData$x, trainData$y, spSeq[j])
#       thetaA[k] = SGDnonk1ChoiceWSAll$thetaChoice
#       iterA[k] = sum(SGDnonk1ChoiceWSAll$sampleSize)
#     }
#     imThetaM[i,j] = mean(thetaA)
#     imIterM[i,j] = sum(iterA)
#     powerHat = exp(testData$x/10)*sin(testData$x)-sqrt(imThetaM[i,j]^2-imThetaM[i,j]+1)*(sin(imThetaM[i,j]*testData$x)+ cos(imThetaM[i,j]*testData$x))
#     imRMSEM[i,j] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
# 
#   }
#   imRMSE[i] = min(imRMSEM[i,])
#   imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
#   imIter[i] = imIterM[i,which(imRMSEM[i,]==imRMSE[i])]
# 
# }
# ntmImperfect = proc.time() - ptm
# 
# ptm <- proc.time()
# #Perform 10 fold cross validation
# for(i in c(1:10)){
#   #Segement your data by fold using the which() function 
#   testIndexes <- which(folds==i,arr.ind=TRUE)
#   testData <- data[testIndexes, ]
#   trainData <- data[-testIndexes, ]
#   trainX = trainData$x
#   trainY = trainData$y
#   for (j in c(1:5)){
#     thetaA <- iterA <-  array(0, dim = 100)
#     for (k in c(1:100)){
#       SGDnonk1ChoiceWSAll = directFindLocal2SGDnonAdap(trainData$x, trainData$y, spSeq[j])
#       thetaA[k] = SGDnonk1ChoiceWSAll$thetaChoice
#       iterA[k] = sum(SGDnonk1ChoiceWSAll$sampleSize)
#     }
#     pThetaM[i,j] = mean(thetaA)
#     pIterM[i,j] = sum(iterA)
#     powerHat1 = exp(testData$x/10)*sin(testData$x)-abs(pThetaM[i,j]+1)*(sin(pThetaM[i,j]*testData$x)+ cos(pThetaM[i,j]*testData$x))
#     pRMSEM[i,j] = sqrt(sum((powerHat1 - testData$y)^2)/length(powerHat))
#   }
# 
#   pRMSE[i] = min(pRMSEM[i,])
#   pTheta[i] = pThetaM[i,which(pRMSEM[i,]==pRMSE[i])]
#   pIter[i] = pIterM[i,which(pRMSEM[i,]==pRMSE[i])]
# }
# ntmPerfect = proc.time() - ptm
# 
# mean(imRMSE)
# mean(imTheta)
# mean(pRMSE)
# mean(pTheta)
# imTheta
# pTheta
# imIter
# pIter
# imIterM
# pIterM
# apply(imIterM,1,sum)
# apply(pIterM,1,sum)
# sum(apply(imIterM,1,sum))
# sum(apply(pIterM,1,sum))
# imThetaM
# pThetaM
# ntmImperfect
# ntmPerfect

testData = read.csv(file = testLoad)



imSample <- pSample <- imIter <- pIter <- imTime <- pTime <- imTheta <- pTheta <- imRMSE <- pRMSE <- array(0, dim = 10)
imSampleM <- pSampleM <- imIterM <- pIterM <- imThetaM <- pThetaM <- imRMSEM <- pRMSEM <- matrix(0, nrow = repetition, ncol = length(spSeq))
for (i in c(1:repetition)){
  ptm <- proc.time()
  for (j in c(1:length(spSeq))){
    SGDnonkChoiceWSAll = directFindLocal1SGDnonAdap(xSeq, ySeq, spSeq[j])
    imIterM[i,j] = SGDnonkChoiceWSAll$iters
    imSampleM[i,j] = sum(SGDnonkChoiceWSAll$sampleSize)
    imThetaM[i,j] = SGDnonkChoiceWSAll$thetaChoice
    powerHat = exp(xSeq/10)*sin(xSeq)-sqrt(imThetaM[i,j]^2-imThetaM[i,j]+1)*(sin(imThetaM[i,j]*xSeq)+ cos(imThetaM[i,j]*xSeq))
    imRMSEM[i,j] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
  }
  ntm = proc.time() - ptm
  imIter[i] = mean(imIterM[i,])
  imSample[i] = sum(imSampleM[i,])
  imTime[i] = as.array(ntm)[3]
  imRMSE[i] = min(imRMSEM[i,])
  imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
  powerHat = exp(testData$x/10)*sin(testData$x)-sqrt(imTheta[i]^2-imTheta[i]+1)*(sin(imTheta[i]*testData$x)+ cos(imTheta[i]*testData$x))
  imRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
  
  ptm <- proc.time()
  for (j in c(1:length(spSeq))){
    SGDnonkChoiceWSAll = directFindLocal2SGDnonAdap(xSeq, ySeq, spSeq[j])
    pIterM[i,j] = SGDnonkChoiceWSAll$iters
    pSampleM[i,j] = sum(SGDnonkChoiceWSAll$sampleSize)
    pThetaM[i,j] = SGDnonkChoiceWSAll$thetaChoice
    powerHat = exp(xSeq/10)*sin(xSeq)-abs(pThetaM[i,j]+1)*(sin(pThetaM[i,j]*xSeq)+ cos(pThetaM[i,j]*xSeq))
    pRMSEM[i,j] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
  }
  ntm = proc.time() - ptm
  pIter[i] = mean(pIterM[i,])
  pSample[i] = sum(pSampleM[i,])
  pTime[i] = as.array(ntm)[3]
  pRMSE[i] = min(pRMSEM[i,])
  pTheta[i] = pThetaM[i,which(pRMSEM[i,]==pRMSE[i])]
  powerHat = exp(testData$x/10)*sin(testData$x)-abs(pTheta[i]+1)*(sin(pTheta[i]*testData$x)+ cos(pTheta[i]*testData$x))
  pRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
}

imResult = data.frame(iter=imIter, sample = imSample, time = imTime, theta = imTheta, RMSE = imRMSE)
pResult = data.frame(iter=pIter, sample = pSample, time = pTime, theta = pTheta, RMSE = pRMSE)

apply(imResult, 2, mean, na.rm = TRUE)
apply(imResult, 2, sd, na.rm = TRUE)
apply(imResult, 2, median, na.rm = TRUE)
apply(pResult, 2, mean, na.rm = TRUE)
apply(pResult, 2, sd, na.rm = TRUE)
apply(pResult, 2, median, na.rm = TRUE)

setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/numExResult')
write.csv(imResult, file = imperfectSave, row.names = F)
write.csv(pResult, file = perfectSave, row.names = F)
