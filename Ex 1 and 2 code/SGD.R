#numerical example SGD only
# numEx from Tuo 2015

library(ggplot2)
library(reshape2)
library(splines)


startPercent = 0.05
spSeq = seq(5,10,length.out = 5)
repetition = 1000

testLoad = "1000TestData.csv"
trainLoad = "1000TrainData.csv"
imperfectSave = "im1000TrainSGDEx1.csv"
perfectSave = "p1000TrainSGDEx1.csv"



# set number of total data and generate data
setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/data')
data = read.csv(file = trainLoad)
numSeq = dim(data)[1]
xSeq = data$x
ySeq = data$y

#####SGD####

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
  Lk = 1*10^-3
  ita = 1.5
  testTheta = 0.9
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = 500
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  # set sample number at each iteration in SGD
  samNum = ceiling(startPercent*length(xSeq))
  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter =j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #get a sample set from data for current step
      samInx = sample(c(1:length(xSeq)), samNum, replace = T)
      samxSeq = xSeq[samInx]
      samySeq = ySeq[samInx]
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = sumPolyLossLocal1SGDnonAdap(theta1, xSeq. = samxSeq, ySeq.=samySeq)
      errorFarm1 = sum(errors1)/length(errors1)
      errors2 = sumPolyLossLocal1SGDnonAdap(theta2, xSeq. = samxSeq, ySeq.=samySeq)
      errorFarm2 = sum(errors2)/length(errors2)
      sumError = (errorFarm1-errorFarm2)/(2*h[k])
      dire = -sumError
      
      # conduct backtracking line search
      ak = 1 + (var(errors1-errors2)/(samNum*sumError^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal1SGDnonAdap(currPara[1], xSeq. = samxSeq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sumError^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
        newLoss = sum(newErrors)/length(newErrors)
      }

      
      # move to next step
      step = 1/currLk
      h = array(0, dim = length(start))
      moving = dire*step
      h[k] = moving
      movingIndex = 0
      currk = currPara[k] + moving
      # give theta a limit range
      while (currk<=-100 || currk>=100){
        moving = moving/2
        h[1] = moving
        currk = currPara[1]+h[1]
        movingIndex  =  movingIndex+1
      }
      currPara[k] = currPara[k] + moving
    }
    paras1[j+1,] = currPara
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]
  newList <- list("thetaChoice" = thetaChoice, "iters" = actualIter, "paras" = paras1)
  return (newList)
}

directFindLocal2SGDnonAdap <- function(xSeq, ySeq, start){
  Lk = 1*10^-3
  ita = 1.5
  testTheta = 0.9
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = 500
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  # set sample number at each iteration in SGD
  samNum = ceiling(startPercent*length(xSeq))
  
  # loop for each iteration
  for (j in c(1:iter)){
    # step size
    step = 1/j
    print(j)
    actualIter = j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #get a sample set from data for current step
      samInx = sample(c(1:length(xSeq)), samNum, replace = T)
      samxSeq = xSeq[samInx]
      samySeq = ySeq[samInx]
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = sumPolyLossLocal2SGDnonAdap(theta1, xSeq. = samxSeq, ySeq.=samySeq)
      errorFarm1 = sum(errors1)/length(errors1)
      errors2 = sumPolyLossLocal2SGDnonAdap(theta2, xSeq. = samxSeq, ySeq.=samySeq)
      errorFarm2 = sum(errors2)/length(errors2)
      sumError = (errorFarm1-errorFarm2)/(2*h[k])
      dire = -sumError
      
      # conduct backtracking line search
      ak = 1 + (var(errors1-errors2)/(samNum*sumError^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      newErrors = sumPolyLossLocal2SGDnonAdap(currPara[1]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal2SGDnonAdap(currPara[1], xSeq. = samxSeq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sumError^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal2SGDnonAdap(currPara[1]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
        newLoss = sum(newErrors)/length(newErrors)
      }

      
      # move to next step
      step = 1/currLk
      h = array(0, dim = length(start))
      moving = dire*step
      h[k] = moving
      movingIndex = 0
      currk = currPara[k] + moving
      # give theta a limit range
      while (currk<=-100 || currk>=100){
        moving = moving/2
        h[1] = moving
        currk = currPara[1]+h[1]
        movingIndex  =  movingIndex+1
      }
      currPara[k] = currPara[k] + moving
    }
    paras1[j+1,] = currPara
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]
  newList <- list("thetaChoice" = thetaChoice, "iters" = actualIter, "paras" = paras1)
  return (newList)
}



testData = read.csv(file = testLoad)


imSample <- pSample <- imIter <- pIter <- imTime <- pTime <- imTheta <- pTheta <- imRMSE <- pRMSE <- array(0, dim = 10)
imSampleM <- pSampleM <- imIterM <- pIterM <- imThetaM <- pThetaM <- imRMSEM <- pRMSEM <- matrix(0, nrow = repetition, ncol = length(spSeq))
for (i in c(1:repetition)){
  ptm <- proc.time()
  for (j in c(1:length(spSeq))){
    SGDnonkChoiceWSAll = directFindLocal1SGDnonAdap(xSeq, ySeq, spSeq[j])
    imIterM[i,j] = SGDnonkChoiceWSAll$iters
    imSampleM[i,j] = imIterM[i,j]*dim(data)[1]*startPercent
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
    pSampleM[i,j] = pIterM[i,j]*dim(data)[1]*startPercent
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
# theta = SGDnonkChoiceWSAll
# powerHat = exp(xSeq/10)*sin(xSeq)-sqrt(theta^2-theta+1)*(sin(theta*xSeq)+ cos(theta*xSeq))
# errorFarm = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
# 
# theta = SGDnonk1ChoiceWSAll
# powerHat1 = exp(xSeq/10)*sin(xSeq)-abs(theta+1)*(sin(theta*xSeq)+ cos(theta*xSeq))
# errorFarm1 = sqrt(sum((powerHat1 - ySeq)^2)/length(powerHat))
# errorFarm
# errorFarm1
# SGDnonkChoiceWSAll
# SGDnonk1ChoiceWSAll
# ntmAdaptive
# ntmAdaptive1
# 
# 
# # cross validation
# 
# #Randomly shuffle the data
# data<-data[sample(nrow(data)),]
# 
# #Create 10 equally size folds
# folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
# 
# imTheta <- pTheta <- imRMSE <- pRMSE <- array(0, dim = 10)
# #Perform 10 fold cross validation
# for(i in 1:10){
#   #Segement your data by fold using the which() function 
#   testIndexes <- which(folds==i,arr.ind=TRUE)
#   testData <- data[testIndexes, ]
#   trainData <- data[-testIndexes, ]
#   trainX = trainData$x
#   trainY = trainData$y
#   SGDnonkChoiceWSAll = directFindLocal1SGDnonAdap(trainData$x, trainData$y, -5)
#   imTheta[i] = SGDnonkChoiceWSAll
#   powerHat = exp(testData$x/10)*sin(testData$x)-sqrt(imTheta[i]^2-imTheta[i]+1)*(sin(imTheta[i]*testData$x)+ cos(imTheta[i]*testData$x))
#   imRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
#   
#   SGDnonk1ChoiceWSAll = directFindLocal2SGDnonAdap(trainData$x, trainData$y, -5)
#   pTheta[i] = SGDnonk1ChoiceWSAll
#   powerHat1 = exp(testData$x/10)*sin(testData$x)-abs(pTheta[i]+1)*(sin(pTheta[i]*testData$x)+ cos(pTheta[i]*testData$x))
#   pRMSE[i] = sqrt(sum((powerHat1 - testData$y)^2)/length(powerHat))
# }
# 
# mean(imRMSE)
# mean(imTheta)
# mean(pRMSE)
# mean(pTheta)
# imTheta
# pTheta

# # cross validation with 5 starting point and the best selected
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
#       iterA[k] =SGDnonk1ChoiceWSAll$iters
#     }
#     imThetaM[i,j] = mean(thetaA)
#     imIterM[i,j] = sum(iterA)
#     powerHat = exp(testData$x/10)*sin(testData$x)-sqrt(imThetaM[i,j]^2-imThetaM[i,j]+1)*(sin(imThetaM[i,j]*testData$x)+ cos(imThetaM[i,j]*testData$x))
#     imRMSEM[i,j] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
#   }
#   imRMSE[i] = min(imRMSEM[i,])
#   imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
#   imIter[i] = imIterM[i,which(imRMSEM[i,]==imRMSE[i])]
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
#       iterA[k] =SGDnonk1ChoiceWSAll$iters
#     }
#     pThetaM[i,j] = mean(thetaA)
#     pIterM[i,j] = sum(iterA)
#     powerHat1 = exp(testData$x/10)*sin(testData$x)-abs(pThetaM[i,j]+1)*(sin(pThetaM[i,j]*testData$x)+ cos(pThetaM[i,j]*testData$x))
#     pRMSEM[i,j] = sqrt(sum((powerHat1 - testData$y)^2)/length(powerHat))
#   }
#   pRMSE[i] = min(pRMSEM[i,])
#   pTheta[i] = pThetaM[i,which(pRMSEM[i,]==pRMSE[i])]
#   pIter[i] = pIterM[i,which(pRMSEM[i,]==pRMSE[i])]
# }
# ntmPerfect = proc.time() - ptm
# 
# 
# t=proc.time()
# proc.time()-t
# samNum = ceiling(0.1*dim(trainData)[1])
# 
# mean(imRMSE)
# mean(imTheta)
# mean(pRMSE)
# mean(pTheta)
# imTheta
# pTheta
# imIter*samNum
# pIter*samNum
# imIterM*samNum
# pIterM*samNum
# apply(imIterM*samNum,1,sum)
# apply(pIterM*samNum,1,sum)
# sum(apply(imIterM*samNum,1,sum))
# sum(apply(pIterM*samNum,1,sum))
# imThetaM
# pThetaM
# ntmImperfect
# ntmPerfect


# RMSESeq1 <- RMSESeq <- kSeq <- seq(from =-1, to = 1, by = 0.01)
# for (i in c(1:length(kSeq))){
#   powerHat =  exp(xSeq/10)*sin(xSeq)-sqrt(kSeq[i]^2-kSeq[i]+1)*(sin(kSeq[i]*xSeq)+ cos(kSeq[i]*xSeq))
#   RMSESeq[i] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
# 
# 
#   powerHat1 = exp(xSeq/10)*sin(xSeq)-abs(kSeq[i]+1)*(sin(kSeq[i]*xSeq)+ cos(kSeq[i]*xSeq))
#   RMSESeq1[i] = sqrt(sum((powerHat1 - ySeq)^2)/length(powerHat))
# }
# plot(kSeq, RMSESeq)
# plot(kSeq, RMSESeq1)
# ploting = cbind(kSeq, RMSESeq, RMSESeq1)
# fploting = as.data.frame(ploting)
# 
# 
# 
# selected = sample(c(1:length(kSeq)), 15)
# 
# 
# p <- ggplot(fploting[selected,], aes(kSeq, RMSESeq)) +
#   geom_point(color = "blue",size =3) + 
#   geom_smooth(method = "loess", se = FALSE)
# 
# p + geom_line(data = fploting, aes(x = kSeq, y = RMSESeq))
# 
# 
# p1 <- ggplot(fploting[selected,], aes(kSeq, RMSESeq1)) +
#   geom_point(color = "blue",size =3) + 
#   geom_smooth(method = "loess", se = FALSE)
# 
# p1 + geom_line(data = fploting, aes(x = kSeq, y = RMSESeq1))
# 
# p <- ggplot(fploting[selected,], aes(kSeq, RMSESeq)) +
#   geom_point(color = "blue",size =3) + 
#   geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE)
# 
# p + geom_line(data = fploting, aes(x = kSeq, y = RMSESeq))
# 
# 
# p1 <- ggplot(fploting[selected,], aes(kSeq, RMSESeq1)) +
#   geom_point(color = "blue",size =3) + 
#   geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE)
# 
# p1 + geom_line(data = fploting, aes(x = kSeq, y = RMSESeq1))
# 
# 
# 
# library(plotly)
# Sys.setenv("plotly_username"="bliu34")
# Sys.setenv("plotly_api_key"="ydOnjbUikQaeAfNbNtal")
# 
# ySeq = exp(xSeq/10)*sin(xSeq)
# RMSEM <- RMSEM1 <-matrix(0, ncol = length(xSeq), nrow = length(kSeq))
# for (i in c(1:length(kSeq))){
#   powerHat =  exp(xSeq/10)*sin(xSeq)-sqrt(kSeq[i]^2-kSeq[i]+1)*(sin(kSeq[i]*xSeq)+ cos(kSeq[i]*xSeq))
#   RMSEM[i,] = abs(powerHat - ySeq)
#   
#   
#   powerHat1 = exp(xSeq/10)*sin(xSeq)-abs(kSeq[i]+1)*(sin(kSeq[i]*xSeq)+ cos(kSeq[i]*xSeq))
#   RMSEM1[i,] = abs(powerHat1 - ySeq)
# }
# 
# p <- plot_ly(x=xSeq, y=kSeq, z=RMSEM, type="surface") 
# p1 <- plot_ly(x=xSeq, y=kSeq, z=RMSEM1, type="surface") 
# 
# library(kernlab)
# library(GPfit)
# library(lhs)
# selected = sample(c(1:length(xSeq)), 10)
# testX = xSeq[selected]
# testY = ySeq[selected]
# testY = testY[order(testX)]
# testX = testX[order(testX)]
# test <- gausspr(x = testX,y= testY)
# GPmodel <- GP_fit(testX/max(testX), testY)
# plot(GPmodel)
# points(xSeq, ySeq, col="green")
# test
# ytest <- predict(test, testX)
# plot(xSeq, ySeq, type ="l")
# points(testX, testY, col="blue")
# lines(testX, ytest, col="red")
# 
# numSam = 50
# selectedX = sample(c(1:length(xSeq)), numSam)
# selectedX = selectedX[order(selectedX)]
# selectedK = sample(c(1:length(kSeq)), numSam)
# selectedK = selectedK[order(selectedK)]
# testX = xSeq[selectedX]
# testK = kSeq[selectedK]
# testRMSE = RMSEM[selectedK,selectedX]
# dfGaussian = as.data.frame(matrix(0,nrow = numSam^2, ncol = 3))
# colnames(dfGaussian) = c("X","K","RMSE")
# count = 1
# for (i in c(1:numSam)){
#   for (j in c(1:numSam)){
#     dfGaussian[count,]$X = testX[i]
#     dfGaussian[count,]$K = testK[j]
#     dfGaussian[count,]$RMSE = testRMSE[j,i]
#     count = count+1
#   }
# }
# 
# test <- gausspr(RMSE~.,data = dfGaussian)
# test
# testRMSE <- predict(test, dfGaussian)
# mtestRMSE = RMSEM[selectedK,selectedX]
# for (i in c(1:numSam)){
#   for (j in c(1:numSam)){
#     mtestRMSE[j,i] = testRMSE[((i-1)*numSam+j)]
#   }
# }
# p <- plot_ly(x=xSeq, y=kSeq, z=RMSEM, type="surface") %>% add_surface(x = testX, y = testK, z = mtestRMSE, colorscale = list(c(0,'#d1d1d1'),c(1,'#000000')))
