#numerical example SGD only
# numEx from Tuo 2015

library(ggplot2)
library(reshape2)
library(splines)


startPercent = 0.1
spSeq = seq(5,10,length.out = 5)
repetition = 1000
setIter = 1000


testLoad = "1000TestDataVariedError2d.csv"
trainLoad = "1000TrainDataVariedError2d.csv"
imperfectSave = "1000Train01SGD2d.csv"


# set number of total data and generate data
setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/data')
data = read.csv(file = trainLoad)
numSeq = dim(data)[1]
x1Seq = data$x1
x2Seq = data$x2
ySeq = data$y

#####SGD####

# the following function calculates the sum of error squre from imperfect model
sumPolyLossLocal1SGDnonAdap<- function(theta, x1Seq., x2Seq., ySeq.){ 
  u1 = 2*theta
  u2 = theta
  powerHat = (1-exp(-1/(2*x2Seq.)))*(u1*1000*x1Seq.^3+1900*x1Seq.^2 +2092*x1Seq.+60)/(100*u2*x1Seq.^3+500*x1Seq.^2 +4*x1Seq.+20)
  errorFarm = (powerHat - ySeq.)^2
  return(errorFarm)
}



# the following function takes data as input and calibrates theta accordingly
directFindLocal1SGDnonAdap <- function(data, start){
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
  samNum = ceiling(startPercent*dim(data)[1])
  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter =j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #get a sample set from data for current step
      samInx = sample(c(1:dim(data)[1]), samNum, replace = T)
      samx1Seq = data$x1[samInx]
      samx2Seq = data$x2[samInx]
      samySeq = data$y[samInx]
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = sumPolyLossLocal1SGDnonAdap(theta1, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      errorFarm1 = mean(errors1)
      errors2 = sumPolyLossLocal1SGDnonAdap(theta2, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      errorFarm2 = sum(errors2)/length(errors2)
      sumError = (errorFarm1-errorFarm2)/(2*h[k])
      dire = -sumError
      
      # conduct backtracking line search
      ak = 1 + (var(errors1-errors2)/(samNum*sumError^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal1SGDnonAdap(currPara[1], x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sumError^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
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
    SGDnonkChoiceWSAll = directFindLocal1SGDnonAdap(data, spSeq[j])
    imIterM[i,j] = SGDnonkChoiceWSAll$iters
    imSampleM[i,j] = imIterM[i,j]*dim(data)[1]*startPercent
    imThetaM[i,j] = SGDnonkChoiceWSAll$thetaChoice
    powerHat = (1-exp(-1/(2*data$x2)))*(2*imThetaM[i,j]*1000*data$x1^3+1900*data$x1^2 +2092*data$x1+60)/(100*imThetaM[i,j]*data$x1^3+500*data$x1^2 +4*data$x1+20)
    imRMSEM[i,j] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
  }
  ntm = proc.time() - ptm
  imIter[i] = mean(imIterM[i,])
  imSample[i] = sum(imSampleM[i,])
  imTime[i] = as.array(ntm)[3]
  imRMSE[i] = min(imRMSEM[i,])
  imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
  powerHat = (1-exp(-1/(2*testData$x2)))*(2*imTheta[i]*1000*testData$x1^3+1900*testData$x1^2 +2092*testData$x1+60)/(100*imTheta[i]*testData$x1^3+500*testData$x1^2 +4*testData$x1+20)
  imRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
  
}

imResult = data.frame(iter=imIter, sample = imSample, time = imTime, theta = imTheta, RMSE = imRMSE)

apply(imResult, 2, mean, na.rm = TRUE)
apply(imResult, 2, sd, na.rm = TRUE)
apply(imResult, 2, median, na.rm = TRUE)


setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/numExResult')
write.csv(imResult, file = imperfectSave, row.names = F)

