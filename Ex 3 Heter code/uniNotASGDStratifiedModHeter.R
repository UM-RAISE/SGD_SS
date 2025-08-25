# numEx from Tuo 2015
# replacement T

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

startPercent = 0.1
spSeq = seq(5,10,length.out = 5)
repetition = 1000

# testLoad = "1000TestData.csv"
# trainLoad = "1000TrainData.csv"
# imperfectSave = "im1000TrainSSGDEx1.csv"
# perfectSave = "p1000TrainSSGDEx1.csv"

testLoad = "1000TestDataVariedError4.csv"
trainLoad = "1000TrainDataVariedError4.csv"
imperfectSave = "im1000TrainSSGD01Ex2Ver4.csv"


# set number of total data and generate data
setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/data')
data = read.csv(file = trainLoad)
numSeq = dim(data)[1]
xSeq = data$x
ySeq = data$y


#####SGD with stratified sample####
# setting in stratified sampling
# a small amount of sample assigned to each stratum
n0= ceiling(0.1*startPercent*numSeq)
# number of stratum
numStra = 4
# set the knot to seperate stratum
knot = seq(0, 4, length.out = numStra+1)[-1]

# assign data to each stratum
samNumSeq = array(0, dim = numStra)
dataNumSeq = array(0, dim = numStra)
data1 = data[which(data$x<knot[1]),]
dataNumSeq[1] = dim(data1)[1]
stratifiedData =list(data1)
for (i in c(1:(numStra-1))){
  stratifiedData[[i+1]] = data[which(data$x>knot[i]&data$x<= knot[i+1]),]
  dataNumSeq[i+1] = dim(stratifiedData[[i+1]])[1]
}
pSeq = dataNumSeq/numSeq
# tentatively assign weight to each stratum based on distribution
wSeq = pSeq


# the following function calculates the sum of error squre from imperfect model
sumPolyLossLocal1SGDnonAdap<- function(theta, xSeq.,ySeq.){ 
  powerHat = -(xSeq.-theta)^2+4
  errorFarm = (powerHat - ySeq.)^2
  return(errorFarm)
}


# the following function takes data as input and calibrates theta accordingly
directFindLocal1SGDnonAdap <- function(data, start){
  # set the limit of sample in SGD 
  sampleLimit = dim(data)[1]
  
  # assign data to each stratum
  samNumSeq = array(0, dim = numStra)
  samNumSeqM = matrix(samNumSeq, nrow = 1, ncol = numStra)
  dataNumSeq = array(0, dim = numStra)
  data = data[order(data$x),]
  data1 = data[which(data$x<knot[1]),]
  dataNumSeq[1] = dim(data1)[1]
  stratifiedData =list(data1)
  for (i in c(1:(numStra-1))){
    stratifiedData[[i+1]] = data[which(data$x>knot[i]&data$x<= knot[i+1]),]
    dataNumSeq[i+1] = dim(stratifiedData[[i+1]])[1]
  }
  pSeq = dataNumSeq/numSeq
  # tentatively assign weight to each stratum based on distribution
  wSeq = pSeq
  
  Lk = 1*10^-3
  ita = 1.5
  testTheta = 0.9
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = 500
  
  # have a matrix to record the number of sample size at each iteration
  sampleSize = array(0,dim = iter)
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  # a matrix to record the change of weight on each strata
  samNumMat = matrix(0, nrow = iter, ncol = numStra)
  # set starting sample number at each iteration in SGD
  startSamNum = ceiling(startPercent*dim(data)[1])
  samNum = startSamNum
  
  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter = j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #get a stratified sample set from data for the first step
      samNumSeq = as.integer(wSeq *(samNum-n0*numStra))+n0
      samNumSeqM = rbind(samNumSeqM,samNumSeq)
      samDataSet = list(stratifiedData[[1]][sample(c(1:dim(stratifiedData[[1]])[1]), samNumSeq[1], replace = T),])
      
      for (ii in c(1:(numStra-1))){
        samDataSet[[ii+1]] = stratifiedData[[ii+1]][sample(c(1:dim(stratifiedData[[ii+1]])[1]), samNumSeq[ii+1], replace = T),]
      }
      
      h = array(0, dim = length(start))
      h[k]=originalH[k]
      
      # calculate the gradient sumError and set -sumError to be direction 
      sumError = 0
      theta1 = currPara[k]+h[k] 
      theta2 = currPara[k]-h[k]
      errors1 = list()
      errors2 = list()
      diff = list()
      errorFarm1 <- errorFarm2 <- sumError <- array(0, dim = numStra)
      productSq <- productKSq <- array(0, dim = numStra)
      for (ii in c(1:(numStra))){
        errors1[[ii]] = sumPolyLossLocal1SGDnonAdap(theta1, xSeq. = samDataSet[[ii]]$x, ySeq.=samDataSet[[ii]]$y)
        errorFarm1[ii] = mean(errors1[[ii]])
        errors2[[ii]] = sumPolyLossLocal1SGDnonAdap(theta2, xSeq. = samDataSet[[ii]]$x, ySeq.=samDataSet[[ii]]$y)
        errorFarm2[ii] = mean(errors2[[ii]])
        diff[[ii]] = (errors1[[ii]]-errors2[[ii]])/(2*h[k])
        sumError[ii] = (errorFarm1[ii]-errorFarm2[ii])/(2*h[k])
      }
      
      for (ii in c(1:(numStra))){
        productSq[ii] = sum((diff[[ii]]*sum(sumError*pSeq)-sum(sumError*pSeq)^2)^2)
        productKSq[ii] = var(diff[[ii]]*sum(sumError*pSeq))
      }
      
      
      dire = -sum(sumError*pSeq)
      
      sumDiffSq = array(0, dim = numStra)
      for (ii in c(1:(numStra))){
        sumDiffSq[ii] = sum(diff[[ii]]^2)
      }
      varHat = sum(sumDiffSq*pSeq^2/(samNumSeq*(samNumSeq-1)))
      
      # conduct backtracking line search
      ak = 1 + (varHat/(sum(sumError*pSeq)^2))
      zetak = max(c(1,2/ak))
      currLk = Lk/zetak
      samxSeq = samDataSet[[1]]$x
      samySeq = samDataSet[[1]]$y
      straIndex = array(0, dim = numStra)
      straIndex[1] = length(samDataSet[[1]]$x)
      for (ii in c(1:(numStra-1))){
        samxSeq = c(samxSeq, samDataSet[[ii+1]]$x)
        samySeq = c(samySeq, samDataSet[[ii+1]]$y)
        straIndex[ii+1] = straIndex[ii]+length(samDataSet[[ii+1]]$x)
      }
      newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, xSeq. = samxSeq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal1SGDnonAdap(currPara[1], xSeq. = samxSeq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sum(sumError*pSeq)^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk,   xSeq. = samxSeq, ySeq.=samySeq)
        newLoss = sum(newErrors)/length(newErrors)
      }
      
      # move to next step
      step = 1/currLk
      h = c(0,0)
      moving = dire*step
      h[k] = moving
      movingIndex = 0
      currk = currPara[k] + moving
      while (currk<=-100 || currk>=100){
        moving = moving/2
        h[1] = moving
        currk = currPara[1]+h[1]
        # currk = (currPara[1]+h[1]) + (currPara[2]+h[2])*currWS + (currPara[3]+h[3])*currTI 
        # currCt = (currPara[5]+h[5]) + (currPara[6]+h[6])*currWS + (currPara[7]+h[7])*currAD
        movingIndex  =  movingIndex+1
        #print(movingIndex)
      }
      currPara[k] = currPara[k] + moving
      
      samSigma = array(0, dim = numStra)
      for (ii in c(1:(numStra))){
        samSigma[ii] = sqrt(productKSq[ii])
      }
      
      wSeq = pSeq*samSigma/sum(pSeq*samSigma)
    }
    paras1[j+1,] = currPara
    
    # record current sample size
    sampleSize[j] = samNum
    samNumMat[j,] = samNumSeq
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]  
  newList <- list("thetaChoice" = thetaChoice, "sampleSize" = sampleSize, "proportion" = samNumMat, "paras" = paras1, "iters" = actualIter, wTrend = samNumSeqM)
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
    imSampleM[i,j] = sum(SGDnonkChoiceWSAll$sampleSize)
    imThetaM[i,j] = SGDnonkChoiceWSAll$thetaChoice
    powerHat = -(xSeq-imThetaM[i,j])^2+4
    imRMSEM[i,j] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
  }
  ntm = proc.time() - ptm
  imIter[i] = mean(imIterM[i,])
  imSample[i] = sum(imSampleM[i,])
  imTime[i] = as.array(ntm)[3]
  imRMSE[i] = min(imRMSEM[i,])
  imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
  powerHat = -(testData$x-imTheta[i])^2+4
  imRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
}

imResult = data.frame(iter=imIter, sample = imSample, time = imTime, theta = imTheta, RMSE = imRMSE)

apply(imResult, 2, mean, na.rm = TRUE)
apply(imResult, 2, sd, na.rm = TRUE)
apply(imResult, 2, median, na.rm = TRUE)

setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/numExResult')
write.csv(imResult, file = imperfectSave, row.names = F)

# losses1 <- thetas <- seq(0,5,length.out = 1000)
# for (i in c(1:length(thetas))){
#   errors1 = sumPolyLossLocal1SGDnonAdap(thetas[i], xSeq. = xSeq, ySeq.=ySeq)
#   losses1[i] = sum(errors1)/length(errors1)
# }
# pdf(file = paste("C:/Users/liubi/Desktop/IOE 366/numEx/SS/graph/Ex2HeterLoss", ".pdf"),
#     width = 4,
#     height = 4)
# plot(thetas,losses1, main= paste('Theta', thetas[1], "to", thetas[1000]) )
# dev.off()
