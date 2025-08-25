# numEx from Tuo 2015
# replacement T

library(ggplot2)
library(reshape2)
library(splines)
library(plyr)
library(rpart)
library(GPLTR)

startPercent = 0.1
addPercent = 1
kappa =0.95
spSeq = seq(5,10,length.out = 5)
repetition = 1000
setIter = 1000

testLoad = "1000TestDataVariedError2d.csv"
trainLoad = "1000TrainDataVariedError2d.csv"
imperfectSave = "SASGDCART01_2dKappa095.csv"


# set number of total data and generate data
setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/data')
data = read.csv(file = trainLoad)
numSeq = dim(data)[1]
x1Seq = data$x1
x2Seq = data$x2
ySeq = data$y

#####SGD with stratified sample####
# setting in stratified sampling
# a small amount of sample assigned to each stratum
n0= ceiling(startPercent/10*numSeq)


# the following function calculates the sum of error squre from imperfect model
# sumPolyLossLocal1SGDnonAdap<- function(theta, x1Seq., x2Seq., ySeq.){ 
#   powerHat = -(x1Seq.-theta/2)^4-(x2Seq.-theta)^2+4
#   errorFarm = (powerHat - ySeq.)^2
#   return(errorFarm)
# }

sumPolyLossLocal1SGDnonAdap<- function(theta, x1Seq., x2Seq., ySeq.){ 
  u1 = 2*theta
  u2 = theta
  powerHat = (1-exp(-1/(2*x2Seq.)))*(u1*1000*x1Seq.^3+1900*x1Seq.^2 +2092*x1Seq.+60)/(100*u2*x1Seq.^3+500*x1Seq.^2 +4*x1Seq.+20)
  errorFarm = (powerHat - ySeq.)^2
  return(errorFarm)
}

# the following function takes data as input and calibrates theta accordingly
directFindLocal1SGDnonAdap <- function(data, start){
  # set the limit of sample in SGD 
  sampleLimit = dim(data)[1]
  
  # number of stratum
  numStra = 4
  
  # assign data to each stratum
  samNumSeq = array(0, dim = numStra)
  samNumSeqM = as.data.frame(matrix(samNumSeq, nrow = 1, ncol = numStra))
  
  Lk = 1*10^-3
  ita = 1.5
  testTheta = kappa
  
  # set the small amount used in SGD
  originalH = c(0.01)
  currPara = start
  # set the limit of iteration number
  iter = setIter
  
  # have a matrix to record the number of sample size at each iteration
  sampleSize = array(0,dim = iter)
  # have a matrix to record the change of theta
  paras1 = matrix(0, nrow = iter+1, ncol = length(start))
  paras1[1,] = start
  
  # a matrix to record the change of weight on each strata
  # samNumMat = matrix(0, nrow = iter, ncol = numStra)
  # set starting sample number at each iteration in SGD
  startSamNum = ceiling(startPercent*dim(data)[1])
  samNum = startSamNum
  
  stratifiedData = list()
  
  sampleData = data[sample(c(1:dim(data)[1]), startSamNum, replace = T),]
  samx1Seq = sampleData$x1
  samx2Seq = sampleData$x2
  samySeq = sampleData$y
  
  # loop for each iteration
  for (j in c(1:iter)){
    # set a default step size
    step = 1/j
    print(j)
    actualIter = j
    
    # loop for each parameter
    for (k in c(1:length(start))){
      #re-stratify data
      h = originalH
      theta1 = currPara+h
      theta2 = currPara-h
      # sampleData = data[sample(c(1:dim(data)[1]), startSamNum, replace = T),]
      errors1 = sumPolyLossLocal1SGDnonAdap(theta1, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      errorFarm1 = mean(errors1)
      errors2 = sumPolyLossLocal1SGDnonAdap(theta2, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      errorFarm2 = mean(errors2)
      diff = (errors1-errors2)/(2*h)
      sumError = (errorFarm1-errorFarm2)/(2*h)
      dataDiff = as.data.frame(cbind(samx1Seq, samx2Seq, samySeq, diff))
      colnames(dataDiff) = c("x1", "x2", "y", "diff")
      
      fit = rpart(diff~.-y, data = dataDiff, control = rpart.control(cp = 0.0001))
      # printcp(fit)
      
      # treeGaussian = best.tree.BIC.AIC(fit, dataDiff, "diff", c("x1","x2"), family = "gaussian")
      # fit.pruned = treeGaussian$tree$BIC
      fit.pruned = prune(fit, fit$cptable[min(max(2,min(8,which.min(fit$cptable[,4]))),dim(fit$cptable)[1]),1])
      
      pred = predict(fit.pruned, data)
      
      splitLevel = fit.pruned$frame$yval[which(fit.pruned$frame$var=="<leaf>")]
      
      numStra = length(splitLevel) 
      dataNumSeq = array()
      for (i in c(1:(numStra))){
        stratifiedData[[i]] = data[which(pred ==splitLevel[i]),]
        dataNumSeq[i] = dim(stratifiedData[[i]])[1]
      }
      pSeq = dataNumSeq/numSeq
      #      if (j ==1){
      wSeq = pSeq
      #      }
      
      #get a stratified sample set from data for the first step
      samNumSeq = as.integer(wSeq *(samNum-n0*numStra))+n0
      samNumSeqM = rbind.fill(samNumSeqM, as.data.frame(t(samNumSeq)))
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
        errors1[[ii]] = sumPolyLossLocal1SGDnonAdap(theta1, x1Seq. = samDataSet[[ii]]$x1, x2Seq. = samDataSet[[ii]]$x2, ySeq.=samDataSet[[ii]]$y)
        errorFarm1[ii] = mean(errors1[[ii]])
        errors2[[ii]] = sumPolyLossLocal1SGDnonAdap(theta2, x1Seq. = samDataSet[[ii]]$x1, x2Seq. = samDataSet[[ii]]$x2, ySeq.=samDataSet[[ii]]$y)
        errorFarm2[ii] = mean(errors2[[ii]])
        diff[[ii]] = (errors1[[ii]]-errors2[[ii]])/(2*h[k])
        sumError[ii] = (errorFarm1[ii]-errorFarm2[ii])/(2*h[k])
      }
      
      for (ii in c(1:(numStra))){
        productSq[ii] = sum((diff[[ii]]*sum(sumError*pSeq)-sum(sumError*pSeq)^2)^2)
        productKSq[ii] = var(diff[[ii]]*sum(sumError*pSeq))
      }
      innerLeft = sum(productKSq*(pSeq^2)/(wSeq*samNumSeq))/samNum
      innerRight = (testTheta^2)*(sum(sumError*pSeq)^4)
      
      # conduct inner product test
      while (innerLeft>innerRight){
        if (samNum >= sampleLimit){
          break()
        }
        # add additional samples
        newSamNum = min(ceiling(startSamNum*addPercent)+samNum, sampleLimit)
        addSamNum = newSamNum- samNum
        addSamNumSeq = ceiling(wSeq *(addSamNum-n0*numStra))+n0
        addSamDataSet = list(stratifiedData[[1]][sample(c(1:dim(stratifiedData[[1]])[1]), max(addSamNumSeq[1],1), replace = T),])
        
        for (ii in c(1:(numStra-1))){
          addSamDataSet[[ii+1]] = stratifiedData[[ii+1]][sample(c(1:dim(stratifiedData[[ii+1]])[1]), max(addSamNumSeq[1],1), replace = T),]
        }
        
        # update
        samNum = newSamNum
        
        for (ii in c(1:(numStra))){
          errors1[[ii]] = sumPolyLossLocal1SGDnonAdap(theta1, x1Seq. = samDataSet[[ii]]$x1, x2Seq. = samDataSet[[ii]]$x2, ySeq.=samDataSet[[ii]]$y)
          errorFarm1[ii] = mean(errors1[[ii]])
          errors2[[ii]] = sumPolyLossLocal1SGDnonAdap(theta2, x1Seq. = samDataSet[[ii]]$x1, x2Seq. = samDataSet[[ii]]$x2, ySeq.=samDataSet[[ii]]$y)
          errorFarm2[ii] = mean(errors2[[ii]])
          diff[[ii]] = (errors1[[ii]]-errors2[[ii]])/(2*h[k])
          sumError[ii] = (errorFarm1[ii]-errorFarm2[ii])/(2*h[k])
        }
        
        for (ii in c(1:(numStra))){
          productSq[ii] = sum((diff[[ii]]*sum(sumError*pSeq)-sum(sumError*pSeq)^2)^2)
          productKSq[ii] = var(diff[[ii]]*sum(sumError*pSeq))
        }
        innerLeft = sum(productKSq*(pSeq^2)/(wSeq*samNumSeq))/samNum
        innerRight = (testTheta^2)*(sum(sumError*pSeq)^4)
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
      samx1Seq = samDataSet[[1]]$x1
      samx2Seq = samDataSet[[1]]$x2
      samySeq = samDataSet[[1]]$y
      straIndex = array(0, dim = numStra)
      straIndex[1] = length(samDataSet[[1]]$x1)
      for (ii in c(1:(numStra-1))){
        samx1Seq = c(samx1Seq, samDataSet[[ii+1]]$x1)
        samx2Seq = c(samx2Seq, samDataSet[[ii+1]]$x2)
        samySeq = c(samySeq, samDataSet[[ii+1]]$y)
        straIndex[ii+1] = straIndex[ii]+length(samDataSet[[ii+1]]$x1)
      }
      newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      newLoss = sum(newErrors)/length(newErrors)
      currErrors = sumPolyLossLocal1SGDnonAdap(currPara[1], x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
      currLoss = sum(currErrors)/length(currErrors)
      
      while (newLoss>currLoss-sum(sumError*pSeq)^2/(2*currLk)){
        currLk = currLk*ita
        newErrors = sumPolyLossLocal1SGDnonAdap(currPara[1]+dire/currLk, x1Seq. = samx1Seq, x2Seq. = samx2Seq, ySeq.=samySeq)
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
    #samNumMat[j,] = samNumSeq
    samNum = startSamNum
    if (abs(paras1[j+1,1]- paras1[j,1])/paras1[j,1] <= 10^-3) {
      print(j)
      paras1[c((j+1):dim(paras1)[1]),1]= paras1[j+1,1]
      break
    }
  }
  thetaChoice = paras1[iter]  
  newList <- list("thetaChoice" = thetaChoice, "sampleSize" = sampleSize, "proportion" = samNumSeqM, "paras" = paras1, "iters" = actualIter)
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
    #powerHat = -(xSeq-imThetaM[i,j])^2+4
    powerHat = (1-exp(-1/(2*data$x2)))*(2*imThetaM[i,j]*1000*data$x1^3+1900*data$x1^2 +2092*data$x1+60)/(100*imThetaM[i,j]*data$x1^3+500*data$x1^2 +4*data$x1+20)
    imRMSEM[i,j] = sqrt(sum((powerHat - ySeq)^2)/length(powerHat))
  }
  ntm = proc.time() - ptm
  imIter[i] = mean(imIterM[i,])
  imSample[i] = sum(imSampleM[i,])
  imTime[i] = as.array(ntm)[3]
  imRMSE[i] = min(imRMSEM[i,])
  imTheta[i] = imThetaM[i,which(imRMSEM[i,]==imRMSE[i])]
  #powerHat = -(testData$x-imTheta[i])^2+4
  powerHat = (1-exp(-1/(2*testData$x2)))*(2*imTheta[i]*1000*testData$x1^3+1900*testData$x1^2 +2092*testData$x1+60)/(100*imTheta[i]*testData$x1^3+500*testData$x1^2 +4*testData$x1+20)
  imRMSE[i] = sqrt(sum((powerHat - testData$y)^2)/length(powerHat))
  
}

imResult = data.frame(iter=imIter, sample = imSample, time = imTime, theta = imTheta, RMSE = imRMSE)

apply(imResult, 2, mean, na.rm = TRUE)
apply(imResult, 2, sd, na.rm = TRUE)
apply(imResult, 2, median, na.rm = TRUE)

setwd('C:/Users/liubi/Desktop/IOE 366/numEx/SS/numExResult')
write.csv(imResult, file = imperfectSave, row.names = F)
