
#' Loads data and prepares datasets for training and testing
#' returns:
#' list of sets (training, validation, test) with auxiliary info
#' input:
#' dge 
#' dataFiltered - normalized reads with ENSEMBL gene names as rownames
#' disease label - name of the disease used in sample info 
#' matrixPath - path in which matrices are stored
#' picWidth
#' picHeight
#' seed - used to split samples
#' 
prepareDataSets = function(dge, dataFiltered, diseaseLabel, matrixPath, picWidth, picHeight,  seed){

  allX =  array(0, dim = c(dim(dataFiltered)[2],picWidth, picHeight))
  for(i in seq(1:dim(dataFiltered)[2]))
  { 
    tmpFrame =  (read.csv(paste(matrixPath,  colnames(dataFiltered)[i], ".txt", sep = ""),
                          header = F, sep = " "))
    sampleMatrix = as.matrix.data.frame(tmpFrame)
    allX[i,,] =  (sampleMatrix)   
  }
  
 
  library(abind)
  # otherSubset = which(  dge$samples$Description == "After resection"
  #                       | dge$samples$Description == "borderline tumor")
  # otherLabels = dge$samples$Description[otherSubset]
  diseaseLabel == "EC"
  ecCtrlSubset =  which(is.element(dge$samples$OriginalGroup, c("HC",
                                                                "EC", "CTRL")))
  # ocRelapseId = which(dge$samples$OriginalGroup == "OC" & is.element(dge$samples$Description, c("Ovarian cancer relapse", "relapse")))
  # dge$samples$OriginalGroup[ocRelapseId] = "Other" 
  # ocResectionId = which(  is.element(dge$samples$Description, c("After resection, before chemotherapy")))
  # dge$samples$OriginalGroup[ocResectionId] = "Other" 
  # otherSubset = c(otherSubset, ocRelapseId)
  # otherSubset = c(otherSubset, ocResectionId)
  # ocCtrlSubset =  which(is.element(dge$samples$OriginalGroup, c("HC", "OC", "CTRL"))
  #                       & dge$samples$Description != "normal tissue"
  #                       & dge$samples$Description != "borderline tumor"
  #                       & dge$samples$Description != "After resection"
  #                       & dge$samples$Description !=  "relapse"
  #                       & dge$samples$Description != "Ovarian cancer relapse"
  #                       & dge$samples$Description != "After resection, before chemotherapy")
  usedSubset = ecCtrlSubset
  #otherX =  allX[otherSubset, , , drop= F]
  allX = allX[usedSubset, , , drop= F]
  library(caret) 
  set.seed(seed)
  # splitFolds = createFolds(dge$samples$OriginalGroup[usedSubset], 5)
  # trainId = (unlist(splitFolds[1:2]))
  # valId = (unlist(splitFolds[5]))
  # 
  # 
  # testId = (unlist(splitFolds[3:4]))
  group = dge$samples$OriginalGroup[usedSubset]
  print(group)
  print(dge$samples$OriginalGroup[usedSubset])
  group[which(group == "CTRL")] = 0
  group[which(group == "HC")] = 0
  group[which(group == diseaseLabel)] = 1
  
  # trainX =  allX[trainId, , , drop = F]
  # testX = allX[testId, , , drop = F]
  # testX = abind(testX, otherX, along = 1)
  # valX =    allX[valId, , , drop = F]
  # trainY = group[trainId] 
  # testY = group[testId] 
  # testY = c(testY, rep(1, length(otherSubset)))
  # #testX = abind(testX, along = 1)
  # #testY = c(testY, rep(1, length(otherSubset)))
  # valY = group[valId]
  # trainY = as.numeric(trainY)
  # testY = as.numeric(testY)
  # valY = as.numeric(valY)
  # trainValY = c(trainY, valY)
  library(abind) 
 #trainValId= c(trainId, valId)
  library(abind)
  #trainValX = abind(trainX, valX, along = 1)
  output = list(#"trainValId" = trainValId,  "trainValX" = trainValX, "trainValY" = trainValY,
                #"testX" = testX, "testY" = testY, 
    "usedSubset" = usedSubset, #"otherSubset" = otherSubset, 
    allX = allX, allY = group)#"testId" = testId)
  print("Datasets prepared")
  return(output)
}
