
spec100 = function(prediction, real, setName,  k)
{
  specificityAt100 = list()
  for(i in 1:k)
  {
    paired = data.frame(prediction[i], real[i])
    colnames(paired) = c("Prediction", "Real")
    paired = paired[order(paired$Prediction),]
    first = 1
    for(j in 1:dim(paired)[1])
    {
      if(paired[j,2] == 1)
      {
        first = j
        break
      }
    }
    print(first)
    tn = length(which(paired[1:first,2] == 0))
    fp = length(which(paired[first:dim(paired)[1],2] == 0))
    specificityAt100[length(specificityAt100)+1] = tn/(tn + fp)
    
  }
  
  testPerformanceMatrix = matrix("0",nrow = 1, ncol = 3)
  rownames(testPerformanceMatrix) = c(paste0("Specificity at 100%")  )
  if(setName != ""){
    
    rownames(testPerformanceMatrix) = c(paste0(setName, " specificity at 100%")  )
  }
  colnames(testPerformanceMatrix) = c("Mean", "SD", "95% CI") 
  meanVal = mean(unlist(specificityAt100))
  sdVal = sd(unlist(specificityAt100))
  nVal = length(specificityAt100)
  error <- qt(0.975,df=nVal-1)*sdVal/sqrt(nVal)
  left <- meanVal-error
  right <- meanVal+error
  testPerformanceMatrix[1, 1] = round( meanVal,4)
  testPerformanceMatrix[1, 2] = round( sdVal,4)
  testPerformanceMatrix[1, 3] = paste0(max(0,round( left, 4)), " - ", min(round( right,4),1))
  
  return(testPerformanceMatrix)
} 


spec100test = function(prediction, real, setName,  k)
{
  specificityAt100 = list()
  
    paired = data.frame(prediction , real )
    colnames(paired) = c("Prediction", "Real")
    paired = paired[order(paired$Prediction),]
    first = 1
    for(j in 1:dim(paired)[1])
    {
      if(paired[j,2] == 1)
      {
        first = j
        break
      }
    }
    print(first)
    tn = length(which(paired[1:first,2] == 0))
    fp = length(which(paired[first:dim(paired)[1],2] == 0))
    specificityAt100[length(specificityAt100)+1] = tn/(tn + fp)
    
  
  testPerformanceMatrix = matrix("0",nrow = 1, ncol = 3)
  rownames(testPerformanceMatrix) = c(paste0("Specificity at 100%")  )
  if(setName != ""){
    
    rownames(testPerformanceMatrix) = c(paste0(setName, " specificity at 100%")  )
  }
  colnames(testPerformanceMatrix) = c("Mean", "SD", "95% CI") 
  meanVal =  unlist (specificityAt100)
  #sdVal = sd(unlist(specificityAt100))
 # nVal = length(specificityAt100)
  #error <- qt(0.975,df=nVal-1)*sdVal/sqrt(nVal)
  #left <- meanVal-error
  #right <- meanVal+error
  testPerformanceMatrix[1, 1] = round( meanVal,4)
  testPerformanceMatrix[1, 2] = "NA"#round( sdVal,4)
  testPerformanceMatrix[1, 3] = "NA"#paste0(max(0,round( left, 4)), " - ", min(round( right,4),1))
  
  return(testPerformanceMatrix)
} 
buildPerformanceMatrix = function(listConf, setName, k)
{
  testSensitivity = vector( length = k)
  testSpecificity = vector(  length = k)
  testAccuracy = vector(  length = k)
  testPrecision= vector(  length = k)
  testRecall = vector(  length = k)
  for(i in 1:k)
  {
    testSensitivity[i] = listConf[[i]]$byClass[2]
    testSpecificity[i] = listConf[[i]]$byClass[1]
    testAccuracy[i] = listConf[[i]]$byClass[11]
    tn = listConf[[i]]$table[1]
    fp = listConf[[i]]$table[2]
    fn = listConf[[i]]$table[3]
    tp = listConf[[i]]$table[4]
    testPrecision[i] = tp/(tp+fp)
    if(is.nan(testPrecision[i]))
    {
      testPrecision[i] = 0
    }
    testRecall[i] = tp/(tp+fn)
  }  
  TestPerformance = list(  testSensitivity, 
                           testSpecificity, testAccuracy,
                           testPrecision, testRecall)
  testPerformanceMatrix = matrix("0",nrow = 5, ncol = 3)
  rownames(testPerformanceMatrix) = c(paste0("Sensitivity"),
                                      paste0("Specificity"),
                                      paste0("Balanced accuracy"),
                                      paste0("Precision"),
                                      paste0("Recall") )
  if(setName != "")
  {
    rownames(testPerformanceMatrix) = c(paste0(setName," sensitivity"),
                                        paste0(setName," specificity"),
                                        paste0(setName," balanced accuracy") )
  } 
  
  colnames(testPerformanceMatrix) = c("Mean", "SD", "95% CI")
  for(i in 1:5)
  {
    meanVal = mean(TestPerformance[[i]])
    sdVal = sd(TestPerformance[[i]])
    nVal = length(TestPerformance[[i]])
    error <- qt(0.975,df=nVal-1)*sdVal/sqrt(nVal)
    left <- meanVal-error
    right <- meanVal+error
    testPerformanceMatrix[i, 1] = round( meanVal,4)
    testPerformanceMatrix[i, 2] = round( sdVal,4)
    testPerformanceMatrix[i, 3] = paste0(max(round( left, 4), 0),  " - ", min(1, round( right,4)))
    
  }
  return(testPerformanceMatrix)
}
buildPerformanceMatrixTest = function(listConf, setName, k)
{
  # testSensitivity = vector( length = k)
  # testSpecificity = vector(  length = k)
  # testAccuracy = vector(  length = k)
  # testPrecision= vector(  length = k)
  # testRecall = vector(  length = k)
 # for(i in 1:k)
 # {
    testSensitivity  = listConf$byClass[2]
    testSpecificity = listConf$byClass[1]
    testAccuracy  = listConf$byClass[11]
    tn = listConf$table[1]
    fp = listConf$table[2]
    fn = listConf$table[3]
    tp = listConf$table[4]
    testPrecision  = tp/(tp+fp)
    if(is.nan(testPrecision ))
    {
      testPrecision  = 0
    }
    testRecall  = tp/(tp+fn)
  #}  
  TestPerformance = list(  testSensitivity, 
                           testSpecificity, testAccuracy,
                           testPrecision, testRecall)
  testPerformanceMatrix = matrix("0",nrow = 5, ncol = 3)
  rownames(testPerformanceMatrix) = c(paste0("Sensitivity"),
                                      paste0("Specificity"),
                                      paste0("Balanced accuracy"),
                                      paste0("Precision"),
                                      paste0("Recall") )
  if(setName != "")
  {
    rownames(testPerformanceMatrix) = c(paste0(setName," sensitivity"),
                                        paste0(setName," specificity"),
                                        paste0(setName," balanced accuracy") )
  } 
  
  colnames(testPerformanceMatrix) = c("Mean", "SD", "95% CI")
  for(i in 1:5)
  {
    meanVal = mean(TestPerformance[[i]])
    sdVal = sd(TestPerformance[[i]])
    nVal = length(TestPerformance[[i]])
    #error <- qt(0.975,df=nVal-1)*sdVal/sqrt(nVal)
    #left <- meanVal-error
    #right <- meanVal+error
    testPerformanceMatrix[i, 1] = round( meanVal,4)
    testPerformanceMatrix[i, 2] = "NA"# round( sdVal,4)
    testPerformanceMatrix[i, 3] = "NA"# paste0(max(round( left, 4), 0),  " - ", min(1, round( right,4)))
    
  }
  return(testPerformanceMatrix)
}

generateResultsReport = function(dge, results, subtitle, reportPath, usedSubset , testId, k, seed)
{  
   
  if(max(unlist(results$trainPredictionResult)) == 2)
  {
    for(i in 1:k)
    {
      results$trainPredictionResult[[i]] = results$trainPredictionResult[[i]] -1
      results$testPredictionResult[[i]] = results$testPredictionResult[[i]] - 1
      results$validationPredictionResult[[i]] = results$validationPredictionResult[[i]] - 1
    }
  }
  
  library(cvAUC)    
  trainAUC = ci.cvAUC(predictions = 
                        (results$trainPredictionScore), labels = (results$trainRealClass), label.ordering = c(0,1))
  validateAUC = ci.cvAUC(predictions = 
                           (results$validationPredictionScore), labels = (results$validationRealClass), label.ordering = c(0,1))
  testAUC = ci.cvAUC(predictions = 
                       (unlist(results$testPredictionScore)), labels = (unlist(results$testRealClass)), label.ordering = c(0,1))# 
  print(testAUC)
  print(trainAUC)
  print(validateAUC)
  results$testPredictionScore = unlist(results$testPredictionScore)
  results$testRealClass = unlist(results$testRealClass)
  results$testPredictionResult = unlist(results$testPredictionResult)
  TestSpecificityAt100 = spec100test(results$testPredictionScore, results$testRealClass, "", k)
  TrainSpecificityAt100 = spec100(results$trainPredictionScore, results$trainRealClass, "",  k)
  ValidateSpecificityAt100 = spec100(results$validationPredictionScore, results$validationRealClass, "", k)
  print(TestSpecificityAt100)
  print(TrainSpecificityAt100)
  print(ValidateSpecificityAt100)
  
  results$testConfusionMatrices  = confusionMatrix(factor(results$testPredictionResult ), factor(results$testRealClass))
  testPerformance = buildPerformanceMatrixTest(listConf =  results$testConfusionMatrices, setName = "" )
  trainPerformance = buildPerformanceMatrix(listConf = results$trainingConfusionMatrices, setName = "", k)
  validatePerformance = buildPerformanceMatrix(listConf = results$validationConfusionMatrices, setName = "", k)
  testPerformance = rbind(testPerformance, TestSpecificityAt100)
  trainPerformance = rbind(trainPerformance, TrainSpecificityAt100)
  validatePerformance = rbind(validatePerformance, ValidateSpecificityAt100)
  print(testPerformance)
  print(trainPerformance)
  print(validatePerformance)
  allPerformance = rbind(c("", "Training set", ""), trainPerformance, 
                         c("", "Validation set", "") , validatePerformance,
                         c("", "Independent test set", "") , testPerformance)
  allPerformance
  allPerformance = rbind(c("", "Training set", ""), trainPerformance,
                         c(  round(trainAUC$cvAUC,4), round(trainAUC$se, 4),
                             paste0(max(round(trainAUC$ci[1],4),0), " - ", min(round(trainAUC$ci[2], 4) ,1))), 
                         c("", "Validation set", "") , validatePerformance, 
                         c(  round(validateAUC$cvAUC, 4), round(validateAUC$se, 4),
                             paste0(max(round(validateAUC$ci[1], 4),0),  " - ", min(round(validateAUC$ci[2], 4),1))),
                         c("", "Independent test set", "") , testPerformance,
                         c(  round(testAUC$cvAUC,4), round(testAUC$se,4),
                             paste0(max(round(testAUC$ci[1],4),0), " - ", min(round(testAUC$ci[2],4),1))))
  
  rownames(allPerformance) = c("",   rownames(trainPerformance),"cvAUC",
                               "",  rownames(validatePerformance),"cvAUC",
                               "",  rownames(testPerformance),"cvAUC")
  allPerformance
  
  write.table(allPerformance, paste0(reportPath, subtitle, Sys.Date(), "_.tsv"), quote = F, sep = "\t" )
  
  
  library(ROCit)
  library(ggplot2)  
  roc_bin_train <- rocit(score = unlist(results$trainPredictionScore), class = unlist(results$trainRealClass))
  pdf(paste0(reportPath, "ROC_", subtitle, Sys.Date(), ".pdf") , width = 21, height = 7)
  par(mfrow = c(1, 3), cex = 1.5)
  plot(roc_bin_train,YIndex = F,  values = F ,   cex=1.5)
  set.seed(200)
  ciROC_bin95_train <- ciROC(roc_bin_train,
                             level = 0.95, nboot = 200)
  lines(ciROC_bin95_train$TPR~ciROC_bin95_train$FPR,
        col = 2, lwd = 2) 
  title(paste0("Training set ROC \n", "Balanced accuracy = ",round(as.numeric(trainPerformance[3,1]) ,2), "   \n",
               "AUC = ", round(trainAUC$cvAUC,4) ))
  
  roc_bin_val <- rocit(score = unlist(results$validationPredictionScore), class = unlist(results$validationRealClass)) 
  plot(roc_bin_val,YIndex = F,  values = F,   cex=1.5)
  set.seed(200)
  ciROC_bin95_val <- ciROC(roc_bin_val,
                           level = 0.95, nboot = 200)
  lines(ciROC_bin95_val$TPR~ciROC_bin95_val$FPR,
        col = 2, lwd = 2) 
  title(paste0("Validation set ROC \n", "Balanced accuracy = ",round(as.numeric(validatePerformance[3,1]) ,2), "   \n",
               "AUC = ", round(validateAUC$cvAUC,4) ))
  roc_bin_test <- rocit(score = unlist(results$testPredictionScore), class = unlist(results$testRealClass))  
  plot(roc_bin_test,YIndex = F,  values = F,   cex=1.5)    
  set.seed(200)
  ciROC_bin95_test <- ciROC(roc_bin_test, 
                            level = 0.95, nboot = 200)  
  lines(ciROC_bin95_test$TPR~ciROC_bin95_test$FPR, 
        col = 2, lwd = 2) 
  title(paste0("Test set ROC \n", "Balanced accuracy = ",round(as.numeric(testPerformance[3,1]) ,2), "  \n",
               "AUC = ", round(testAUC$cvAUC,4) ))
  dev.off() 
  
  par(mfrow = c(1, 3), cex = 1.5)  
  pdf(paste0(reportPath, "Boxplot_", subtitle, "_", Sys.Date(), ".pdf"))
  
  reference = as.character(results$originalClass)
  if(diseaseLabel == "EC")
    reference[which(reference == "EC")] = "Endometrial cancer"
  
  reference[which(reference == "CTRL")] = "Benign control"
  reference[which(reference == "HC")] = "Healthy donor"
  prediction = data.frame(results$testPredictionScore)
  prediction = as.matrix(prediction)
  #prediction = sapply(1:dim(prediction)[1], function(i) {mean(prediction[i, ])})  
  if(diseaseLabel == "EC")
    reference = factor(reference, levels = c("Endometrial cancer", "Benign control", "Healthy donor"))
  #colors[which(colors == "CTRL")] = "green" 
  colorDescription = c()
  fillColor = c()
  colors = c((dge$samples$Description ))#, dge$samples$Description[otherSubset])
  stages  = c((dge$samples$Stage ))#,  dge$samples$Stage[otherSubset])
  colors[which(colors == "Healthy control")] = "blue"
  if(length(which(colors == "relapse")) > 0)
    colors[which(colors == "relapse")] = "deeppink3"
  if(length(which(colors == "Relapse")) > 0)
    colors[which(colors == "Relapse")] = "deeppink3"
  if(length(which(colors == "Ovarian cancer relapse")) > 0)
  { 
    colorDescription = c(colorDescription, "Ovarian cancer relapse")
    fillColor = c(fillColor,"deeppink3" )
    colors[which(colors == "Ovarian cancer relapse")] = "deeppink3"
  }
  colors[which(is.element(stages, c("IA", "IB", "IC")))] = "orange"
  colors[which(is.element(stages, c("II", "IIA", "IIB","IIC")))] = "pink"
  colors[which(is.element(stages, c("IIIA", "IIIB", "IIIC", "IIIC1", "IIIC2")))] = "red"
  colors[which(is.element(stages, c("IV", "IVB")))] = "black"
  if(length(which(is.element(stages, c("Missing"))))>0)
  { 
    colorDescription = c(colorDescription, "Missing")
    fillColor = c(fillColor,"lightpink" )
    colors[which(is.element(stages, c("Missing")))] = "lightpink"
    if(length(which(is.element(colors, c("NULL"))))>0)
       colors[which(is.element(colors, c("NULL")))] = "lightpink"
    
  }
  if(length(which(colors == "Serrous ovarian cancer")) > 0)
    colors[which(colors == "Serrous ovarian cancer")] = "lightpink" 
  if(length(which(colors == "borderline tumor")) > 0)
  {
    colorDescription = c(colorDescription, "Borderline tumor")
    fillColor = c(fillColor,"yellow" )
    colors[which(colors == "borderline tumor")] = "yellow" 
  }
  if(length(which(colors == "After resection, before chemotherapy")) > 0)
  {
    colorDescription = c(colorDescription, "After resection, before chemotherapy")
    fillColor = c(fillColor,"green" )
    colors[which(colors == "After resection, before chemotherapy")] = "green" 
  }
  if(length(which(colors == "Brenner tumor")) > 0) 
  { 
    colorDescription = c(colorDescription, "Brenner tumor")
    fillColor = c(fillColor,"darksalmon" )
    colors[which(colors == "Brenner tumor")] = "darksalmon" 
  }
  if(length(which(colors == "adenomyosis")) > 0)
  { 
    colorDescription = c(colorDescription, "Adenomyosis")
    fillColor = c(fillColor,"darkolivegreen1" )
    colors[which(colors == "adenomyosis")] = "darkolivegreen1"
  }
  if(length(which(colors == "Endometriosis")) > 0)
  { 
    colorDescription = c(colorDescription, "Endometriosis")
    fillColor = c(fillColor,"darkolivegreen2" )
    colors[which(colors == "Endometriosis")] = "darkolivegreen2"
  }
  if(length(which(colors == "myoma uteri")) > 0)
  { 
    colorDescription = c(colorDescription, "Myoma uteri")
    fillColor = c(fillColor,"violet" )
    colors[which(colors == "myoma uteri")] = "violet"
  } 
  if(length(which(colors == "cyst")) > 0)
  { 
    colorDescription = c(colorDescription, "Cyst")
    fillColor = c(fillColor,"lightgreen" )
    colors[which(colors == "cyst")] = "lightgreen"
  }
  if(length(which(colors == "malignant peritoneal pseudomyxoma high grade ")) > 0)
    colors[which(colors == "malignant peritoneal pseudomyxoma high grade ")] = "cyan"
  if(length(which(colors == "mature teratoma ")) > 0)
  { 
    colorDescription = c(colorDescription, "Mature teratoma")
    fillColor = c(fillColor,"brown" )
    colors[which(colors == "mature teratoma ")] = "brown"
  }
  if(length(which(colors == "turned out healthy")) > 0)
    colors[which(colors == "turned out healthy")] = "lightblue"
  if(length(which(colors == "complex atypical endometrial hyperplasia")) > 0)
  { 
    colorDescription = c(colorDescription, "Complex atypical endometrial hyperplasia")
    fillColor = c(fillColor,"mistyrose" )
    colors[which(colors == "complex atypical endometrial hyperplasia")] = "mistyrose"
  }
  if(length(which(colors == "polypus adenoleiomyomatosus")) > 0)
  { 
    colorDescription = c(colorDescription, "Polypus adenoleiomyomatosus")
    fillColor = c(fillColor,"lavender" )
    colors[which(colors == "polypus adenoleiomyomatosus")] = "lavender"
  }
  if(length(which(colors == "endometrial hyperplastic polyp")) > 0)
  { 
    colorDescription = c(colorDescription, "Endometrial hyperplastic polyp")
    fillColor = c(fillColor,"peachpuff4" )
    colors[which(colors == "endometrial hyperplastic polyp")] = "peachpuff4"
  }
  if(length(which(colors == "leiomyoma ")) > 0)
  { 
    colorDescription = c(colorDescription, "Leiomyoma")
    fillColor = c(fillColor,"thistle1" )
    colors[which(colors == "leiomyoma ")] = "thistle1"
  }
  if(length(which(colors == "normal tissue")) > 0)
  { 
    colorDescription = c(colorDescription, "normal tissue")
    fillColor = c(fillColor,"grey" )
    colors[which(colors == "normal tissue")] = "grey"
  }
  if(length(which(colors == "pseudomyxoma peritonei")) > 0)
  { 
    colorDescription = c(colorDescription, "Peudomyxoma peritonei")
    fillColor = c(fillColor,"gray71" )
    colors[which(colors == "pseudomyxoma peritonei")] = "gray71"
  }
  if(length(colors[which(colors == "serous adenocarcinoma")]) > 0)
    colors[which(colors == "serous adenocarcinoma")] = "lightpink"
  if(length(colors[which(colors == "serous carcinoma")]) > 0)
    colors[which(colors == "serous carcinoma")] = "lightpink"
  if(diseaseLabel == "EC")
    colorDescription =  c("Healthy donor","Endometrial cancer stage I","Endometrial cancer stage II",
                          "Endometrial cancer stage III", "Endometrial cancer stage IV",
                          "Endometrial cancer stage unknown",   colorDescription) 
  fillColor = c("blue", "yellow", "orange", "pink", "red", "black", fillColor) 
  
  if(diseaseLabel == "EC")
    boxplot(prediction[which(reference == "Endometrial cancer")]  
            ~ reference[which(reference == "Endometrial cancer")] , outpch=NA, main = "Independent test set", ylab = "TEPS score", xlab = "Group" , #list(paste(paste(colnames(boxplot_array)[k], sep = ""), sep = "\n")
            ylim = c(0,1), cex = 1, lex.order = T) 
  boxplot(prediction[which(reference == "Benign control")]  
          ~ reference[which(reference == "Benign control")] , add = T, outpch=NA, main = "Independent test set", ylab = "TEPS score", xlab = "Group" , #list(paste(paste(colnames(boxplot_array)[k], sep = ""), sep = "\n")
          ylim = c(0,1), cex = 1, lex.order = T) 
  boxplot(prediction[which(reference == "Healthy donor")]    
          ~ reference[which(reference == "Healthy donor")]   , add = T,outpch=NA, main = "Independent test set", ylab = "TEPS score", xlab = "Group" , #list(paste(paste(colnames(boxplot_array)[k], sep = ""), sep = "\n")
          ylim = c(0,1), cex = 1, lex.order = T) 
  
  
  
  if(diseaseLabel == "EC")
    stripchart(prediction[which(reference == "Endometrial cancer")]     ~reference[which(reference == "Endometrial cancer")] , vertical = TRUE, method = "jitter",
               pch = 21,  col =colors[which(reference == "Endometrial cancer")] , bg = colors[which(reference == "Endometrial cancer")] , 
               ylim = c(0,1), add = TRUE, #at = 1  , 
               ylab = "Prediction score", xlab = toPlot , cex = 1)
  stripchart(prediction[which(reference == "Benign control")]    ~reference[which(reference == "Benign control")], vertical = TRUE, method = "jitter",
             pch = 21,  col =colors[which(reference == "Benign control")], bg = colors[which(reference == "Benign control")], 
             ylim = c(0,1), add = TRUE,  ylab = "Prediction score", xlab = toPlot , cex = 1)
  stripchart(prediction[which(reference == "Healthy donor")]    ~reference[which(reference == "Healthy donor")], vertical = TRUE, method = "jitter",
             pch = 21,  col =colors[which(reference == "Healthy donor")], bg = colors[which(reference == "Healthy donor")], 
             ylim = c(0,1), add = TRUE,    ylab = "Prediction score", xlab = toPlot , cex = 1)
  par(xpd=TRUE)
  legend("topright", #inset=.32, 
         title="Group",
         colorDescription, 
         fill=fillColor, horiz=F, cex=0.4)
  
  abline(h = 0.5, lty = 2)
  dev.off()
  
  
  upperLeft = c()
  upperRight = c()
  lowerLeft = c()
  lowerRight = c()
  ConfusionToSave = c("", "Median (range)", "")
  lconf = results$trainingConfusionMatrices
  for(i in 1:k)
  {
    upperLeft = c(upperLeft, lconf[[i]]$table[4])
    upperRight = c(upperRight, lconf[[i]]$table[2])
    lowerLeft = c(lowerLeft, lconf[[i]]$table[3])
    lowerRight = c(lowerRight, lconf[[i]]$table[1])
  }
  paste0(median(upperLeft), " (", min(upperLeft), "-", max(upperLeft), ")")
  paste0(median(upperRight), " (", min(upperRight), "-", max(upperRight), ")")
  paste0(median(lowerLeft), " (", min(lowerLeft), "-", max(lowerLeft), ")")
  paste0(median(lowerRight), " (", min(lowerRight), "-", max(lowerRight), ")")
  
  
  ConfusionToSave = rbind(ConfusionToSave, c("", "Training set", ""),
                          c("Prediction", diseaseLabel, "Control"),
                          c(diseaseLabel, paste0(median(upperLeft), " (", min(upperLeft), "-", max(upperLeft), ")"),
                            paste0(median(upperRight), " (", min(upperRight), "-", max(upperRight), ")")),
                          c("Control", 
                            paste0(median(lowerLeft), " (", min(lowerLeft), "-", max(lowerLeft), ")"),
                            paste0(median(lowerRight), " (", min(lowerRight), "-", max(lowerRight), ")")))
  
  upperLeft = c()
  upperRight = c()
  lowerLeft = c()
  lowerRight = c()
  lconf = results$validationConfusionMatrices
  for(i in 1:k)
  {
    upperLeft = c(upperLeft, lconf[[i]]$table[4])
    upperRight = c(upperRight, lconf[[i]]$table[2])
    lowerLeft = c(lowerLeft, lconf[[i]]$table[3])
    lowerRight = c(lowerRight, lconf[[i]]$table[1])
  }
  ConfusionToSave = rbind(ConfusionToSave, c("", "Validation set", ""),
                          c("Prediction", diseaseLabel, "Control"),
                          c(diseaseLabel, paste0(median(upperLeft), " (", min(upperLeft), "-", max(upperLeft), ")"),
                            paste0(median(upperRight), " (", min(upperRight), "-", max(upperRight), ")")),
                          c("Control", 
                            paste0(median(lowerLeft), " (", min(lowerLeft), "-", max(lowerLeft), ")"),
                            paste0(median(lowerRight), " (", min(lowerRight), "-", max(lowerRight), ")")))
  
  upperLeft = c()
  upperRight = c()
  lowerLeft = c()
  lowerRight = c()
  
  lconf = results$testConfusionMatrices 
    upperLeft = c(upperLeft, lconf$table[4])
    upperRight = c(upperRight, lconf$table[2])
    lowerLeft = c(lowerLeft, lconf$table[3])
    lowerRight = c(lowerRight, lconf$table[1]) 
  ConfusionToSave = rbind(ConfusionToSave, c("", "Independent test set", ""),
                          c("Prediction", diseaseLabel, "Control"),
                          c(diseaseLabel, paste0(median(upperLeft), " (", min(upperLeft), "-", max(upperLeft), ")"),
                            paste0(median(upperRight), " (", min(upperRight), "-", max(upperRight), ")")),
                          c("Control", 
                            paste0(median(lowerLeft), " (", min(lowerLeft), "-", max(lowerLeft), ")"),
                            paste0(median(lowerRight), " (", min(lowerRight), "-", max(lowerRight), ")")))
  
  write.table(ConfusionToSave, paste0(reportPath, "Confusion_matrix_", subtitle, "_fold_",i, "_", Sys.Date(), ".tsv"), sep = "\t", col.names = F, row.names = F)  
  
 print("Report generation finished")
}
