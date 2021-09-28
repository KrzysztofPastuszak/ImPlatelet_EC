

#' Performs k-fold cross-valdidation
#' returns:
#' result - list of performance metrics, classification outputs and real classes
#' input:
#' model - keras CNN model
#' dge
#' dataSets - prepared with function prepareDataSets
#' pathToWeights - path to initial weights file
#' iters - number of iteriations
#' k - k used in k-fold cv
#' seed - used for generating k-fold splits

runCrossValidation = function(model, dge, dataSets, pathToWeights, iters, k, seed)
{
  # print(dim(usedSubset$testX))
  # print(length(usedSubset$testY))
  results = list()
  results$trainingConfusionMatrices = list() 
  results$testConfusionMatrices = list()
  results$testScore = list()
  results$validationConfusionMatrices  = list()
  lscoreValidate2 = list()
  
  results$trainPredictionResult =  list()
  results$trainPredictionScore =  list()
  results$trainRealClass  =  list()
  results$validationRealClass  =  list()
  results$testPredictionResult  =  list()
  results$testPredictionScore  =  list() 
  results$testRealClass  =  list()
  results$validationPredictionResult =     list()
  results$validationPredictionScore =     list()
  dataSets$allY = as.numeric(dataSets$allY)
  results$originalClass = dge$samples$OriginalGroup[usedSubset]
  for(i in 1:length(dataSets$usedSubset))
  {
    k = 5
    set.seed(seed)
    folds = createFolds(y = as.factor(dataSets$allY), k = k)
    trainFold = unlist(folds[3:5])
    valFold = (unlist(folds[1:2,drop = F]))
    trainFoldY = dataSets$allY[trainFold]
    trainFoldX = dataSets$allX[trainFold, , , drop = F]
    valFoldY = dataSets$allY[valFold,drop = F]
    valFoldX = dataSets$allX[valFold, , , drop = F]
    j = i
    testY = dataSets$allY[j,drop = F]
    testX = dataSets$allX[j, , ,drop = F]
    if(is.element(j, trainFold))
      trainFold = trainFold[!trainFold==j]
    if(is.element(j, valFold))
      valFold = valFold[!valFold==j]
    model %>% load_model_weights_hdf5(filepath = pathToWeights)
    trainWeights = trainFoldY
    trainWeights[which(trainFoldY == 1)] =  1.06*length(which(trainFoldY == 0.0))/length(trainFoldY)
    trainWeights[which(trainFoldY == 0)] =  length(which(trainFoldY == 1.0))/length(trainFoldY) 
    model %>% fit( trainFoldX, trainFoldY, validation_data = list(valFoldX, valFoldY) , epochs =  iters,
                    batch_size = 15, shuffle =  T,  sample_weight = trainWeights 
    )
    
    
     
    trainPredictionClass = model %>% predict_classes(trainFoldX)
    trainPredictionScore = model %>% predict(trainFoldX) 
    trainConfusionMatrix = confusionMatrix(factor(trainPredictionClass[, 1]), factor(trainFoldY)) 
    testScore <- model %>% evaluate( testX, testY)
    testPredictionClass = model %>% predict_classes( testX)
    testPredictionScore = model %>% predict(testX) 
    #testConfusionMatrix = confusionMatrix(factor(testPredictionClass[, 1]), factor(dataSets$testY))  
    validationPredictionClass = model %>% predict_classes(valFoldX)
    validationPredictionScore = model %>% predict(valFoldX) 
    validationConfusionMatrix = confusionMatrix(factor(validationPredictionClass[, 1]), factor(valFoldY))
    print(trainConfusionMatrix$table )
   # print(testConfusionMatrix$table )
    
    
    results$trainingConfusionMatrices[[length(results$trainingConfusionMatrices)+1]] =   trainConfusionMatrix 
   # results$testConfusionMatrices[[length(results$testConfusionMatrices)+1]] =   testConfusionMatrix
    results$testScore[[length(results$testScore)+1]] =   testScore
    results$validationConfusionMatrices[[length(results$validationConfusionMatrices)+1]] =   validationConfusionMatrix 
    results$trainPredictionResult[[length(results$trainPredictionResult)+1]] =   trainPredictionClass
    results$trainPredictionScore[[length(results$trainPredictionScore)+1]] =   trainPredictionScore
    results$trainRealClass[[length(results$trainRealClass)+1]] =   trainFoldY
    results$testPredictionResult[[length(results$testPredictionResult)+1]] =   testPredictionClass
    results$testPredictionScore[[length(results$testPredictionScore)+1]] =   testPredictionScore 
    results$testRealClass[[length(results$testRealClass)+1]] =   testY
    results$validationPredictionResult[[length(results$validationPredictionResult)+1]] =   validationPredictionClass
    results$validationPredictionScore[[length(results$validationPredictionScore)+1]] =   validationPredictionScore
    results$validationRealClass[[length(results$validationRealClass)+1]] =   valFoldY
    print(paste("Finished cross-validation for fold: ", i))
  }
  return(results)
}
