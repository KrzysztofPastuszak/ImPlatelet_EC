 

#' Prepares keras CNN model
#' returns:
#' model
#' input:
#' path  - path to the directory in which matrices and images should beheld
#' weightsPath - path in which weights should be stored
#' weightsFile - optional, name of the file with initial weights, default modelInitialWeights.hdf5
#' picWidth
#' picHeight

 
buildModel = function(path, weightsPath = "weights/", weightsFile = "modelInitialWeights.hdf5", picWidth = 345, picHeight = 243)
{
  library(caret)
  library(keras)
  gpu <- tensorflow::tf$config$experimental$get_visible_devices('GPU')[[1]]
  tensorflow::tf$config$experimental$set_memory_growth(device = gpu, enable = TRUE)
  picWidth = 345
  picHeight = 243
  #set.seed(seed)
  model <- keras_model_sequential() 
  model %>%
    layer_reshape(c( picWidth ,picHeight, 1), input_shape = c(picWidth,picHeight))%>%
    layer_conv_2d(kernel_size = c(3,3), 
                  activation = "tanh", filters = 4, padding = "same") %>%
    layer_conv_2d(kernel_size = c(3,3), 
                  activation = "tanh", filters = 4, padding = "same") %>% 
    layer_flatten() %>%
    
    layer_dense(units = 256, activation = 'relu') %>%
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate = 0.15) %>%
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate = 0.15) %>% 
    layer_dense(units = 64, activation = 'relu') %>%
    
    layer_dense(units = 1, activation = 'sigmoid')
   
  model %>% compile(
    optimizer = 'AdaDelta', 
    loss = loss_binary_crossentropy,
    metrics = metric_binary_accuracy
  )
  
  model %>% save_model_weights_hdf5(filepath = paste0(path, weightsPath, weightsFile))
  print("Model prepared")
  return(model)
}