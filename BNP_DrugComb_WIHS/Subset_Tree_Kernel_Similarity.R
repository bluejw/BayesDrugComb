RegimenReconstruct <- function(x, data){
  
  # reconstruct the regimen 
  # divide drugs into different classes and sort them alphabetically
  # 
  # args: x: regimen used by individual i at visit j
  #       data: the preprocessed data from WIHS dataset 
  # returns: regimen: regimen with drugs divided into different classes and sorted alphabetically in each classes  
  
  # construct a new regimen
  x <- as.vector(na.omit(x))
  regimen <- matrix(NA, nrow=data$num_classes, ncol=length(x))
  index <- rep(0, data$num_classes)
  
  # divide drugs into different classes
  for (i in 1:length(x)){
    the_class <- as.vector(data$abbr_group[x[i]][1,]) 
    # only if the class is in "Classes" vector
    if (the_class %in% data$Classes){
      the_class_index <- which(data$Classes == the_class) 
      index[the_class_index] <- index[the_class_index]+1 
      regimen[the_class_index,index[the_class_index]] <- x[i]
    }
  }
  
  # sort drugs in each class alphabetically
  for (j in 1:data$num_classes){
    if (index[j] != 0){
      regimen[j,1:index[j]] <- sort(regimen[j,1:index[j]]) 
    }
  }
  return(regimen)
}



ThirdLayerSimilarity <- function(x, y, data){
  
  # similarity in the third layer 
  # in this layer, their child nodes are leaves
  # we only need to count the number of same drugs
  # 
  # args: x, y: two regimens, 
  #       data: the preprocessed data from WIHS dataset 
  # returns: rho_third: similarities for each nodes
  #          simi_third: total similarities in third layer     
  
  rho_third <- rep(0, data$num_classes) 
  for (i in 1:data$num_classes){
    x_i <- as.vector(na.omit(x[i,]))
    y_i <- as.vector(na.omit(y[i,]))
    
    # count the number of same drugs in x_i and y_i
    if (length(x_i)>0 & length(y_i)>0){
      rho_third[i] <- sum(!is.na(match(x_i, y_i)))
    }
  }
  
  # similarity in third layer
  simi_third <- data$eta*sum(rho_third)
  return_list <- list(rho_third=rho_third, simi_third=simi_third)
  return(return_list)
}



SecondLayerSimilarity <- function(x, y, rho_third, data){
  
  # similarity in the second layer 
  # if the same parts of two trees contain the same number of elements
  # the corresponding nodes of two trees on the second layer share the same subset tree structure
  # 
  # args: x, y: two regimens
  #       rho_third: similarties for each nodes in third layer
  #       data: the preprocessed data from WIHS dataset 
  # returns: rho_second: similarities for each nodes
  #          simi_second: total similarities in second layer   
  
  rho_second <- rep(0, data$num_classes) 
  for (i in 1:data$num_classes){
    x_i <- as.vector(na.omit(x[i,]))
    y_i <- as.vector(na.omit(y[i,]))
    
    # only if both x_i and y_i contain the same drugs
    if (length(x_i)>0 & length(x_i)==length(y_i)){
      rho_second[i] <- data$eta*(1+data$eta)^(rho_third[i])
    }
  }
  
  simi_second <- sum(rho_second)
  return_list <- list(rho_second=rho_second, simi_second=simi_second)
  return(return_list)
}



FirstLayerSimilarity <- function(x, y, rho_second, data){
  
  # similarity in the first layer
  # for the root, only need to check whether the drugs contained in these two regimens belong to the same classes
  # 
  # args: x, y: two regimens, 
  #       rho_second: similarties for each nodes in second layer
  #       data: the preprocessed data from WIHS dataset 
  # returns: simi_first: similarities in first layer   
  
  simi_first <- 0
  index_x <- which(!is.na(x[,1]))
  index_y <- which(!is.na(y[,1]))
  # only if x and y contain the same classes
  if (identical(index_x, index_y)){
    simi_first <- data$eta
    for (i in index_x){
      simi_first <- simi_first*(1+rho_second[i])
    }
  }
  return(simi_first)
}



ZeroLayerSimilarity <- function(regimens_x, regimens_y, J_x, J_y, index_regimens_x, index_regimens_y, data){
  
  # similarity in the "zero" layer
  # need to check whether two subjects have the same regimen sequence
  # and then sum all the similarities between pairs of regimens in sequences
  # 
  # args: x, y: index for two subjects, 
  #       data: the preprocessed data from WIHS dataset 
  # returns: simi_zero: similarities in "zero" layer   
  
  simi_zero <- 0
  
  # only if subject x and y share the same regime sequence 
  if (identical(regimens_x, regimens_y)){
    simi_zero <- data$eta
    for (j in 1:J_x){
      # the j-th regimen in subject x's regimen sequence
      regimen <- RegimenReconstruct(regimens_x[j,], data) 
      rho_third <- ThirdLayerSimilarity(regimen, regimen, data)$rho_third
      rho_second <- SecondLayerSimilarity(regimen, regimen, rho_third, data)$rho_second
      simi_zero <- simi_zero*(1+FirstLayerSimilarity(regimen, regimen, rho_second, data))
    }
  }
  
  # sum all the similarities between pairs of regimens in sequences
  for (jx in 1:J_x){
    for (jy in 1:J_y){
      # similarities between regimens used by subject (x, y) at visit (jx, jy)
      simi_zero <- simi_zero+data$kappa[index_regimens_x[jx], index_regimens_y[jy]]
    }
  }
  return(simi_zero)
}