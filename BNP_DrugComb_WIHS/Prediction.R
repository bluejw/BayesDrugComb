Similarity <- function(J_pred, regimen_pred, index_regimen_pred, data){
  
  # similarity matrix for a hypothetical individual 
  # 
  # args: J_pred: number of visits for the hypothetical individual 
  #       regimen_pred, index_regimen_pred: the ART regimens and their index used by the hypothetical individual 
  #       data: the preprocessed data from WIHS dataset 
  # returns: Kappa: 
  
  Kappa <- array(NA, dim=c(J_pred, data$n+1, data$n+1))
  
  for (j in 1:J_pred){
    if (j==1) { 
      regimen_pred_j <- matrix(regimen_pred[1:j,], j, byrow=T) 
    }else { regimen_pred_j <- regimen_pred[1:j,] }
    index_regimen_pred_j <- index_regimen_pred[1:j]
    
    Kappa_pred <- rep(NA, data$n+1) # the similarity between the hypothetical individual and others 
    Kappa_pred[data$n+1] <- ZeroLayerSimilarity(regimen_pred_j, regimen_pred_j, j, j, 
                                                index_regimen_pred_j, index_regimen_pred_j, data) # self-similarity 
    for (i in 1:data$n){ # (standardized) similarity between the hypothetical individual and others 
      Kappa_pred[i] <- ZeroLayerSimilarity(regimen_pred_j, data$Z[i,1:data$J[i],], j, data$J[i], index_regimen_pred_j, 
                                           data$index_regimens[i,], data) / sqrt(Kappa_pred[data$n+1] * data$Kappa[i,i]) }
    Kappa[j,1:data$n,1:data$n] <- data$Kappa # n*n similarity matrix is the same 
    Kappa[j,data$n+1,] <- Kappa_pred; Kappa[j,,data$n+1] <- Kappa_pred # update the n-th row and column of the matrix
  }
  return(Kappa)
}



Prediction <- function(Nit, J_pred, J_data, data, result, post, Kappa_pred, X_pred, H_pred){
  
  # prediction on depression scores for an individual 
  # 
  # args: Nit: number of itertaions for posterior samples 
  #       J_pred: number of visits for the individual in total
  #       J_data: number of visits as data
  #       data: preprocessed data 
  #       result: MCMC posterior samples 
  #       post: posterior estimation of partition
  #       Kappa_pred: similarity matrix for n+1 individuals at the visit j = 1,2,...,J 
  #       X_pred, H_pred: the design matrices for covariates and kernel weights for the individual 
  # returns: Y_pred: the predictive depression scores for an individual 
  
  Y_pred <- array(0, dim=c(Nit, J_pred, data$Q)) # predicitve depression scores 
  
  Sk <- matrix(NA, nrow=post$r_n, ncol=data$n) # individuals in each cluster
  for (k in 1:post$r_n){ Sk[k,1:length(which(post$pi_n == k))] <- which(post$pi_n == k) } 
  X_tilde_pred <- cbind(X_pred, H_pred) # combined covariates 
  
  for (nit in 1:Nit){
    
    print(nit)
    
    Y_pred[nit,1:J_data,1:data$Q] <- data$Y_pred 
    
    e_tilde <- matrix(NA, nrow=data$Q, ncol=data$S+data$D)
    B_tilde_inv <- array(NA, dim=c(data$S+data$D, data$S+data$D, (post$r_n+1)*data$Q))
    for (q in 1:data$Q){
      e_tilde[q,] <- c(result$e[nit,q,], result$f[nit,q,])
      for (k in 1:(post$r_n+1)){
        B_tilde_inv[,,data$B_index[k,]] <- chol2inv(chol(as.matrix(bdiag(result$B[nit,q,,], result$Lambda[nit,q,,]))))
      }
    }
    
    for (j in (J_data+1):J_pred){
      
      logpmf <- rep(NA, post$r_n+1) # log probability mass function for the partition 
      pi_n_pred <- rep(NA, data$n+1); pi_n_pred[1:data$n] <- post$pi_n
      
      for (k in 1:(post$r_n+1)){
     
        pi_n_pred[data$n+1] <- k # assign the individual to cluster 1,2,...r_n+1
        logpmf[k] <- logpmf_pi_n_rcpp(pi_n_pred, result$m_0[nit], c(result$sigma[nit,], data$n+1), Kappa_pred[j,,])
        if (k <= post$r_n){
          # if the predictive cluster exists
          logpmf[k] <- logpmf[k] + loglikelihood_rcpp(Sk[k,][!is.na(Sk[k,])]-1, t(e_tilde), B_tilde_inv[,,data$B_index[k,]], 
                                   X_tilde_pred[1:(j-1),], Y_pred[nit,1:(j-1),], data$X_tilde_sum_trans, data$Xy_tilde_sum_trans, 
                                   result$sigma2[nit], max(data$J), j-1, data$Q, data$S, data$D)
        }else{
          # if the predictive cluster does not exist 
          logpmf[k] <- logpmf[k] + loglikelihood0_rcpp(t(e_tilde), B_tilde_inv[,,data$B_index[k,]], X_tilde_pred[1:(j-1),], 
                                    Y_pred[nit,1:(j-1),], result$sigma2[nit], max(data$J), j-1, data$Q, data$S, data$D)
        }  
      }
      logpmf <- logpmf - max(logpmf)
      pmf <- exp(logpmf) / sum(exp(logpmf)) # probability mass function for partition
      
      for (k in 1:(post$r_n+1)){
        if (k<=post$r_n){ # if the predictive cluster exists
          Y_pred[nit,j,] <- Y_pred[nit,j,] + pmf[k]*(result$beta[nit,k,,] %*% X_pred[j,] + result$gamma[nit,k,,] %*% H_pred[j,])
        }else{ # if the predictive cluster does not exist 
          Y_pred[nit,j,] <- Y_pred[nit,j,] + pmf[k]*(result$e[nit,,] %*% X_pred[j,] + result$f[nit,,] %*% H_pred[j,]) 
        }
      }
    }
  }
  
  return(Y_pred)
}