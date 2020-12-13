BNP_DrugComb_MCMC <- function(Nit, burn.in, thin.fac, r_n_max){
  
  # MCMC main function for BNP_DrugComb model 
  # 
  # args: Nit: number of intertaions 
  #       burn.in: the burn-in iterations
  #       thin.fac: thinning factor for post burn-in samples 
  #       r_n_max: maximum number of clusters 
  # returns: MCMC posterior samples after burn-in with thinning factor 
  
  mcmc <- NULL
  mcmc$pi_n <- matrix(NA, nrow=Nit, ncol=n)
  mcmc$r_n <- rep(NA, Nit)
  mcmc$beta <- array(0, dim=c(Nit, r_n_max, Q, S))
  mcmc$gamma <- array(0, dim=c(Nit, r_n_max, Q, D))
  mcmc$Sigma_omega <- array(NA, dim=c(Nit, Q, Q))
  
  mcmc$m_0 <- rep(NA, Nit)
  mcmc$sigma <- matrix(NA, nrow=Nit, ncol=n)
  mcmc$e <- array(NA, dim=c(Nit, Q, S))
  mcmc$B <- array(NA, dim=c(Nit, Q, S, S))
  mcmc$f <- array(NA, dim=c(Nit, Q, D))
  mcmc$Lambda <- array(NA, dim=c(Nit, Q, D, D))
  mcmc$omega <- array(NA, dim=c(Nit, n, max(J), Q))
  mcmc$sigma2 <- rep(NA, Nit)
  
  # Initialize
  initial <- Init()
  mcmc$pi_n[1,] <- initial$pi_n
  mcmc$r_n[1] <- length(unique(mcmc$pi_n[1,]))
  mcmc$beta[1,1:mcmc$r_n[1],,] <- initial$beta
  mcmc$gamma[1,1:mcmc$r_n[1],,] <- initial$gamma
  mcmc$Sigma_omega[1,,] <- initial$Sigma_omega
  
  mcmc$m_0[1] <- initial$m_0
  mcmc$sigma[1,] <- initial$sigma
  mcmc$e[1,,] <- initial$e
  mcmc$B[1,,,] <- initial$B
  mcmc$f[1,,] <- initial$f
  mcmc$Lambda[1,,,] <- initial$Lambda
  mcmc$omega[1,,,] <- initial$omega
  mcmc$sigma2[1] <- initial$sigma2
  
  # Start of the chain
  start.time = proc.time()
  
  for (nit in 2:Nit){
    
    print(nit)
    
    # Update partition pi_n
    pi_n_update <- Update_pi_n(r_n_max, mcmc$pi_n[nit-1,], mcmc$m_0[nit-1], mcmc$sigma[nit-1,], Kappa, mcmc$omega[nit-1,,,], 
                      mcmc$e[nit-1,,], mcmc$B[nit-1,,,], mcmc$f[nit-1,,], mcmc$Lambda[nit-1,,,], mcmc$sigma2[nit-1])
    r_n_update <- max(pi_n_update)
    mcmc$pi_n[nit,] <- PartitionReconstruct(pi_n_update)$pi_n
    mcmc$r_n[nit] <- length(unique(mcmc$pi_n[nit,])) # number of classes in current partition pi_n
    
    print(sort(table(mcmc$pi_n[nit,]),decreasing = T)) # print current partition 
    
    # Update label switching of parameters
    if (r_n_update > mcmc$r_n[nit]){
      old_labels <- PartitionReconstruct(pi_n_update)$old_labels
      new_labels <- PartitionReconstruct(pi_n_update)$new_labels
      mcmc$beta[nit-1,new_labels,,] <- mcmc$beta[nit-1,old_labels,,]
      mcmc$gamma[nit-1,new_labels,,] <- mcmc$gamma[nit-1,old_labels,,]
    }
    
    # Update mass parameter m_0 and permutation sigma
    mcmc$m_0[nit] <- Update_m_0(mcmc$pi_n[nit,], mcmc$m_0[nit-1])
    mcmc$sigma[nit,] <- Update_sigma(mcmc$pi_n[nit,], mcmc$m_0[nit], mcmc$sigma[nit-1,], Kappa)
    
    # Update unique values of parameters beta 
    mcmc$beta[nit,1:mcmc$r_n[nit],,] <- Update_beta(mcmc$pi_n[nit,], mcmc$m_0[nit], mcmc$sigma[nit,], Kappa, 
         mcmc$gamma[nit-1,1:mcmc$r_n[nit],,], mcmc$e[nit-1,,], mcmc$B[nit-1,,,], mcmc$omega[nit-1,,,], mcmc$sigma2[nit-1])
    mcmc$e[nit,,] <- Update_e(mcmc$r_n[nit], mcmc$beta[nit,1:mcmc$r_n[nit],,], mcmc$B[nit-1,,,])
    mcmc$B[nit,,,] <- Update_B(mcmc$r_n[nit], mcmc$beta[nit,1:mcmc$r_n[nit],,], mcmc$e[nit,,])
   
    # Update unique values of parameters gamma
    mcmc$gamma[nit,1:mcmc$r_n[nit],,] <- Update_gamma(mcmc$pi_n[nit,], mcmc$m_0[nit], mcmc$sigma[nit,], Kappa, 
       mcmc$beta[nit,1:mcmc$r_n[nit],,], mcmc$f[nit-1,,], mcmc$Lambda[nit-1,,,], mcmc$omega[nit-1,,,], mcmc$sigma2[nit-1])
    mcmc$f[nit,,] <- Update_f(mcmc$r_n[nit], mcmc$gamma[nit,1:mcmc$r_n[nit],,], mcmc$Lambda[nit-1,,,])
    mcmc$Lambda[nit,,,] <- Update_Lambda(mcmc$r_n[nit], mcmc$gamma[nit,1:mcmc$r_n[nit],,], mcmc$f[nit,,])
    
    # Update normal correlation term omega and correlation matrix Sigma_omega
    mcmc$omega[nit,,,] <- Update_omega(mcmc$pi_n[nit,], mcmc$beta[nit,1:mcmc$r_n[nit],,], 
                          mcmc$gamma[nit,1:mcmc$r_n[nit],,], mcmc$Sigma_omega[nit-1,,], mcmc$sigma2[nit-1])
    mcmc$Sigma_omega[nit,,] <- Update_Sigma_omega(mcmc$omega[nit,,,], mcmc$Sigma_omega[nit-1,,], mcmc$sigma2[nit-1])
    
    # Update variance of i.i.d normal error 
    mcmc$sigma2[nit] <- Update_sigma2(mcmc$pi_n[nit,], mcmc$beta[nit,1:mcmc$r_n[nit],,], 
                                      mcmc$gamma[nit,1:mcmc$r_n[nit],,], mcmc$omega[nit,,,])
  }
  
  duration = proc.time()-start.time
  print(duration) # print the running time 
  
  # MCMC results 
  result <- NULL
  sample_index <- seq(burn.in+1, Nit, by=thin.fac)
  
  result$pi_n <- mcmc$pi_n[sample_index,] 
  result$r_n <- mcmc$r_n[sample_index]
  result$beta <- mcmc$beta[sample_index,,,]
  result$gamma <- mcmc$gamma[sample_index,,,]
  result$Sigma_omega <- mcmc$Sigma_omega[sample_index,,]
  
  result$m_0 <- mcmc$m_0[sample_index]
  result$sigma <- mcmc$sigma[sample_index,]
  result$e <- mcmc$e[sample_index,,]
  result$B <- mcmc$B[sample_index,,,]
  result$f <- mcmc$f[sample_index,,]
  result$Lambda <- mcmc$Lambda[sample_index,,,]
  result$omega <- mcmc$omega[sample_index,,,]
  result$sigma2 <- mcmc$sigma2[sample_index]
  
  return(result)
}