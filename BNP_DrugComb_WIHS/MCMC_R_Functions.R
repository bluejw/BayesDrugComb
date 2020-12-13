Init <- function(){
  
  # MCMC initialization 
  # initialize parameters need to be estimated in MCMC
  # returns: init_list: pi_n, m_0, sigma, beta, gamma, e, B, f, Lambda, sigma2, omega, Sigma_omega
  
  # pi_n_init <- rep(1, n)
  pi_n_init <- sample(1:2, n, replace = T)
  m_0_init <- rgamma(1, c, d)
  sigma_init <- seq(1, n, by=1)
  
  beta_init <- rmvn_rcpp(Q, rep(0,S), diag(1,S))
  gamma_init <- rmvn_rcpp(Q, rep(0,D), diag(1,D))
  
  e_init <- matrix(0, nrow=Q, ncol=S)
  B_init <- array(NA, dim=c(Q, S, S))
  # for (q in 1:Q){ B_init[q,,] <- diag(1, S) }
  for (q in 1:Q){ B_init[q,,] <- diag(100, S) }
  
  f_init <- matrix(0, nrow=Q, ncol=D)
  Lambda_init <- array(NA, dim=c(Q, D, D))
  # for (q in 1:Q){ Lambda_init[q,,] <- diag(1, D) }
  for (q in 1:Q){ Lambda_init[q,,] <- diag(100, D) }
  
  sigma2_init <- 1
  Sigma_omega_init <- diag(1.001, Q)
  omega_init <- array(0, dim=c(n, max(J), Q))
  for (i in 1:n){ for (j in 1:J[i]){ omega_init[i,j,] <- rmvn_rcpp(1, rep(0, Q), sigma2_init * Sigma_omega_init) }}
  
  # return list of initialized parameters
  init_list <- list(pi_n=pi_n_init, m_0=m_0_init, sigma=sigma_init, beta=beta_init, gamma=gamma_init, 
                    e=e_init, B=B_init, f=f_init, Lambda=Lambda_init, sigma2=sigma2_init,
                    omega=omega_init, Sigma_omega=Sigma_omega_init)
  return(init_list)
}



Pi_n <- function(n, m_0, sigma, Kappa){
  
  # generate i.i.d partition sample pi_n from ddCRP prior
  # 
  # args: n: number of subjects, pi_n: partition, m_0: mass parameter,  sigma: permutation,
  #       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  # returns: pi_n: i.i.d partition sample pi_n from ddCRP prior
  
  pi_n <- rep(NA, n); pi_n[sigma[1]] <- 1
  
  # pmf of partition pi is defined iteratively for t > 1
  for (t in 2:n){
    
    previous_subjects <- sigma[1:(t-1)] # permutation at t-1: sigma_{1},...sigma_{t-1} 
    current_subject <- sigma[t] # current subject index at t-1: sigma_{t}
    previous_classes <- pi_n[previous_subjects] # previous partition at t-1: pi_n(sigma_{1},...sigma_{t-1})
    
    r_n <- length(unique(previous_classes))
    pmf <- rep(1, r_n+1); log_pmf <- rep(0, r_n+1)
    
    for (k in 1:(r_n+1)){
      
      current_class <- k # current class number at t: pi_n(sigma_{t})
      
      # if the class number of t-th subject 
      # exists in previous t-1 subjects' class numbers
      if (current_class %in% previous_classes){
        
        # previous t-1 subjects in the same class of t-th subject 
        previous_subjects_same_class <- previous_subjects[previous_classes==current_class] 
        # similarities between current subject and 
        # previous subjects in the same class or all previous subjects
        similarity_sum_same_class <- sum(Kappa[current_subject, previous_subjects_same_class])
        similarity_sum <- sum(Kappa[current_subject, previous_subjects])
        
        if (similarity_sum == 0){
          # if the similarities bewteen t-th subject and all previous subjects are 0
          pmf[k] <- pmf[k]*(t-1)/(m_0+t-1)*1/length(unique(previous_classes))
          log_pmf[k] <- log_pmf[k]+log(t-1)-log(m_0+t-1)-log(length(unique(previous_classes)))
        }else{
          if (similarity_sum_same_class == 0){
            # if the similarities bewteen t-th subject and all previous subjects in the same class are 0
            eps <- 1e-300 # avoid log_pmf being exactly -Inf
            similarity_sum_same_class <- similarity_sum_same_class+eps 
          }
          pmf[k] <- pmf[k]*(t-1)/(m_0+t-1)*similarity_sum_same_class/similarity_sum
          log_pmf[k] <- log_pmf[k]+log(t-1)-log(m_0+t-1)+log(similarity_sum_same_class)-log(similarity_sum)
        }
      }else{
        
        # if the class number of t-th subject is new 
        # i.e., it does not exist in previous t-1 subjects' class numbers
        pmf[k] <- pmf[k]*m_0/(m_0+t-1)
        log_pmf[k] <- log_pmf[k]+log(m_0)-log(m_0+t-1)
      }
    }
    
    log_prob_classes <- log_pmf - max(log_pmf)
    prob_classes <- exp(log_prob_classes) / sum(exp(log_prob_classes))   
    # class_number <- sample(1:(r_n+1), 1, replace = TRUE, prob = pmf) 
    class_number <- sample(1:(r_n+1), 1, replace = TRUE, prob = prob_classes) 
    pi_n[current_subject] <- class_number
  }
  return(pi_n)
}



PartitionReconstruct <- function(pi_n){
  
  # reconstruct the partition pi_n to a standard partition
  # for example, if the partition is (6,1,3,3,4)
  # then it will return a standard partition (4,1,2,2,3)
  # 
  # args: pi_n: any partition
  # returns: pi_standard: standard partition without "jumps" 
  #                 i.e., return (1,2,2,3) but not (1,3,3,6)
  #          old_labels, new_labels: previous and current labels
  
  pi_standard <- rep(NA, length(pi_n)) # standard partition
  classes <- unique(sort(pi_n)) # vector of unique class numbers 
  classes_change <- rank(classes) # "guidebook" to change unique class numbers
  
  for (i in 1:length(pi_n)){
    index_classes <- which(classes == pi_n[i]) # find position of class number of i-th subject 
    pi_standard[i] <- classes_change[index_classes] # assign new class number according to "guidebook"
  }
  returnlist = list(pi_n = pi_standard, old_labels = classes, new_labels = classes_change)
  return(returnlist)
}



Update_pi_n <- function(r_n_max, pi_n, m_0, sigma, Kappa, omega, e, B, f, Lambda, sigma2){
  
  # update the partition pi_n 
  # 
  # args: r_n_max: maximum number of clusters 
  #       pi_n: current partition, m_0: mass parameter, sigma: permutation,
  #       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  #       omega: (n*J)*Q-dimensional normal correlation vector
  #       e, B, f, Lambda: hyper-parameters for base measure G_0
  #       sigma2: variance of i.i.d normal error term
  # returns:  pi_n_update: new partition after update pi_n, r_n: number of clusters of new partition
  
  index_start <- 1
  B_index <- matrix(NA, nrow=r_n_max, ncol=Q)
  for (k in 1:r_n_max){
    index_end <- index_start + Q - 1
    B_index[k,] <- seq(index_start, index_end, by=1)
    index_start <- index_end + 1
  }
  
  e_tilde <- matrix(NA, nrow=Q, ncol=S+D)
  B_tilde_inv <- array(NA, dim=c(S+D, S+D, r_n_max*Q))
  
  for (q in 1:Q){
    e_tilde[q,] <- c(e[q,], f[q,])
    for (k in 1:r_n_max){
      B_tilde_inv[,,B_index[k,]] <- chol2inv(chol(as.matrix(bdiag(B[q,,], Lambda[q,,]))))
    }
  }
  
  Y_tilde <- Y - omega
  Xy_tilde_sum <- array(0, dim=c(n, Q, S+D))
  for (i in 1:n){ Xy_tilde_sum[i,,] <- t(Y_tilde[i,,]) %*% X_tilde[i,,] }
  Y_tilde_trans <- aperm(Y_tilde, c(2,3,1))
  Xy_tilde_sum_trans <- aperm(Xy_tilde_sum, c(2,3,1))
  
  # update the class number of each subjects
  pi_n_update <- update_pi_n_rcpp(pi_n, m_0, sigma, Kappa, t(e_tilde), B_tilde_inv, t(B_index)-1, X_tilde_trans, 
                                  Y_tilde_trans, X_tilde_sum_trans, Xy_tilde_sum_trans, n, r_n_max, sigma2, J, Q, S, D)
  
  return(pi_n_update)
}



Update_m_0 <- function(pi_n, m_0){
  
  # update the mass parameter m_0 in ddCRP
  # 
  # args: pi_n: current partition
  #       m_0: previous mass parameter
  # returns: m_0_update: new mass parameter m_0
  
  num <- length(pi_n) # number of data in partition pi_n
  r_n <- length(unique(pi_n)) # number of classes in partition pi_n
  
  tau_0 <- rbeta(1, m_0+1, num) # auxiliary variable tau_0
  lambda_0 <- (c+r_n-1)/(c+r_n-1+num*(d-log(tau_0))) # coefficient of mixture of two Gamma distributions
  
  # update m_0 by sampling from a mixture of two Gamma distributions
  if (runif(1) < lambda_0){
    m_0_update <- rgamma(1, c+r_n, d-log(tau_0))
  }else{ m_0_update <- rgamma(1, c+r_n-1, d-log(tau_0)) }
  return(m_0_update)
}



Update_sigma <- function(pi_n, m_0, sigma, Kappa){
  
  # update the permutation sigma
  # 
  # args: pi_n: current partition, m_0: mass parameter, sigma: permutation,
  #       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  # returns:  sigma_update: new permutation after update
  
  num <- length(pi_n) # number of data in partition pi_n 
  num_shuffle <- 10 # number of data to be shuffled 
  index_shuffle <- sample(1:num, num_shuffle, replace = F) # index of data to be shuffled 
  
  sigma_star <- sigma # sigma_star: proposed new permutation 
  sigma_star[index_shuffle] <- sample(sigma[index_shuffle], num_shuffle, replace = F) 
  
  # decide accept the new proposed permutation or not
  ratio <- logpmf_pi_n_rcpp(pi_n, m_0, sigma_star, Kappa) - logpmf_pi_n_rcpp(pi_n, m_0, sigma, Kappa)
  if (log(runif(1)) < ratio){
    sigma_update <- sigma_star
  }else{ sigma_update <- sigma }   
  return(sigma_update)
}



Update_beta <- function(pi_n, m_0, sigma, Kappa, gamma, e, B, omega, sigma2){
  
  # update values of parameters beta
  # 
  # args: pi_n: current partition, m_0: mass parameter, sigma: permutation,
  #       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  #       gamma:  unique values of parameters
  #       e, B: hyper-parameters for base measure G_0
  #       omega: (n*J)*Q-dimensional normal correlation vector
  #       sigma2: variance of i.i.d normal error term
  # returns:  beta_update: new beta after update
  
  r_n <- length(unique(pi_n)) # number of classes in partition pi_n
  if (length(dim(gamma))!=3) { gamma <- array(gamma, dim=c(1, Q, D))}
  beta_update <- array(NA, dim=c(r_n, Q, S))
  
  for (k in 1:r_n){
    index_class_k <- which(pi_n == k) # index of subjects in class k
    for (q in 1:Q){
      # data from prior distribution G_0   
      e_kq <- e[q,]; B_kq <- B[q,,]; B_kq_inv <- chol2inv(chol(B_kq))      
      # Update beta_{kq}, k = 1,2,...,r_n, q = 1,2,...,Q
      beta_update[k,q,] <- update_beta_rcpp(k, q, index_class_k, gamma, e_kq, B_kq_inv,
                                            omega, sigma2, Y, X, H, J, Q, S, D)
    }
  }
  return(beta_update)
}



Update_e <- function(r_n, beta, B){
  
  # update the mean of base measuere G_0 for beta
  # 
  # args: r_n: number of classes in partition pi_n
  #       beta: unique values of parameters
  #       B: covaraince matrix of base measuere G_0
  # returns: e_update: new e after update
  
  e_update <- matrix(NA, nrow=Q, ncol=S)
  if (length(dim(beta))!=3) { beta <- array(beta, dim=c(1, Q, S))}
  
  for (q in 1:Q){
    beta_sum <- rep(0, S)
    for (k in 1:r_n){ beta_sum <- beta_sum + beta[k,q,] }
    B_q_inv <- chol2inv(chol(B[q,,])) # inverse of covariance matrix B_q
    
    # Gaussian posterior distribution 
    V_n <- chol2inv(chol(E_0_inv + r_n * B_q_inv))
    mu_n <- V_n %*% (e_0_star + B_q_inv %*% beta_sum)
    e_update[q,] <- rmvn_rcpp(1, mu_n, V_n)
  }
  
  return(e_update)
}



Update_B <- function(r_n, beta, e){
  
  # update the covaraince matrix of base measuere G_0
  # 
  # args: r_n: number of classes in partition pi_n
  #       beta: unique values of parameters
  #       e: mean of base measuere G_0
  # returns: B_update: new B after update
  
  B_update <- array(NA, dim=c(Q, S, S))
  if (length(dim(beta))!=3) { beta <- array(beta, dim=c(1, Q, S))}
  
  for (q in 1:Q){
    beta_sum <- matrix(0, nrow=S, ncol=S)
    for (k in 1:r_n){ beta_sum <- beta_sum + (beta[k,q,]-e[q,]) %*% t(beta[k,q,]-e[q,]) }
    
    # Inverse Wishart posterior distribution 
    b_0_star <- b_0 + r_n; B_0_star <- B_0 + beta_sum
    B_update[q,,] <- riwish_rcpp(b_0_star, B_0_star) 
  }
  
  return(B_update)
}



Update_gamma <- function(pi_n, m_0, sigma, Kappa, beta, f, Lambda, omega, sigma2){
  
  # update values of parameters gamma
  # 
  # args: pi_n: current partition, m_0: mass parameter, sigma: permutation,
  #       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  #       beta:  unique values of parameters
  #       f, Lambda: hyper-parameters for base measure G_0
  #       omega: (n*J)*Q-dimensional normal correlation vector
  #       sigma2: variance of i.i.d normal error term
  # returns:  gamma_update: new gamma after update
  
  r_n <- length(unique(pi_n)) # number of classes in partition pi_n
  if (length(dim(beta))!=3) { beta <- array(beta, dim=c(1, Q, S))}
  gamma_update <- array(NA, dim=c(r_n, Q, D))
  
  for (k in 1:r_n){
    index_class_k <- which(pi_n == k) # index of subjects in class k
    for (q in 1:Q){
      # inverse of diagonal matrix Lambda 
      Lambda_inv <- chol2inv(chol(Lambda[q,,]))
      # data from likelihood 
      index_start <- 1; n_k <- sum(J[index_class_k])
      H_tilde <- matrix(NA, nrow=n_k, ncol=D) # sum(J)*D matrix used for covariance matrix V_n
      Y_tilde <- rep(NA, n_k) # sum(J)-dimension vector used for mean mu_n
      for (i in index_class_k){
        H_tilde[index_start:(index_start+J[i]-1),] <- H[i,1:J[i],]
        Y_tilde[index_start:(index_start+J[i]-1)] <- Y[i,1:J[i],q] - X[i,1:J[i],] %*% beta[k,q,] - omega[i,1:J[i],q] 
        index_start <- index_start + J[i]
      }
      
      # Update gamma_{kq}, k = 1,2,...,r_n, q = 1,2,...,Q
      # V_n <- chol2inv(chol(t(H_tilde) %*% H_tilde / sigma2 + Lambda_inv))
      V_n <- chol2inv(chol(mat_mult_rcpp(t(H_tilde), H_tilde) / sigma2 + Lambda_inv))
      # mu_n <- as.vector(V_n %*% (t(H_tilde) %*% Y_tilde / sigma2))
      mu_n <- as.vector(V_n %*% (crossprod(H_tilde, Y_tilde) / sigma2 + Lambda_inv %*% f[q,]))
      gamma_update[k,q,] <- rmvn_rcpp(1, mu_n, V_n)
    }
  }
  return(gamma_update)
}



Update_f <- function(r_n, gamma, Lambda){
  
  # update the mean of base measuere G_0 for gamma
  # 
  # args: r_n: number of classes in partition pi_n
  #       gamma: unique values of parameters
  #       Lambda: covaraince matrix of base measuere G_0
  # returns: f_update: new f after update
  
  f_update <- matrix(NA, nrow=Q, ncol=D)
  if (length(dim(gamma))!=3) { gamma <- array(gamma, dim=c(1, Q, D))}
  
  for (q in 1:Q){
    gamma_sum <- rep(0, D)
    for (k in 1:r_n){ gamma_sum <- gamma_sum + gamma[k,q,] }
    Lambda_q_inv <- chol2inv(chol(Lambda[q,,])) # inverse of covariance matrix  Lambda_q
    
    # Gaussian posterior distribution 
    V_n <- chol2inv(chol(F_0_inv + r_n *  Lambda_q_inv))
    mu_n <- V_n %*% (f_0_star +  Lambda_q_inv %*% gamma_sum)
    f_update[q,] <- rmvn_rcpp(1, mu_n, V_n)
  }
  
  return(f_update)
}



Update_Lambda <- function(r_n, gamma, f){
  
  # update the covaraince matrix of base measuere G_0
  # 
  # args: r_n: number of classes in partition pi_n
  #       gamma: unique values of parameters
  #       f: mean of base measuere G_0
  # returns:  Lambda_update: new Lambda after update
  
  Lambda_update <- array(NA, dim=c(Q, D, D))
  if (length(dim(gamma))!=3) { gamma <- array(gamma, dim=c(1, Q, D))}
  
  for (q in 1:Q){
    gamma_sum <- matrix(0, nrow=D, ncol=D)
    for (k in 1:r_n){ gamma_sum <- gamma_sum + (gamma[k,q,]-f[q,]) %*% t(gamma[k,q,]-f[q,]) }
    
    # Inverse Wishart posterior distribution 
    lambda_0_star <- lambda_0 + r_n; Lambda_0_star <- Lambda_0 + gamma_sum
    Lambda_update[q,,] <- riwish_rcpp(lambda_0_star, Lambda_0_star) 
  }
  
  return(Lambda_update)
}



Update_omega <- function(pi_n, beta, gamma, Sigma_omega, sigma2){
  
  # update normal correlation term omega
  # 
  # args: pi_n: current partition
  #       beta, gamma: unique values of parameters
  #       Sigma_omega: covariance matrix of omega_{ij}
  #       sigma2: variance of i.i.d normal error term
  # returns: omega_update: new omega after update
  
  omega_update <- array(NA, dim=c(n, max(J), Q))
  
  Sigma_omega_inv <- chol2inv(chol(sigma2 * Sigma_omega)) 
  if (length(dim(beta))!=3) { beta <- array(beta, dim=c(1, Q, S))}
  if (length(dim(gamma))!=3) { gamma <- array(gamma, dim=c(1, Q, D))}
  
  omega_update <- update_omega_rcpp(Y, pi_n, beta, gamma, Sigma_omega_inv, sigma2, n, X, H, J, Q, S, D)
  return(omega_update) 
}



Update_Sigma_omega <- function(omega, Sigma_omega, sigma2){
  
  # update covariance matrix of omega
  # 
  # args: omega: normal correlation term
  #       Sigma_omega: correlation matrix of omega
  #       sigma2: variance of i.i.d normal error term
  # returns: Sigma_omega_update: new Sigma_omega after update
  
  Sigma_omega_update <- update_sigma_omega_rcpp(omega, Sigma_omega, n, J, Q, sigma2)
  return(Sigma_omega_update)
}



Update_sigma2 <- function(pi_n, beta, gamma, omega){
  
  # update variance of i.i.d normal error term
  # 
  # args: pi_n: current partition
  #       beta, gamma: unique values of parameters
  #       omega: normal correlation term
  # returns: sigma2_update: new sigma2 after update
  
  sigma2_update <- rep(NA, 1)
  
  Y_tilde <- array(NA, dim=c(n, max(J), Q))
  if (length(dim(beta))!=3) { beta <- array(beta, dim=c(1, Q, S))}
  if (length(dim(gamma))!=3) { gamma <- array(gamma, dim=c(1, Q, D))}
  for (i in 1:n){
    class_index <- pi_n[i]
    beta_i <- beta[class_index,,]; gamma_i <- gamma[class_index,,]
    for (j in 1:J[i]){
      Y_tilde[i,j,] <- (Y[i,j,] - beta_i%*%X[i,j,] - gamma_i%*%H[i,j,] - omega[i,j,])^2
    }
  }
  
  # Inverse Gamma posterior distribution 
  g_1_star <- g_1 + Q*sum(J)/2
  g_2_star <- g_2 + sum(Y_tilde[!is.na(Y_tilde)])/2
  sigma2_update <- rinvgamma(1, g_1_star, g_2_star)
  return(sigma2_update)
}