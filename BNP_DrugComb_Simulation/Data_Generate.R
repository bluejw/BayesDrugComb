library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(Matrix)
library(LaplacesDemon)

###################################### Load the Preprocessed Data ##################################

# load(file = "Simu.Data.Preprocess.Rdata") # preprocessed data from WIHS dataset
load(file = "/Users/BlueJ/Desktop/Simu.Data.Preprocess.Rdata")

n <- data_preprocess$n # number of patients
J <- data_preprocess$J # number of visits for each patient
Z <- data_preprocess$Z # HAART regimens for each patient at each visit
D <- data_preprocess$D # number of "knots" HAART regimens for kernel regression
kappa <- data_preprocess$kappa # similarity between HAART regimens
index_regimens <- data_preprocess$index_regimens # index of HAART regimens in similarity matrix 
index_kernel <- data_preprocess$index_kernel # index of "knots" HAART regimens in similarity matrix 
Kappa <- data_preprocess$Kappa # patient level similarity matrix 

###################################### Start to Generate Simulated Data ##################################

set.seed(123) # set seed for generating simulated data

Q <- 3 # number of depression items
S <- 3 # dimension of covariates X_{ij}

# generate covariates from normal distribution
X <- array(0, dim=c(n, max(J), S))
for (i in 1:n){
  X[i,1:J[i],1] <- rep(1, J[i]) # intercept term
  X[i,1:J[i],2] <- rep(rnorm(1,0,1), J[i]) # time-invariant covariate
  X[i,1:J[i],3] <- sort(rnorm(J[i],0,1)) # time-varying covariate
}

# generate parition pi_n from ddCRP
set.seed(123); pi_n_test <- Pi_n(n, 3.6, seq(1, n, by=1), Kappa) 
index_class_1 <- which(pi_n_test!=2 & pi_n_test!=3)
index_class_2 <- which(pi_n_test==2)
index_class_3 <- which(pi_n_test==3)
pi_n <- rep(1,n); pi_n[index_class_2] <- 2; pi_n[index_class_3] <- 3 
r_n <- length(unique(pi_n)) # number of classes r_n for n subjects 

# generate beta from normal distribution
beta_index_class_1 <- t(rbind(rmvn_rcpp(S,rep(0,Q),diag(1,Q))))
beta_index_class_2 <- t(rbind(rmvn_rcpp(S,rep(0,Q),diag(1,Q))))
beta_index_class_3 <- t(rbind(rmvn_rcpp(S,rep(0,Q),diag(1,Q))))
beta <- array(NA, dim=c(n, Q, S)) # beta_{i}: Q*S coefficient matrix
for (i in 1:n){
  if (i %in% index_class_1){ beta[i,,] <- beta_index_class_1
  }else if (i %in% index_class_2){ beta[i,,] <- beta_index_class_2
  }else { beta[i,,] <- beta_index_class_3 }
}

# generate gamma from normal distribution 
gamma_index_class_1 <- t(rbind(rmvn_rcpp(D,rep(0,Q),diag(1,Q))))
gamma_index_class_2 <- t(rbind(rmvn_rcpp(D,rep(0,Q),diag(1,Q))))
gamma_index_class_3 <- t(rbind(rmvn_rcpp(D,rep(0,Q),diag(1,Q))))
gamma <- array(0, dim=c(n, Q, D)) # gamma_{i}: Q*D coefficient matrix
for (i in 1:n){
  if (i %in% index_class_1){ gamma[i,,] <- gamma_index_class_1
  }else if (i %in% index_class_2){ gamma[i,,] <- gamma_index_class_2
  }else { gamma[i,,] <- gamma_index_class_3 }
}

# calculate the kernel weight matrix 
H <- array(NA, dim=c(n, max(J), D)) # H_{ij}: D-dimensional vector, whose d-th element k(Z_{ij},z_d)/sum_{d}(k(Z_{ij},z_d))
for (i in 1:n){ 
  for (j in 1:J[i]){
    current_regimen <- as.vector(na.omit(Z[i,j,]))
    if (length(current_regimen)!=0 & sum(kappa[index_regimens[i,j],index_kernel])!=0){
      # H_{ij} = k(Z_{ij},z_d)/sum_{d}(k(Z_{ij},z_d))
      H[i,j,1:D] <- kappa[index_regimens[i,j],index_kernel] / sum(kappa[index_regimens[i,j],index_kernel]) 
    }else{ 
      # H_{ij} = 0 for empty regimens and zero-similarity regimens
      H[i,j,1:D] <- 0 
    }
  }
}

# normalize covariates H for each dimension d
for (d in 1:D){
  H_hat <- H[,,d][!is.na(H[,,d])]
  H[,,d] <- (H[,,d] - mean(H_hat)) / sd(H_hat)
}

# calculate combination effects
h <- array(NA, dim=c(n, max(J), Q)) 
for (i in 1:n){
  for (j in 1:J[i]){
    # h_{ij} = gamma_{i}H_{ij} is a Q-dimensional vector
    h[i,j,] <- gamma[i,,] %*% H[i,j,] 
  }
}

# generate correlation terms, correlation matrix and variace of i.i.d errors 
sigma2 <- 1 # variance of i.i.d normal error
Sigma_omega <- matrix(1, nrow=Q, ncol=Q) # correlation matrix 
Sigma_omega[1,2] <- 0.25; Sigma_omega[2,1] <- 0.25
Sigma_omega[1,3] <- 0.5; Sigma_omega[3,1] <- 0.5
Sigma_omega[2,3] <- 0.75; Sigma_omega[3,2] <- 0.75
omega <- array(0, dim=c(n, max(J), Q)) # omega_{ij}: Q-dimensional correlation term
for (i in 1:n){ for (j in 1:J[i]){ omega[i,j,] <- rmvn_rcpp(1, rep(0, Q), sigma2*Sigma_omega) }}

# genereate depression scores
Y <- array(0, dim=c(n, max(J), Q)) # Y_{ij}: Q-dimensional continuous depression scores
for (i in 1:n){
  for (j in 1:J[i]){
    Y[i,j,] <- rmvn_rcpp(1, as.vector(beta[i,,]%*%X[i,j,]+h[i,j,]+omega[i,j,]), diag(sigma2, Q)) 
  }
}

###################################### Principal Component Analysis ####################################

index_start <- 1; index_end <- sum(J)
H_tilde <- matrix(NA, nrow=index_end, ncol=D) # sum(J)*D matrix 
for (i in 1:n){
  H_tilde[index_start:(index_start+J[i]-1),] <- H[i,1:J[i],]
  index_start <- index_start + J[i]
}

# singular value decomposition 
H_svd <- svd(H_tilde)
eigenval <- H_svd$d^2 # eigenvalues 
# sum(eigenval[1:39])/sum(eigenval); plot(eigenval)
eigenvec <- H_svd$v # eigenvectors
D_pca <- 39 # number of principle components 
H_tilde_proj <- H_tilde %*% eigenvec[,1:D_pca] 

index_start <- 1; index_end <- sum(J)
H_proj <- array(NA, dim=c(n, max(J), D_pca))
for (i in 1:n){
  H_proj[i,1:J[i],] <- H_tilde_proj[index_start:(index_start+J[i]-1),]
  index_start <- index_start + J[i]
}

# standardization 
D <- D_pca; H <- H_proj
for (d in 1:D){
  H_hat <- H[,,d][!is.na(H[,,d])]
  H[,,d] <- (H[,,d] - mean(H_hat)) / sd(H_hat)
}

###################################### Hyper-parameters and Pre-calculated data ##################################

# hyper-parameters
c <- 1; d <- 1 # c, d: m_0 \sim Gamma(c,d)
e_0 <- rep(0, S); E_0 <- diag(100, S) # e_0, E_0: e \sim N(e_0, E_0)
f_0 <- rep(0, D); F_0 <- diag(100, D) # f_0, F_0: f \sim N(f_0, F_0)
E_0_inv <- diag(0.01, S); e_0_star <- as.vector(E_0_inv %*% e_0) 
F_0_inv <- diag(0.01, D); f_0_star <- as.vector(F_0_inv %*% f_0) 
b_0 <- S+1; B_0 <- diag(100, S) # b_0, B_0: B \sim Inverse-Wishart(b_0, B_0^{-1})
lambda_0 <- D+1; Lambda_0 <- diag(100, D) # lambda_0, Lambda_0: Lambda \sim Inverse-Wishart(lambda_0, Lambda_0^{-1})
g_1 <- 1; g_2 <- 1 # sigma2 \sim Inverse-Gamma(g_1, g_2)

# pre-calculated data matrix for ease of computation
X_tilde <- array(0, dim=c(n, max(J), S+D))
for (i in 1:n){ X_tilde[i,1:J[i],] <- cbind(X[i,1:J[i],], H[i,1:J[i],]) }
X_tilde_sum <- array(0, dim=c(n, S+D, S+D))
for (i in 1:n){ X_tilde_sum[i,,] <- X_tilde_sum[i,,] + t(X_tilde[i,,]) %*% X_tilde[i,,] }
X_tilde_trans <- aperm(X_tilde, c(2,3,1))
X_tilde_sum_trans <- aperm(X_tilde_sum, c(2,3,1))


######################################## Save Simulation Truths ######################################

data <- NULL # simulation truths

# parameters 
data$pi_n <- pi_n
data$r_n <- r_n
data$index_class_1 <- index_class_1
data$index_class_2 <- index_class_2
data$index_class_3 <- index_class_3
data$beta1 <- beta_index_class_1
data$beta2 <- beta_index_class_2
data$beta3 <- beta_index_class_3
data$gamma1 <- gamma_index_class_1
data$gamma2 <- gamma_index_class_2
data$gamma3 <- gamma_index_class_3
data$Sigma_omega <- Sigma_omega
data$sigma2 <- sigma2

# data
data$Y <- Y
data$X <- X
data$H <- H
data$h <- h
data$omega <- omega
data$n <- n
data$Q <- Q
data$S <- S
data$D <- D
data$J <- J
data$Z <- Z
data$kappa <- kappa
data$index_kernel <- index_kernel
data$index_regimens <- index_regimens
data$Kappa <- Kappa

# others 
data$D_pca <- c(26,33,39,46,51) # number of principal components used in sensitivity analysis

# save(data, file="Simu.Truth.Rdata")
save(data, file="/Users/BlueJ/Desktop/Simu.Truths.Rdata")