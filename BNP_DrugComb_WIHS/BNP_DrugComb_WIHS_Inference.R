##########################################################################################
#                      Posterior Inference for WIHS Dataset                              #
##########################################################################################

library(Rcpp)
library(Matrix)
library(ggplot2)
library(gridExtra)

source("MCMC_R_Functions.R")
sourceCpp("MCMC_Rcpp_Functions.cpp")
source("BNP_DrugComb_MCMC.R")

load(file = "WIHS.Data.Preprocess.Rdata") # data: the preprocessed data from WIHS dataset

##########################################################################################
#                           MCMC Posterior Inference                                     #
##########################################################################################

load(file="WIHS.MCMC.Results.ddCRP.ST.Rdata") # result: MCMC posterior samples for WIHS dataset analysis

Clustering_Summary <- function(posterior_samples){
  nit <- dim(posterior_samples)[1] # number of itertaions
  n <- dim(posterior_samples)[2] # number of patients
  cl_num <- max(posterior_samples) # maximum number of clusters
  cl_prob <- matrix(0, nrow=n, ncol=cl_num) # posterior probability in each cluster
  for (cl in 1:cl_num){ cl_prob[,cl] <- colSums(posterior_samples == cl)/nit }
  pi_n_est <- rep(NA, n) # point estimation of posterior paritions 
  for (i in 1:n) { pi_n_est[i] <- which.max(cl_prob[i,]) }
  return(pi_n_est)
}

post <- NULL # posterior clustering results
post$pi_n <- Clustering_Summary(result$pi_n)
post$r_n <- length(unique(post$pi_n))

mean <- NULL # posterior means of parameters 
mean$beta <- colMeans(result$beta[,1:post$r_n,,],dims=1)
mean$gamma <- colMeans(result$gamma[,1:post$r_n,,],dims=1)

Cred_Interval <- function(post_samples){
  # 95% credible intervals for posterior samples of matrix 
  num_cluster <- dim(post_samples)[2]; num_row <- dim(post_samples)[3]; num_col <- dim(post_samples)[4]
  cred_interval <- array(NA, dim=c(num_cluster, num_row, num_col, 2)) # credible intervals for matrix 
  for (k in 1:num_cluster){ for (i in 1:num_row){ for (j in 1:num_col){ 
    cred_interval[k,i,j,] <- quantile(post_samples[,k,i,j], c(0.025, 0.975)) }}}
  return(cred_interval)
}

ci <- NULL # 95% credible intervals of parameters 
ci$beta <- Cred_Interval(result$beta[,1:post$r_n,,])
ci$gamma <- Cred_Interval(result$gamma[,1:post$r_n,,])

##########################################################################################
#         Figure 5: Posterior Means and 95% CIs for the Estimated Coefficients           #
##########################################################################################

s <- 2 # Age
df1 <- data.frame(beta = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$beta[,,s]),
                  lower = as.vector(ci$beta[,,s,1]), upper = as.vector(ci$beta[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p1 <- ggplot(df1, aes(beta, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("coef_esti_a.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p1
dev.off()

s <- 6 # CD4 count
df2 <- data.frame(beta = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$beta[,,s]),
                  lower = as.vector(ci$beta[,,s,1]), upper = as.vector(ci$beta[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p2 <- ggplot(df2, aes(beta, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("coef_esti_b.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p2
dev.off()

s <- 7 # Viral load
df3 <- data.frame(beta = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$beta[,,s]),
                  lower = as.vector(ci$beta[,,s,1]), upper = as.vector(ci$beta[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p3 <- ggplot(df3, aes(beta, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("coef_esti_c.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p3
dev.off()

s <- 10 # Substance use
df4 <- data.frame(beta = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$beta[,,s]),
                  lower = as.vector(ci$beta[,,s,1]), upper = as.vector(ci$beta[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p4 <- ggplot(df4, aes(beta, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("coef_esti_d.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p4
dev.off()

##########################################################################################
#              Figure 6: Estimated Coefficients for Combination Effects                  #
##########################################################################################

s <- 1 # PC1
df1 <- data.frame(gamma = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$gamma[,,s]),
                  lower = as.vector(ci$gamma[,,s,1]), upper = as.vector(ci$gamma[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p1 <- ggplot(df1, aes(gamma, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("drug_effect_esti_a.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p1
dev.off()

s <- 2 # PC2
df2 <- data.frame(gamma = rep(c("Cluster 1", "Cluster 2", "Cluster 3"), 4), value = as.vector(mean$gamma[,,s]),
                  lower = as.vector(ci$gamma[,,s,1]), upper = as.vector(ci$gamma[,,s,2]),
                  depression = factor(c(rep("Somatic",3), rep("Negative",3), rep("Positive",3), rep("Interpersonal",3)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p2 <- ggplot(df2, aes(gamma, value, colour = depression)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "Item") +
  theme(legend.title=element_text(size = 40, face='bold')) + theme(legend.text=element_text(size = 40)) + 
  theme(axis.title.x = element_text(size = 40, face='bold')) + theme(axis.title.y = element_text(size = 40, face='bold')) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) 
name <- paste("drug_effect_esti_b.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p2
dev.off()

##########################################################################################
#        Table 1: Top Five Positively and Negatively Related ART Regimens                #
#                    for the First Two Principal Components                              #
##########################################################################################

index <- c(87,86,85,84,83,1,2,3,4,5) # the index of top five positively and negatively related

# 1st principal components       
print(sort(data$eigenvec[,1])[index])
print(order(data$eigenvec[,1])[index])
print(data$z[data$index_kernel[order(data$eigenvec[,1])[index]],1:8])

# 2nd principal components       
print(sort(data$eigenvec[,2])[index])
print(order(data$eigenvec[,2])[index])
print(data$z[data$index_kernel[order(data$eigenvec[,2])[index]],1:8])

##########################################################################################
#  Figure 7: Predictive Depression Scores for an Individual with Hypothetical Scenarios  #
##########################################################################################

source("Subset_Tree_Kernel_Similarity.R")
source("Prediction.R")

Nit <- 500 # number of MCMC posterior iterations
J_pred <- 7 # number of visits in total
J_data <- 6 # number of visits used as data for prediction

# Design matrix for the covariates after standardization
X_pred <- matrix(NA, nrow=J_pred, ncol=data$S)
X_pred <- data$X_pred

############################################################
#                      Scenario #1                         #
############################################################

# ART regimens for scenario #1
regimen_pred <- matrix(NA, nrow=J_pred, ncol=8)
regimen_pred[1:J_data,] <- data$regimen_pred
index_regimen_pred <- rep(NA, J_pred)
index_regimen_pred[1:J_data] <- data$index_regimen_pred

regimen_pred[7,1:3] <- c("AZT", "LAM", "LPV")
index_regimen_pred[7] <- 79 # index of the regimen in all the ART regimens 

# Kernel weight matrix after performing principal component analysis
H_pred <- matrix(NA,nrow=J_pred, ncol=data$D)
H_pred[1:J_data,] <- data$H_pred
H_pred[7,] <- data$H[60,6,] # index_regimens[60,6] = 79

# Individual similarity matrix 
Kappa_pred <- Similarity(J_pred, regimen_pred, index_regimen_pred, data)
# Predictive depression scores for the hypothetical individual 
Y_pred <- Prediction(Nit, J_pred, J_data, data, result, post, Kappa_pred, X_pred, H_pred)
Y_pred_est <- colMeans(Y_pred, dims=1)

Cred_Interval <- function(post_samples){
  # 95% credible intervals for posterior samples of matrix 
  num_row <- dim(post_samples)[2]; num_col <- dim(post_samples)[3]
  cred_interval <- array(NA, dim=c(num_row, num_col, 2)) # credible intervals for matrix 
  for (i in 1:num_row){ for (j in 1:num_col){ 
    cred_interval[i,j,] <- quantile(post_samples[,i,j], c(0.025, 0.975)) }}
  return(cred_interval)
}
Y_pred_ci <- Cred_Interval(Y_pred); Y_pred_ci[Y_pred_ci<0] <- 0

# Plot ART regimens and prediction results
df11a <- data.frame(visit = factor(rep(1:J_pred, 5), levels=c("1", "2", "3", "4", "5", "6", "7")),
                    Class = rep(c("NRTI", "NNRTI", "PI", "EI", "INSTI"), each=J_pred),
                    num = c(c(1, 2, 2, 2, 2, 2, 2), rep(0, J_pred), c(0, 0, 0, 1, 1, 1, 1), rep(0, J_pred), rep(0, J_pred)))
df11b <- data.frame(drug_name1 = c(c("AZT", " ", " "," "), c("AZT", "LAM", " "," "),  c("AZT", "LAM", " ", " "),
                                   c("SQV", "AZT"," LAM"," "), c("SQV", "AZT"," LAM"," "), c("SQV", "AZT"," LAM"," "),
                                   c("LPV", "AZT"," LAM"," ")),
                    xpos1 = c(c(1, 0, 0, 0), c(2, 2, 0, 0), c(3, 3, 0, 0), c(4, 4, 4, 0),
                              c(5, 5, 5, 0), c(6, 6, 6, 0), c(7, 7, 7, 0)),
                    ypos1 = c(c(1, 0, 0, 0), c(1, 2, 0, 0), c(1, 2, 0, 0), c(1, 2, 3, 0), 
                              c(1, 2, 3, 0), c(1, 2, 3, 0), c(1, 2, 3, 0)))
p11 <- ggplot(data=df11a, aes(x=visit, y=num, fill=Class)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Paired") + 
  geom_text(data=df11b, aes(y=ypos1, x=xpos1, label=drug_name1), vjust=1.5, color="white", size=10, inherit.aes = FALSE) +
  scale_y_continuous(limits=c(0,4)) + xlab(" ") + ylab("Number of drugs used") + 
  theme(legend.title=element_text(size=30, face='bold')) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size = 30, face='bold')) + theme(axis.title.y = element_text(size = 30, face='bold')) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30))            
name <- paste("dep_pred_a.pdf")
pdf(name, width = 12, height = 6, onefile = TRUE)
p11
dev.off()

df12 <- data.frame(time = rep(1:J_pred,4), score = c(Y_pred_est[,1],Y_pred_est[,2],Y_pred_est[,3],Y_pred_est[,4]),
                  lower = c(Y_pred_ci[,1,1],Y_pred_ci[,2,1],Y_pred_ci[,3,1], Y_pred_ci[,4,1]),
                  upper = c(Y_pred_ci[,1,2],Y_pred_ci[,2,2],Y_pred_ci[,3,2], Y_pred_ci[,4,2]),
                  depression = factor(c(rep("Somatic",J_pred), rep("Negative",J_pred), rep("Positive",J_pred), rep("Interpersonal",J_pred)),
                                      levels = c("Somatic","Negative","Positive","Interpersonal")))
p12 <- ggplot(df12, aes(time, score, group=depression, colour=depression, shape=depression)) +
  # geom_ribbon(aes(time, ymin=lower, ymax=upper, group=depression), fill="darkgray", colour="darkgray") +
  geom_line(size = 1.5) + geom_point(size = 4) + 
  geom_line(aes(time, lower, group=depression, colour=depression), size = 1, linetype = "dashed") +
  geom_line(aes(time, upper, group=depression, colour=depression), size = 1, linetype = "dashed") +
  theme(axis.title.x = element_text(size = 30, face='bold')) + theme(axis.title.y = element_text(size = 30, face='bold')) + 
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) + 
  theme(legend.title=element_text(size = 30, face='bold')) + theme(legend.text=element_text(size = 30)) + 
  xlab(" ") + ylab("Depression score") + labs(colour="Item", shape="Item") + scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(labels = as.character(1:J_pred), breaks = 1:J_pred)
name <- paste("dep_pred_c.pdf")
pdf(name, width = 13, height = 6.5, onefile = TRUE)
p12
dev.off()

############################################################
#                      Scenario #2                         #
############################################################

# ART regimens for scenario #2
regimen_pred <- matrix(NA, nrow=J_pred, ncol=8)
regimen_pred[1:J_data,] <- data$regimen_pred
index_regimen_pred <- rep(NA, J_pred)
index_regimen_pred[1:J_data] <- data$index_regimen_pred

regimen_pred[7,1:3] <- c("EFV","FTC","TDF")
index_regimen_pred[7] <- 11 # index of the regimen in all the ART regimens 

# Kernel weight matrix after performing principal component analysis
H_pred <- matrix(NA,nrow=J_pred, ncol=data$D)
H_pred[1:J_data,] <- data$H_pred
H_pred[7,] <- data$H[3,4,] # index_regimens[3,4] = 11

# Individual similarity matrix 
Kappa_pred <- Similarity(J_pred, regimen_pred, index_regimen_pred, data)
# Predictive depression scores for the hypothetical individual 
Y_pred <- Prediction(Nit, J_pred, J_data, data, result, post, Kappa_pred, X_pred, H_pred)
Y_pred_est <- colMeans(Y_pred, dims=1)

Cred_Interval <- function(post_samples){
  # 95% credible intervals for posterior samples of matrix 
  num_row <- dim(post_samples)[2]; num_col <- dim(post_samples)[3]
  cred_interval <- array(NA, dim=c(num_row, num_col, 2)) # credible intervals for matrix 
  for (i in 1:num_row){ for (j in 1:num_col){ 
    cred_interval[i,j,] <- quantile(post_samples[,i,j], c(0.025, 0.975)) }}
  return(cred_interval)
}
Y_pred_ci <- Cred_Interval(Y_pred); Y_pred_ci[Y_pred_ci<0] <- 0


# Plot ART regimens and prediction results
df21a <- data.frame(visit = factor(rep(1:J_pred, 5), levels=c("1", "2", "3", "4", "5", "6", "7")),
                    Class = rep(c("NRTI", "NNRTI", "PI", "EI", "INSTI"), each=J_pred),
                    num = c(c(1, 2, 2, 2, 2, 2, 2), c(0, 0, 0, 0, 0, 0, 1), c(0, 0, 0, 1, 1, 1, 0), rep(0, J_pred), rep(0, J_pred)))
df21b <- data.frame(drug_name1 = c(c("AZT", " ", " "," "), c("AZT", "LAM", " "," "),  c("AZT", "LAM", " ", " "),
                                   c("SQV", "AZT"," LAM"," "), c("SQV", "AZT"," LAM"," "), c("SQV", "AZT"," LAM"," "), c("FTC", "TDF", "EFV", " ")),
                    xpos1 = c(c(1, 0, 0, 0), c(2, 2, 0, 0), c(3, 3, 0, 0), c(4, 4, 4, 0),
                              c(5, 5, 5, 0), c(6, 6, 6, 0), c(7, 7, 7, 0)),
                    ypos1 = c(c(1, 0, 0, 0), c(1, 2, 0, 0), c(1, 2, 0, 0), c(1, 2, 3, 0), 
                              c(1, 2, 3, 0), c(1, 2, 3, 0), c(1, 2, 3, 0)))
p21 <- ggplot(data=df21a, aes(x=visit, y=num, fill=Class)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Paired") + 
  geom_text(data=df21b, aes(y=ypos1, x=xpos1, label=drug_name1), vjust=1.5, color="white", size=10, inherit.aes = FALSE) +
  scale_y_continuous(limits=c(0,4)) + xlab(" ") + ylab("Number of drugs used") + 
  theme(legend.title=element_text(size=30, face='bold')) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size = 30, face='bold')) + theme(axis.title.y = element_text(size = 30, face='bold')) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30))            
name <- paste("dep_pred_b.pdf")
pdf(name, width = 12, height = 6, onefile = TRUE)
p21
dev.off()

df22 <- data.frame(time = rep(1:J_pred,4), score = c(Y_pred_est[,1],Y_pred_est[,2],Y_pred_est[,3],Y_pred_est[,4]),
                   lower = c(Y_pred_ci[,1,1],Y_pred_ci[,2,1],Y_pred_ci[,3,1], Y_pred_ci[,4,1]),
                   upper = c(Y_pred_ci[,1,2],Y_pred_ci[,2,2],Y_pred_ci[,3,2], Y_pred_ci[,4,2]),
                   depression = factor(c(rep("Somatic",J_pred), rep("Negative",J_pred), rep("Positive",J_pred), rep("Interpersonal",J_pred)),
                                       levels = c("Somatic","Negative","Positive","Interpersonal")))
p22 <- ggplot(df22, aes(time, score, group=depression, colour=depression, shape=depression)) +
  # geom_ribbon(aes(time, ymin=lower, ymax=upper, group=depression), fill="darkgray", colour="darkgray") +
  geom_line(size = 1.5) + geom_point(size = 4) + 
  geom_line(aes(time, lower, group=depression, colour=depression), size = 1, linetype = "dashed") +
  geom_line(aes(time, upper, group=depression, colour=depression), size = 1, linetype = "dashed") +
  theme(axis.title.x = element_text(size = 30, face='bold')) + theme(axis.title.y = element_text(size = 30, face='bold')) + 
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) + 
  theme(legend.title=element_text(size = 30, face='bold')) + theme(legend.text=element_text(size = 30)) + 
  xlab(" ") + ylab("Depression score") + labs(colour="Item", shape="Item") + scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(labels = as.character(1:J_pred), breaks = 1:J_pred)
name <- paste("dep_pred_d.pdf")
pdf(name, width = 13, height = 6.5, onefile = TRUE)
p22
dev.off()

##########################################################################################
#         Table S11: WAIC scores for five different methods in WIHS data analysis        #
##########################################################################################

Log_Likelihood_ST <- function(post, post_num){
  
  # calcluate the log-likelihood of the model
  # args: post: posterior samples
  # returns: logll: a n*post_num matrix of log-likelihood for n data points and post_num posterior samples
  
  logll <- matrix(0, nrow=n, ncol=post_num)
  
  # calculate the log-likelihood for each individual 
  for (nit in 1:post_num){
    print(nit)
    for (i in 1:n){
      cl <- post$pi_n[nit,i] # clustering membership for individual i
      for (j in 1:J[i]){
        # the mean and covariance matrix for the data
        mean <- post$beta[nit,cl,,] %*% data$X[i,j,] + post$gamma[nit,cl,,] %*% data$H[i,j,]
        var <- post$Sigma_omega[nit,,] + diag(post$sigma2[nit],Q,Q)
        logll[i,nit] <- logll[i,nit] + dmvn_rcpp(data$Y[i,j,], mean, var, logd=TRUE)
      }
    }
  }
  return(logll)
}

# calculate the WAIC for ddCRP+ST
load(file="WIHS.MCMC.Results.ddCRP.ST.Rdata")
logll_ddcrp_st <- Log_Likelihood_ST(result, post_num=500)
waic_ddcrp_st <- WAIC(logll_ddcrp_st)

# calculate the WAIC for Normal+ST
load(file="WIHS.MCMC.Results.Normal.ST.Rdata")
logll_normal_st <- Log_Likelihood_ST(result, post_num=500)
waic_normal_st <- WAIC(logll_normal_st)

# calculate the WAIC for DP+ST
load(file="WIHS.MCMC.Results.DP.ST.Rdata")
logll_dp_st <- Log_Likelihood_ST(result, post_num=500)
waic_dp_st <- WAIC(logll_dp_st)

Log_Likelihood_Linear <- function(post, post_num){
  
  # calcluate the log-likelihood of the model
  # args: post: posterior samples
  # returns: logll: a n*post_num matrix of log-likelihood for n data points and post_num posterior samples
  
  logll <- matrix(0, nrow=n, ncol=post_num)
  
  # calculate the log-likelihood for each individual 
  for (nit in 1:post_num){
    print(nit)
    for (i in 1:n){
      cl <- post$pi_n[nit,i] # clustering membership for individual i
      for (j in 1:J[i]){
        # the mean and covariance matrix for the data
        mean <- post$beta[nit,cl,,] %*% data$X[i,j,] + post$gamma[nit,cl,,] %*% data$H_linear[i,j,]
        var <- post$Sigma_omega[nit,,] + diag(post$sigma2[nit],Q,Q)
        logll[i,nit] <- logll[i,nit] + dmvn_rcpp(data$Y[i,j,], mean, var, logd=TRUE)
      }
    }
  }
  return(logll)
}

# calculate the WAIC for Normal+Lineaer
load(file="WIHS.MCMC.Results.Normal.Linear.Rdata")
logll_normal_linear <- Log_Likelihood_Linear(result, post_num=500)
waic_normal_linear <- WAIC(logll_normal_linear)

# calculate the WAIC for DP+Linear
load(file="WIHS.MCMC.Results.DP.Linear.Rdata")
logll_dp_linear <- Log_Likelihood_Linear(result, post_num=500)
waic_dp_linear <- WAIC(logll_dp_linear)

# WAIC Score Summary
waic_ddcrp_st$WAIC # 90378.26
waic_normal_st$WAIC # 100704.2
waic_dp_st$WAIC # 118164.1
waic_normal_linear$WAIC # 100258.4
waic_dp_linear$WAIC # 219565.6