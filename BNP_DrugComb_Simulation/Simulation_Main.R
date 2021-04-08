##########################################################################################
#                  Posterior Inference for Simulation Study                              #
##########################################################################################

library(Rcpp)
library(ggplot2)
library(gridExtra)
library(lattice)

source("MCMC_R_Functions.R")
sourceCpp("MCMC_Rcpp_Functions.cpp") # warnings may be generated when include RcppEigen, but it will not affect the usage
source("Data_Generate.R") # data: the generated simulation truths

##########################################################################################
#                           Run MCMC on Simulated Data                                   #
##########################################################################################

source("BNP_DrugComb_MCMC.R")

Nit <- 10000 # number of total MCMC intertaions 
burn.in <- 5000 # the burn-in iterations
thin.fac <- 10 # thinning factor for post burn-in samples 
r_n_max <- 10 # maximum number of clusters 
# Run MCMC (it takes about 1.5 hours for 10000 iterations)
MCMC_Result <- BNP_DrugComb_MCMC(Nit, burn.in, thin.fac, r_n_max)

##########################################################################################
#                           MCMC Posterior Inference                                     #
##########################################################################################

load(file="Simu_MCMC_Results.Rdata") # result: MCMC posterior samples for one simulated dataset

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
#                     Figure 4: Combination Effects Estimation                           #
##########################################################################################

load(file="Simu.Combination.Effects.Rdata") # estimations of combination effects under five different methods

i <- 2; q <- 2 # index of first patient and depression item
visit <- which(data$index_regimens[i,]!=1466) # delete the empty regimens 
df1 <- data.frame(x = rep(1:length(visit),6), 
                  h = c(data$h[i,visit,q], 
                        colMeans(h_est$h_est_ddcrp_st[,i,visit,q],dims=1), 
                        colMeans(h_est$h_est_normal_st[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_normal_linear[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_dp_st[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_dp_linear[,i,visit,q],dims=1)),
                  lower = rep(apply(h_est$h_est_ddcrp_st[,i,visit,q], 2, quantile, prob=c(0.025)), 6),
                  upper = rep(apply(h_est$h_est_ddcrp_st[,i,visit,q], 2, quantile, prob=c(0.975)), 6),
                  est = factor(c(rep("Truths", length(visit)), 
                                 rep("ddCRP+ST", length(visit)),
                                 rep("Normal+ST", length(visit)),
                                 rep("Normal+Linear", length(visit)),
                                 rep("DP+ST", length(visit)),
                                 rep("DP+Linear",length(visit))), 
                               levels = c("Truths", "ddCRP+ST", "Normal+ST", "Normal+Linear", "DP+ST", "DP+Linear")))
p1 <- ggplot(data=df1, aes(x, h, colour=est, shape=est)) + 
  geom_ribbon(aes(x, ymin=lower, ymax=upper), fill="darkgray", colour="darkgray") + geom_line(size=1.25) + geom_point(size=6) +
  theme(axis.title.x = element_text(size = 50, face="bold")) + theme(axis.title.y = element_text(size = 50, face="bold")) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) + 
  theme(legend.title=element_text(size = 42.5)) + theme(legend.text=element_text(size = 42.5)) + 
  xlab("Visit") + ylab("Combination effect") + labs(col=" ", shape=" ") +  
  scale_x_continuous(labels = as.character(1:length(visit)), breaks = 1:length(visit)) +
  theme(legend.position = "top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(breaks = c("Truths", "ddCRP+ST", "Normal+ST", "Normal+Linear", "DP+ST", "DP+Linear"), 
                     values = c("black", "red", "orange", "green", "purple", "blue"))
name <- paste("comb_effect_a.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p1
dev.off()

i <- 29; q <- 2 # index of second patient and depression item
visit <- which(data$index_regimens[i,]!=1466) # delete the empty regimens 
df2 <- data.frame(x = rep(1:length(visit),6), 
                  h = c(data$h[i,visit,q], 
                        colMeans(h_est$h_est_ddcrp_st[,i,visit,q],dims=1), 
                        colMeans(h_est$h_est_normal_st[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_normal_linear[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_dp_st[,i,visit,q],dims=1),
                        colMeans(h_est$h_est_dp_linear[,i,visit,q],dims=1)),
                  lower = rep(apply(h_est$h_est_ddcrp_st[,i,visit,q], 2, quantile, prob=c(0.025)), 6),
                  upper = rep(apply(h_est$h_est_ddcrp_st[,i,visit,q], 2, quantile, prob=c(0.975)), 6),
                  est = factor(c(rep("Truths", length(visit)), 
                                 rep("ddCRP+ST", length(visit)),
                                 rep("Normal+ST", length(visit)),
                                 rep("Normal+Linear", length(visit)),
                                 rep("DP+ST", length(visit)),
                                 rep("DP+Linear",length(visit))), 
                               levels = c("Truths", "ddCRP+ST", "Normal+ST", "Normal+Linear", "DP+ST", "DP+Linear")))
p2 <- ggplot(data=df2, aes(x, h, colour=est, shape=est)) + 
  geom_ribbon(aes(x, ymin=lower, ymax=upper), fill="darkgray", colour="darkgray") + geom_line(size=1.25) + geom_point(size=6) +
  theme(axis.title.x = element_text(size = 50, face="bold")) + theme(axis.title.y = element_text(size = 50, face="bold")) +
  theme(axis.text.x = element_text(size = 40)) + theme(axis.text.y = element_text(size = 40)) + 
  theme(legend.title=element_text(size = 42.5)) + theme(legend.text=element_text(size = 42.5)) + 
  xlab("Visit") + ylab("Combination effect") + labs(col=" ", shape=" ") +  
  scale_x_continuous(labels = as.character(1:length(visit)), breaks = 1:length(visit)) +
  theme(legend.position = "top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(breaks = c("Truths", "ddCRP+ST", "Normal+ST", "Normal+Linear", "DP+ST", "DP+Linear"), 
                     values = c("black", "red", "orange", "green", "purple", "blue"))
name <- paste("comb_effect_b.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p2
dev.off()

##########################################################################################
#         Figure S1: Heatmaps for True Clustering Scheme and Patient Co-clustering       #
##########################################################################################

load(file="Simu.Co.Clustering.Rdata") # Csub: posterior patient co-clustering based on 100 repeated simulations

palette <- rep(NA, 100) # heatmap color palette
for(i in 1:100){ palette[i] = paste("gray", as.character(101-i), sep = "") }
x.scale <- list(cex=2, alternating=1); y.scale <- list(cex=2, alternating=1)

Ord <- order(data$pi_n) 
pi_n_0 <- data$pi_n[Ord] # reorder the patients by the true clustering scheme
Csub_0 <- matrix(0, nrow=data$n, ncol=data$n) # true clustering scheme
for(i in 1:data$n){
  Csub_0[i,] <- Csub_0[i,] + (pi_n_0==pi_n_0[i])
}

p1 <- levelplot(Csub_0, col.regions = palette, scales=list(x=x.scale, y=y.scale),
                xlab=list("Individual number", cex=2.25), ylab=list("Individual number", cex=2.25),
                colorkey=list(labels=list(cex=2))) 
name <- paste("co_clustering_a.pdf") # heatmaps of true clustering scheme 
pdf(name, width = 10, height = 10, onefile = TRUE)
p1
dev.off()

p2 <- levelplot(Csub, col.regions = palette, scales=list(x=x.scale, y=y.scale),
                xlab=list("Individual number", cex=2.25), ylab=list("Individual number", cex=2.25),
                colorkey=list(labels=list(cex=2)))
name <- paste("co_clustering_b.pdf") # heatmaps of posterior patient co-clustering
pdf(name, width = 10, height = 10, onefile = TRUE)
p2
dev.off()

##########################################################################################
#                        Table S2: Simulation Truths                                     #
##########################################################################################

print(data$beta1); print(data$beta2); print(data$beta3)

##########################################################################################
#                 Figure S3: Autocorrelation Plots for Parameters                        #
##########################################################################################

name <- paste("acf.pdf")
pdf(name, width = 18, height = 18, onefile = TRUE)
par(mfrow=c(3, 3), ps = 20, cex = 1)
p11 <- acf(result$beta[,1,1,1], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['1,1,1']))
p12 <- acf(result$beta[,1,2,2], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['1,2,2']))
p13 <- acf(result$beta[,1,3,3], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['1,3,3']))
p21 <- acf(result$beta[,2,1,1], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['2,1,1']))
p22 <- acf(result$beta[,2,2,2], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['2,2,2']))
p23 <- acf(result$beta[,2,3,3], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['2,3,3']))
p31 <- acf(result$beta[,3,1,1], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['3,1,1']))
p32 <- acf(result$beta[,3,2,2], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['3,2,2']))
p33 <- acf(result$beta[,3,3,3], lag.max = 30, ylab = "Autocorrelation", 
           main = expression(bold("Autocorrelation plot for ") ~ beta['3,3,3']))
dev.off()

##########################################################################################
#                       Figure S4: Trace Plots for Parameters                            #
##########################################################################################

theme_update(plot.title = element_text(hjust = 0.5, size = 24))

df11 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,1,1,1])
p11 <- ggplot(df11, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta1[1,1]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['1,1,1'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df12 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,1,2,2])
p12 <- ggplot(df12, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta1[2,2]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['1,2,2'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df13 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,1,3,3])
p13 <- ggplot(df13, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta1[3,3]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['1,3,3'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df21 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,2,1,1])
p21 <- ggplot(df21, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta2[1,1]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['2,1,1'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df22 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,2,2,2])
p22 <- ggplot(df22, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta2[2,2]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['2,2,2'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df23 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,2,3,3])
p23 <- ggplot(df23, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta2[3,3]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['2,3,3'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df31 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,3,1,1])
p31 <- ggplot(df31, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta3[1,1]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['3,1,1'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df32 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,3,2,2])
p32 <- ggplot(df32, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta3[2,2]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['3,2,2'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 
df33 <- data.frame(itertion = seq(1, 500, by=1), value = result$beta[,3,3,3])
p33 <- ggplot(df33, aes(iteration, value)) + geom_line(colour = "blue") +  
  geom_hline(aes(yintercept = data$beta3[3,3]), linetype = "dashed", colour = "red") +
  labs(title=expression(bold("Trace plot for ") ~ beta['3,3,3'])) + xlab("Iteration") + ylab(" ") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

name <- paste("trace_plot.pdf")
pdf(name, width = 18, height = 18, onefile = TRUE)
grid.arrange(p11, p12, p13, p21, p22, p23, p31, p32, p33,  widths = c(1, 1, 1),
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
dev.off()

##########################################################################################
#                     Figure S5: Historgram of Estimated Number of Clusters              #
##########################################################################################

load(file="Simu.Cluster.Number.Rdata") # estimated number of clusters in 100 repeated simulations

df <- data.frame(r_n = c(1,2,3,4,5), 
                 prob = c(sum(R_n==1)/length(R_n), sum(R_n==2)/length(R_n), sum(R_n==3)/length(R_n),
                          sum(R_n==4)/length(R_n), sum(R_n==5)/length(R_n)))
p <- ggplot(data=df, aes(x = r_n, y = prob)) + geom_bar(stat="identity") + 
  scale_x_continuous(breaks = seq(1, 5, by = 1)) + scale_y_continuous(limits=c(0,1)) + 
  xlab("Number of clusters") + ylab("Posterior probability") + 
  theme(legend.title=element_text(size = 24, face='bold')) + 
  theme(legend.text=element_text(size = 24)) + 
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(axis.text.x = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 24))
name <- paste("cluster_number.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#             Figure S6: Histogram and Density Plots for Combination Effects             #
##########################################################################################

h_est <- array(NA, dim=c(data$n, max(data$J), data$Q)) # posterior expected values of combination effects
for (i in 1:data$n){ for (j in 1:data$J[i]){ for (q in 1:data$Q){
  h_est[i,j,q] <- colMeans(result$gamma[,post$pi_n[i],q,] %*% data$H[i,j,]) }}}
h_lq <- array(NA, dim=c(data$n, max(data$J), data$Q)) # 2.5% lower quantile of estimated combination effects
for (i in 1:data$n){ for (j in 1:data$J[i]){ for (q in 1:data$Q){
  h_lq[i,j,q] <- quantile(result$gamma[,post$pi_n[i],q,] %*% data$H[i,j,], 0.025) }}}
h_uq <- array(NA, dim=c(n, max(J), Q)) # 97.5% upper quantile of estimated combination effects
for (i in 1:data$n){ for (j in 1:data$J[i]){ for (q in 1:data$Q){
  h_uq[i,j,q] <- quantile(result$gamma[,post$pi_n[i],q,] %*% data$H[i,j,], 0.975) }}}

name <- paste("comb_effect_density.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
hist(data$h[!is.na(data$h)], freq=F, main=" ", xlab="Combination effect", ylab="Density", cex.lab=1.5, cex.axis=1.5)
lines(density(h_est[!is.na(h_est)], bw=1.5), col='red', lty=1, lwd=2)
text <- c("Truth", "Estimated")
legend("topleft", legend=text, text.width = strwidth(text)[1], col=c("black", "red"), lty=rep(1,2), cex=1.5, box.lty=0)
dev.off()

##########################################################################################
#                        Table S7: Mean Squared Error                                    #
##########################################################################################

load(file="Simu.Mean.Squared.Error.Rdata") # mse: the mean squared error of beta accross 100 repeated simulations

formatC(colMeans(mse$beta1, dims=1), format = "e", digits = 3)
formatC(colMeans(mse$beta2, dims=1), format = "e", digits = 3)
formatC(colMeans(mse$beta3, dims=1), format = "e", digits = 3)

##########################################################################################
#                  Figure S8: 95% Credible Intervals for Parameters                      #
##########################################################################################

df1 <- data.frame(beta = seq(1, 9, by=1), value = c(as.vector(t(data$beta1))),
                  lower = c(as.vector(t(ci$beta[1,,,1]))), upper = c(as.vector(t(ci$beta[1,,,2]))))
p1 <- ggplot(df1, aes(beta, value)) + geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, size = 1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 1.5) + xlab(" ") + ylab(" ") +
  scale_x_discrete(limit = c(1:9),
                   labels = c(expression(beta['1,1,1']), expression(beta['1,1,2']), expression(beta['1,1,3']), 
                              expression(beta['1,2,1']), expression(beta['1,2,2']), expression(beta['1,2,3']),  
                              expression(beta['1,3,1']), expression(beta['1,3,2']), expression(beta['1,3,3']))) + 
  theme(axis.title.x = element_text(size = 45, face='bold')) + theme(axis.title.y = element_text(size = 50, face='bold')) +
  theme(axis.text.x = element_text(size = 45, hjust = 0)) + theme(axis.text.y = element_text(size = 50)) 
name <- paste("simu_results_a.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p1
dev.off()

df2 <- data.frame(beta = seq(1, 9, by=1), value = c(as.vector(t(data$beta2))),
                  lower = c(as.vector(t(ci$beta[2,,,1]))), upper = c(as.vector(t(ci$beta[2,,,2]))))
p2 <- ggplot(df2, aes(beta, value)) + geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, size = 1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 1.5) + xlab(" ") + ylab(" ") +
  scale_x_discrete(limit = c(1:9),
                   labels = c(expression(beta['2,1,1']), expression(beta['2,1,2']), expression(beta['2,1,3']), 
                              expression(beta['2,2,1']), expression(beta['2,2,2']), expression(beta['2,2,3']),  
                              expression(beta['2,3,1']), expression(beta['2,3,2']), expression(beta['2,3,3']))) + 
  theme(axis.title.x = element_text(size = 45, face='bold')) + theme(axis.title.y = element_text(size = 50, face='bold')) +
  theme(axis.text.x = element_text(size = 45, hjust = 0)) + theme(axis.text.y = element_text(size = 50)) 
name <- paste("simu_results_b.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p2
dev.off()

df3 <- data.frame(beta = seq(1, 9, by=1), value = c(as.vector(t(data$beta3))),
                  lower = c(as.vector(t(ci$beta[3,,,1]))), upper = c(as.vector(t(ci$beta[3,,,2]))))
p3 <- ggplot(df3, aes(beta, value)) + geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, size = 1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 1.5) + xlab(" ") + ylab(" ") +
  scale_x_discrete(limit = c(1:9),
                   labels = c(expression(beta['3,1,1']), expression(beta['3,1,2']), expression(beta['3,1,3']), 
                              expression(beta['3,2,1']), expression(beta['3,2,2']), expression(beta['3,2,3']),  
                              expression(beta['3,3,1']), expression(beta['3,3,2']), expression(beta['3,3,3']))) + 
  theme(axis.title.x = element_text(size = 45, face='bold')) + theme(axis.title.y = element_text(size = 50, face='bold')) +
  theme(axis.text.x = element_text(size = 45, hjust = 0)) + theme(axis.text.y = element_text(size = 50)) 
name <- paste("simu_results_c.pdf")
pdf(name, width = 18, height = 12, onefile = TRUE)
p3
dev.off()

##########################################################################################
#                     Table S9: Mean Squared Error of Combination Effects Estimation     #
##########################################################################################

# MSE mean of combination effect estimations
mean(h_est$h_mse_ddcrp_st) # 0.1821469
mean(h_est$h_mse_dp_st) # 0.1823971
mean(h_est$h_mse_normal_st) # 50.532
mean(h_est$h_mse_dp_linear) # 12.2538
mean(h_est$h_mse_normal_linear) # 56.94874

# MSE sd of combination effect estimations
sd(h_est$h_mse_ddcrp_st) # 0.005638027
sd(h_est$h_mse_dp_st) # 0.005636313
sd(h_est$h_mse_normal_st) # 0.0106698
sd(h_est$h_mse_dp_linear) # 4.934936
sd(h_est$h_mse_normal_linear) # 0.008427978

##########################################################################################
#              Figure S14: Sensitivity Analysis  (eta = 0.1, 0.3, 0.5, 0.8, 1)           #
##########################################################################################

load(file="Sensi.MCMC.Results.Rdata") # sensi_result: MCMC posterior samples for sensitivity analysis

post <- NULL # posterior clustering results
post$pi_n1 <- Clustering_Summary(sensi_result$pi_n1)
post$r_n1 <- length(unique(post$pi_n1))
post$pi_n2 <- Clustering_Summary(sensi_result$pi_n2)
post$r_n2 <- length(unique(post$pi_n2))
post$pi_n3 <- Clustering_Summary(sensi_result$pi_n3)
post$r_n3 <- length(unique(post$pi_n3))
post$pi_n4 <- Clustering_Summary(sensi_result$pi_n4)
post$r_n4 <- length(unique(post$pi_n4))
post$pi_n5 <- Clustering_Summary(sensi_result$pi_n5)
post$r_n5 <- length(unique(post$pi_n5))

ci <- NULL # 95% credible intervals of parameters 
ci$beta1 <- Cred_Interval(sensi_result$beta1[,1:post$r_n1,,])
ci$gamma1 <- Cred_Interval(sensi_result$gamma1[,1:post$r_n1,,])
ci$beta2 <- Cred_Interval(sensi_result$beta2[,1:post$r_n2,,])
ci$gamma2 <- Cred_Interval(sensi_result$gamma2[,1:post$r_n2,,])
ci$beta3 <- Cred_Interval(sensi_result$beta3[,1:post$r_n3,,])
ci$gamma3 <- Cred_Interval(sensi_result$gamma3[,1:post$r_n3,,])
ci$beta4 <- Cred_Interval(sensi_result$beta4[,1:post$r_n4,,])
ci$gamma4 <- Cred_Interval(sensi_result$gamma4[,1:post$r_n4,,])
ci$beta5 <- Cred_Interval(sensi_result$beta5[,1:post$r_n5,,])
ci$gamma5 <- Cred_Interval(sensi_result$gamma5[,1:post$r_n5,,])


index <- c(1,5,9) # index of selected parameters 
df1 <- data.frame(beta = rep(c(1:3), 5),
                  value = rep(c(as.vector(t(data$beta1))[index]),5),
                  lower = c(as.vector(t(ci$beta1[1,,,1]))[index], as.vector(t(ci$beta2[1,,,1]))[index], as.vector(t(ci$beta3[1,,,1]))[index], 
                            as.vector(t(ci$beta4[1,,,1]))[index], as.vector(t(ci$beta5[1,,,1]))[index]), 
                  upper = c(as.vector(t(ci$beta1[1,,,2]))[index], as.vector(t(ci$beta2[1,,,2]))[index], as.vector(t(ci$beta3[1,,,2]))[index], 
                            as.vector(t(ci$beta4[1,,,2]))[index], as.vector(t(ci$beta5[1,,,2]))[index]), 
                  eta = factor(c(rep("0.1",3), rep("0.3",3), rep("0.5",3), rep("0.8",3), rep("1",3))))
p1 <- ggplot(df1, aes(beta, value, colour = eta)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, position=position_dodge(0.8), size=1.75) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge(0.8), size=2) + 
  scale_x_discrete(limit = c(1:3), labels = c(expression(beta['1,1,1']), expression(beta['1,2,2']), expression(beta['1,3,3']))) + 
  xlab(" ") + ylab(" ") + labs(col = expression(eta)) + 
  theme(legend.title=element_text(size = 65, face='bold')) + theme(legend.text=element_text(size = 65)) + 
  theme(axis.title.x = element_text(size = 65, face='bold')) + theme(axis.title.y = element_text(size = 55, face='bold')) +
  theme(axis.text.x = element_text(size = 65)) + theme(axis.text.y = element_text(size = 55))
name <- paste("sensi_results_a.pdf")
pdf(name, width = 18, height = 16, onefile = TRUE)
p1
dev.off()

df2 <- data.frame(beta = rep(c(1:3), 5),
                  value = rep(c(as.vector(t(data$beta2))[index]),5),
                  lower = c(as.vector(t(ci$beta1[2,,,1]))[index], as.vector(t(ci$beta2[2,,,1]))[index], as.vector(t(ci$beta3[2,,,1]))[index], 
                            as.vector(t(ci$beta4[2,,,1]))[index], as.vector(t(ci$beta5[2,,,1]))[index]), 
                  upper = c(as.vector(t(ci$beta1[2,,,2]))[index], as.vector(t(ci$beta2[2,,,2]))[index], as.vector(t(ci$beta3[2,,,2]))[index], 
                            as.vector(t(ci$beta4[2,,,2]))[index], as.vector(t(ci$beta5[2,,,2]))[index]), 
                  eta = factor(c(rep("0.1",3), rep("0.3",3), rep("0.5",3), rep("0.8",3), rep("1",3))))
p2 <- ggplot(df2, aes(beta, value, colour = eta)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, position=position_dodge(0.8), size=1.75) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge(0.8), size=2) + 
  scale_x_discrete(limit = c(1:3), labels = c(expression(beta['2,1,1']), expression(beta['2,2,2']), expression(beta['2,3,3']))) + 
  xlab(" ") + ylab(" ") + labs(col = expression(eta)) + 
  theme(legend.title=element_text(size = 65, face='bold')) + theme(legend.text=element_text(size = 65)) + 
  theme(axis.title.x = element_text(size = 65, face='bold')) + theme(axis.title.y = element_text(size = 55, face='bold')) +
  theme(axis.text.x = element_text(size = 65)) + theme(axis.text.y = element_text(size = 55))
name <- paste("sensi_results_b.pdf")
pdf(name, width = 18, height = 16, onefile = TRUE)
p2
dev.off()

df3 <- data.frame(beta = rep(c(1:3), 5),
                  value = rep(c(as.vector(t(data$beta3))[index]),5),
                  lower = c(as.vector(t(ci$beta1[3,,,1]))[index], as.vector(t(ci$beta2[3,,,1]))[index], as.vector(t(ci$beta3[3,,,1]))[index], 
                            as.vector(t(ci$beta4[3,,,1]))[index], as.vector(t(ci$beta5[3,,,1]))[index]), 
                  upper = c(as.vector(t(ci$beta1[3,,,2]))[index], as.vector(t(ci$beta2[3,,,2]))[index], as.vector(t(ci$beta3[3,,,2]))[index], 
                            as.vector(t(ci$beta4[3,,,2]))[index], as.vector(t(ci$beta5[3,,,2]))[index]), 
                  eta = factor(c(rep("0.1",3), rep("0.3",3), rep("0.5",3), rep("0.8",3), rep("1",3))))
p3 <- ggplot(df3, aes(beta, value, colour = eta)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), shape=17, position=position_dodge(0.8), size=1.75) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge(0.8), size=2) + 
  scale_x_discrete(limit = c(1:3), labels = c(expression(beta['3,1,1']), expression(beta['3,2,2']), expression(beta['3,3,3']))) + 
  xlab(" ") + ylab(" ") + labs(col = expression(eta)) + 
  theme(legend.title=element_text(size = 65, face='bold')) + theme(legend.text=element_text(size = 65)) + 
  theme(axis.title.x = element_text(size = 65, face='bold')) + theme(axis.title.y = element_text(size = 55, face='bold')) +
  theme(axis.text.x = element_text(size = 65)) + theme(axis.text.y = element_text(size = 55))
name <- paste("sensi_results_c.pdf")
pdf(name, width = 18, height = 16, onefile = TRUE)
p3
dev.off()
