# BayesDrugComb
R data and code for paper "A Bayesian Nonparametric Approach for Inferring Drug Combination Effects on Mental Health in People with HIV".

## Data

The data analyzed in Section 4 “Application: WIHS Data Analysis” of the paper are from The Women's Interagency HIV Study (WIHS), which is a multisite, longitudinal cohort study of the natural and treated history of women living with HIV and women at-risk for HIV in the United States.
The data are publicly available. However, one need to fill in a request form for access. Full details of the data are available at https://statepi.jhsph.edu/wihs/wordpress/. R data for the simulation study and WIHS data analysis of the paper are available at https://drive.google.com/open?id=1FB8o0cHx0lVq-PdEGZciVoknB8nCUXCI.

## Code 

The R scripts in the folder “BNP_DrugComb_Simulation” are for Section 3 “Simulation Study”, and the R scripts in the folder “BNP_DrugComb_WIHS” are for Section 4 “Application: WIHS Data Analysis”. 

Libraries and Version Numbers: R 3.5.1, Rcpp 1.0.0, RcppArmadillo 0.9.700.2.0, RcppEigen 0.3.3.5.0, Matrix 1.2-14, LaplacesDemon 16.1.1, ggplot2 3.1.0, gridExtra 2.3, lattice 0.20-35.

### Instructions for Use

In the folder “BNP_DrugComb_Simulation”:

* The R script “Simulation_Main.R” reproduces Figure 4, Figure 5, Figure 6 in the manuscript, and Table S1, Figure S2, Figure S3, Figure S4, Table S5, Figure S6 in the Supplementary Material;
    
* The R data file “Simu.Data.Preprocess.Rdata” contains the preprocessed data from the WIHS dataset used for generating simulation truths. The R script “Data_Generate.R” generates the simulated dataset, which is saved in the R data file “Simu.Truths.Rdata”;  

* The R script “MCMC_R_Functions.R” provides R functions used for MCMC, the Rcpp script “MCMC_Rcpp_Functions.cpp” provides Rcpp functions used for MCMC, and the R
script “BNP_DrugComb_MCMC.R” provides the MCMC main function;

* We save the MCMC posterior samples for one randomly selected simulated dataset in the simulation study in the R data file “Simu.MCMC.Results.Rdata", and the MCMC posterior samples for the same simulated dataset in the sensitivity analysis in the R data file “Sensi.MCMC.Results.Rdata";

* We save the posterior co-clustering results based on 100 repeated simulations in the R data file “Simu.Co.Clustering.Rdata”, and the mean squared errors for all 100 repeated simulations in the R data file “Simu.Mean.Squared.Error.Rdata”;

* We save the estimated combination effects under two methods Normal+Linear and DP+Linear in the R data file “Simu.Combination.Effects.Linear.Kernel.Rdata”.

In the folder “BNP_DrugComb_WIHS”:

* The R script “BNP_DrugComb_WIHS_Inference.R” reproduces Figure 7, Figure 8, Table 1, and Figure 9 in the manuscript;

* The R script “MCMC_R_Functions.R” provides R functions used for MCMC, the Rcpp script “MCMC_Rcpp_Functions.cpp” provides Rcpp functions used for MCMC, and the R script “BNP_DrugComb_MCMC.R” provides the MCMC main function;

* The R script “Subset_Tree_Kernel_Similarity.R” provides functions for calculating the similarity matrix induced by subset-tree kernel, and the R script “Prediction.R” provides functions for predictions；

* We save the preprocessed data from the WIHS dataset for inference in the R data file “WIHS.Data.Preprocess.Rdata”, and the MCMC posterior samples for the WIHS data analysis in the R data file “WIHS.MCMC.Results.Rdata”.







