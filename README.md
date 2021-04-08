# BayesDrugComb
R data and code for the paper:

"A Bayesian Nonparametric Approach for Inferring Drug Combination Effects on Mental Health in People with HIV".

## Data

The data analyzed in Section 4 “Application: WIHS Data Analysis” of the paper are from The Women's Interagency HIV Study (WIHS), which is a multisite, longitudinal cohort study of the natural and treated history of women living with HIV and women at-risk for HIV in the United States.
The data are publicly available. However, one need to fill in a request form for access. Full details of the data are available at https://statepi.jhsph.edu/wihs/wordpress/. R data for the simulation study and the WIHS data analysis of the paper are available at the following link: https://drive.google.com/file/d/1U2lZK9UwS60ABzhO5AtOV7ohPsnG49l1/view?usp=sharing.

## Code 

The R scripts in the folder “BNP_DrugComb_Simulation” are for Section 3 “Simulation Study”, and the R scripts in the folder “BNP_DrugComb_WIHS” are for Section 4 “Application: WIHS Data Analysis”. 

Libraries and Version Numbers: R 3.5.1, Rcpp 1.0.0, RcppArmadillo 0.9.700.2.0, RcppEigen 0.3.3.5.0, Matrix 1.2-14, LaplacesDemon 16.1.1, ggplot2 3.1.0, gridExtra 2.3, lattice 0.20-35.

### Instructions for Use

In the folder “BNP_DrugComb_Simulation”:

* The R script “Simulation_Main.R” reproduces Figure 4 in the manuscript, and Figure/Table S1-S9 and S14 in the Supplementary Material;
    
* The R data file “Simu.Data.Preprocess.Rdata” contains the preprocessed data from the WIHS dataset used for generating simulation truths. The R script “Data_Generate.R” generates the simulated dataset, which is saved in the R data file “Simu.Truths.Rdata”;  

* The R script “MCMC_R_Functions.R” provides R functions used for MCMC, the Rcpp script “MCMC_Rcpp_Functions.cpp” provides Rcpp functions used for MCMC, and the R script “BNP_DrugComb_MCMC.R” provides the MCMC main function;

* We save the MCMC posterior samples for one randomly selected simulated dataset in the simulation study in the R data file “Simu.MCMC.Results.Rdata", and the MCMC posterior samples for the same simulated dataset in the sensitivity analysis in the R data file “Sensi.MCMC.Results.Rdata";

* We save the posterior co-clustering results based on 100 repeated simulations in the R data file “Simu.Co.Clustering.Rdata”, the posterior estimated number of clusters for all the 100 repeated simulations in the R data file "Simu.Cluster.Number.Rdata", and the mean squared errors for all 100 repeated simulations in the R data file “Simu.Mean.Squared.Error.Rdata”;

* We save the estimated combination effects and the MSE in the 100 repeated simulations under our proposed method ddCRP+ST, and four alternative methods Normal+ST, Normal+Linear, DP+ST, and DP+Linear in the R data file “Simu.Combination.Effects.Rdata”.

In the folder “BNP_DrugComb_WIHS”:

* The R script “BNP_DrugComb_WIHS_Inference.R” reproduces Figure 5-7, Table 1 in the manuscript, and Table S10 in the Supplementary Material;

* The R script “MCMC_R_Functions.R” provides R functions used for MCMC, the Rcpp script “MCMC_Rcpp_Functions.cpp” provides Rcpp functions used for MCMC, and the R script “BNP_DrugComb_MCMC.R” provides the MCMC main function;

* The R script “Subset_Tree_Kernel_Similarity.R” provides functions for calculating the similarity matrix induced by subset-tree kernel, and the R script “Prediction.R” provides functions for predictions；

* We save the preprocessed data from the WIHS dataset for inference in the R data file “WIHS.Data.Preprocess.Rdata”, and the MCMC posterior samples for the WIHS data analysis under our proposed method ddCRP+ST and four altervative methods Normal+ST, Normal+Linear, DP+ST, and DP+Linear in R data files “WIHS.MCMC.Results.ddCRP.ST.Rdata”, “WIHS.MCMC.Results.Normal.ST.Rdata”, “WIHS.MCMC.Results.Normal.Linear.Rdata”, “WIHS.MCMC.Results.DP.ST.Rdata”, and “WIHS.MCMC.Results.DP.Linear.Rdata", respectively.
