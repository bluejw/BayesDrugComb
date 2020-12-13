# BayesDrugComb
R data and code for the paper "A Bayesian Nonparametric Approach for Inferring Drug Combination Effects on Mental Health in People with HIV".

## Data

The data analyzed in Section 4 “Application: WIHS Data Analysis” of the paper are from The Women's Interagency HIV Study (WIHS), which is a multisite, longitudinal cohort study of the natural and treated history of women living with HIV and women at-risk for HIV in the United States.
The data are publicly available. However, one need to fill in a request form for access. Full details of the data are available at https://statepi.jhsph.edu/wihs/wordpress/. R data for the simulation study and WIHS data analysis of the paper are available at https://drive.google.com/open?id=1FB8o0cHx0lVq-PdEGZciVoknB8nCUXCI.

## Code 

The R scripts in the folder “BNP_DrugComb_Simulation” are for Section 3 “Simulation Study”, and the R scripts in the folder “BNP_DrugComb_WIHS” are for Section 4 “Application: WIHS Data Analysis”. 

Libraries and Version Numbers: R 3.5.1, Rcpp 1.0.0, RcppArmadillo 0.9.700.2.0, RcppEigen 0.3.3.5.0, Matrix 1.2-14, LaplacesDemon 16.1.1, ggplot2 3.1.0, gridExtra 2.3, lattice 0.20-35.

### Instructions for Use

In the folder “BNP_DrugComb_Simulation”:

* The R script “Simulation_Main.R” reproduces Figure 4, Figure 5, Figure 6 in the manuscript, and Table S1, Figure S2, Figure S3, Figure S4, Table S5, Figure S6 in the Supplementary Material;
    
    
