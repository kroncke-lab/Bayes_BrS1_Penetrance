# A Bayesian method to estimate disease penetrance from genetic variant properties
Here we propose a method to estimate the penetrance of the heart arrhythmia Brugada syndrome attributable to variants in the cardiac sodium channel gene *SCN5A*. All code and data used and referenced in the manuscript (not yet accepted) are included here. An overview of the method and results are provided in "SCN5A-report.html". The R markdown file used to generate this report is "SCN5A-BrS1-penetrance-report.Rmd", which describes the methods and code used in the resulting manuscript. The raw data used are contained in "VariantSCN5A-second-revision.db" and in the "covariates" folder. The "distance_file" contains raw distances between centroid atoms in the three-dimensional structure of NaV1.5 (protein product of *SCN5A*); these distance data are used to generate the BrS penetrance density in the chunk starting line 379 using the "func_dist_seq.R" script. "scn5a_dataset.csv" is an output csv file which contains raw observation data, covariates, and penetrance priors. 

## Folders
The "covariates" folder contains all the predictive covariates/variant specific features used in the regression model. All other folders are the output of "SCN5A-BrS1-penetrance-report.Rmd" for the resulting "SCN5A-report.html" report.
