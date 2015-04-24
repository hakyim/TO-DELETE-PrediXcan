Kaanan Shah
kaanan@uchicago.edu
April 8, 2015
Running PrediXcan

The set of scripts in this folder are set up to run prediXcan on tarbell. 

PrediXcan_with_DGN_WTCCC.pl: master script that will run all processes. Each step must be run separately by using the 0/1 switches at the top of the master script. The steps are outlined below with their individual script

1. WTCCC_imputation_QCcheck.pl: creates simple histograms for imputation quality (r2) and MAF to make sure things generally look as expected. This step reads in vcf info files downloaded from the Michigan imputation server (this step is optional). 

2. vcf_to_dose_hapmap2.pl: this step converts vcf files outputted from imputation using the Michigan imputation server and converts them to the dos format needed to run prediXcan. This step only retains hapmap2 snps because our PrediXcan predictors are restricted to this set of snps. This was done primarily to keep files smaller and speed computational time. 

3. runPrediXcan3.pl: this step uses a set of predictors (beta file created by Heather) to predict gene expression values in the cohort of interest. cohort

4. prediXcanAssocation_jointimpute.r: this script is R code to test for association between case/control status and the predicted gene expression values using logistic regression with no covariates. This script assumes the cases and controls are in the same file (jointly imputed).

5. runGWAS.pl: this optional final script runs a standard GWAS in plink for the imputed data. again with no covariates and assuming cases and controls are imputed in the same file. DO NOT run this part of the script at the same time as others. Plink requires unzipped files while all other scripts above read zipped files so running them at the same time will mess things up. 

WTCCC_Summaries.R  is the code used to created figures and tables for this part of the PrediXcan paper. 

