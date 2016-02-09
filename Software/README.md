PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

##Reference
Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK*. (2015) A gene-based association method for mapping traits using reference transcriptome data. Nat Genet. doi:10.1038/ng.3367.

† equal contribution

[An open access preprint can be found on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164)

##Instructions

To run PrediXcan you will need

Software Requirements:

- Linux or Mac OS
- Python 2.7
    - numpy package
- R

Scripts:

- PrediXcan.py
- PrediXcanAssociation.R

Input Files: 

- genotype file 
- sample file
- phenotype file
- filter file - Specifies a subset of rows on which to perform association tests (optional)
- transcriptome prediction model (sqlite db to be downloaded from [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB-unscaled_0.5.db "DGN-WB-EN-unscaled_0.5") (A 33MB file will be downloaded if you click the link).

All the scripts used to develop the prediction models can be found [here](https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts/Heather/DGN-calc-weights "Prediction Model Pipeline")

###Imputing Expression

To predict the transcriptome from a given genotype file, include the `--predict` flag when running PrediXcan.py and specify the following arguments:

1. genelist: list of genes. Optional. By default it will use all available genes in model database
2. dosages: imputed genotype file path. Default value: 'data/dosages/'
3. dosage_prefix: prefix of dosage file. Default value: 'chr' 
4. weights: full name of database. Default value: 'weight.db'
5. output_dir: path to the desired output directory.  Default value: 'output/'

This will produce a file in the specified output directory called `predicted_expression.txt`, which contains all of the predicted expression levels.

####Dosage File Format
- Columns are snpid rsid position allele1 allele2 MAF id1 ..... idn.
- Dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2).
- It is expected that there will be one file per chromosome.
- There must be a file of the individuals with id #'s listed in the same order as the genotype columns. This should be in the format of a PLINK .fam file.

####Usage
> ./PrediXcan.py  --predict --dosages dosagefile_path  --dosages_prefix chr --weights prediction_db --output_dir output

####Prediction Example
- Download and untar this file [Prediction Example tar file](https://s3.amazonaws.com/imlab-open/Data/PredictDB/predixcan_predict_example.tar)
    - To untar, go the folder and run `tar -xvf predixcan_predict_example.tar`
- Go to folder and run the following

> ./PrediXcan.py --predict --dosages dosages --dosages_prefix chr --weights DGN-WB_0.5.db --output_dir output

###Running Association with Phenotype

To perform an association test between the predicted expression levels and phenotype, include the `--assoc` flag when running PrediXcan.py and specify the following arguments:

1. pred_exp: predicted transcriptome from previous run of PrediXcan.  Default value: 'output/predicted_expression.txt'.
2. pheno: phenotype file.  No default value.  See below for file format.
3. filter: filter file to specify which rows to include in test and a number to filter on.  Optional. See below for details.
4. linear or logistic: specify one of these to perform a linear or logistic regression between the expression levels of each gene and phenotype.  Default is logistic.
5. output_dir: path to the desired output directory. Default value: 'output'

This will produce a file in the specified output directory called `association.txt`, with summary statistics on the association between each gene and the phenotype.

####Phenotype File Format

Phenotype files are expected to be in a format similar to the format required for PLINK.  Most commonly, the phenotype file is tab delimited, and preferably has a header.  By default, PrediXcan will assume the first column is the Family ID, the second column is the Individual ID, and the *last* column is the phenotype column.

**Note**: If the phenotype file has a header line, which preferably it will, the first two columns **must** be labeled FID and IID, respectively.  If there are multiple phenotype columns, you can specify which column to perform the association on with the `--pheno_name` flag.

If there is more than one phenotype column in the file, you can specify which phenotype to perform the association on with the `--mpheno` option.  For example `--mpheno 1` will do the association with the 3rd column in the phenotype file, as columns 1 and 2 are ID numbers, `--mpheno 2` does the association on 4th, etc. This option will mainly be used for when there is no header line, and may behave unexpectedly if the user does not specify options carefully.

Unlike PLINK, for logistic tests on qualititative traits, by default the trait is assumed to be encoded as 0 for unaffected and 1 for affected.  0 is NOT a missing value.

By default, NA specifies a missing phenotype value.  To specify a missing phenotype value that is encoded numerically, say -9 for example, include `--missing_phenotype -9'.

####Filter File Format

Filter files can specify a subset rows in the pheno file to perform the association on.  It is a tab delimited file with the first 2 columns identical to the pheno file.  The third column holds numerical values on which to filter.  If the filter file is called filter.txt, with filter values 1 and 2, including `--filter filter.txt 2` will perform the association test only on individuals marked 2 in the filter file.

####Usage
> ./PrediXcan.py --assoc --pheno phenotype_file --pred_exp predicted_expression_file --linear --filter filter_file filter_val --output_dir output_directory

####Association Example
- Download and untar this file [Association Example tar file](https://s3.amazonaws.com/imlab-open/Data/PredictDB/predixcan_association_example.tar)
- Go to folder and run the following

> ./PrediXcan.py --assoc --pheno igrowth_phenotype.txt --filter igrowth_phasefilter.txt 3 --pred_exp predicted_expression.txt --linear

** to speed up the process the dosage files can be filtered to SNPs in HapMapSnpsCEU.list.gz.

hapmapSnpsCEU.list.gz: List of SNPs used to develop the models is [here](https://app.box.com/s/6ftz3lr5h6detnf2iwzc7soyo5szrrej "HapMap2 SNP set") downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz "HapMap2 UCSC"). 
SNPs not included on this list is not used to make predictions.

#### Helper Scripts
Conversion from Plink to Dosage (provided by scottritchie73 via pull request, thank you!) [link](https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py)
