PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

##Reference
Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK*. (2015) A gene-based association method for mapping traits using reference transcriptome data. Nat Genet. doi:10.1038/ng.3367.

† equal contribution

[An open access preprint can be found on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164)

##Instructions

These instructions are for generating predicted expression levels. The association has to be performed separately. 

A beta version of the software that does both prediction and association can be found here [HOWTO-beta.md](https://github.com/hakyimlab/PrediXcan/blob/master/Software/HOWTO-beta.md)

To run PrediXcan you will need 

Input: 

- genotype file 
- phenotype file
- transcriptome prediction model (sqlite db to be downloaded from [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.db "DGN-WB-EN0.50") (A 28MB file will be downloded if you click the link).

 	- tissue: Whole Blood (default)
	- source: DGN (default)
	- model: Elastic Net 0.50 (default)

download [this](https://github.com/hakyimlab/PrediXmod/blob/master/PrediXcan/predict_gene_expression.py "Prediction Script") python script to compute predicted expression levels. Python 2.7 is needed.

The cross validated performance measures for each gene can be found in each db.

The script predict\_gene\_expression.py predicts gene expression levels using prediction models (stored in sqlite db such as DGN-WB_0.5.db) and whole genome variation data.

For now the association with phenotype needs to be performed manually. We are currently working on an R package that will do both the prediction of expression levels and the association with the phenotype.

All the scripts used to develop the prediction models can be found [here](https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts/Heather/DGN-calc-weights "Prediction Model Pipeline")


Supported operating systems:
Linux or Mac Os.


###Transcriptome Prediction


The following arguments are allowed, with default values as follows

1. genelist: list of genes. By default it will use all available genes in model database
2. dosages: imputed genotype file path. Default value: 'data/dosages/'
3. dosage_prefix: prefix of dosage file. Default value: 'chr' 
4. weights: full name of database. Default value: 'weight.db'
5. output: output file name 'output'

####dosage file format
- columns are snpid rsid position allele1 allele2 MAF dosage_1 ..... dosage_n 
- dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2)
- it is expected that there will be one file per chromosome

####USAGE
> ./predict_gene_expression.py  --dosages dosagefile_path  --dosages_prefix chr --weights prediction_db --output output

####Example
- Download and untar this file [Working Example tar file](https://s3.amazonaws.com/imlab-open/Data/PredictDB/predixcan-working-example.tar)
- Go to folder and run the following

> ./predict_gene_expression.py  --dosages dosages  --dosages_prefix chr --weights DGN-WB_0.5.db --output output

** to speed up the process the dosage files can be filtered to SNPs in HapMapSnpsCEU.list.gz.

hapmapSnpsCEU.list.gz: List of SNPs used to develop the models is [here](https://app.box.com/s/6ftz3lr5h6detnf2iwzc7soyo5szrrej "HapMap2 SNP set") downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz "HapMap2 UCSC"). 
SNPs not included on this list is not used to make predictions.

#### Helper Scripts
Conversion from Plink to Dosage (provided by scottritchie73 via pull request, thank you!) [link](https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py)
