PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

PredictDB
=========
PredictDB hosts genetic prediction models of transcriptome levels to be used with PrediXcan. 

The following models are available for download. 

- DGN Whole Blood Elastic Net SQTLite db [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.db "DGN-WB-EN0.50.db") (A ~28MB file will be downloaded)
- DGN Whole Blood Elastic Net text [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.txt "DGN-WB-EN0.50.txt") (A ~18MB file will be downloaded)
- DGN Whole Blood LASSO SQTLite db [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_1.db "DGN-WB-EN1.db") (A ~18MB file will be downloaded)

##Reference
Eric R. Gamazon†, Heather E. Wheeler†, Sahar Mozaffari†, Kaanan P. Shah, Keston Aquino-Michaels, Robert J. Carroll, Anne E. Eyler, Joshua C. Denny, GTEx Consortium, Dan L. Nicolae, Nancy J. Cox, and Hae Kyung Im* **PrediXcan: Trait Mapping Using Human Transcriptome Regulation** (2015 Under revision)

##Instructions

To run PrediXcan you will need 

Input: 

- genotype file 
- phenotype file
- transcriptome prediction model (sqlite db to be downloaded from [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.db "DGN-WB-EN0.50") (A 28MB file will be downloded if you click the link).

 	- tissue: Whole Blood (default)
	- source: DGN (default)
	- model: Elastic Net 0.50 (default)

download [this](https://github.com/hakyimlab/PrediXmod/blob/master/PrediXcan/predict_gene_expression.py "Prediction Script") python script to compute predicted expression levels. Python 2.7 is needed.

The cross validated performance measures for each gene will be added to the db. (- [ ] TODO)

The script predict\_gene\_expression.py predicts gene expression levels using prediction models (stored in sqlite db such as DGN-WB_0.5.db) and whole genome variation data.

For now the association with phenotype needs to be performed manually. We are currently working on an R package that will do both the prediction of expression levels and the association with the phenotype.

All the scripts used to develop the prediction models can be found [here](https://github.com/hwheeler01/PrediXmod "Prediction Model Pipeline")


Supported operating systems:
Linux or Mac Os.

A Perl version is also available 

- Download scripts [here](https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts/Kaanan)
- Download prediction model in text format [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.txt "DGN-WB-EN0.50.txt") (A ~18MB file will be downloaded)

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
> ./predict_gene_expression.py  --dosages dosages  --dosages_prefix chr --weights DGN-WB_0.5.db --output output

** to speed up the process the dosage files can be filtered to SNPs in HapMapSnpsCEU.list.gz.

hapmapSnpsCEU.list.gz: List of SNPs used to develop the models is [here](https://app.box.com/s/6ftz3lr5h6detnf2iwzc7soyo5szrrej "HapMap2 SNP set") downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz "HapMap2 UCSC")

##Acknowledgements

### DGN RNA-seq data

Data downloaded from [NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M., McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et al. (2014). Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.


