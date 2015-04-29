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

G2Pdb
=========
G2Pdb, Gene to Phenotype database, hosts the results of PrediXcan applied to a variety of phenotypes. [Link to prototype](http://52.10.85.136/imlab1/results.php). Currently, the prototype hosts the results of PrediXcan applied to WTCCC (Wellcome Trust Case Control Consoritium) diseases using DGN whole blood prediction models.


##Reference
Eric R. Gamazon†, Heather E. Wheeler†, Sahar Mozaffari†, Kaanan P. Shah, Keston Aquino-Michaels, Robert J. Carroll, Anne E. Eyler, Joshua C. Denny, GTEx Consortium, Dan L. Nicolae, Nancy J. Cox, and Hae Kyung Im* **PrediXcan: Trait Mapping Using Human Transcriptome Regulation** (2015 Under revision) [download link](https://s3.amazonaws.com/imlab-open/Webdata/Papers/2015/PrediXcan/PrediXcan.pdf)


##Software

### Python version

- Gene expression prediction script [link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

### Perl version 

- Download scripts [here](https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts/Kaanan)
- Download prediction model in text format [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.txt "DGN-WB-EN0.50.txt") (A ~18MB file will be downloaded)

### R version

- Gene expression prediction script [link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

##Acknowledgements

### DGN RNA-seq data

Data downloaded from [NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M., McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et al. (2014). Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.


