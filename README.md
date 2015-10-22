PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 


##Reference
Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK*. (2015) A gene-based association method for mapping traits using reference transcriptome data. Nat Genet. doi:10.1038/ng.3367. [Preprint on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164)

†:equal contribution *:correspondence


##Software

### Python version

- Download software from this [link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

PredictDB
=========
PredictDB hosts genetic prediction models of transcriptome levels to be used with PrediXcan. 

The following models are available for download. 

- DGN Whole Blood Elastic Net SQTLite db [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.db "DGN-WB-EN0.50.db") (A ~28MB file will be downloaded)
- DGN Whole Blood Elastic Net text [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.txt "DGN-WB-EN0.50.txt") (A ~18MB file will be downloaded)
- DGN Whole Blood LASSO SQTLite db [link](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_1.db "DGN-WB-EN1.db") (A ~18MB file will be downloaded)

G2Pdb
=========
G2Pdb, Gene to Phenotype database, hosts the results of PrediXcan applied to a variety of phenotypes. [Link to prototype](http://www.gene2pheno.org/). Currently, the prototype hosts the results of PrediXcan applied to WTCCC (Wellcome Trust Case Control Consoritium) diseases using DGN whole blood prediction models.


##Acknowledgements

### DGN RNA-seq data

Data downloaded from [NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M., McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et al. (2014). Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.


