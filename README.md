PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 


##Reference
Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK*. (2015) A gene-based association method for mapping traits using reference transcriptome data. Nat Genet. doi:10.1038/ng.3367. 

[Link to paper](http://www.nature.com/ng/journal/v47/n9/full/ng.3367.html)

[Preprint on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164)

†:equal contribution *:correspondence haky at uchicago dot edu


##Software

### Python version

- Download software from this [link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

PredictDB
=========
PredictDB hosts genetic prediction models of transcriptome levels to be used with PrediXcan. 

The following models are available for download. 

- DGN Whole Blood Elastic Net SQLite db [Download](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.db "DGN-WB-EN0.50.db") (A ~28MB file will be downloaded)
- DGN Whole Blood Elastic Net text [Download](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_0.5.txt "DGN-WB-EN0.50.txt") (A ~18MB file will be downloaded)
- DGN Whole Blood LASSO SQLite db [Download](https://s3.amazonaws.com/imlab-open/Data/PredictDB/DGN-WB_1.db "DGN-WB-EN1.db") (A ~18MB file will be downloaded)

- GTEx models (internal release 6/13/2014, models ran on 11/13/2015 by Kaanan) [Box Link](https://app.box.com/s/gujt4m6njqjfqqc9tu0oqgtjvtz9860w).
  - these models have now ENSid's instead of gene symbols in the gene column. The gene symbol is in the 'extra' table in the db under the genename columns.

G2Pdb
=========
G2Pdb, Gene to Phenotype database, hosts the results of PrediXcan applied to a variety of phenotypes. [Link to prototype](http://www.gene2pheno.org/). Currently, the prototype hosts the results of PrediXcan applied to WTCCC (Wellcome Trust Case Control Consoritium) diseases using DGN whole blood prediction models.


##Acknowledgements

### DGN RNA-seq data

Data downloaded from [NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M., McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et al. (2014). Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.


