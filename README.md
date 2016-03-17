PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

You have only summary results? Try this new extension of PrediXcan that uses only summary statistics, no individual level data [link](https://github.com/hakyimlab/MetaXcan)


##Reference
- Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K, Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK. (2015) **A gene-based association method for mapping traits using reference transcriptome data**. Nat Genet. doi:10.1038/ng.3367. ([Link to paper](http://www.nature.com/ng/journal/v47/n9/full/ng.3367.html), [Link to Preprint on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164) )

  †:equal contribution 
  
  *:correspondence haky at uchicago dot edu

##Software

### Python version

- Download software from this [link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

PredictDB
=========
PredictDB hosts genetic prediction models of transcriptome levels to be used with PrediXcan [PredictDB](http://predictdb.org). 

Below are older links

- GTEx and DGN models (internal release 6/13/2014, models ran on 11/13/2015 by Kaanan, db and cov processed by Alvao) 
  [Box Link](https://app.box.com/s/gujt4m6njqjfqqc9tu0oqgtjvtz9860w).
  - these models have now ENSid's instead of gene symbols in the gene column. The gene symbol is in the 'extra' table in the db under the genename columns.

G2Pdb
=========
G2Pdb, Gene to Phenotype database, hosts the results of PrediXcan applied to a variety of phenotypes. [Link to prototype](http://www.gene2pheno.org/). Currently, the prototype hosts the results of PrediXcan applied to WTCCC (Wellcome Trust Case Control Consoritium) diseases using DGN whole blood prediction models.

Genetic Architecture of Gene Expression Traits
=======
- Heather E Wheeler, Kaanan P Shah, Jonathon Brenner, Tzintzuni Garcia, Keston Aquino-Michaels, GTEx Consortium, Nancy J Cox, Dan L Nicolae, and Hae Kyung Im (2016) **Survey of the Heritability and Sparsity of Gene Expression Traits Across Human Tissues** [Link to Preprint](http://dx.doi.org/10.1101/043653), ; correspondence hwheeler at luc dot edu and haky at uchicago dot edu
- Database of heritability estimates [link](https://s3.amazonaws.com/imlab-open/Webdata/Paper-Links/h2r2-2016-03-17-no-TS.db) or [link](https://github.com/WheelerLab/GenArchDB)

##Acknowledgements

### DGN RNA-seq data

Data downloaded from [NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M., McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et al. (2014). Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.


