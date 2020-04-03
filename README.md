# Deprecation notice

This repository contains the original reference implementation of 
the [PrediXcan](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552594/) method.
It is now considered deprecated and exists only for reference purposes.

Active development is now conducted at the [MetaXcan](https://github.com/hakyimlab/MetaXcan) repository. Tutorial for this new version is [here](https://github.com/hakyimlab/MetaXcan/wiki/Individual-level-PrediXcan:-introduction,-tutorials-and-manual)

## PrediXcan

PrediXcan is a gene-based association test that prioritizes genes that
are likely to be causal for the phenotype. 

Do you have only summary results? Try
[MetaXcan](https://github.com/hakyimlab/MetaXcan), a new extension of
PrediXcan that uses only summary statistics. No individual level data
necessary.

## Mailing List

Please join this [Google Group](https://groups.google.com/forum/#!forum/predixcanmetaxcan)
for news on releases, features, etc. For support and feature requests,
you can use this repository's issue tracker.


## Reference
- Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K,
Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK. (2015)
**A gene-based association method for mapping traits using reference
transcriptome data**. Nat Genet. doi:10.1038/ng.3367.
([Link to paper](http://www.nature.com/ng/journal/v47/n9/full/ng.3367.html),
[Link to Preprint on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164))

  †:equal contribution 
  
  *:correspondence haky at uchicago dot edu

- Alvaro Barbeira, Kaanan P Shah, Jason M Torres, Heather E Wheeler,
Eric S Torstenson, Todd Edwards, Tzintzuni Garcia, Graeme I Bell,
Dan Nicolae, Nancy J Cox, Hae Kyung Im. (2016) **MetaXcan: Summary
Statistics Based Gene-Level Association Method Infers Accurate PrediXcan
Results** [link to preprint](http://dx.doi.org/10.1101/045260)

- Heather E Wheeler, Kaanan P Shah, Jonathon Brenner, Tzintzuni Garcia,
Keston Aquino-Michaels, GTEx Consortium, Nancy J Cox, Dan L Nicolae, Hae
Kyung Im. (2016) **Survey of the Heritability and Sparsity of Gene
Expression Traits Across Human Tissues**.
[link to preprint](http://dx.doi.org/10.1101/043653)

## Software

### Python version

- Download software from this
[link](https://github.com/hakyimlab/PrediXcan/tree/master/Software)

## PredictDB

[PredictDB](http://predictdb.org/) hosts genetic prediction
models of transcriptome levels to be used with PrediXcan. See
[our wiki](https://github.com/hakyimlab/PrediXcan/wiki/PredictDB-Update:-Aug-18,-2016)
for a report of a recent update of the prediction models.

## Gene2Pheno database of results

G2Pdb, Gene to Phenotype database, hosts the results of PrediXcan
applied to a variety of phenotypes. [Link to prototype](http://gene2pheno.org/). 

# Genetic Architecture of Gene Expression Traits

- Heather E Wheeler, Kaanan P Shah, Jonathon Brenner, Tzintzuni Garcia,
Keston Aquino-Michaels, GTEx Consortium, Nancy J Cox, Dan L Nicolae, and
Hae Kyung Im (2016) **Survey of the Heritability and Sparsity of Gene
Expression Traits Across Human Tissues**
[Link to Preprint](http://dx.doi.org/10.1101/043653); correspondence
hwheeler at luc dot edu and haky at uchicago dot edu
- Database of heritability estimates
[link](https://github.com/jlbren/GenArchDB)
[older link](https://s3.amazonaws.com/imlab-open/Webdata/Paper-Links/h2r2-2016-03-17-no-TS.db)
or [older link](https://github.com/WheelerLab/GenArchDB)

## Acknowledgements

### GTEx data

Data downloaded from dbGaP [link](http://www.gtexportal.org/)

### DGN RNA-seq data

Data downloaded from
[NIMH Repository and Genomics Resource](https://www.nimhgenetics.org )

Battle, A., Mostafavi, S., Zhu, X., Potash, J.B., Weissman, M.M.,
McCormick, C., Haudenschild, C.D., Beckman, K.B., Shi, J., Mei, R., et
al. (2014). Characterizing the genetic basis of transcriptome diversity
through RNA-sequencing of 922 individuals. Genome Research 24, 14–24.

___
[![Analytics](https://ga-beacon.appspot.com/UA-61894206-3/PrediXcan-Readme-Github?useReferrer)](https://github.com/hakyim/PrediXcan)
