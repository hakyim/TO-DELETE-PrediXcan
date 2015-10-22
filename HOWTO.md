# Introduction

*gene_association.py* is a command line tool to generate a gene expression
matrix and run association tests between the expression matrix and a phenotype file.
 
The following need to be installed in order to run properly:

- Python 2.7
- R
- The python packages numpy and rpy2
- The R package dplyr

In addition, the script *PrediXcanAssociation.R* needs to be located in the same directory as
**gene_association.py* to work correctly.

# Input Data

This script needs the following data:

- genotype biallelic SNP dosage data
- a Transcriptome Prediction Model (see sample files below)
- a phenotype file
- a filter file

It expects this input in a folder structure (at the same level the script is run) such as:

```bash
data/dosages     	# folder containing gzipped dosages, such as "chr1.dosage.txt.gz"
data/weights.db  	# sqlite database with weights from transcriptome model
data/GWAS.fam   	# phenotype file corresponding to the dosages file
data/GWAS.filter 	# filter for which individuals to use
```

If the input files are named exactly as in the previous example, you can just run:

```bash
$ ./gene_association.py
```

And it will output a txt file with the transcription matrix and association results.

## Parameters

The script supports some command line arguments. Type:
```bash
$ ./gene_association.py --help
```
for reference. 

For example:

```bash
$ ./gene_association.py \
    --weights data/DGN-WB_0.5.db \
    --dosages another_folder \
    --pheno data/phenotype.fam \
    --filter data/phenotype_filter.filter \
    --output results
```

will cause the dosages at *another_folder* and the database *data/DGN-WB_0.5.db* to be read,
and the result will be saved in *results/predicted_expression.txt*.  Then, the phenotype
and filter files will be read from *data/phenotype.fam* and *data/phenotype_filter.filter* 
respectively. The results from the association between the phenotype and predicted expression
will be saved in *results/association.txt*.

## Sample data

### Transcriptome models

You can download sqlite database models from [here](https://app.box.com/s/5nejbvzgsis77wtrxt8xokyn7pcdgnnp)

### Dosage data

These files are expected to be compressed files in the following format:
```
# chr1.dosage.txt.gz sample content, made up of the following columns:
# CHROMOSOME RSID POSITION ALLELE_1 ALLELE_2 ALLELE_2_FREQUENCY_AVERAGE (...)
# where (...) means a list of allele frequency values for different persons
...
chr1 rs10399749 55299 C T 0.164302600472813 0 0 1 ...
...
```

If you have access to Hae Kyung Im's public data repository, and have Amazon Web Services command line tools,
you can get dosage by executing:

``` bash
mkdir data
mkdir data/dosages
cd data/dosages
aws s3 cp s3://imlab-open/Data/1000Genomes/Transcriptome/GEUVADIS/dosagefiles-hapmap2/ . --recursive
```

### Phenotype and Filter Files

The phenotype file is expected to be formatted with the following columns:
```
FAMILY INDIVIDUAL PATERNAL MATERNAL SEX PHENOTYPE
```

The filter file should be formatted similarly with the following columns:
```
FAMILY INDIVIDUAL FILTER
```

##Association Test

The association test works by performing a linear regression on the
phenotype and the expression level of each individual gene.  The results of each
regression are saved as a row in file *association.txt*  with the following
entries:
```
gene beta z-stat p-val se(beta)
```