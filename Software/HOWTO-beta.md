# Introduction

*gene_association.py* is a command line tool to generate a gene expression
matrix and run association tests between the expression matrix and a phenotype file.
 
The following need to be installed in order to run properly:

- Python 2.7
- R
- The python packages numpy and rpy2
- The R package dplyr

In addition, the script *PrediXcanAssociation.R* needs to be located in the same directory as
*gene_association.py* to work correctly.

# Input Data

This script needs the following data:

- genotype biallelic SNP dosage data
- a Transcriptome Prediction Model (see sample files below)
- a phenotype file (see below)
- a filter file (see below)

It expects this input in a folder structure (at the same level the script is run) such as:

```bash
data/dosages     	# folder containing gzipped dosages, such as "chr1.dosage.txt.gz"
data/weights.db  	# sqlite database with weights from transcriptome model
data/GWAS.fam   	# phenotype file corresponding to the dosages file
data/GWASfilter.txt	# filter for which individuals to use
```

If the input files are named exactly as in the previous example, you can just run:

```bash
$ ./gene_association.py
```

And it will output txt files with the transcription matrix and association results.

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
    --filter data/phenotype_filter.txt \
    --output results
```

will cause the dosages at *another_folder* and the database *data/DGN-WB_0.5.db* to be read,
and the result will be saved in *results/predicted_expression.txt*.  Then, the phenotype
and filter files will be read from *data/phenotype.fam* and *data/phenotype_filter.txt* 
respectively. The results from the association between the phenotype and predicted expression
will be saved in *results/association.txt*.

## Sample data

### Transcriptome models

You can download sqlite database models from here: [PredictDB Link](https://github.com/hakyimlab/PrediXcan/blob/master/README.md#predictdb)

### Dosage File Format

These files are expected to be compressed files in the following format:
```
# chr1.dosage.txt.gz sample content, made up of the following columns:
# snpid rsid position Allele1 tAllele2 tAlleleFrequency tdosage1 ... tdosage_n
#
...
chr1 rs10399749 55299 C T 0.164302600472813 0 0 1 ... 0
...
```

NB:
- dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2).
- it is expected that there will be one file per chromosome.


If you have access to Hae Kyung Im's public data repository, and have Amazon Web Services command line tools,
you can get dosage by executing:

``` bash
mkdir data
mkdir data/dosages
cd data/dosages
aws s3 cp s3://imlab-open/Data/1000Genomes/Transcriptome/GEUVADIS/dosagefiles-hapmap2/ . --recursive
```

#### Converting from other file formats to the dosage file format

The script [convert_plink_to_dosage.py](convert_plink_to_dosage.py) will convert data in the [binary PED format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) used by the popular [plink](http://pngu.mgh.harvard.edu/~purcell/plink/) program to the dosage file format.

To use [convert_plink_to_dosage.py](convert_plink_to_dosage.py) you will need to have the latest version of plink installed: [plink 1.90 beta](https://www.cog-genomics.org/plink2/). By default, the script expects this available as `plink2`, but the path to the binary can be specified at runtime. 

If your data is in another format, you can use plink to convert your file into the plink binary PED format first.
Plink can be used to convert data in the [oxford format](https://www.cog-genomics.org/plink2/input#oxford) output by the imputation program [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html), data in the [23andMe format](https://www.cog-genomics.org/plink2/input#23file), and [VCF or BCF files](https://www.cog-genomics.org/plink2/input#vcf). The [older version of plink (1.70)](http://pngu.mgh.harvard.edu/~purcell/plink/) can be used to convert data in the [dosage format](http://pngu.mgh.harvard.edu/~purcell/plink/dosage.shtml#format) output by the imputation programs [BEAGLE](http://faculty.washington.edu/browning/beagle/beagle.html) and [MACH](www.sph.umich.edu/csg/abecasis/mach/).

### Phenotype and Filter Files

The phenotype file is expected to be formatted with the following columns:

- Family ID 
- Individual ID
- Paternal ID
- Maternal ID
- Sex
- Phenotype


Do NOT include the column names in the header.

Example Rows:
```
NA12341 NA12341 0 0 0 NA
NA12761 NA12761 0 0 0 61720
```

The filter file tells the software which rows of data in the phenotype file to use while doing the association.

It is expected to have the following column names:

- Family ID
- Individual ID
- Filter (where 1 means include, 0 exclude.)

Example Row:
```
NA12341 NA12341 1
```

Again, do not include the column names in the header.  Currently, a filter file is required to run the software.

##Association Test

The association test works by performing a linear regression on the
phenotype and the expression level of each individual gene.  The results of each
regression are saved as a row in file *association.txt*  with the following
entries:
```
gene beta z-stat p-val se(beta)
```

##Example

A working example can be download from [here](https://s3.amazonaws.com/imlab-open/Data/PredictDB/association_working_example.tar).  Note: a 250 MB file will be downloaded by clicking on the link.

Untar this file, go to the directory, and run the following in the commandline:

>./gene_association.py --dosages dosages --dosages_prefix chr --weights DGN-WB_0.5.db --pheno GWAS.fam --filter GWASfilter.txt 

** to speed up the process the dosage files can be filtered to SNPs in HapMapSnpsCEU.list.gz.

hapmapSnpsCEU.list.gz: List of SNPs used to develop the models is [here](https://app.box.com/s/6ftz3lr5h6detnf2iwzc7soyo5szrrej "HapMap2 SNP set") downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz "HapMap2 UCSC"). 
SNPs not included on this list is not used to make predictions.
