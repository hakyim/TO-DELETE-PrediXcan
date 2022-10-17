#!/bin/bash

###Note on 20150126
## I originally ran the pre-GCTA QC on genegate interactively. I have attempted to update the paths to what is now on tarbell
## See '/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/run_01_gcta_QC.sh' for genegate paths

##load software used on tarbell
module load R
module load plink/1.09 ##plink2
module load vcftools
module load tabix/0.2.6

$DIR=/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/

####pull hapmap2 SNPs
plink --bfile /group/im-lab/hwheeler/DGN-genotypes/GenRED.II.autosomal.Final --extract /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/gEUVADIS_genotypes/hapmap2.CEU.SNP.list --make-bed --out ${DIR}DGN.hapmap2

####Check sex
plink --bfile /group/im-lab/hwheeler/DGN-genotypes/GenRED.II.X.final --check-sex --out ${DIR}GenRED.II.X.final
####results
## no problems detected

####keep only individuals with expression data
perl ${DIR}match_ids.pl
plink --bfile ${DIR}DGN.hapmap2 --keep ${DIR}DGN-WB.exp.FID.IID --make-bed --out ${DIR}DGN.hapmap2

####Check call rates and rm poorly called SNPs
plink --bfile ${DIR}DGN.hapmap2 --missing --out ${DIR}DGN.hapmap2
plink --bfile ${DIR}DGN.hapmap2 --geno 0.01 --make-bed --out ${DIR}tmp.geno
plink --bfile ${DIR}tmp.geno --mind 0.01 --make-bed --out ${DIR}tmp.geno.mind
## 4739 variants removed due to missing genotype data (--geno).
##  0 people removed due to missing genotype data (--mind).

####Calculate HWE statistics and rm SNPs with P < 0.05, strict for GCTA
plink --bfile ${DIR}tmp.geno.mind --hwe 0.05 --make-bed --out ${DIR}tmp.hwe
##  --hwe: 27940 variants removed due to Hardy-Weinberg exact test.

####Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F < or > 3 sd from mean)
plink --bfile ${DIR}tmp.hwe --het --out ${DIR}DGN.hapmap2
R --vanilla < pull_het_outliers.r --args ${DIR}DGN.hapmap2.het
##  DGN.hapmap2.het.pdf plot looks good, don't remove anyone


####Relationship check -- rm one of pair with pi_hat > 0.05, rm the one with no expression if possible
plink --bfile ${DIR}tmp.hwe --genome --min 0.05 --maf 0.01 --out ${DIR}DGN.hapmap2.pihat.0.05
##  no relateds to remove in DGN.hapmap2.pihat.0.05.genome

####include snps with MAF>0.01
plink --bfile ${DIR}tmp.hwe --maf 0.01 --make-bed --out ${DIR}tmp.maf
##  2225 variants removed due to MAF threshold

####Remove ambiguous A/T or C/G SNPs
perl /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/8_pull_unamb_SNPs.pl ${DIR}tmp.maf.bim > ${DIR}tmp.snplist
plink --bfile ${DIR}tmp.maf --extract ${DIR}tmp.snplist --make-bed --out ${DIR}DGN.hapmap2.chr1-22.QC
##  649515 variants and 922 people pass filters and QC.
##  shorten FID to match exp data
cut -d'_' -f 3- ${DIR}DGN.hapmap2.chr1-22.QC.fam >o
mv o ${DIR}DGN.hapmap2.chr1-22.QC.fam

rm ${DIR}tmp*

##Final pre-impute QC plink files: ${DIR}DGN.hapmap2.chr1-22.QC.*
