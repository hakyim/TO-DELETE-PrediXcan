## Some QC on samples. 
## Eric R. Gamazon
##

## Retrive the samples from the genotype files (formatted in MACH)
##
##   Generated from the following: 
##
##	for i in `seq 1 22`;
##	do
##       	zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr$i.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr$i.samples
##               echo "zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr$i.mldose.gz | cut -d \" \" -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . \"\n\" ' > tmp/chr$i.samples"
##     	done
##

zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr22.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr22.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr21.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr21.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr20.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr20.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr19.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr19.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr18.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr18.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr17.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr17.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr16.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr16.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr15.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr15.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr14.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr14.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr13.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr13.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr12.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr12.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr11.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr11.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr10.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr10.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr9.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr9.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr8.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr8.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr7.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr7.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr6.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr6.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr5.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr5.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr4.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr4.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr3.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr3.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr2.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr2.samples 
zcat DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr1.mldose.gz | cut -d " " -f 1 | perl -ane '@data = split /->/, $F[0]; print $data[0] . "\n" ' > tmp/chr1.samples 


## Compare samples from expression file and genotype files
##
for i in `seq 1 22`;
do
	diff DGN-WB.exp.ID.list tmp/chr$i.samples >> tmp/diff_all_chroms 
	echo "DGN-WB.exp.ID.list tmp/chr$i.samples"
done    
 
wc -l tmp/diff_all_chroms       



