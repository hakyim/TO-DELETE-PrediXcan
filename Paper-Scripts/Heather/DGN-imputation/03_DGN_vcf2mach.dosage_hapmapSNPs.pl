#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw[min max];

####This perl script takes the DGN imputed vcf files (MAF>0.05) as input, removes ambiguous-strand SNPs (A/T and C/G), removes non-HapMap2 SNPs,
#### removes SNPs with R2<0.8, finds the rsID for each SNP, and makes several output files for each autosome for future parallel computing:
#### .mlinfo.gz and .mldose.gz MACH files for GCTA
#### .SNPxID matrix for quick scanning into R
#### .bim plink bim file with SNP pos info in .SNPxID
#### .ID.list colnames of matrix
#### .SNP.list rownames of matrix


my $dir = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/UMich-imputation-results/results/";
my $refdir = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/";


open(HAP, "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list"); ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz

my %hapmapsnps;
while(<HAP>){
    chomp;
    my ($snp) = split(/\n/);
    $hapmapsnps{$snp} = 1;
}

#parse by chr
for(my $i = 1; $i <= 22; $i++){
    print "$i\n";
    my $snpxidhandle = "SNPxID" . $i;
    my $mlinfohandle = "MLINFO" . $i;
    my $bimhandle = "BIM" . $i;
    open($snpxidhandle, ">DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.SNPxID");
    open($bimhandle, ">DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.bim");
    open($mlinfohandle, ">DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.mlinfo");
    open(ID, ">DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.ID.list");
    open(SNP, ">DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.SNP.list");

    open(REF, "${refdir}ALL.chr${i}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.all.noSingleton.RECORDS");

    my %rsid;
    while(<REF>){
	chomp;
	my ($chrom, $pos, $rs) = split(/\t/);
	$rsid{$pos} = $rs;
    }

    open(VCF, "${dir}chr${i}.dose_maf0.05_rm.indel.recode.vcf");
    open(INTRO, ">intro");
    while(<VCF>){
	chomp;
	my ($first) = split(/\t/);
	if($first =~ m/##/){ ##skip vcf info lines
	    next;
	}
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
	my ($expfreq, $impr2, $imper2) = split(/;/,$info);
	my ($b, $rsq) = split(/=/,$impr2);
	my $quality = $rsq; ##only one quality score per SNP in minimac, output twice in mach files 
	my $rs = $rsid{$pos}; ##get rsid from .RECORDS hash

	if($chr eq "#CHROM"){
	    my $mlinfohandle = "MLINFO" . $i;
	    print $mlinfohandle "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n";
	
	    foreach my $id (@genos){
		my ($shortid,$e) = split(/_/,$id);
		print ID "$shortid\n";
		print INTRO "$shortid->$shortid MLDOSE\n";
	    }
	}
	if($ref eq "A" && $alt eq "T"){ ##rm potentially ambiguous strand SNPs
	    next;
	}elsif($ref eq "T" && $alt eq "A"){
	    next;
	}elsif($ref eq "C" && $alt eq "G"){
	    next;
	}elsif($ref eq "G" && $alt eq "C"){
	    next;
	}elsif(defined($hapmapsnps{$rs}) && $rsq >= 0.8 && $pos =~ m/\d+/){ ###only pull hapmap SNPs with R2>0.8 & don't print header rows
	    my $snpxidhandle = "SNPxID" . $chr;		
	    my $bimhandle = "BIM" . $chr;
	    my $mlinfohandle = "MLINFO" . $chr;
	    my $sum = 0; #for $freqalt calc
	    foreach my $geno (@genos){
		my ($gt, $dos) = split(/:/,$geno);
		print $snpxidhandle "$dos\t";
		$sum = $sum + $dos; #add up the dosages to calc ALT allele freq
	    }
	    my $n = @genos; #n samples
	    my $freqalt = $sum/($n*2); #calc ALT allele freq (I found that ALT is not always the minor allele)
	    my $freqref = 1 - $freqalt;
	    my $maf = min($freqref,$freqalt);
	    $maf = sprintf("%.5f", $maf); #round to 5 places
	    $freqref = sprintf("%.5f", $freqref);
	    print SNP "$rs\n";
	    print $bimhandle "$chr\t$rs\t0\t$pos\t$ref\t$alt\n";
            print $mlinfohandle "$rs\t$ref\t$alt\t$freqref\t$maf\t$quality\t$rsq\n";
	    print $snpxidhandle "\n";
	}
    }
}


for(my $i = 1; $i <= 22; $i++){
    print "$i\n";
    open(R, ">runR.R") or die "cant mae runR.R\n";
    print R "dat<-scan(\"DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.SNPxID\")\n";
    print R "gtidlist<-scan(\"DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.ID.list\",\"character\")\n";
    print R "dat<-matrix(dat, ncol=length(gtidlist), byrow=T)\n";
    print R "dat<-t(dat)\n";

    print R "write.table(dat,\"t.dos.chr${i}\",col.names=F,row.names=F,quote=F)\n";
    close(R);
    system("R --vanilla < runR.R");
    system("paste -d\' \' intro t.dos.chr${i} > DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.mldose");
}

system("gzip DGN-imputed-for-PrediXcan/*.mldose");
system("gzip DGN-imputed-for-PrediXcan/*.mlinfo");
system("rm intro t.dos.chr* runR.R");

