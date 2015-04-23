#!/usr/bin/perl
use strict;
use warnings;

############################################################
# script to pull out meta analysis data for predictor snps #
# from predixcan and search for top genes                  #
############################################################

my $betafile = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/DGN-WB_elasticNet_alpha0.5_weights_all_chr1-22_2015-02-02.txt";
my %BETASNPS = ();
open (BETA, "$betafile") or die "cant open $betafile\n";
while (my $line = <BETA>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    $BETASNPS{$tmp[1]} = 1;
}
close(BETA);

## BD
## PGC BPD meta analysis
my $outfilebd = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/BDmeta_beta_overlap.txt";
open (OUT, ">$outfilebd") or die "cant open $outfilebd\n";
my $bdfile = "/group/im-lab/nas40t2/haky/Signatures/rawdata/PGC/pgc.bip.full.2012-04.txt";
open (BD, "$bdfile") or die "cant open $bdfile\n";
while (my $line = <BD>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    #    if ($tmp[0] eq "SNPID") {print OUT "$line\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$line\n";}
}
close(BD);
close(OUT);

## HT
#rsid,chr.hg18,pos.hg18,pval.GC.SBP,pval.GC.DBP
my $outfileht = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/HTmeta_beta_overlap.txt";
open (OUT, ">$outfileht") or die "cant open $outfileht\n";
my $htfile = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ICBP-summary-Nature.csv";
open (HT, "$htfile") or die "cant open $htfile\n";
while (my $line = <HT>) {
    chomp($line);
    my @tmp = split(/,/,$line);
    if ($tmp[0] eq "rsid") {print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\n";}
}
close(HT);
close(OUT);

## RA
my $outfilera = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/RAmeta_beta_overlap.txt";
open (OUT, ">$outfilera") or die "cant open $outfilera\n";
my $rafile = "/group/im-lab/nas40t2/haky/Signatures/rawdata/RA-2014/RA_GWASmeta_European_v2.txt";
open (RA, "$rafile") or die "cant open $rafile\n";
while (my $line = <RA>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "SNPID") {print OUT "$line\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$line\n";}
}
close(RA);
close(OUT);

## CD
my $outfilecd = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/CDmeta_beta_overlap.txt";
open (OUT, ">$outfilecd") or die "cant open $outfilecd\n";
my $cdfile = "/group/im-lab/nas40t2/haky/Signatures/rawdata/IBD/cd-meta.txt";
system("cp -p ${cdfile}.gz /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/");
system("gunzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/cd-meta.txt.gz");
open (CD, "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/cd-meta.txt") or die "cant open $cdfile\n";
while (my $line = <CD>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "SNP") {print OUT "$line\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$line\n";}
}
close(CD);
close(OUT);
system("rm -rf /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/cd-meta.txt");

## T1D
my $outfilet1d = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/T1Dmeta_beta_overlap.txt";
open (OUT, ">$outfilet1d") or die "cant open $outfilet1d\n";
my $t1dfile1 = "/group/im-lab/nas40t2/haky/Signatures/rawdata/T1D-CHOP/T1D-CHOP-Illumina-Plink.txt";
open (T1D, "$t1dfile1") or die "cant open $t1dfile1\n";
while (my $line = <T1D>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[2] eq "SNP") {print OUT "Illumina\t$line\n";}
    if (exists($BETASNPS{$tmp[2]})) {print OUT "Illumina\t$line\n";}
}
close(T1D);
my $t1dfile2 = "/group/im-lab/nas40t2/haky/Signatures/rawdata/T1D-CHOP/T1D-CHOP-Affy-Plink.txt";
open (T1D, "$t1dfile2") or die "cant open $t1dfile2\n";
while (my $line = <T1D>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if (exists($BETASNPS{$tmp[2]})) {print OUT "Affy\t$line\n";}
}
close(T1D);
close(OUT);

## T2D
my $outfilet2d = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/T2Dmeta_beta_overlap.txt";
open (OUT, ">$outfilet2d") or die "cant open $outfilet2d\n";
my $t2dfile = "/group/im-lab/nas40t2/haky/Signatures/rawdata/DIAGRAM/DIAGRAMv3.2012DEC17.txt";
open (T2D, "$t2dfile") or die "cant open $t2dfile\n";
while (my $line = <T2D>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "SNP") {print OUT "$line\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$line\n";}
}
close(T2D);
close(OUT);

## CAD
my $outfilecad = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/CADmeta_beta_overlap.txt";
open (OUT, ">$outfilecad") or die "cant open $outfilecad\n";
my $cadfile = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/CARDIoGRAM_GWAS_RESULTS.txt";
open (CAD, "$cadfile") or die "cant open $cadfile\n";
while (my $line = <CAD>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "SNP") {print OUT "$line\n";}
    if (exists($BETASNPS{$tmp[0]})) {print OUT "$line\n";}
}
close(CAD);
close(OUT);
