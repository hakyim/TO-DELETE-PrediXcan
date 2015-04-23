#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;

####This perl script takes imputed vcf file as input, removes ambiguous-strand SNPs (A/T and C/G)
#### and makes several output files for each autosome for future parallel computing:

### VCF HEADER ###
##fileformat=VCFv4.1
##filedate=2014.11.1
##source=Minimac3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+P(1/1)]">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
#cols:
#0 CHROM
#1 POS
#2 ID
#3 REF
#4 ALT
#5 QUAL
#6 FILTER
#7 INFO
#8 FORMAT
#9 FAM_BD_WTCCC170629 ...

if (scalar(@ARGV)< 2) {print "I need a cohort and chromosome to run"; die;}
my $cohort = $ARGV[0];
my $chr = $ARGV[1];

my $rthresh = 0.8; # only keep snps with Rsq > 0.8

#get list of only hapmap2 snps to keep
my $hapmapsnplist = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list";
my %HAPMAP = ();
open (HAP, "$hapmapsnplist") or die "cant open $hapmapsnplist\n";
while (my $line = <HAP>) {
    chomp($line);
    $HAPMAP{$line} = 1;
}
close(HAP);

#get position for each hapmap SNP
my %RS = ();  #$RS{chr}{pos}
my $positionfile = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/ALL.chr${chr}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.all.noSingleton.RECORDS";
open (POSS, "$positionfile") or die "cant open $positionfile\n";
while (my $line = <POSS>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if (!exists($HAPMAP{$tmp[2]})) {next;}
    $RS{$tmp[0]}{$tmp[1]} = $tmp[2];
}
close(POSS);

# pull out 1kgp SNP names from 1kgp phase 1 variants downloaded from minimac help pages
#my %RS = ();  #$RS{chr}{pos}
#foreach my $chr (1 .. 22) {
#my $snpfile = "/group/im-lab/nas40t2/kaanan/PrediXcan/1kgp/chr${chr}.map";
#system("gunzip ${snpfile}.gz");
#open (KGP, "$snpfile") or die "cant open $snpfile\n";
#while (my $line = <KGP>) {
#    chomp($line);
#    my @tmp = split(/\s+/,$line);
#    if ($tmp[0] eq "CHROM") {next;} #skip header line
#    $RS{$tmp[0]}{$tmp[1]} = $tmp[2];
#}
#close(KGP);
#system("gzip $snpfile");
#}

# convert vcfs for all wtccc datasets
#foreach my $cohort (@ARGV) {
my $vcfdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/032415_imputation/results/";
my $outdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/";
system("mkdir $outdir");
# my @variants = ();
#my %DOSE = ();
#my %FILTER = ();
#my %CHR = ();
#my %POS = ();
#my %MAF = ();


#foreach my $chr (22) {
#if ($cohort ne "T2D") {
#    system("cp -p /group/im-lab/nas40t2/jung/${cohort}/Imputation/results/chr${chr}.dose.vcf.gz $vcfdir");
#} else {
#    system("cp -p /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/RA/032415_imputation/results/chr${chr}.dose.vcf.gz $vcfdir");
#}
if ($chr == 22) {open (SAMPLES, ">${outdir}${cohort}samples.txt") or die "cant open ${outdir}${cohort}samples.txt";}
my @cols = ();
my $infile = "${vcfdir}chr${chr}.dose.vcf.gz";
my $outfile = "${outdir}${cohort}_chr${chr}.hapmap2.dos";
#system("gunzip ${infile}.gz");
#open (IN, "$infile") or die "cant open $infile \n";
open (OUT, ">$outfile") or die "cant open $outfile\n";
my $gz = gzopen($infile, "rb") or die "Error reading $infile\n";
while ($gz->gzreadline(my $line) > 0) {
    #while (my $line = <IN>) {
    chomp($line);
    if ($line =~ /^##/) {next;} # skip vcf header
        my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "#CHROM") {
        foreach my $ind (0 .. $#tmp) {
            push(@cols,$ind);
            if ($chr == 22 && $ind>8) {print SAMPLES "$tmp[$ind]\n";}
        }
        next;
    }
    #push(@variants,$tmp[2]);
    ### check rsq for variant, and potentially skip
    my @tmpinfo = split(/;/,$tmp[7]);
    my @tmpinfo2 = split(/=/,$tmpinfo[1]);
    if ($tmpinfo2[1] <= 0.8) {next;}
    my $snp = $tmp[2];
    if (!exists($RS{$tmp[0]}{$tmp[1]})) {next;}
    $snp = $RS{$tmp[0]}{$tmp[1]};
    #    }
    #$CHR{$tmp[2]} = $tmp[0];
    #$POS{$tmp[2]} = $tmp[1];
    #$REF{$tmp[2]} = $tmp[3];
    #$ALT{$tmp[2]} = $tmp[4];
    #$FILTER{$tmp[2]} = $tmp[6];
    my @tmpmaf1 = split(/;/,$tmp[7]);
    my @tmpmaf2 = split(/=/,$tmpmaf1[0]);
    #$MAF{$tmp[2]} = $tmpmaf2[1];
    print OUT "$snp $snp $tmp[1] $tmp[3] $tmp[4] $tmpmaf2[1]";
    foreach my $ind (9 .. $#tmp) {
        my @dosetmp = split(/:/,$tmp[$ind]);
        #$DOSE{$tmp[2]}{$cols[$ind]} = $dosetmp[1];
        print OUT " $dosetmp[1]";
    }
    print OUT "\n";
}
if ($chr == 22) {close(SAMPLES)};
close(OUT);
#close(IN);
$gz->gzclose() ;
#system("gzip ${infile}");
#system("rm -rf $infile");
system("gzip $outfile");
# }
#}

