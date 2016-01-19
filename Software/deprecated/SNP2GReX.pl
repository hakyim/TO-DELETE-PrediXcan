#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;

####################################
# Kaanan Shah
# 2/10/14
# Cox Lab
# This script will  predict gene expression values using betas from a weight file and snp dosages for specific cohort of people.
#
# The script requires 4 arguments:
# 1. The weight file in the following format:
# gene	SNP	refAllele	effectAllele	beta
# A2LD1	rs10508048	C	T	 2.467627e-03
# A2LD1	rs1055705	G	A	 3.480362e-03
# The columns must be in the exact order!
# the effect allele must be listed 2nd
#
# 2. The director of the dosages. Each gzipped dosage file should follow the naming structure: [cohort]_chr[chr#].dos.gz and have the following format:
# columns are snp snp position allele1 allele2 MAF dosage.....dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2).
# in the example below, the 4th individual is likely homozygous TT for both SNPs.
# rs56934710 rs56934710 96736 T C 0.10693 0.933 0.934 0.023 0.024 ...
# rs72710816 rs72710816 97562 T C 0.10576 0.936 0.936 0.021 0.022 ...
# There should be 1 file for each chromosome (1-22) and a sample file ([cohort]sample.txt that has the individual IDs 1 per line (no header) in the order they appear in the dosage files.
# The columns must be in the exact order!
# The dosage allele must be listed 2nd.
#
# 3. the prefix for the cohort being tested
#
# 4. the name of the output file which will include predicted epxression for all genes and all indiviausl : 1 row per gene, 1 column per individual.
####################################


if (scalar(@ARGV) != 4) {print "I need more information to run\n runPrediXcan.pl [betafile] [genotypedirectory] [cohort] [outfile]\nThe script requires 4 arguments:\n1. The weight file in the following format:\n     gene	SNP	refAllele	effectAllele	beta\n     A2LD1	rs10508048	C	T	 2.467627e-03\n     A2LD1	rs1055705	G	A	 3.480362e-03\n     The columns must be in the exact order!\n     the effect allele must be listed 2nd\n\n    2. The director of the dosages. Each gzipped dosage file should follow the naming structure: [cohort]_chr[chr#].dos.gz and have the following format:\n     columns are snp snp position allele1 allele2 MAF dosage.....dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2).\n     in the example below, the 4th individual is likely homozygous TT for both SNPs.\n     rs56934710 rs56934710 96736 T C 0.10693 0.933 0.934 0.023 0.024 ...\n     rs72710816 rs72710816 97562 T C 0.10576 0.936 0.936 0.021 0.022 ...\n     There should be 1 file for each chromosome (1-22) and a sample file ([cohort]sample.txt that has the individual IDs 1 per line (no header) in the order they appear in the dosage files.\n     The columns must be in the exact order!\n     The dosage allele must be listed 2nd.\n\n\     3. the prefix for the cohort being tested\n\n     4. the name of the output file which will include predicted epxression for all genes and all indiviausl : 1 row per gene, 1 column per individual.\n"; die;}

my $betafile = $ARGV[0];
#"/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/DGN-WB_elasticNet_alpha0.5_weights_all_chr1-22_2015-02-02.txt";
# gene	SNP	refAllele	effectAllele	beta
# A2LD1	rs10508048	C	T	 2.467627e-03
# A2LD1	rs1055705	G	A	 3.480362e-03
# A2LD1	rs11069419	T	C	-4.681995e-04
# A2LD1	rs11842969	C	T	-2.138703e-02

my $genodir = $ARGV[1];
#"/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/BD/";
# rs56934710 rs56934710 96736 T C 0.10693 0.933 0.934 0.023 0.024 ...
# rs72710816 rs72710816 97562 T C 0.10576 0.936 0.936 0.021 0.022 ...

my $cohort = $ARGV[2];
my $outfile = $ARGV[3];

my $samplefile = "${genodir}${cohort}samples.txt"; # 1 sample ID per line, no header

## read in beta file information
my %SNPS = (); # $SNP{rs#} = 1 if snp is used for prediction
my %BETAS = (); # $BETAS{rs#}{gene} = beta value
my %REF = (); # $REF{rs#}{gene} = reference allele
my %EFF = (); # $EFF{rs#}{gene} = effect allele
open (BETA, "$betafile") or die "cant open $betafile\n";
while (my $line = <BETA>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    if ($tmp[0] eq "gene") {
        if ($tmp[3] ne "effectAllele") {print "make sure the 4th column contains the effectAllele\n"; die;}
        if ($tmp[2] ne "refAllele") {print "make sure the 3rd column contains the refAllele\n"; die;}
        next;
    }
    $SNPS{$tmp[1]} = 1;
    $BETAS{$tmp[1]}{$tmp[0]} = $tmp[4];
    $REF{$tmp[1]}{$tmp[0]} = $tmp[2];
    $EFF{$tmp[1]}{$tmp[0]} = $tmp[3];
}
close(BETA);

## sample file: store the order of sample ids by expected index in genotype dosage files
my %INDEX = ();
my $index = 5; #starting index from genotype files (first columns are snp information cols)
open (SAMPLE, "$samplefile") or die "cant open $samplefile\n";
while (my $line = <SAMPLE>) {
    chomp($line);
    $index++;
    $INDEX{$index} = $line;
    
}
close(SAMPLE);

## read dosage files and do prediction!
my %PREDICT = ();
my %predictedgenes = ();
foreach my $chr (1 .. 22) {
    my $dosefile = "${genodir}${cohort}_chr${chr}.dos.gz";
    #system("gunzip ${dosefile}.gz");
    my $gz = gzopen($dosefile, "rb") or die "Error reading $dosefile\n";
    while ($gz->gzreadline(my $line) > 0) {
        
        #    gzopen (DOS, "$dosefile") or die "cant open $dosefile\n";
        #while (my $line = <DOS>) {
        chomp($line);
        my @tmp = split(/\s+/,$line);
        if ($#tmp != $index) {die "ERROR: $dosefile: number of individuals in sample file does not match dosage file\n";}
        if (!exists($SNPS{$tmp[1]})) {next;} # skip this snp, it is not used in any predictors
        foreach my $gene (keys(%{$BETAS{$tmp[1]}})) { #go through all genes that use that snp
            if (!exists($predictedgenes{$gene})) {$predictedgenes{$gene}=1;}
            if ($tmp[4] eq $EFF{$tmp[1]}{$gene}) { #if the effect snp is the same as the alternate (dosage) allele
                foreach my $i (6 .. $index) { #for each individual, i
                    if ($tmp[$i]<0) {die "ERROR: dosages should be between 0-2. Check $INDEX{$i} for $tmp[1]\n";}
                    if ($tmp[$i]>2) {die "ERROR: dosages should be between 0-2. Check $INDEX{$i} for $tmp[1]\n";}
                    if (exists($PREDICT{$gene}{$i})) {
                        $PREDICT{$gene}{$i} = $PREDICT{$gene}{$i} + $tmp[$i]*$BETAS{$tmp[1]}{$gene};
                    } else {
                        $PREDICT{$gene}{$i} = $tmp[$i]*$BETAS{$tmp[1]}{$gene};
                    }
                }
                next;
            }
            if ($tmp[3] eq $EFF{$tmp[1]}{$gene}) {#if the effect snp is the same as the reference (non-dosage) allele... dosage needs to be flipped to match alleles.
                foreach my $i (6 .. $index) { #for each individual, i
                    if (exists($PREDICT{$gene}{$i})) {
                        $PREDICT{$gene}{$i} = $PREDICT{$gene}{$i} + (2-$tmp[$i])*$BETAS{$tmp[1]}{$gene};
                    } else {
                        $PREDICT{$gene}{$i} = (2-$tmp[$i])*$BETAS{$tmp[1]}{$gene};
                    }
                }
                next;
            }
            print "ERROR: $tmp[1] effect SNP allele does not match either allele in dosage file, this SNP was removed from prediction model. Please check!\n";
        }
    }
    #close(DOS);
    $gz->gzclose() ;
    
    #system("gzip ${dosefile}");
}
open (OUT, ">$outfile") or die "cant make $outfile\n";
#print header row
print OUT "gene";
foreach my $i (6 .. $index) {
    print OUT "\t$INDEX{$i}"
}
print OUT "\n";
#print out predicted gene expression, 1 gene per row
foreach my $gene (keys(%predictedgenes)) {
    print  OUT "$gene";
    foreach my $i (6 .. $index) {
        print OUT "\t$PREDICT{$gene}{$i}"
    }
    print OUT "\n";
}
close(OUT);
print "Prediction complete for $cohort\n";


