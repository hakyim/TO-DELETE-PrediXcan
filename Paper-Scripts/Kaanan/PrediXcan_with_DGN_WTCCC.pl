#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;

####################################
# Kaanan Shah
# 2/10/14
# Cox Lab
# This script sets up all the jobs and submits to tarbell to run predixcan on DNG whole blood for wtccc data. It will create all the prediction pipelines and run prediction pipelines.
####################################
my $checkimputation = 1;
my $makeconvert = 0;
my $runconvert = 0;
my $makescripts = 0;
my $runPrediction = 0;
my $makeassoc = 0;
my $runAssociation = 0;
my $makegwas = 0; #dont run gwas at same time as predixcan...
my $rungwas = 0;

my @tissues = ("DGNWholeBlood");

my @disease = ("BD","RA","T1D","T2D","HT","CAD","CD");
my @alphas = (1,0.5);


#### check R2 quality distribution of imputed data:
foreach my $co (@disease) {
    my $scriptdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/";
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/");
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/logs/");
    if ($checkimputation == 1) {
        open (OUT, ">${scriptdir}QCpipeline${co}.txt") or die "cant open ${scriptdir}QCpipeline${co}.txt\n";
        print OUT "#!/bin/bash\n";
        print OUT "#PBS -N QC_${co}\n";
        print OUT "#PBS -S /bin/bash\n";
        print OUT "#PBS -l mem=10gb\n";
        print OUT "#PBS -o ${scriptdir}logs/QC_${co}.out\n";
        print OUT "#PBS -e ${scriptdir}logs/QC_${co}.err\n";
        print OUT "#PBS -l walltime=100:00:00\n";
        print OUT "module load R/3.1.0\n";
        print OUT "perl /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/WTCCC_imputation_QCcheck.pl ${co}\n";
        close(OUT);
        system("qsub ${scriptdir}QCpipeline${co}.txt");
    }
    
}

#### convert vcf files to dos files
foreach my $co (@disease) {
    my $scriptdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/";
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/");
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/logs/");
    foreach my $chr (1 .. 22) {
        if ($makeconvert == 1) {
            open (OUT, ">${scriptdir}convertVCFpipeline${co}${chr}.txt") or die "cant open ${scriptdir}convertVCFpipeline${co}${chr}.txt\n";
            print OUT "#!/bin/bash\n";
            print OUT "#PBS -N convert${co}_${chr}\n";
            print OUT "#PBS -S /bin/bash\n";
            print OUT "#PBS -l mem=20gb\n";
            print OUT "#PBS -o ${scriptdir}logs/convert${co}_${chr}.out\n";
            print OUT "#PBS -e ${scriptdir}logs/convert${co}_${chr}.err\n";
            print OUT "#PBS -l walltime=100:00:00\n";
            print OUT "perl /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/vcf_to_dose_hapmap2.pl ${co} ${chr}\n";
            close(OUT);
        }
        if ($runconvert == 1) {
            system("qsub ${scriptdir}convertVCFpipeline${co}${chr}.txt");
            
        }
    }
}


#### run gwas
foreach my $co (@disease) {
    my $scriptdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/";
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/");
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/logs/");
    if ($makegwas == 1) {
        #cohort
        open (OUT, ">${scriptdir}GWAS_${co}.txt") or die "cant open ${scriptdir}GWAS_${co}.txt\n";
        print OUT "#!/bin/bash\n";
        print OUT "#PBS -N GWAS_${co}\n";
        print OUT "#PBS -S /bin/bash\n";
        print OUT "#PBS -l mem=10gb\n";
        print OUT "#PBS -o ${scriptdir}logs/GWAS_${co}.out\n";
        print OUT "#PBS -e ${scriptdir}logs/GWAS_${co}.err\n";
        print OUT "#PBS -l walltime=100:00:00\n";
        print OUT "module load plink/1.07\n";
        print OUT "module load R/3.1.0\n";
        print OUT "perl /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/runGWAS2.pl ${co}\n";
        close(OUT);
    }
    if ($rungwas == 1) {
        system("qsub ${scriptdir}GWAS_${co}.txt");
    }
}

#### run predixcan
foreach my $co (@disease) {
    my $scriptdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/";
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/");
    system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/${co}/logs/");
    foreach my $tiss (@tissues) {
        foreach my $a (@alphas) {
            if ($makescripts == 1) {
                #cohort
                open (OUT, ">${scriptdir}Prediction_${co}_${tiss}_EN${a}.txt") or die "cant open ${scriptdir}Prediction_${co}_${tiss}_EN${a}.txt\n";
                print OUT "#!/bin/bash\n";
                print OUT "#PBS -N ${co}_${tiss}_${a}\n";
                print OUT "#PBS -S /bin/bash\n";
                print OUT "#PBS -l mem=20gb\n";
                print OUT "#PBS -o ${scriptdir}logs/${co}_${tiss}_${a}_Predict.out\n";
                print OUT "#PBS -e ${scriptdir}logs/${co}_${tiss}_${a}_Predict.err\n";
                print OUT "#PBS -l walltime=100:00:00\n";
                print OUT "perl /group/im-lab/nas40t2/kaanan/PrediXcan/scripts/runPrediXcan3.pl /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/DGN-WB_elasticNet_alpha${a}_weights_all_chr1-22_2015-02-02.txt /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/ ${co} /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/PredictedExpression_${co}_${tiss}_EN${a}.txt";
                close(OUT);
            }
            if ($makeassoc == 1) {
                #Association
                open (OUT, ">${scriptdir}Association_${co}_${tiss}_EN${a}.R") or die "cant open ${scriptdir}Association_${co}_${tiss}_EN${a}.R\n";
                print OUT "source(\"/group/im-lab/nas40t2/kaanan/PrediXcan/scripts/prediXcanAssociation_jointimpute.r\")\n";
                print OUT "predixcanassociation(\"${co}\",\"${tiss}\",\"${a}\")\n";
                close(OUT);
                
                open (OUT, ">${scriptdir}Association_${co}_${tiss}_EN${a}.txt") or die "cant open ${scriptdir}Association_${co}_${tiss}_EN${a}.txt\n";
                print OUT "#!/bin/bash\n";
                print OUT "#PBS -N ${co}_${tiss}_${a}_Assoc\n";
                print OUT "#PBS -S /bin/bash\n";
                print OUT "#PBS -l mem=20gb\n";
                print OUT "#PBS -o ${scriptdir}logs/Assoc_${co}_${tiss}_${a}.out\n";
                print OUT "#PBS -e ${scriptdir}logs/Assoc_${co}_${tiss}_${a}.err\n";
                print OUT "#PBS -l walltime=100:00:00\n";
                print OUT "module load R/3.1.0\n";
                print OUT "R --no-save < ${scriptdir}Association_${co}_${tiss}_EN${a}.R\n";
                # print OUT "gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/GTExTissues/PredictedExpression_${co}_${tiss}_EN${a}.txt\n";
                #print OUT "gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/GTExTissues/PredictedExpression_58C_${tiss}_EN${a}.txt\n";
                #print OUT "gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/GTExTissues/PredictedExpression_NBS_${tiss}_EN${a}.txt\n";
                #print OUT "gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/GTExTissues/PrediXcan_${co}_${tiss}_EN${a}.txt\n";
                close(OUT);
            }
            if ($runPrediction == 1) {
                system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${co}/PrediXcan/");
                system("qsub ${scriptdir}Prediction_${co}_${tiss}_EN${a}.txt");
                #  system("qsub ${scriptdir}Prediction_58C_${tiss}_EN${a}.txt");
                #system("qsub ${scriptdir}Prediction_NBS_${tiss}_EN${a}.txt");
            }
            if ($runAssociation ==1) {
                system("qsub ${scriptdir}Association_${co}_${tiss}_EN${a}.txt");
            }
        }
    }
}
