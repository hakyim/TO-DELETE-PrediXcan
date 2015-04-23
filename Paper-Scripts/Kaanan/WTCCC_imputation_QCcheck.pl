#!/usr/bin/perl
use strict;
use warnings;


if (scalar(@ARGV)< 1) {print "I need a cohort to run"; die;}

foreach my $cohort (@ARGV) {
    my $vcfdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/032415_imputation/results/";
    my $outdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/QCcheck/";
    system("mkdir $outdir");
    foreach my $chr (1 .. 22) {
        #system("cp -p /group/im-lab/nas40t2/jung/${cohort}/results/chr${chr}.info.gz $vcfdir");
        system("gunzip ${vcfdir}chr${chr}.info.gz");
        open (IN, "${vcfdir}chr${chr}.info") or die "cant open ${vcfdir}chr${chr}.info\n";
        open (OUT, ">${vcfdir}chr${chr}.info.cleaned") or die "cant open ${vcfdir}chr${chr}.info.cleaned\n";
        while (my $line = <IN>) {
            chomp($line);
            my @tmp = split(/\s+/,$line);
            if (scalar(@tmp) != 14) {next;}
            print OUT "$line\n";
        }
        close(IN); close(OUT);
        #system("rm -rf ${vcfdir}chr${chr}.info");
    }
    
    
    open (R, ">${cohort}runR.R") or die "cnat make runR.R\n";
    # print R "alldat<-NULL\n";
    print R "for (i in 1:22) {\n";
    print R "dat<-read.table(paste(\"${vcfdir}chr\",i,\".info.cleaned\",sep=\"\"),header=T)\n";
    print R "file = paste(\"${outdir}${cohort}\",i,\"Rsq_hists.pdf\",sep=\"\")\n";
    print R "pdf(file)\n";
    print R "par(mfrow=c(2,2))\n";
    print R "hist(dat\$Rsq,main=\"All SNPs\",xlab=\"Rsq\")\n";
    print R "abline(v=0.8,col=\"red\")\n";
    print R "hist(dat\$Rsq[dat\$MAF<0.05],main=\"MAF < 0.05\",xlab=\"Rsq\")\n";
    print R "abline(v=0.8,col=\"red\")\n";
    print R "hist(dat\$Rsq[dat\$MAF>=0.05],main=\"MAF >= 0.05\",xlab=\"Rsq\")\n";
    print R "abline(v=0.8,col=\"red\")\n";
    print R "hist(dat\$Rsq[dat\$MAF>=0.2],main=\"MAF >= 0.2\",xlab=\"Rsq\")\n";
    print R "abline(v=0.8,col=\"red\")\n";
    print R "dev.off()\n";
    #print R "alldat<-rbind(alldat,dat)\n";
    print R "}\n";
    #print R "file = paste(\"${outdir}${cohort}allsnpsRsq_hists.pdf\",sep=\"\")\n";
    #print R "pdf(file)\n";
    #print R "par(mfrow=c(2,2))\n";
    #print R "hist(alldat\$Rsq,main=\"All SNPs\",xlab=\"Rsq\")\n";
    #print R "abline(v=0.8,col=\"red\")\n";
    #print R "hist(alldat\$Rsq[alldat\$MAF<0.05],main=\"MAF < 0.05\",xlab=\"Rsq\")\n";
    #print R "abline(v=0.8,col=\"red\")\n";
    #print R "hist(alldat\$Rsq[alldat\$MAF>=0.05],main=\"MAF >= 0.05\",xlab=\"Rsq\")\n";
    #print R "abline(v=0.8,col=\"red\")\n";
    #print R "hist(alldat\$Rsq[alldat\$MAF>=0.2],main=\"MAF >= 0.2\",xlab=\"Rsq\")\n";
    #print R "abline(v=0.8,col=\"red\")\n";
    #print R "dev.off()\n";
    close(R);
    system("R --no-save < ${cohort}runR.R > ${cohort}runR.out");
    system("rm -rf ${cohort}runR.R");
}
