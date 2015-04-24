#!/usr/bin/perl
use strict;
use warnings;


if (scalar(@ARGV) != 1) {print "I need a single cohort to run"; die;}
#my @junk = split(/\s+/,@ARGV);
my $cohort = $ARGV[0];
my $dogwas = 1;
my $makeplot = 1;
system("mkdir /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/GWAS/");
my $outdir = "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/GWAS/";

if ($dogwas == 1) {
    open (LST, ">${outdir}${cohort}.lst") or die "cant open ${outdir}${cohort}.lst\n";
    open (FAM, ">${outdir}${cohort}.fam") or die "cant open ${outdir}${cohort}.fam\n";
    foreach my $set ("${cohort}") { #,"NBS","58C") {
        open (IN, "/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}samples.txt") or die "cant open ${cohort}samples.txt\n";
        while (my $line = <IN>) {
            chomp($line);
            print FAM "$line\t$line\t0\t0\t0\t";
            if ($line =~ /^FAM_${cohort}/)  {print FAM "1\n";} else {print FAM "0\n";}
            print LST "$line\t$line\n";
        }
        close(IN);
    }
    close(FAM);
    close(LST);
   
    ## make list of individuals to remove
    open (OUT, ">/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/GWAS/removeindividuals.txt") or die "cant open removeindividuals.txt\n";
    #cohort
    open (IN, "/group/im-lab/nas40t2/haky/Signatures/data/cohorts/WTCCC1/exclusion-list-05-02-2007-${cohort}.txt") or die "cant open exclusion-list-05-02-2007-${cohort}.txt\n";
    while (my $line = <IN>) {
        chomp($line);
        if ($line =~ /^#/) {next;}
            my @tmp = split(/\s+/,$line);
        print OUT "FAM_${cohort}_$tmp[2]\tFAM_${cohort}_$tmp[2]\n";
    }
    close(IN);
    #58C
    open (IN, "/group/im-lab/nas40t2/haky/Signatures/data/cohorts/WTCCC1/exclusion-list-05-02-2007-58C.txt") or die "cant open exclusion-list-05-02-2007-58C.txt\n";
    while (my $line = <IN>) {
        chomp($line);
        if ($line =~ /^#/) {next;}
            my @tmp = split(/\s+/,$line);
        print OUT "FAM_58C_$tmp[2]\tFAM_58C_$tmp[2]\n";
    }
    close(IN);
    #NBS
    open (IN, "/group/im-lab/nas40t2/haky/Signatures/data/cohorts/WTCCC1/exclusion-list-NBS.txt") or die "cant open exclusion-list-NBS.txt\n";
    while (my $line = <IN>) {
        chomp($line);
        if ($line =~ /^#/) {next;}
            my @tmp = split(/\s+/,$line);
        print OUT "FAM_NBS_$tmp[2]\tFAM_NBS_$tmp[2]\n";
    }
    close(IN);
    close(OUT);
    
    foreach my $chr (1 .. 22) { #1 .. 22
        system("gunzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}_chr${chr}.hapmap2.dos.gz");
        
        ## make list files
        open (LIST, ">${outdir}${cohort}files${chr}.list") or die "cant make ${outdir}${cohort}files.list\n";
        print LIST "1 /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}_chr${chr}.hapmap2.dos ${outdir}${cohort}.lst\n";
        close(LIST);
        system("plink --fam ${outdir}${cohort}.fam --dosage ${outdir}${cohort}files${chr}.list list sepheader format=1 skip1=2 skip2=1 --out ${outdir}${cohort}_chr${chr} --allow-no-sex --1 --remove /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/GWAS/removeindividuals.txt");
        system("gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}_chr${chr}.hapmap2.dos.gz");
        
    }
}

if ($makeplot == 1) {
    open (R, ">${outdir}runR.R") or die "cant make runR.R\n";
    print R "source(\"/home/kshah2/qqplots.R\")\n";
    print R "bad<-read.table(\"/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/58C/GWAS/badsnps_controls\",header=F)\n";
    print R "bad\$bad<-1\n";
    print R "bad<-bad[,c(1,9)]\n";
    
    print R "dat<-NULL\n";
    foreach my $i (1 .. 22) {
        system("zcat /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}_chr${i}.hapmap2.dos.gz | cut -f 2-3 -d \" \" > ${outdir}tmp.txt");
        system("awk '{print \"$i \" \$0;}' ${outdir}tmp.txt >> ${outdir}${cohort}snps.map");
        system("gunzip ${outdir}${cohort}_chr${i}.assoc.dosage.gz");
        print R "tmp<-read.table(\"${outdir}${cohort}_chr${i}.assoc.dosage\", header=T)\n";
        print R "tmp<-tmp[!is.na(tmp\$P),]\n";
        print R "tmp<-merge(tmp,bad,by.x=1,by.y=1,all.x=T,all.y=F)\n";
        print R "tmp<-tmp[is.na(tmp\$bad),]\n";
        print R "dat<-rbind(dat,tmp[,c(1,8)])\n";
        #print R "write.table(tmp[tmp\$P<=10^-5,],\"${outdir}badsnps_controls\",col.names=F,append=T,row.names=F,quote=F)\n";
    }
    print R "map<-read.table(\"${outdir}${cohort}snps.map\",header=F)\n";
    print R "colnames(map)<-c(\"chr\",\"SNP\",\"pos\")\n";
    print R "dat<-merge(dat,map,by.x=1,by.y=2,all.x=T,all.y=F)\n";
    print R "maxs<-NULL\n";
    print R "currentadd<-0\n";
    print R "jpeg(\"${outdir}${cohort}_manhanttanplot.jpeg\")\n";
    print R "plot(0,0,xlim=c(0,3*10^9),ylim=c(0,10),main=\"${cohort} vs 58C+NBS\",ylab=\"-log10(p-value)\", xlab=\"Genomic Position\",type=\"n\",xaxt=\"n\")\n";
    print R "for (i in 1:22) {\n";
    print R "tmp<-dat[dat\$chr == i,]\n";
    print R "tmp\$pos<-tmp\$pos+currentadd\n";
    print R "maxs<-c(maxs,max(tmp\$pos))\n";
    print R "currentadd<-max(tmp\$pos)\n";
    print R "if (i %% 2 == 0) {color<-\"grey65\"} else {color<-\"black\"}\n";
    print R "points(tmp\$pos,-log10(tmp\$P),pch=20,cex=0.5,col=color)\n";
    print R "}\n";
    print R "axis(side=1,at=maxs,labels=1:22)\n";
    print R "abline(h=5,col=\"blue\",lty=\"dashed\")\n";
    print R "abline(h=-log10(5*10^-8),col=\"red\")\n";
    print R "dev.off()\n";
    print R "jpeg(\"${outdir}${cohort}_qqlot.jpeg\")\n";
    print R "qqunif(dat\$P,plot=T)\n";
    print R "dev.off()\n";
    close(R);
    system("R --no-save < ${outdir}runR.R > ${outdir}runR.out");
    foreach my $i (1 .. 22) {
        system("gzip ${outdir}${cohort}_chr${i}.assoc.dosage");
        system("gzip /group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/${cohort}/${cohort}_chr${i}.hapmap2.dos.gz");
    }
}