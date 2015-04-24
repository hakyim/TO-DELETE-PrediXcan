### Kaanan Shah
### Cox Lab
### April 8, 2015
### Figures and tables to summarize WTCCC results for Predixcan Paper
#### all summaries are for wtccc imputed data done by Eric Gamazon and DGN whole blood expression elastic net predictors with alpha = 0.5 developed by Heather Wheeler
################################################################
#local computer
setwd("Desktop/predixcanreviewstuff/FiguresandTables/")
qqunif = function(p,BH=T,CI=T,FDRthres=0.05,plot=F,...)
{
    p<-p[!is.na(p)]
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    dat<-cbind(sort(p),1:nn)
    q<-(nn*dat[,1])/dat[,2] # calculate q-values from p-values
    dat<-cbind(dat,q)
    if (sum(dat[,3]<FDRthres) == 0) {
        nsnps=0
    } else {
        nsnps<-round((sum(p<=max(dat[dat[,3]<FDRthres,1]))/nn)*100,4)
    }
    if (plot) {
        plot( xx,  -sort(log10(p)),
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        pch=20,cex.lab=1.3,cex.axis=1,las=1,cex=0.5,mgp=c(2,1,0),bty="n")
        if(CI)
        {
            ## create the confidence intervals
            c95 <- rep(0,nn)
            c05 <- rep(0,nn)
            ## the jth order statistic from a
            ## uniform(0,1) sample
            ## has a beta(j,n-j+1) distribution
            ## (Casella & Berger, 2002,
            ## 2nd edition, pg 230, Duxbury)
            ## this portion was posted by anonymous on
            ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
            for(i in 1:nn)
            {
                c95[i] <- qbeta(0.95,i,nn-i+1)
                c05[i] <- qbeta(0.05,i,nn-i+1)
            }
            polygon(c(xx, rev(xx)), c(-log10(c95), rev(-log10(c05))),col = "grey", border = NA)
            #    lines(xx,-log10(c95),col='gray')
            #    lines(xx,-log10(c05),col='gray')
        }
        points(xx,  -sort(log10(p)),pch=20,cex=0.5)
        #points(xx[dat[,3]<0.05],  -log10(dat[dat[,3]<0.05,1]),pch=20,col="pink")
 
        y<-max(-log10(p))
        #text(0,y,paste(nsnps,"% of genes have a q-value <= ",FDRthres,sep=""),pos=4)
        #  print(paste(nsnps,"% of SNPs have a q-value <= ",FDRthres,sep=""))
        abline(0,1,col='red')
        if(BH)
        {
            abline(-log10(0.05),1, col='black',lty=2)
            abline(-log10(0.10),1, col='black',lty=3)
            abline(-log10(0.25),1, col='black',lty=4)
            # legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),col=c('black','black','black'),lty=2:4, cex=1)
            abline(h=-log10(0.05/nn),col="blue") ## bonferroni
        }
    } else {
        return(nsnps)
    }
}

#### random help code ###
# example of layout
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#hist(wt)
#hist(mpg)
#hist(disp)

######################################################################
### make 1 figure for T1D with manhattan , qq, and gwas enrichment  ##
######################################################################

cohorts<-c("T1D") #"RA","BD","CD","CAD","HT","RA","T2D")
r2thres<-0.01
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)

for (co in cohorts) {
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
	file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
    dat<-dat[dat$R2>=r2thres,]
    bf<-0.05/dim(dat)[1]
#    figurefile<-paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/",co,"_PrediXcanResults_EN0.5.ps",sep="")
    figurefile<-paste(co,"_PrediXcanResults_EN0.5.ps",sep="")
    postscript(figurefile,height=6, width= 6)
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    par(mar=c(5,4,1.5,1)+0.1)
    ######################################
    ## Manhattan plots for each disease ##
    ######################################
    allgenes<-read.table("gencode.v12.V1.summary.protein.nodup.genenames",header=F)
#    allgenes<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/gencode.v12.V1.summary.protein.nodup.genenames",header=F)
    annot<-merge(allgenes,qual,by.x=6,by.y=1,all.x=F,all.y=T)
    annot$mid<-(annot$V4+annot$V3)/2
    annot$pos<-0
    start<-0
    ends<-NULL
    starts<-NULL
    annot$color<-"grey65"
    ticks<-NULL
    for (chr in 1:22) {
        starts<-c(starts,start)
        ch<-paste("chr",chr,sep="")
        if (chr%%2==1) {annot$color[annot$V1==ch]<-"black"}
        annot$pos[annot$V1==ch]<-annot$mid[annot$V1==ch]+start
        start<-max(annot$pos)+1
        ends<-c(ends,start)
        ticks<-c(ticks,mean(annot$pos[annot$V1==ch]))
    }
    data<-merge(dat,annot,by.x=1,by.y=1,all.x=T,all.y=F)
    plot(data$pos,-log10(data$p.val),pch=20,ylab=expression(-log[10](italic(p))),xlab="Chromosome",xaxt="n",main="",cex.lab=1.3,cex.axis=1,col=data$color,las=1,cex=0.5,mgp=c(2,1,0),bty="n")
    axis(side=1,at=ticks,labels=1:22,cex.lab=1)
    outer = FALSE
    line = 0.25
    cex = 2
    adj  = 0.025
    title(outer=outer,adj=adj,main="a",cex.main=cex,col="black",font=2,line=line)
    #plot(rnorm(100),col="blue")
    #title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
    #plot(rnorm(100),col="green")
    #title(outer=outer,adj=adj,main="C",cex.main=cex,col="black",font=2,line=line)
    #MakeLetter("c")
    abline(h=-log10(bf),col="blue")
    tmp<-data[data$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]
    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>3) {
            tmp<-tmp[1:3,]
        }
        text(tmp$pos,-log10(tmp$p.val),tmp$gene,pos=4,cex=0.7)
    }
    #################################################
    ##qqplots with gene names labled for top genes ##
    #################################################
    qqunif(dat$p.val,plot=T)
    nn = length(dat$p.val)
    dat<-dat[order(dat$p.val),]
    dat$xx =  -log10((1:nn)/(nn+1))
    
    # title(paste("PrediXcan Results\n",co,", EN0.5",sep=""))
    tmp<-dat[dat$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]
    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>3) {
            tmp<-tmp[1:3,]
        }
        text(tmp$xx,-log10(tmp$p.val),tmp$gene,pos=2,cex=0.7)
    }
    outer = FALSE
        line = 0.25
        cex = 2
        adj  = 0.025
        title(outer=outer,adj=adj,main="b",cex.main=cex,col="black",font=2,line=line)

    #################################################################
    ## ENRICHMENT OF KNOWN DISEASE GENES WITHIN PREDIXCAN RESULTS ###
    #################################################################
    #known<-read.table(paste("/group/im-lab/nas40t2/haky/main/Data/NHGRI-GWAS/",co,"-genes.txt",sep=""),header=F)
    known<-read.table(paste(co,"-genes.txt",sep=""),header=F)
    known$ragene<-1
    data<-merge(dat,known,by.x=1,by.y=1,all.x=T,all.y=F)
    data$ragene[is.na(data$ragene)]<-0
    a<-glm(data$ragene~data$p.val+data$n.snps+data$R2,family="binomial")
    #plot(-log10(dat$p.val),dat$R2)
    print(summary(a))
    pthres<-0.01
    nsig<-NULL
    nreps<-10000
    set.seed(110285)
    for (i in 1:nreps) {
        sum(data$ragene)
        rows<-sample(1:dim(data)[1],sum(data$ragene),replace=F)
        nsig<-c(nsig,sum(data$p.val[rows]<=pthres))
    }
    truesig<-sum(data$p.val[data$ragene==1]<=pthres)
    hist(nsig,xlim=c(0,max(c(nsig,truesig))),xlab=paste("Expected genes with p<",pthres,sep=""),main="",ylab="Frequency",cex.lab=1.3,cex.axis=1,las=1,bty="n")
    points(truesig,0,col="black",pch=20, cex=1.5)
    outer = FALSE
    line = 0.25
    cex = 2
    adj  = 0.025
    title(outer=outer,adj=adj,main="c",cex.main=cex,col="black",font=2,line=line)
    penrich<-sum(nsig>=truesig)/nreps
    print(penrich)
    dev.off()
}


    #################################################################
    ## ENRICHMENT OF KNOWN DISEASE GENES WITHIN PREDIXCAN RESULTS ###
    #################################################################
cohorts<-c("RA","CD","BD","CAD","HT","T2D")
r2thres<-0.01
ps<-c(0.05,0.01,0.001)
#set.seed(110285)
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
for (pthres in ps) {
    figurefile<-paste("Enrichment_pthres",pthres,".ps",sep="")
    postscript(figurefile,height=9, width= 9)
    layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
   # par(mar=c(5,4,1.5,1)+0.1)
letters<-c("a","b","c","d","e","f")
#par(mfrow=c(3,3))
for (j in 1:6) {
	co<-cohorts[j]
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
	file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
    dat<-dat[dat$R2>=r2thres,]
    bf<-0.05/dim(dat)[1]
    #known<-read.table(paste("/group/im-lab/nas40t2/haky/main/Data/NHGRI-GWAS/",co,"-genes.txt",sep=""),header=F)
    known<-read.table(paste(co,"-genes.txt",sep=""),header=F)
    known$ragene<-1
    data<-merge(dat,known,by.x=1,by.y=1,all.x=T,all.y=F)
    data$ragene[is.na(data$ragene)]<-0
    a<-glm(data$ragene~data$p.val+data$n.snps+data$R2,family="binomial")
    #plot(-log10(dat$p.val),dat$R2)
    #print(summary(a))
    nsig<-NULL
    nreps<-10000
    set.seed(110285)
	for (i in 1:nreps) {
        sum(data$ragene)
        rows<-sample(1:dim(data)[1],sum(data$ragene),replace=F)
        nsig<-c(nsig,sum(data$p.val[rows]<pthres))
    }
    truesig<-sum(data$p.val[data$ragene==1]<pthres)
    hist(nsig,xlim=c(0,max(c(nsig,truesig))),xlab=paste("Expected genes with p<",pthres,sep=""),main=co,ylab="Frequency",cex.lab=1.3,cex.axis=1,las=1,bty="n")
    points(truesig,0,col="black",pch=20, cex=1.5)
    outer = FALSE
    line = 2
    cex = 2
    adj  = 0.025
    letter<-letters[j]
    title(outer=outer,adj=adj,main=letter,cex.main=cex,col="black",font=2,line=line)
    penrich<-sum(nsig>=truesig)/nreps
    print(penrich)
}
   dev.off()
}


    ####################################################################
    ##qqplots with gene names labled for top genes for other diseases ##
    ####################################################################
cohorts<-c("RA","CD","BD","CAD","HT","T2D")
r2thres<-0.01
#set.seed(110285)
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
    figurefile<-paste("6disease_qqplot.ps",sep="")
    postscript(figurefile,height=9, width= 9)
    layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
   # par(mar=c(5,4,1.5,1)+0.1)
letters<-c("a","b","c","d","e","f")
#par(mfrow=c(3,3))
for (j in 1:6) {
	co<-cohorts[j]
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
	file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
    dat<-dat[dat$R2>=r2thres,]
    bf<-0.05/dim(dat)[1]
     qqunif(dat$p.val,plot=T)
    nn = length(dat$p.val)
    dat<-dat[order(dat$p.val),]
    dat$xx =  -log10((1:nn)/(nn+1))
    # title(paste("PrediXcan Results\n",co,", EN0.5",sep=""))
    tmp<-dat[dat$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]
    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>3) {
            tmp<-tmp[1:3,]
        }
        text(tmp$xx,-log10(tmp$p.val),tmp$gene,pos=2,cex=0.7)
    }
      outer = FALSE
    line = 2
    cex = 2
    adj  = 0.025
    letter<-letters[j]
    title(outer=outer,adj=adj,main=letter,cex.main=cex,col="black",font=2,line=line)
	title(co)  
}
dev.off()

    ###################################
    ## Manhattan plots for 6 disease ##
    ###################################
cohorts<-c("RA","CD","BD","CAD","HT","T2D")
r2thres<-0.01
#set.seed(110285)
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
    figurefile<-paste("6disease_manhattan.ps",sep="")
    postscript(figurefile,height=11, width= 8)
    layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
    par(mar=c(5,4,1.5,1)+0.1)
letters<-c("a","b","c","d","e","f")
#par(mfrow=c(3,3))
for (j in 1:6) {
	co<-cohorts[j]
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
	file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
    dat<-dat[dat$R2>=r2thres,]
    bf<-0.05/dim(dat)[1]
    allgenes<-read.table("gencode.v12.V1.summary.protein.nodup.genenames",header=F)
#    allgenes<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/gencode.v12.V1.summary.protein.nodup.genenames",header=F)
    annot<-merge(allgenes,qual,by.x=6,by.y=1,all.x=F,all.y=T)
    annot$mid<-(annot$V4+annot$V3)/2
    annot$pos<-0
    start<-0
    ends<-NULL
    starts<-NULL
    annot$color<-"grey65"
    ticks<-NULL
    for (chr in 1:22) {
        starts<-c(starts,start)
        ch<-paste("chr",chr,sep="")
        if (chr%%2==1) {annot$color[annot$V1==ch]<-"black"}
        annot$pos[annot$V1==ch]<-annot$mid[annot$V1==ch]+start
        start<-max(annot$pos)+1
        ends<-c(ends,start)
        ticks<-c(ticks,mean(annot$pos[annot$V1==ch]))
    }
    data<-merge(dat,annot,by.x=1,by.y=1,all.x=T,all.y=F)
    plot(data$pos,-log10(data$p.val),pch=20,ylab=expression(-log[10](italic(p))),xlab="Chromosome",xaxt="n",main="",cex.lab=1.3,cex.axis=1,col=data$color,las=1,cex=0.5,mgp=c(2,1,0),bty="n",ylim=c(0,max(c(-log10(bf),-log10(data$p.val)))))
    axis(side=1,at=ticks,labels=1:22,cex.lab=1)
      outer = FALSE
    line = 0.4
    cex = 2
    adj  = 0.025
    letter<-letters[j]
    title(outer=outer,adj=adj,main=letter,cex.main=cex,col="black",font=2,line=line)
	title(co)  
    #plot(rnorm(100),col="blue")
    #title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
    #plot(rnorm(100),col="green")
    #title(outer=outer,adj=adj,main="C",cex.main=cex,col="black",font=2,line=line)
    #MakeLetter("c")
    abline(h=-log10(bf),col="blue")
    tmp<-data[data$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]
    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>3) {
            tmp<-tmp[1:3,]
        }
        text(tmp$pos,-log10(tmp$p.val),tmp$gene,pos=4,cex=0.7)
    }
}
dev.off() 
 

############################################################################
### table of top result- bonferroni corrected and "moderately" associated ##
############################################################################

cohorts<-c("RA","BD","CD","CAD","HT","T2D","T1D")
r2thres<-0.01
BFTab<-NULL
moderate<-10^-4
MODTab<-NULL
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
for (co in cohorts) {
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
 	file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    dat$disease<-co
    dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
    dat<-dat[dat$R2>=r2thres,]
    print(dim(dat))
    bf<-0.05/dim(dat)[1]
    tmp<-dat[dat$p.val<=bf,]
    if (dim(tmp)[1]>0) {
        #tmp$disease<-co
        # tmp<-merge(tmp,qual,by.x=1,by.y=1,all.x=T,all.y=F)
        tmp<-tmp[,c(5,1,3,4,10,11)]
        tmp<-tmp[order(tmp$p.val),]
        BFTab<-rbind(BFTab,tmp)
        
    }
    tmp<-dat[dat$p.val<=moderate,]
    if (dim(tmp)[1]>0) {
        # tmp$disease<-co
        #tmp<-merge(tmp,qual,by.x=1,by.y=1,all.x=T,all.y=F)
        tmp<-tmp[,c(5,1,3,4,10,11)]
        tmp<-tmp[order(tmp$p.val),]
        MODTab<-rbind(MODTab,tmp)
    }
}
colnames(BFTab)<-c("disease","gene","PrediXcan.z.stat","PrediXcan.p.val","n.snps.DGN.predictor","DGN.CV.R2")
colnames(MODTab)<-c("disease","gene","PrediXcan.z.stat","PrediXcan.p.val","n.snps.DGN.predictor","DGN.CV.R2")
#write.table(BFTab,paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/summarytables/WTCCC_PrediXcanResults_EN0.5_BonferroniSignificant_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
#write.table(MODTab,paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/summarytables/WTCCC_PrediXcanResults_EN0.5_",moderate,"_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(BFTab,paste("WTCCC_PrediXcanResults_EN0.5_BonferroniSignificant_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(MODTab,paste("WTCCC_PrediXcanResults_EN0.5_",moderate,"_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

#############################################
### Predictor SNP enrichment in gwas metas ##
#############################################
cohorts<-c("RA","CD","BD","CAD","HT","T2D","T1D")
r2thres<-0.01
#set.seed(110285)
#qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
qual<-read.table("DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
#betas<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/DGN-WB_elasticNet_alpha0.5_weights_all_chr1-22_2015-02-02.txt", header=T)
betas<-read.table("DGN-WB_elasticNet_alpha0.5_weights_all_chr1-22_2015-02-02.txt", header=T)

### RA
co<-"RA"
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
dat<-read.table(file,header=T)
dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
dat<-dat[dat$R2>=r2thres,]
bf<-0.05/dim(dat)[1]
genes<-dat$gene[dat$p.val<=bf]
keepbetas<-merge(genes,betas,all.x=F,all.y=F,by.x=1,by.y=1)
metafile<-paste(co,"meta_beta_overlap.txt",sep="")
meta<-read.table(metafile,header=T)
keepmeta<-merge(keepbetas,meta,all.x=F,all.y=F,by.x=2,by.y=1)
keepmeta<-keepmeta[,c(1,2,3,4,5,13)]
colnames(keepmeta)<-c("SNP","gene","predix_refallele","predix_effectallele","predix_beta","meta_pvalue")
write.table(keepmeta, "RA_topgenes_meta_pvalues.txt",col.names=T,row.names=F,quote=F)

### CD
co<-"CD"
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
dat<-read.table(file,header=T)
dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
dat<-dat[dat$R2>=r2thres,]
bf<-0.05/dim(dat)[1]
genes<-dat$gene[dat$p.val<=bf]
keepbetas<-merge(genes,betas,all.x=F,all.y=F,by.x=1,by.y=1)
metafile<-paste(co,"meta_beta_overlap.txt",sep="")
meta<-read.table(metafile,header=T)
meta<-meta[,c(1,5)]
keepmeta<-merge(keepbetas,meta,all.x=F,all.y=F,by.x=2,by.y=1)
colnames(keepmeta)<-c("SNP","gene","predix_refallele","predix_effectallele","predix_beta","meta_pvalue")
write.table(keepmeta, "CD_topgenes_meta_pvalues.txt",col.names=T,row.names=F,quote=F)

### BD
co<-"BD"
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
dat<-read.table(file,header=T)
dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
dat<-dat[dat$R2>=r2thres,]
bf<-0.05/dim(dat)[1]
genes<-dat$gene[dat$p.val<=bf]
keepbetas<-merge(genes,betas,all.x=F,all.y=F,by.x=1,by.y=1)
metafile<-paste(co,"meta_beta_overlap.txt",sep="")
meta<-read.table(metafile,header=F)
meta<-meta[,c(1,8)]
keepmeta<-merge(keepbetas,meta,all.x=F,all.y=F,by.x=2,by.y=1)
colnames(keepmeta)<-c("SNP","gene","predix_refallele","predix_effectallele","predix_beta","meta_pvalue")
write.table(keepmeta, "BD_topgenes_meta_pvalues.txt",col.names=T,row.names=F,quote=F)

### CAD-- no genome wide sig results

### HT
co<-"HT"
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
dat<-read.table(file,header=T)
dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
dat<-dat[dat$R2>=r2thres,]
bf<-0.05/dim(dat)[1]
genes<-dat$gene[dat$p.val<=bf]
keepbetas<-merge(genes,betas,all.x=F,all.y=F,by.x=1,by.y=1)
metafile<-paste(co,"meta_beta_overlap.txt",sep="")
meta<-read.table(metafile,header=T)
meta<-meta[,c(1,4,5)]
keepmeta<-merge(keepbetas,meta,all.x=F,all.y=F,by.x=2,by.y=1)
colnames(keepmeta)<-c("SNP","gene","predix_refallele","predix_effectallele","predix_beta","meta_pvalue_SBP","meta_pvalue_DBP")
write.table(keepmeta, "HT_topgenes_meta_pvalues.txt",col.names=T,row.names=F,quote=F)

### T1D
co<-"T1D"
#    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
file = paste("PrediXcan_",co,"_DGNWholeBlood_EN0.5.txt",sep="")
dat<-read.table(file,header=T)
dat<-merge(dat,qual,by.x=1,by.y=1,all.x=T,all.y=F)
dat<-dat[dat$R2>=r2thres,]
bf<-0.05/dim(dat)[1]
genes<-dat$gene[dat$p.val<=bf]
keepbetas<-merge(genes,betas,all.x=F,all.y=F,by.x=1,by.y=1)
metafile<-paste(co,"meta_beta_overlap.txt",sep="")
meta<-read.table(metafile,header=T)
meta<-meta[,c(4,9,1)]
keepmeta<-merge(keepbetas,meta,all.x=F,all.y=F,by.x=2,by.y=1)
colnames(keepmeta)<-c("SNP","gene","predix_refallele","predix_effectallele","predix_beta","meta_pvalue","platform_filename")
write.table(keepmeta, "T1D_topgenes_meta_pvalues.txt",col.names=T,row.names=F,quote=F)

### T2D # no significnat genes



