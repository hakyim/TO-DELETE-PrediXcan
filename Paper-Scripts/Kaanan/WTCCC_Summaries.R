#### all summaries are for elastic net results with alpha = 0.5 based on the wtccc imputed data done by jung
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
        cex.lab=1,mgp=c(2,1,0),pch=20 )
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
        points(xx,  -sort(log10(p)),pch=20)
        #points(xx[dat[,3]<0.05],  -log10(dat[dat[,3]<0.05,1]),pch=20,col="pink")
        
        
        y<-max(-log10(p))
        text(0,y,paste(nsnps,"% of genes have a q-value <= ",FDRthres,sep=""),pos=4)
        
        #  print(paste(nsnps,"% of SNPs have a q-value <= ",FDRthres,sep=""))
        abline(0,1,col='red')
        if(BH)
        {
            abline(-log10(0.05),1, col='black',lty=2)
            abline(-log10(0.10),1, col='black',lty=3)
            abline(-log10(0.25),1, col='black',lty=4)
            legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
            col=c('black','black','black'),lty=2:4, cex=1)
            abline(h=-log10(0.05/nn),col="blue") ## bonferroni
        }
    } else {
        return(nsnps)
    }
}


### table of top result- bonferroni corrected and "moderately" associated
cohorts<-c("BD","CD","CAD","HT","RA","T1D","T2D")
BFTab<-NULL
moderate<-10^-4
r2thres<-0
MODTab<-NULL
qual<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_exp_10-foldCV_elasticNet_alpha0.5_imputedSNPs_chr1-22_2015-02-02.txt",header=T)
for (co in cohorts) {
    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/",co,"_PrediXcan_EN0.5.txt",sep="")
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
write.table(BFTab,paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/summarytables/WTCCC_PrediXcanResults_EN0.5_BonferroniSignificant_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(MODTab,paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/summarytables/WTCCC_PrediXcanResults_EN0.5_",moderate,"_r2",r2thres,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

##qqplots with gene names labled for top genes
#source("/home/kshah2/qqplots.R")
for (co in cohorts) {
    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/",co,"_PrediXcan_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    bf<-0.05/dim(dat)[1]
    jpeg(paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/qqplots/",co,"_PrediXcanResults_EN0.5_genenames_qq.jpeg",sep=""))
    qqunif(dat$p.val,plot=T)
    nn = length(dat$p.val)
    dat<-dat[order(dat$p.val),]
    dat$xx =  -log10((1:nn)/(nn+1))
    
    title(paste("PrediXcan Results\n",co,", EN0.5",sep=""))
    tmp<-dat[dat$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]
    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>10) {
            tmp<-tmp[1:5,]
        }
        text(tmp$xx,-log10(tmp$p.val),tmp$gene,pos=2,cex=0.7)
    }
    dev.off()
}


#Manhattan plots for each disease
allgenes<-read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/gencode.v12.V1.summary.protein.nodup.genenames",header=F)
annot<-merge(allgenes,qual,by.x=6,by.y=1,all.x=F,all.y=T)
annot$mid<-(annot$V4+annot$V3)/2
annot$pos<-0
start<-0
ends<-NULL
starts<-NULL
annot$color<-"grey25"
ticks<-NULL
for (chr in 1:22) {
    starts<-c(starts,start)
    ch<-paste("chr",chr,sep="")
    if (chr%%2==1) {annot$color[annot$V1==ch]<-"grey85"}
    annot$pos[annot$V1==ch]<-annot$mid[annot$V1==ch]+start
    start<-max(annot$pos)+1
    ends<-c(ends,start)
    ticks<-c(ticks,median(annot$pos[annot$V1==ch]))
    
}
for (co in cohorts) {
    file = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/",co,"_PrediXcan_EN0.5.txt",sep="")
    dat<-read.table(file,header=T)
    bf<-0.05/dim(dat)[1]
    jpeg(paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/ResultsFiguresAndTables/manhattanplots/",co,"_PrediXcanResults_EN0.5_genenames_manhattan.jpeg",sep=""))
    dat<-merge(dat,annot,by.x=1,by.y=1,all.x=T,all.y=F)
    plot(dat$pos,-log10(dat$p.val),pch=20,ylab="-log10(pvalue)",xlab="Genomic Position",xaxt="n",main=paste(co,"\nDGN, EN 0.5",sep=""),col=dat$color)
    axis(side=1,at=ticks,labels=1:22)
    abline(h=-log10(bf),col="blue",lty="dashed")
    tmp<-dat[dat$p.val<=bf,]
    tmp<-tmp[order(tmp$p.val),]

    if (dim(tmp)[1]>0) {
        if (dim(tmp)[1]>10) {
            tmp<-tmp[1:5,]
        }
        text(tmp$pos,-log10(tmp$p.val),tmp$gene,pos=3,cex=0.7)
    }
    dev.off()
}
