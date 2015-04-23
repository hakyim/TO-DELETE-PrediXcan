####by Heather E. Wheeler 20150202####
##see runscripts/run_01_imputedDGN-WB_CV_elasticNet_chr*sh and qsub.txt for tarbell job submission scripts
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','1')
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables
#pre <- "/Users/heather/Dropbox/elasticNet_testing"
pre <- ""

my.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
ct.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
gt.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/DGN-imputed-for-PrediXcan/"
en.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/"

#snpset <- ".hapmapSnpsCEU" ##2015-02-02 results
#snpset <- ".wtcccGenotypedSNPs" ##2015-03-12 results
snpset <- "_1000G"

k <- 10 ### k-fold CV
tis <- "DGN-WB"  
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

alpha <- as.numeric(args[2]) #alpha to test in CV

################################################
### Functions & Libraries

library(glmnet)
#library(doMC) ##slower on tarbell than using 1 core, not sure why
#registerDoMC(10)
#getDoParWorkers()

################################################
rpkmid <- ct.dir %&% tis %&% ".exp.ID.list"
expid <- scan(rpkmid,"character")
rpkmgene <- ct.dir %&% tis %&% ".exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- ct.dir %&% tis %&% ".exp.IDxGENE"
expdata <- scan(rpkmfile) 
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

t.expdata <- expdata #don't need to transpose DGN

gencodefile <- my.dir %&% 'gencode.v12.V1.summary.protein.nodup.genenames'
gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,6]
gencode <- gencode[gencode[,1]==chrname,] ##pull genes on chr of interest
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(t.expdata) ###samples with exp data###

#bimfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".bim" ###get SNP position information###
#bim <- read.table(bimfile)
bimfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".bim.rds" ###get SNP position information###
bim <- readRDS(bimfile)
rownames(bim) <- bim$V2
                
famfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".ID.list" ###samples with imputed gt data###
fam <- scan(famfile,"character")
samplelist <- intersect(fam,expsamplelist)
                        
exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

#gtfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".SNPxID"
#gtX <- scan(gtfile)
gtfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".SNPxID.rds" #smaller file format, faster readin
gtX <-readRDS(gtfile)
gtX <- matrix(gtX, ncol = length(fam), byrow=TRUE)
colnames(gtX) <- fam
rownames(gtX) <- bim$V2
X <- gtX[,samplelist]

X <- t(X) #transpose to match code below (one ID per row)

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- "working_" %&% tis %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene","SNP","refAllele","effectAllele","beta")
workingweight <- en.dir %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncol=5,sep="\t")

set.seed(1001) ##forgot to include in 2/2/15 run, should I re-run?

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  geneinfo <- gencode[gene,]
  chr <- geneinfo[1]
  c <- substr(chr$V1,4,5)
  start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
  cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
  if(is.null(dim(cisgenos))){
    bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps]
    cisgenos <- scale(cisgenos, center=T, scale=T)
    cisgenos[is.na(cisgenos)] <- 0
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{

      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)
  
      ##run Cross-Validation over alphalist
      fit <- cv.glmnet(cisgenos,exppheno,nfolds=k,alpha=alpha,keep=T,parallel=F) ##parallel=T is slower on tarbell, not sure why

      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
      best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
      cvm.best = best.lam[,1]
      lambda.best = best.lam[,2]
      nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
      
      ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
      ret[ret == 0.0] <- NA
      bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
      names(bestbetas) = rownames(ret)[which(!is.na(ret))]

      pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda

    }
  }
  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    genename <- as.character(gencode[gene,6])
    rsq <- res$r.squared
    pval <- res$coef[2,4]

    resultsarray[gene,] <- c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)

    
    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- bim[bestbetalist,]
    betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
    betafile<-cbind(genename,betatable[,2],betatable[,5],betatable[,6],betatable[,7]) ##output "gene","SNP","refAllele","effectAllele","beta"
    write(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function

  }else{
    genename <- as.character(gencode[gene,6])
    resultsarray[gene,1] <- genename
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)

  }
  write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
}


write.table(resultsarray,file=tis %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
