####by Heather E. Wheeler 20150209####
args <- commandArgs(trailingOnly=T)
#args <- c('21','0.0001')
date <- Sys.Date() 
"%&%" = function(a,b) paste(a,b,sep="")
###############################################
### Directories & Variables
pre <- "/group/im-lab/nas40t2/hwheeler/"

my.dir <- pre %&% "PrediXcan_CV/"
ct.dir <- pre %&% "PrediXcan_CV/cis.v.trans.prediction/"
gt.dir <- pre %&% "PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/DGN-imputed-for-PrediXcan/"
en.dir <- pre %&% "PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-calc-weights/DGN-WB_weights/"

tis <- "DGN-WB"
k <- 10 ### k-fold CV
chrom <- as.numeric(args[1])
chrname <- "chr" %&% chrom
pthresh <- as.numeric(args[2]) ##p-vaule threshold to test in CV
 
################################################
### Functions & Libraries
 
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)

## fit betas by Vasa edited by Heather to get pvals
fit.betas <- function(G, Y) {
require(RcppGSL)
betas <- rep(NA, ncol(G))
stderrs <- rep(NA, ncol(G))
dfs <- rep(NA, ncol(G))
for (i in 1:ncol(G)) {
res <- fastLmPure(y=Y, X=matrix(G[,i]))
betas[i] <- res$coefficients
stderrs[i] <- res$stderr
dfs[i] <- res$df
}
t <- betas/stderrs
pvals <- 2*pt(-abs(t),df=dfs)
output <- cbind(betas,pvals)
return(output)
}

 
################################################
###input adjusted expression data###
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

bimfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% chrom %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

famfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.ID.list" ###samples with imputed gt data###
fam <- scan(famfile,"character")
samplelist <- intersect(fam,expsamplelist)

exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

gtfile <- gt.dir %&% 'DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr' %&% chrom %&% '.SNPxID'
gtX <- scan(gtfile)
gtX <- matrix(gtX, ncol = length(fam), byrow=TRUE)
colnames(gtX) <- fam
rownames(gtX) <- bim$V2
X <- gtX[,samplelist]

X <- t(X) #transpose to match code below (one ID per row)
 
###create results array
resultsarray <- array(NA,c(length(explist),5))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","Pthresh","n.snps","CV.R2","CV.pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- "working_" %&% tis %&% "_exp_" %&% k %&% "-foldCV_polyscore_Pthresh" %&% pthresh %&% "_imputedSNPs_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(resultscol,file=workingbest,ncolumns=5,sep="\t")

weightcol = c("gene","SNP","refAllele","effectAllele","beta")
workingweight <- en.dir %&% tis %&% "_polyscore_Pthresh" %&% pthresh %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncol=5,sep="\t")

###run polyscore CV

set.seed(1001)

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
    ###effectively skips genes with <2 cis-SNPs
    resinfo <- c(gene,pthresh,"no.cis.snps",NA,NA)
    resultsarray[gene,] <- resinfo
  }else{
    minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps]
    cisgenos <- scale(cisgenos, center=T, scale=T)
    cisgenos[is.na(cisgenos)] <- 0
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      ###effectively skips genes with <2 cis-SNPs
      resinfo <- c(gene,pthresh,"no.cis.snps",NA,NA)
      resultsarray[gene,] <- resinfo   
    }else{
      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)

      ### run k-fold CV

      ###randomize samples into CV groups
      samplelength <- length(exppheno)
      g <- 1:k ##k-fold CV
      groupid <- sample(g,samplelength,replace=T)
      newiddata <- data.frame(groupid,samplelist)

      ### create empty matrix to fill in with polyscore vals
      cvset <- data.frame(samplelist,matrix(NA,samplelength,1),exppheno)
      rownames(cvset) <- samplelist
      names(cvset) <- c("ID","Pred","Exp")
      for(idsubset in 1:k){
        trainset <- data.frame(newiddata$samplelist[newiddata$groupid != idsubset])
        testset <- data.frame(newiddata$samplelist[newiddata$groupid == idsubset])
        ### pull trainset genos and phenos
        traingenos <- cisgenos[intersect(trainset[,1],rownames(cisgenos)),]
        trainexp <- exppheno[match(trainset[,1],rownames(exppheno))]  ###if scaled, use rownames
        ### calculate betas and pvals from training set
        betares <- fit.betas(traingenos,trainexp)
        rownames(betares) <- colnames(traingenos)
        topres <- subset(betares,betares[,2] < pthresh) ###pull snps with P < pthresh
        betas <- topres[,1]
        if(length(betas) <= 1){
          names(betas) <- rownames(topres)
        }
        ### polyscore
        ### multiply betas by test genotypes and take sum
        testgenos <- cisgenos[intersect(testset[,1],rownames(cisgenos)),intersect(names(betas),colnames(cisgenos))]
        if(is.array(testgenos)=='FALSE'){ ### calcs polyscore if only one individual in test set, happens when expression sample size is low (<100)
          if(dim(testset)[1] < 2){ ### calcs polyscore if only one individual in test set, happens when expression sample size is low (<100)
            pred.polyscore <- sum(testgenos * betas)
            names(pred.polyscore)<-testset[,1]
          }else{ ### calcs polyscore if only one beta with p<pthresh
            pred.polyscore <- testgenos * betas
          }
          cvset$Pred[match(names(pred.polyscore),rownames(cvset))] <- pred.polyscore
        }else{
       	  testsweep <- sweep(testgenos, 2, betas, "*")
          pred.polyscore <- as.vector(rowSums(testsweep))
          names(pred.polyscore) <- rownames(testsweep)
          ### add predictions to cvset
          cvset$Pred[match(names(pred.polyscore),rownames(cvset))] <- pred.polyscore
        }
      }


      ### output bestbetas for PrediXcan
      allbetares <- fit.betas(cisgenos,exppheno) ###~30sec/gene
      rownames(allbetares) <- colnames(cisgenos)
      alltopres <- subset(allbetares,allbetares[,2] < pthresh) ###pull snps with P < pthresh
      allbetas <- alltopres[,1]
      nallbetas <- length(allbetas)
      if(nallbetas <= 1){
        names(allbetas) <- rownames(alltopres)
      }
      if(nallbetas==0){
        ### calculate correlation between predicted and observed expression
        res<-summary(lm(cvset$Pred~cvset$Exp))
        resinfo <- c(gene,pthresh,nallbetas,res$r.squared,res$coef[2,4])
        resultsarray[gene,] <- resinfo
      }else{
        allbetalist <- names(allbetas)
        allbetainfo <- bim[allbetalist,]
        betatable <- as.matrix(cbind(allbetainfo,alltopres))
        betafile<-cbind(gene,betatable[,2],betatable[,5],betatable[,6],betatable[,7]) ##output "gene","SNP","refAllele","effectAllele","beta"
        write(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function

        ### calculate correlation between predicted and observed expression
        res<-summary(lm(cvset$Pred~cvset$Exp))
        resinfo <- c(gene,pthresh,nallbetas,res$r.squared,res$coef[2,4])
        resultsarray[gene,] <- resinfo  
      }
    }
  }
  write(resultsarray[gene,],file=workingbest,ncolumns=5,append=T,sep="\t")
}

write.table(resultsarray,file=tis %&% "_exp_" %&% k %&% "-foldCV_polyscore_Pthresh" %&% pthresh %&% "_imputedSNPs_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")



