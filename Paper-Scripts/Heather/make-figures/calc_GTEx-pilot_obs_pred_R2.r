###Calc R2 between DGN-predicted and observed expression in GTEx pilot 9 tissues
date <- Sys.Date()
"%&%" = function(a,b) paste(a,b,sep="")

#predfile<-'/Volumes/im-lab/nas40t2/smozaffari/Elastic_Net/GTEx_pilot_predicted' ##on tarbell
predfile<-'/Users/heather/Dropbox/PrediXcan/GTEx/GTEx_pilot_predicted_from_DGN_elasticNet_alpha0.5' ##local copy

#obsdir<-'/Volumes/im-lab/nas40t2/nknoblauch/predixmod_Inputs/peer_factors/' ##on tarbell
obsdir<-'/Users/heather/Dropbox/PrediXcan/GTEx/gtex-exp/' ##local copy

my.dir <- '/Volumes/im-lab/nas40t2/hwheeler/PrediXcan-Paper/scripts/Heather/make-figures/'

pilotfile<-my.dir %&% 'pilot-GTEx.SubjectDS.v10.1.DRAFT.txt'
pilot<-read.table(pilotfile,header=T)
pilotlist<-pilot$SUBJID

pred<-read.table(predfile,header=T)
rownames(pred)<-pred[,1]
tpred<-t(pred[,-1])
predid<-scan(my.dir %&% 'GTEx_pilot_predicted_from_DGN_elasticNet_alpha0.5_ID.list','character')
rownames(tpred)<-predid

tissuelist<-scan(my.dir %&% 'GTEx_nine_tissues','character')

nmatrix<-matrix(0,ncol=2,nrow=length(tissuelist))

for(i in 1:length(tissuelist)){
  tis <- tissuelist[i]
  obs<-readRDS(obsdir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS")
  genefile <- obsdir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENENAME.list"
  idfile <- obsdir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list.fixed"
  genelist <- scan(genefile,"character")
  idlist <- scan(idfile,"character")
  colnames(obs) <- genelist
  rownames(obs) <- idlist
  
  commonIDs<-intersect(rownames(tpred),rownames(obs))
  commonIDs<-intersect(commonIDs,pilotlist)
  newpred<-tpred[commonIDs,]
  newobs<-obs[commonIDs,]
  n<-length(commonIDs)
  nmatrix[i,]<-c(tis,n)
  
  commonGenes<-intersect(colnames(newpred),colnames(obs))
  
  resultsmat <- matrix(0,ncol=2,nrow=length(commonGenes))
  colnames(resultsmat) <- c('R2','p')
  rownames(resultsmat) <- commonGenes
  
  for(j in 1:length(commonGenes)){
    gene<-commonGenes[j]
    res<-cor.test(newobs[,gene],newpred[,gene])
    info<-c(round(res$estimate^2,5),signif(res$p.value,5))
    resultsmat[j,]<-info
  }
  
  finalres<-resultsmat[order(resultsmat[,1],decreasing=T),]
  write.table(finalres,file="train.DGN-WB_test.GTEX_" %&% tis %&% "_R2.txt",quote=F)
}

write.table(nmatrix,file="GTEx_pilot_sample_sizes.txt",quote=F,row.names=F,col.names=F)
