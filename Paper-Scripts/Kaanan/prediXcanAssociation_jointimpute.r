
source("/home/kshah2/qqplots.R")
predixcanassociation<-function(co,tis,alpha) {
    #cohorts<-c("T2D")
    #alpha<-1
    #tissues<-c("Thyroid","Adipose-Visceral-Omentum") #"Adipose-Subcutaneous","WholeBlood","Brain-Cerebellum","Brain-Hypothalamus","Liver","Muscle-Skeletal","Pancreas","Pituitary","SmallIntestine-TerminalIleum","Stomach","Colon-Sigmoid","Nerve-Tibial","Colon-Transverse",)
    #for (tis in tissues) {
    #for (co in cohorts) {
    #control1<-read.table(paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/GTExTissues/PredictedExpression_","58C","_",tis,"_EN",alpha,".txt",sep=""),header=T)
    #control1[(dim(control1)[1]+1),]<-c(NA,rep(0,(dim(control1)[2]-1)))
    #control2<-read.table(paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/GTExTissues/PredictedExpression_","NBS","_",tis,"_EN",alpha,".txt",sep=""),header=T)
    #control2[(dim(control2)[1]+1),]<-c(NA,rep(0,(dim(control2)[2]-1)))
    #control<-merge(control1,control2,by.x=1,by.y=1,all.x=F,all.y=F)
    scorefile = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PredictedExpression_",co,"_",tis,"_EN",alpha,".txt",sep="")
    preds<-read.table(scorefile,header=T)
    phenorow<-dim(preds)[1]+1
    preds[phenorow,]<-c(NA,rep(1,(dim(preds)[2]-1)))
    preds[phenorow,grep("58C",colnames(preds))]<-0
    preds[phenorow,grep("NBS",colnames(preds))]<-0
    #dat<-merge(preds,control,by.x=1,by.y=1,all.x=F,all.y=F)
    dat<-preds
    OUT<-NULL
    for (i in 1:(dim(dat)[1]-1)) {
        testpheno<-as.numeric(as.vector(t(dat[phenorow,-1])))
        geneexp<-as.numeric(as.vector(t(dat[i,-1])))
        tmp = coef(summary(glm(testpheno~geneexp,family="binomial",maxit=10)))[c(2,6,8)]
        tmp<-c(as.character(dat$gene[i]),tmp)
        OUT<-rbind(OUT,tmp)
    }
    outfile = paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_",tis,"_EN",alpha,".txt",sep="")
    colnames(OUT)<-c("gene","beta","z-stat","p-val")
    write.table(OUT,outfile,col.names=T,row.names=F,quote=F)
    jpeg(paste("/group/im-lab/nas40t2/kaanan/PrediXcan/WTCCC/",co,"/PrediXcan/PrediXcan_",co,"_",tis,"_EN",alpha,".jpeg",sep=""))
    qqunif(OUT[,4],plot=T)
    dev.off()
}
#}
#}

