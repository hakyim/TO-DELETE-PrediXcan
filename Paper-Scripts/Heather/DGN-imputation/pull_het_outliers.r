args <- commandArgs(trailingOnly=T)
"%&%" <- function(a, b) paste(a, b, sep="")

hetfile <- args[1]


HET=read.table(hetfile,header=T,as.is=T)
H = (HET$N.NM.-HET$O.HOM.)/HET$N.NM.
pdf(file=hetfile %&% ".pdf")
oldpar=par(mfrow=c(1,2))
hist(H,50)
hist(HET$F,50)
par(oldpar)
dev.off()
HET[order(HET$F),]

outliers<-data.frame()

for(i in 1:length(HET$F)){
	if(HET[i,6] > (mean(HET$F)+3*sd(HET$F))){
		outliers <- rbind(outliers,HET[i,])
	}
	if(HET[i,6] < (mean(HET$F)-3*sd(HET$F))){
		outliers <- rbind(outliers,HET[i,])
	}
}

write.table(outliers,file=hetfile %&% ".outlier.list",quote=F,col.names=F,row.names=F)
