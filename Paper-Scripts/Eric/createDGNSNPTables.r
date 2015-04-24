# Eric R. Gamazon
# Create genotype table for eQTL mapping

"%&%" = function(a,b) paste(a,b,sep="")
mydir = "/nas40t0/egamazon/VANDY/PREDIXCAN/"; 

for (i in (1:22))
{

	a <- read.table(gzfile(mydir %&% 'DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr' %&% i %&% '.mldose.gz'), header=F)
	a <- a[,-1]
	a <- a[,-1]
	a.t <- t(a)
	dim(a.t)
	snps <- read.table(gzfile(mydir %&% 'DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr' %&% i %&% '.mlinfo.gz'), header=T)
	row.names(a.t) = t(snps[,1])
	write.table(a.t, file=mydir %&% "DGN.imputed" %&% ".SNPxID" %&% i, sep="\t", row.names=T, quote=F, col.names=F)

}

