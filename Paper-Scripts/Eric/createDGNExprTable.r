# Eric R. Gamazon
# Create required Expr Table

a <- read.table('/nas40t0/egamazon/VANDY/PREDIXCAN/DGN-WB.exp.IDxGENE', header=F)
a.t <- t(a)
dim(a.t)
genes <- read.table('/nas40t0/egamazon/VANDY/PREDIXCAN/DGN-WB.exp.GENE.list', header=F)
row.names(a.t) = t(genes)
write.table(a.t, file="DGN-WB.forMatrixeQTL.GENExID", sep="\t", row.names=T, quote=F, col.names=F)


