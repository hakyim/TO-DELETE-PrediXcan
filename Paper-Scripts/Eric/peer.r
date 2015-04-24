
# INSTALL IMPUTE
# Eric Gamazon:  PEER calculations for expression data

source("http://bioconductor.org/biocLite.R")
biocLite("impute")

library(impute)

a <- read.table('extracted_expr/DGN-WB_expression', header=T)
#b <- a[,-1]
#b <- b[,-1]
bmat <- as.matrix(b)

#bimp <- impute.knn(bmat, rowmax=1, colmax=1)


# log transform!
cmat <- log2(bmat+1)

#c <- t(bmat)
c <-t(cmat)

.libPaths("/userhome/genegateRPackages/lib.2.14.0/")
library(peer)

dim(c)

model = PEER()
PEER_setPhenoMean(model,as.matrix(c))
dim(PEER_getPhenoMean(model))
PEER_setNk(model,10)
PEER_getNk(model)

# this doesn't seem to work even though it's in NAture Protocol paper.. Got it to work ... too many NAs before
#PEER_setNMax_iterations(model, 10000)

# this was causing a crash in PEER
PEER_update(model)

factors = PEER_getX(model)

dim(factors)

weights = PEER_getW(model)

dim(weights)

precision = PEER_getAlpha(model)

dim(precision)

residuals = PEER_getResiduals(model)
dim(residuals)

jpeg('PEER/diagnostics_peer.whole_blood.jpg')
PEER_plotModel(model)
dev.off()

