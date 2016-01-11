
read_pheno <- function(pheno_file) {
    pheno <- read.table(pheno_file, as.is=TRUE)
    return(pheno)
}

read_filter <- function(filter_file, filter_column) {
    fil <- read.table(filter_file, as.is=TRUE)
    # Only keep columns of filter file necessary
    fil <- fil[c(1, 2, filter_column)]
    names(fil) <- c("FID", "IID", "FIL_VAL")
    return(fil)
}
    
read_predicted <- function(pred_exp_file) {
    pred_exp <- read.delim(pred_exp_file)
    return(pred_exp)
}

# Filter out rows of the phenotype dataframe
filter_pheno <- function(pheno, fil, filter_on) {
    # merge the pheno df and fil by the Family ID and Individul ID columns
    pheno <- merge(pheno, fil, by=c(1,2), sort=FALSE)
    pheno <- filter(pheno, FIL_VAL==filter_on)
    pheno$FIL_VAL <- NULL
    return(pheno)
}

association <- function(pheno, pred_exp, fil, filter_on) {
    genes <- colnames(pred_exp)
    merged <- merge(cbind(pheno, pred_exp), fil, by=c(1,2), sort=FALSE)
    filtered <- filter(merged, FIL_VAL=filter_on)

    pheno <- filtered$PHENO
    filtered <- subset(filtered, select = genes)

    OUT<-NULL
    for (gene in genes) {
      data <- filtered[[gene]]
      results <- coef(summary(lm(pheno ~ data)))[c(2,6,8,4)]
      line <- c(gene,results)
      OUT <- rbind(OUT,line)
    }
    colnames(OUT) <- c("gene","beta","z-stat","p-val", "se(beta)")
    return(OUT)
}

write_association <- function(OUT, output_file) {
    write.table(OUT,output_file,col.names=T,row.names=F,quote=F)
}
