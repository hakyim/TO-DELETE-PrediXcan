library(dplyr, warn.conflicts=FALSE)

read_pheno <- function(fam_file) {
    fam = read.table(fam_file)
    colnames(fam) <- c("FAMILY","INDIVIDUAL", "PATERNAL", "MATERNAL", "SEX","PHENOTYPE")
    return(fam)
}

read_filter <- function(filter_file) {
    filter = read.table(filter_file)
    colnames(filter) <- c("FAMILY", "INDIVIDUAL", "FILTER")
    return(filter)
}
    
read_predicted <- function(predicted_file) {
    predicted = read.delim(predicted_file)
    return(predicted)
}

association <- function(fam, filter, predicted) {
    genes = colnames(predicted)
    merged_people = merge(fam, filter)
    merged <- cbind(predicted, merged_people)
    filtered <- merged %>% filter(FILTER==1)

    pheno <- filtered$PHENO
    filtered <- subset(filtered, select = genes)

    OUT<-NULL
    for (gene in genes) {
      data <- filtered[[gene]]
      results = coef(summary(lm(pheno ~ data)))[c(2,6,8,4)]
      line = c(gene,results)
      OUT<-rbind(OUT,line)
    }
    colnames(OUT)<-c("gene","beta","z-stat","p-val", "se(beta)")
    return(OUT)
}

write_association <- function(OUT, output_file) {
    write.table(OUT,output_file,col.names=T,row.names=F,quote=F)
}