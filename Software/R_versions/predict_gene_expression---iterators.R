#!/usr/bin/Rscript --vanilla
library(RSQLite)
library(data.table)

Rprof()

D <- NULL
genelist_vector <- NULL

get_all_dosages <-function(directory, prefix) {
  chrfiles <- list.files(directory)
  chrfiles <- chrfiles[grepl(paste('^', prefix, sep=''), chrfiles)]
  chrfile <- NA
  chrindex <- 0
  con <- NA
  con.skip <- 0
  lineBuffer <- NA
  lineBuffer.skip <- 0
  
  parsedLine <- function(line) {  
    arr <- unlist(strsplit(line, '\t'))
    rsid <- arr[2]
    refallele <- arr[5]        
    dosages <- as.numeric(arr[7:length(arr)])
    return(list(rsid=rsid, refallele=refallele, dosages=dosages))
  }
  
  getLine <- function() {
    if(is.na(con)) {return(NA)}
    if (is.na(lineBuffer) || lineBuffer.skip > length(lineBuffer)) # Renew the buffer
     {
      lineBuffer <<- scan(con, what="character", sep='\n', skip=con.skip, nmax=100000)  
      lineBuffer.skip <<- 0
      con.skip <<- length(lineBuffer)
    }
    if (length(lineBuffer) == 0) # Nothing more in the file
     {return(NA)}
    else {
      lineBuffer.skip <<- lineBuffer.skip + 1 # Iterate over the buffer
      return(parsedLine(lineBuffer[lineBuffer.skip]))}
    }    

  function() {
    line <- getLine()
    if (any(is.na(line))){ # First invocation, or finished the file.
      chrindex <<- chrindex + 1
      if (chrindex > length(chrfiles)) {return(NA) } # No more files
      # Otherwise, let's process the next file
      chrfile <<- chrfiles[chrindex]  
      print(paste(Sys.time(), "Processing", chrfile))    
      con <<- gzfile(file.path(directory, chrfile), open="r")
      return(getLine()) # return the first line of the new file
    } else {
    return(line)
    }        
  }
}

get_all_applications_of <- function(rsid) {
  query <- paste("SELECT gene, weight, eff_allele FROM weights WHERE rsid=\"", rsid, "\"", sep='')
  apps <- dbGetQuery(weights.db.connection, query)
  idx <- 0
  function () {
    if (idx >= nrow(apps)) {return(NA)}
    idx <<- idx + 1
    return(as.list(apps[idx,]))
  }}
                                          
transcription_matrix_update <- function(genelist, gene, weight, ref_allele, allele, dosage.row) {
  if (is.null(D)) {
    genelist_vector <<- read_genelist_vector(genelist)
    D <<- matrix(rep(0, length(genelist_vector)*length(dosage.row)), 
                nrow=length(genelist_vector),
                ncol=length(dosage.row))
    rownames(D) <<- genelist_vector
  }
  
  if (gene %in% genelist_vector) {
    if (ref_allele == allele) {multiplier = 1} else {multiplier = -1}
    D[gene, ] <<- D[gene,] + dosage.row * weight * multiplier
  }
}

read_genelist_vector <- function(filename){
  if (is.null(filename))
  {
    dbGetQuery(weights.db.connection, "SELECT DISTINCT gene FROM weights ORDER BY gene")$gene
  } else{
  lines <- readLines(file(filename))
  x<-unlist(regmatches(lines, gregexpr( "^.*\t"  ,lines), invert=T))
  x[grepl(pattern = ".", x)]  
}
}

save_transcription_matrix <- function () {  
  write.table(as.data.frame(t(D)), file="output.txt", sep="\t", col.names=T, row.names=F, quote=F)
}

predict_gene_expression <- function(genelist, dosages, prefix, weights, output) {
  weights.db.connection <<- dbConnect(RSQLite::SQLite(), dbname=weights)
  get_next_dosage <- get_all_dosages(dosages, prefix)    
  while(!any(is.na(ret <- get_next_dosage()))){   
   rsid <- ret$rsid
   allele <- ret$refallele
   dosage.row <- ret$dosages
   get_next_application <- get_all_applications_of(rsid=rsid)
   while(!any(is.na(ret <- get_next_application()))) {     
     gene <- ret$gene
     weight <- ret$weight
     ref_allele <- ret$eff_allele
     transcription_matrix_update(genelist, gene, weight, ref_allele, allele, dosage.row)
   }}
  save_transcription_matrix()
}
  
predict_gene_expression(NULL, 'data/dosages/', 'chr', 'data/weights.db', 'output.txt')
