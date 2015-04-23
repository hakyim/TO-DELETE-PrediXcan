# this is a file that will download NHGRI catalog and perform some annotations

# functions
'%&%' <- function(a, b) paste(a, b, sep="") 

# grab the catalog 04/15/2015
wget_string = "wget -qO- http://www.genome.gov/admin/gwascatalog.txt"
syst_string = system(wget_string, intern=TRUE)

# set locale
Sys.setlocale(locale="C")

# WTCCC 7 traits
traits = c("Type 2 diabetes", 
           "Type 1 diabetes",
           "Rheumatoid arthritis",
           "Crohn's disease",
           "Hypertension",
           "Bipolar disorder",
           "Coronary heart disease")

# parse catalog
parse.df = data.frame()
for (i in 2:(length(syst_string)-1)) {
  cat(i, "\n")
  temp = strsplit(syst_string[i], "\t")
  temp = t(as.data.frame(temp[[1]]))
  colnames(temp) = strsplit(syst_string[1], "\t")[[1]]
  parse.df = rbind(parse.df, temp)
}

# identifying which papers had WTCCC data
out.vec = vector()
out.met = vector()
for (i in 1:nrow(parse.df)) {
  cat(i, "\n")
  pull.1 = grep("WTCCC",
                as.character(as.matrix(parse.df[i,])))
  if (length(pull.1) > 0) {
    out.vec[i] = 1
  } else {
    out.vec[i] = 0
  }
  
  pull.2 = grep("Wellcome",
                as.character(as.matrix(parse.df[i,])))
  if (length(pull.2) > 0) {
    out.vec[i] = out.vec[i] + 1
  } else { 
    out.vec[i] = out.vec[i] + 0
  }

  pull.3 = grep("Meta-analy", 
                as.character(as.matrix(parse.df[i,])),
                ignore.case=TRUE)
  if (length(pull.3) > 0) {
    out.met[i] = 1
  } else {
    out.met[i] = 0
  }
}

# attaching this new info
parse.df$WTCCC = out.vec
parse.df$METAA = out.met

# selecting papers not involved in WT
noWT = parse.df[parse.df[,35] == 0,]

filter = data.frame(PMID=noWT[, 2], 
                    trai=noWT[, 8],
                    gene=noWT[,14],
                    snps=noWT[,22],
                    chid=noWT[,12],
                    cpos=noWT[,13],
                    pval=noWT[,28],
                    meta=noWT[,36])

# subset by trait, should be ordered by traits vec
sep.traits = list(t2d=filter[filter$trai == traits[1],],
                  t1d=filter[filter$trai == traits[2],],
                  rad=filter[filter$trai == traits[3],],
                  crd=filter[filter$trai == traits[4],],
                  hyt=filter[filter$trai == traits[5],],
                  bpd=filter[filter$trai == traits[6],],
                  chd=filter[filter$trai == traits[7],])

# write to file
for (i in 1:length(sep.traits)) {
  temp = strsplit(traits[i], " ")
  temp = paste(temp[[1]][1:length(temp[[1]])], collapse = "_")
  write.table(sep.traits[[i]],
              file = "NHGRI_"%&% temp %&%"_noWTCCC.txt",
              quote=F,
              col.names=T,
              row.names=F,
              sep = "\t")
}

#
