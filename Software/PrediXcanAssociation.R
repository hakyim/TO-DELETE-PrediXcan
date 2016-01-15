#!/usr/bin/env Rscript

read_pheno <- function(pheno_file, pheno_column = NULL, pheno_name = NULL) {
  pheno <- read.table(pheno_file, header = F, as.is = T)
  # Fix dataframe if there is a header row.
  if (pheno[1,1] == "FID" & pheno[1,2] == "IID") {
    names(pheno) <- pheno[1,]
    pheno <- pheno[-1,]
  }
  # Only keep first 2 columns and the phenotype column we care about.
  # If user does not specify pheno_column, phenotype vals are in last column.
  if (pheno_column != NULL) {
    pheno <- pheno[c(1, 2, pheno_column)]
  } else if (pheno_name != NULL) {
    pheno <- cbind(pheno[c(1, 2)], pheno[pheno_name], 
  } else {
    pheno <- pheno[c(1, 2, ncol(pheno))]
  }
  names(pheno) <- c("fid", "iid", "phenotype")
  return(pheno)
}

read_filter <- function(filter_file, filter_column = 3) {
  fil <- read.table(filter_file, as.is = T)
  # Only keep columns of filter file necessary
  fil <- fil[c(1, 2, filter_column)]
  names(fil) <- c("fid", "iid", "fil_val")
  return(fil)
}
    
read_predicted <- function(pred_exp_file) {
  pred_exp <- read.delim(pred_exp_file)
  return(pred_exp)
}

merge_and_filter <- function(pheno, pred_exp, fil = NULL, filter_val = 1) {
  # Column binds the phenotype and predicted expression dataframes and
  # filters them based on the optionally provided filter file (fil)
  # 
  # IMPORTANT: Each row in the pheno and pred_exp dfs represent people so
  # for the cbind to make sense, the row numbers in each data frame must
  # correspond to the same person.
  if (fil != NULL) {
    merged <- merge(cbind(pheno, pred_exp), fil, by = c(1, 2), sort = F)
    merged <- filter(merged, fil_val = filter_val)
  } else {
    merged <- cbind(pheno, pred_exp)
  }
  return(merged)
}

association <- function(merged, genes, test_type = "logistic") {
  assoc_df <- NULL # Init association dataframe

  # Perform test between each pred_gene_exp column and phenotype
  for (gene in genes) {
    pred_gene_exp <- merged[[gene]]
    if (test_type == "logistic") { 
      model <- glm(phenotype ~ pred_gene_exp, data = merged, family = binomial)
    } else if (test_type == "linear") {
      model <- lm(phenotype ~ pred_gene_exp, data = merged)
    } else if (test_type == "survival") {
      # TODO: survival analysis
      model <- NULL
    }
    results <- coef(summary(model))[c(2,6,8,4)]
    line <- c(gene,results)
    assoc_df <- rbind(assoc_df,line)
  }

  # Specify column names of assoc_df
  if (test_type == "logistic") {
    colnames(assoc_df) <- c("gene", "beta", "z", "p", "se(beta)")
  } else if (test_type == "linear") {
    colnames(assoc_df) <- c("gene", "beta", "t", "p", "se(beta)")
  } else if (test_type == "survival") {
    # TODO
  }
  return(assoc_df)
}

write_association <- function(assoc_df, output_file) {
  write.table(assoc_df, output_file, col.names = T, row.names = F, quote = F)
}

# Get Arguments----------------------------------------------------------------
argv <- commandArgs(trailingOnly = T)
names <- seq(1, length(argv), 2)
vals <- seq(2, length(argv), 2)
argv <- as.data.frame(
    t(as.data.frame(argv[vals], row.names = argv[names])),
    stringsAsFactors = F
)

# Set default values for arguments and set to correct data types---------------
if (argv$PHENO_FILE == NULL) {
  cat("Error: User must supply a phenotype file to for association test.\n")
  stop()
}
if (argv$PHENO_COLUMN != NULL) {
  argv$PHENO_COLUMN <- as.numeric(argv$PHENO_COLUMN)
}
if (argv$FILTER_COLUMN == NULL) {
  argv$FILTER_COLUMN <- 3
} else {
  argv$FILTER_COLUMN <- as.numeric(argv$FILTER_COLUMN)
}
if (argv$FILTER_VAL == NULL) {
  argv$FILTER_VAL <- 1
} else {
  argv$FILTER_VAL <- as.numeric(argv$FILTER_VAL)
}
if (argv$TEST_TYPE == NULL) {
  argv$TEST_TYPE <- "logistic"
}
if (argv$ONE_FLAG == NULL) {
  argv$ONE_FLAG <- FALSE
} else {
  argv$ONE_FLAG <- as.logical(argv$ONE_FLAG)
}
if (argv$MISSING_PHENOTYPE == NULL) {
  argv$MISSING_PHENOTYPE <- -9
} else {
  argv$MISSING_PHENOTYPE <- as.numeric(argv$MISSING_PHENOTYPE)
}

# Run functions----------------------------------------------------------------
# Read pheno----
cat(c(as.character(Sys.time()), "Reading phenotype file...\n"))
pheno <- read_pheno(
    argv$PHENO_FILE,
    pheno_column = argv$PHENO_COLUMN,
    pheno_name = argv$PHENO_NAME
)
# Read filter file if given----
if (argv$FILTER_FILE == NULL) {
  fil_df <- NULL
} else {
  cat(c(as.character(Sys.time()), "Reading filter file...\n"))
  fil_df <- read_filter(argv$FILTER_FILE, filter_column = argv$FILTER_COLUMN)
}
# Read Transcription----
cat(c(as.character(Sys.time()), "Reading transcription file...\n"))
pred_exp <- read_predicted(argv$PRED_EXP_FILE)

genes <- colnames(pred_exp)

cat(c(as.character(Sys.time()), "Processing data...\n"))
merged <- merge_and_filter(pheno, pred_exp, fil = fil_df, filter_val = argv$FILTER_VAL)
# Remove rows with missing phenotype data, and if doing a logistic regression,
# Make sure affected == 1 and unaffected == 0.
if (argv$TEST_TYPE == "logistic" & argv$ONE_FLAG == FALSE) {
  merged <- subset(merged, phenotype != argv$MISSING_PHENOTYPE | phenotype != 0)
  # Normal input for unaffected and affected is 1 and 2. Change to 0 and 1. 
  merged$phenotype <- merged$phenotype - 1
} else {
  merged <- subset(merged, phenotype != argv$MISSING_PHENOTYPE)
} 
cat(c(as.character(Sys.time()), "Performing association test..."))
assoc_df <- association(merged, genes, test_type = argv$TEST_TYPE)
cat(c(as.character(Sys.time()), "Writing results to file..."))
write_association(assoc_df, argv$OUT)
cat(c(as.character(Sys.time()), "Done. Results saved in", argv$OUT))
