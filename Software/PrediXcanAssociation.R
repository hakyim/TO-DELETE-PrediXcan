#!/usr/bin/env Rscript
options(digits.secs=6)

read_pheno <- function(pheno_file) {
  pheno <- read.table(pheno_file, header = F, as.is = T)
  # Fix dataframe if there is a header row. Returned df will have apropos column names
  if (pheno[1,1] == "FID" & pheno[1,2] == "IID") {
    names(pheno) <- pheno[1,]
    pheno <- pheno[-1,]
  }
  return(pheno)
}

reduce_pheno <- function(pheno, pheno_name = colnames(pheno)[ncol(pheno)]) {
  # Only keep first 2 columns and the phenotype column we care about.
  # If user does not specify pheno column, phenotype vals are in last column.
  if (!is.null(pheno_name)) {
    pheno <- cbind(pheno[c(1, 2)], pheno[pheno_name])
  } else {
    pheno <- pheno[c(1, 2, ncol(pheno))]
  }
  names(pheno) <- c("fid", "iid", "phenotype")
  pheno$phenotype <- as.numeric(pheno$phenotype)
  return(pheno)
}

read_filter <- function(filter_file, filter_column = 3) {
  fil <- read.table(filter_file, header = F, as.is = T)
  # Fix dataframe if there is a header row.
  if (fil[1,1] == "FID" & fil[1,2] == "IID") {
    names(fil) <- fil[1,]
    fil <- fil[-1,]
  }
  # Only keep columns of filter file necessary
  fil <- fil[c(1, 2, filter_column)]
  names(fil) <- c("fid", "iid", "fil_val")
  return(fil)
}
    
read_predicted <- function(pred_exp_file) {
  pred_exp <- read.table(pred_exp_file, header = T, as.is = T)
  return(pred_exp)
}

merge_and_filter <- function(pheno, pred_exp, fil = NULL, filter_val = 1) {
  # Column binds the phenotype and predicted expression dataframes and
  # filters them based on the optionally provided filter file (fil)
  # 
  # IMPORTANT: Each row in the pheno and pred_exp dfs represent people so
  # for the cbind to make sense, the row numbers in each data frame must
  # correspond to the same person.
  npheno <- nrow(pheno)
  npred_exp <- nrow(pred_exp)
  cat(c(as.character(Sys.time()), npheno, "individuals found in phenotype file.\n"))
  cat(c(as.character(Sys.time()), npred_exp, "individuals found in predicted expression file.\n"))
  merged <- merge(pheno, pred_exp, by = c(1, 2), sort = F)
  nmerged <- nrow(merged)
  cat(c(as.character(Sys.time()), nmerged, "individuals common to both.  Performing assoc on intersection.\n"))
  if (!is.null(fil)) {
    merged <- merge(merged, fil, by = c(1, 2), sort = F)
    merged <- merged[merged$fil_val == filter_val, ]
    nfilter <- nrow(merged)
    cat(c(as.character(Sys.time()), nfilter, "individuals remain for association after filtering.\n"))
  }
  return(merged)
}

association <- function(merged, genes, test_type = "linear") {
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
  return(as.data.frame(assoc_df))
}

write_association <- function(assoc_df, output_file) {
  write.table(assoc_df, output_file, col.names = T, row.names = F, quote = F)
}

# Get Arguments----------------------------------------------------------------
argv <- commandArgs(trailingOnly = T)
# Make data.frame of arguments by selecting odd numbered argvs as options
# and even numbered argvs as values of options
names <- seq(1, length(argv), 2)
vals <- seq(2, length(argv), 2)
argv <- as.data.frame(
    t(as.data.frame(argv[vals], row.names = argv[names])),
    stringsAsFactors = F
)

# Set default values for arguments and set to correct data types---------------
if (is.null(argv$PHENO_FILE)) {
  cat(c(as.character(Sys.time()), "Error: User must supply a phenotype file to for association test.\n"))
  stop()
}
# Default PHENO_COLUMN: NULL
if (!is.null(argv$PHENO_COLUMN)) {
  if (argv$PHENO_COLUMN != 'None') {
    argv$PHENO_COLUMN <- suppressWarnings(as.numeric(argv$PHENO_COLUMN))
  } else {
    argv$PHENO_COLUMN <- NULL
  }
}
# Default PHENO_NAME: NULL
if (!is.null(argv$PHENO_NAME)) {
  if (argv$PHENO_NAME == 'None') {
    argv$PHENO_NAME <- NULL
  }
}
# Default FILTER_COLUMN: 3
if (is.null(argv$FILTER_COLUMN) | argv$FILTER_COLUMN == 'None') {
  argv$FILTER_COLUMN <- 3
} else {
  argv$FILTER_COLUMN <- suppressWarnings(as.numeric(argv$FILTER_COLUMN) + 2)
}
# Default FILTER_VAL: 1
if (is.null(argv$FILTER_VAL)) {
  argv$FILTER_VAL <- 1
} else {
  argv$FILTER_VAL <- suppressWarnings(as.numeric(argv$FILTER_VAL))
}
# Default TEST_TYPE: linear
if (is.null(argv$TEST_TYPE)) {
  argv$TEST_TYPE <- "linear"
}
# Default MISSING_PHENOTYPE: NA
if (is.null(argv$MISSING_PHENOTYPE)) {
  argv$MISSING_PHENOTYPE <- NA
} else {
  argv$MISSING_PHENOTYPE <- suppressWarnings(as.numeric(argv$MISSING_PHENOTYPE))
}

# Run functions----------------------------------------------------------------
# Read pheno-------------------------------------------------------------------
cat(c(as.character(Sys.time()), "Reading phenotype file...\n"))
pheno <- read_pheno(argv$PHENO_FILE)

# Parse info on pheno df-------------------------------------------------------
if (is.null(argv$PHENO_COLUMN) & is.null(argv$PHENO_NAME)) {
  argv$PHENO_NAME <- colnames(pheno)[ncol(pheno)]
} else if (is.null(argv$PHENO_NAME)) {
  argv$PHENO_NAME <- colnames(pheno)[argv$PHENO_COLUMN]
} else if (!is.null(argv$PHENO_NAME) & !is.null(argv$PHENO_COLUMN)) {
  # Error checking
  if (argv$PHENO_NAME != colnames(pheno)[argv$PHENO_COLUMN]) {
    cat(c(as.character(Sys.time()), "Mismatch between phenotype name and column number.\n"))
    stop()
  }
}
pheno <- reduce_pheno(pheno, pheno_name = argv$PHENO_NAME)

# Read filter file if given----------------------------------------------------
if (argv$FILTER_FILE == 'None' | is.null(argv$FILTER_FILE)) {
  fil_df <- NULL
} else {
  cat(c(as.character(Sys.time()), "Reading filter file...\n"))
  fil_df <- read_filter(argv$FILTER_FILE, filter_column = argv$FILTER_COLUMN)
}

# Read Transcription-----------------------------------------------------------
cat(c(as.character(Sys.time()), "Reading transcription file...\n"))
pred_exp <- read_predicted(argv$PRED_EXP_FILE)

genes <- colnames(pred_exp)[c(-1, -2)] # First 2 cols are FID and IID

cat(c(as.character(Sys.time()), "Processing data...\n"))
merged <- merge_and_filter(pheno, pred_exp, fil = fil_df, filter_val = argv$FILTER_VAL)
# Remove rows with missing phenotype data
if (is.na(argv$MISSING_PHENOTYPE)) {
  merged <- merged[which(!is.na(merged$phenotype)), ]
} else {
  merged <- merged[which(merged$phenotype != argv$MISSING_PHENOTYPE), ]
}

# If all rows filtered out, throw error.
if (nrow(merged) == 0) {
  cat(c(as.character(Sys.time()), "ERROR: Filtered out all rows of phenotype\n"))
  stop()
}

# If logistic, and more than 2 values for phenotype, throw error.
if (argv$TEST_TYPE == 'logistic' & length(table(merged$phenotype)) > 2) {
  cat(c(as.character(Sys.time()), "ERROR: For logistic tests, phenotype column can only have 2 values: 0 - unaffected, 1 - affected\n"))
  stop()
}

# Association Tests------------------------------------------------------------
cat(c(as.character(Sys.time()), "Performing association test...\n"))
assoc_df <- association(merged, genes, test_type = argv$TEST_TYPE)
cat(c(as.character(Sys.time()), "Writing results to file...\n"))
write_association(assoc_df, argv$OUT)
cat(c(as.character(Sys.time()), "Done. Results saved in", argv$OUT, "\n"))
