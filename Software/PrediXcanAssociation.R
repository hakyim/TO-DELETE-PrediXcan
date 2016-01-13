
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
    merged <- filter(merged, FIL_VAL = filter_val)
  } else {
    merged <- cbind(pheno, pred_exp)
  }
  return(merged)
}

association <- function(merged, genes, test_type = "logistic") {
  assoc_df <- NULL # Init association dataframe
  # Perform test between each pred_gene_exp column and phenotype----
  for (gene in genes) {
    pred_gene_exp <- merged[[gene]]
    if (test_type == "logistic") { 
      model <- glm(phenotype ~ pred_gene_exp, data = merged, family = binomial)
    } else if (test_type == "linear") {
      model <- lm(phenotype ~ pred_gene_exp, data = merged)
    } else if (test_type == "survival") {
      # TODO
      model <- NULL
    }
    results <- coef(summary(model))[c(2,6,8,4)]
    line <- c(gene,results)
    assoc_df <- rbind(assoc_df,line)
  }

  # Specify column names----
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

# Get Arguments
argv <- commandArgs(trailingOnly = T)

# TODO: Parse Arguments/Assign to appropriate variables.
# If variables missing, assign default values here.

# Run functions.
pheno <- read_pheno(PHENO_FILE, pheno_column = PHENO_COLUMN, pheno_name = PHENO_NAME)
if (FILTER_FILE == NULL) {
  fil_df <- NULL
} else {
  fil_df <- read_filter(FILTER_FILE, filter_column = FILTER_COLUMN)
}
pred_exp <- read_predicted(PRED_EXP_FILE)
genes <- colnames(pred_exp)
merged <- merge_and_filter(pheno, pred_exp, fil = fil_df, filter_val = FILTER_VAL)
# Remove rows with missing phenotype data, and if doing a logistic regression,
# Make sure affected = 1 and unaffected = 0.
if (TEST_TYPE == "logistic" & ONE_FLAG == FALSE) {
  merged <- subset(merged, phenotype != MISSING_PHENOTYPE | phenotype != 0)
  merged$phenotype <- merged$phenotype - 1
} else {
  merged <- subset(merged, phenotype != MISSING_PHENOTYPE)
} 
assoc_df <- association(merged, genes, test_type = TEST_TYPE)
write_association(assoc_df, OUT)
