#install.packages("/nas40t2/egamazon/matrixeqtl/MatrixEQTL_1.6.2.tar.gz", repos = NULL, type="source")
library(MatrixEQTL)


"%&%" = function(a,b) paste(a,b,sep="")

args <- commandArgs(TRUE);

## Settings
chr = as.numeric(args[1]);
 
#mydir = "/nas40t2/egamazon/PREDIXCAN/"; 
mydir = ""; 

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR ; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = mydir %&% "DGN.imputed.SNPxID" %&% chr;
snps_location_file_name = mydir %&% "snp_location_info.MOD";

# Gene expression file name
expression_file_name = mydir %&% "DGN-WB.forMatrixeQTL.GENExID";
gene_location_file_name = mydir %&% "gene_location_info";

# Covariates file name
# Set to character() for no covariates
# covariates_file_name = "covariates_new/lung.txt";
covariates_file_name = character();

# Output file name
output_file_name_cis = "OUTPUT/eQTL_results_R_cis.DGN-WB." %&% chr  %&% ".txt";
output_file_name_tra = "OUTPUT/eQTL_results_R_tra.DGN-WB." %&% chr  %&% ".txt";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 0; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 0; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 0; # one row of column labels
cvrt$fileSkipColumns = 0; # one column of row labels
cvrt$fileSliceSize = 2000; # read file in one piece
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = FALSE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = FALSE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot");

## Plot the Q-Q plot of local and distant p-values
pdf('DGN-WB.' %&% chr  %&% '.pdf')
plot(me)
dev.off()

