###!/opt/python/bin/python -O
## Hae Kyung Im
## This script go through the imputed files, compute the polygenic scores for all genes on the list and save a matrix of predicted gene expression levels for individuals in the phenotype file and genes in the gene list.

## This needs imputed files with dosages
## format: chr, rsid, A1, A2, MAF, n columns with dosages for each individuals

import argparse,gzip,os,sys,datetime,StringIO,math
import numpy as np
import MySQLdb as db 


## INPUT PARAMETERS/FILENAMES


## define paths
bios = False
if bios:
     prebios = '/glusterfs/users/im_lab/' 
else:
     prebios = '/'

data_dir = prebios + 'nas40t2/haky/Signatures/data/'
annot_dir = prebios + 'nas40t2/haky/main/Annotations/'

disease_name = 'CD'

## PHENOFILENAMEfirewall
phenofilename = data_dir + 'cohorts/WTCCC1/'+ disease_name + '/imputed/' + 'Affx_sample_' + disease_name + '.txt' 

## GENE LIST FILENAME / GENEDATAFILENAME
genedatafilename = prebios + 'nas40t2/haky/main/Data/landmark/Landmark-Genes-n978-2.txt'

## GENOTYPE/DOSAGE FILE FORM: header + chr + tail
genoheader = data_dir + 'cohorts/WTCCC1/'+ disease_name + '/imputed/'+ disease_name + '_chr'
genotail =  '_aligned_imputed' + '.dos.gz'

## BETA FILE NAMES FORM: betaheader + gene_name + betatail
source = 'GTEx'
betaheader = data_dir + 'betas/transcriptome/'
betatail = '-WB.' + source  + '.txt'

## PREDICTED FILE NAMES / OUTFILENAME
outfilename = data_dir + 'cohorts/WTCCC1/'+ disease_name  + '/score/GENEX.GTEx.ver.txt.profile'

## EXCLUDE SNPsoffice space
excludeSNPfilename = data_dir + 'cohorts/WTCCC1/' + 'exclusion-list-snps-26_04_2007.txt'



## ----------------------------------------
## END OF PARAMETERS THAT NEED TO BE MODIFIED
## ----------------------------------------

## ----------------------------------------
# FUNCTION DEFINITIONS
## ----------------------------------------
"""Parses beta file names for the gene, tissue, and study it records"""
def parse_title(filename):
    step1 = filename.rsplit("-",1)
    gene = step1[0]
    step2 = step1[1].split(".")
    tissue = step2[0]
    study = step2[1]
    return (gene, study, tissue)


def readArray(fname,delim=None):
    print(fname)
    part = fname.split('.')
    isgz = part[-1] == 'gz'
    if isgz:
        fin=gzip.open(fname,'rb')
    else:
        fin = open(fname,'rb')
    ret=[]
    ret=[x.strip().split(delim) for x in fin if(not x.startswith(r'#'))]
    fin.close()
    return ret

def readHeader(fname,delim=None):
    fin = open(fname,'rb')
    ret=fin.readline().strip().split(delim)
    while(ret[0][0]=='#'):
        ret[0] = ret[0][1:]
    fin.close()
    return ret

def writeArray(fname,arr,delim='\t'):
    part = fname.split('.')
    print(part)
    isgz = part[-1] == 'gz'
    if isgz:
        fout=gzip.open(fname,'wb')
        print 'opening gzip'
    else:
        print 'opening regular file'
        fout = open(fname,'wb')
    for x in arr:
        fout.write(delim.join(x)+'\n')
    fout.close()
    return True

def convertToNPArray(tup):
     newlist = list(tup)
     newarray = np.asarray(newlist)
     assert(newarray[2] == tup[2])
     return newarray
## ----------------------------------------
## END OF FUNCTION DEFINITIONS
## ----------------------------------------

## ----------------------------------------
## MAIN CODE STARTS HERE 
## -------------------- 
##Parse command line arguments for file paths. Defaultsfirewallfirewall defined above

parser = argparse.ArgumentParser(description="Parse input/output files.")
parser.add_argument("-pfn", default=phenofilename, help="Pheno File path")
parser.add_argument("-gdfn",default=genedatafilename, help = "Gene list file path")
parser.add_argument("--genoheader", default=genoheader, help="header of genotype/dosage file") 
parser.add_argument("--genotail", default=genotail, help="tail of genotype/dosage file") 
parser.add_argument("--betaheader", default=betaheader, help="header of beta file") 
parser.add_argument("--betatail", default=betatail, help="tail of beta file")
parser.add_argument("--outfile", default=outfilename, help="Outfile path")
parser.add_argument("--excludeSNPfilename", default=excludeSNPfilename, help="ExSNP file path")
parser.add_argument("--db", default=0, help="Get beta from DB or from local files")
parser.add_argument("--dbUser", required=True, help="username for DB")
parser.add_argument("--dbPass", required=True, help="password for DB")

args = parser.parse_args()

phenofilename = args.pfn
genedatafilename = args.gdfn 
genoheader = args.genoheader
genotail = args.genotail
betaheader = args.betaheader
betatail = args.betatail
outfilename = args.outfile
excludeSNPfilename = args.excludeSNPfilename
useDB = args.db
usern = args.dbUser
password = args.dbPass

## read phenotype file
phenoheader = readHeader(phenofilename)
phenodata = readArray(phenofilename,'\t')
nsamp = len(phenodata)

## read gene list to compute predictions
genedata = readArray(genedatafilename,'\t')
genedata = np.asarray(genedata[1:])
genelist = genedata[:,1]

## list of SNPs to be excluded
## claimsclaimsTODO: WOULD A DICTIONARY WORK FASTER THAN LIST/ARRAY?
excludeSNPdata = readArray(excludeSNPfilename)
excludeSNPdata = np.asarray(excludeSNPdata)
excludeSNPdata = excludeSNPdata[1:] ## exclude title
excludeSNPlist = np.unique(excludeSNPdata[:,2])
exsnpindex = {}
for ss in range(len(excludeSNPlist)):
    SNP = excludeSNPlist[ss]
    exsnpindex[SNP] = ''

## INDEX ALL BETA FILES IN GENELIST
indexindex = {}

database = db.connect(host="192.170.232.66", # your host 
                     user=usern, # your username
                      passwd=password, # your password
                      db="mysql",port=3306) # name of the data base
cur = database.cursor()
statement = """SELECT * FROM SNPs where genename = %s AND study_name = %s AND tissue = %s;"""

## READ BETA FILE
if useDB == 'db':
     print "Using DB method"
     for gg in genelist:
          betafilename =  betaheader + gg + betatail
          gene, study, tissue = parse_title(gg + betatail)
          
          cur.execute(statement, (gene,study,tissue))
          betarray = cur.fetchall()
          betaindex = {}
          for beta in betarray:    
               conv = convertToNPArray(beta)
               if(conv[2] == beta[2]):
                    pass 
               else:
                    print "bad conversion"
                    print conv[2]
                    print beta[2]
                    raw_input("Continue? Press ctrl-C to kill.  ")
               rsid = conv[0]
               betaindex[rsid] = conv
          indexindex[gg] = betaindex
else:
     print "Using local file method"
     for gg in genelist:
          betafilename =  betaheader + gg + betatail
          if(os.path.isfile(betafilename)):       #should be changed to query for gene  + study + tissue
               betalistdata = readArray(betafilename)
               betarray = np.asarray(betalistdata[1:]) ## exclude title
               nsnps = len(betarray)
        ## INDEX BETA FILE
               betaindex = {}
               
               for rr in range(nsnps):
                    rsid = betarray[rr,0]
                    betaindex[rsid] = betarray[rr,:]
               indexindex[gg] = betaindex
               
## NEW GENELIST, ONLY THOSE THAT HAVE PREDICTIVE MODELS
print('old ngen ' + str(len(genelist)))
genelist = indexindex.keys()
ngen = len(genelist) ## genes in gendata that also have predictive models

## PREDARRAY(ngen,nsamp)  MATRIX OF PREDICTIONS/POLYSCORES, will transpose at the end
predarray = np.zeros((ngen,nsamp))
## LOOP OVER CHR
for cc in range(1,23):
    chr = str(cc).zfill(2)
    infilename = genoheader + chr + genotail
    print(infilename)
    ## READ IMPUTED DOSAGES, GO THROUGH ROWS AND COMPUTE CONTRIBUTION TO POLYSCORE
    dosagefile = gzip.open(infilename)
    for line in dosagefile:
         part = line.split(None,6)
         rsid = part[1]
         # IF RSID IS NOT IN EXCLUSION LIST
         if not(rsid in exsnpindex):  
              refalele = part[3]
              ## CHECK FOR ALL GENES WHETHER RSID IS EQTL
              numarrayed = False
              dosagerow = []
              cont = 0
              ## LOOP OVER GENES IN GENELIST
              for gg in genelist:
                   betaindex = indexindex[gg]
                   if rsid in betaindex:
                        betaA1 = betaindex[rsid][1]
                        print(betaA1 + ' ' + refalele)
                        if refalele == betaA1:
                             beta = float(betaindex[rsid][2])
                        else:
                             beta = -float(betaindex[rsid][2])
                    ## dosagerow
                        if(np.logical_not(numarrayed)): 
                             dosagerow = np.asarray(part[6:][0].split(),float)
                             numarrayed = True
                        predarray[cont,] += beta * dosagerow
                   cont += 1
    dosagefile.close()

## SAVE PREDICTIONS TO FILE
outarray = predarray.transpose()
outarray = np.vstack((genelist,outarray))
phenoarray = np.vstack((phenoheader,phenodata)) ## np.vstack turns lists into arrays
outarray = np.hstack((phenoarray,outarray))

## TODO use the same format as other scores FID, IID, etc.
writeArray(outfilename,outarray)

print('FINISHED SUCCESSFULLY ' + disease_name)
