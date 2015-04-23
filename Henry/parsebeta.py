import argparse,gzip,os,sys,datetime,StringIO,math
import numpy as np
import MySQLdb as db 
import pandas as pd

#ABC11-48400900C8.1-WB.GTEx.txt
def parse_title(filename,type='regular'):
        if type == 'regular':
                step1 = filename.rsplit("-",1)
                #print step1
                gene = step1[0]
                step2 = step1[1].split(".")
                #print step2
                tissue = step2[0]
                study = step2[1]
                return (gene, tissue, study)
        elif type == 'LASSO':
                step1 = filename.rsplit('.',1)[0].rsplit('-',2)
                gene = step1[0]
                study = step1[1]
                tissue = step1[2]
                """
                step1 = filename.split('-')
                gene = step1[0]
                study = step1[1]
                tissue = step1[2].split('.')[0] #parameter gore 
                """
                return (gene, tissue, study) 

## MAIN ## 


parser = argparse.ArgumentParser(description="Parse input/output files.")
parser.add_argument("--betapath", help="path for betafiles")
parser.add_argument("--filelist", help="list of beta-files to process", default="trfilelist.txt")
parser.add_argument("--filetype", help="format of beta files: regular vs LASSO",default='regular')
args = parser.parse_args()

fl = args.filelist 
pathname = args.betapath
filetype = args.filetype 

filelist = open(fl)

database = db.connect(host="localhost", # your host 
                     user="root", # your username
                      passwd="mathtype5", # your password
                      db="mysql") # name of the data base
cur = database.cursor()

print filetype
for fname in filelist.readlines():
	fname = fname.strip('\n')
	fullpath = pathname + fname
	gene,tissue,study = parse_title(fname,filetype)
        #print (gene, tissue, study)
        #raw_input("good?")
        if tissue == "GTEx":
                print fname 
                raw_input("bad parsing, Proceed?")

	try:
		snpframe = pd.read_table(fullpath)
	except:
		print "error in read_table", sys.exc_info()[0]
		print fullpath
		continue       

        if filetype == "regular" and (snpframe.cis.isnull()[0]):
                snpframe.cis = snpframe.N
                snpframe.N = snpframe.N.map(lambda x: np.nan if x == True else x)
                snpframe = snpframe.where(pd.notnull(snpframe), None)
                nullparts = True
        else:
                nullparts = False

	for row in snpframe.iterrows():

                if row[1]["beta"] == 'Inf': #skip rows with inf in beta  
                        continue
                        
                if nullparts == True:
                        numparts = "Null"
                elif filetype=='regular':
                        numparts = row[1]["N"]

                statement = """INSERT INTO SNPs (rsnum, eff_allele, beta, p_value, N, cis, genename, tissue, study_name, version) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);"""
                
                if filetype=='regular':
                        try:
                                cur.execute(statement,(row[1]["SNP"],row[1]["eff.allele"],row[1]["beta"],row[1]["p.value"],numparts,row[1]["cis"],gene,tissue,study,1) )                             
                        except:
                                err = sys.exc_info()[0]
                                print "Error: %s" % err
                                print "On file: %s" % fname 
                                continue
                elif filetype=='LASSO':
                        try:
                                print gene, tissue, study,row[1]["SNP"],row[1]["eff.allele"],row[1]["beta"] 
                                cur.execute(statement,(row[1]["SNP"],row[1]["eff.allele"],row[1]["beta"],"Null","Null","Null",gene,tissue,study, 2) )
                                 
                        except:
                                err = sys.exc_info()[1]
        
                                raw_input('whatif?')
                                print "Error: %s" % err
                                print "On file: %s" % fname 
                                continue
                        
               
                """
                cur.execute("SELECT * FROM SNPs where genename = %s;", gene)
                print cur.fetchall()
                raw_input("continue2?")
                """
database.commit()
database.close()
