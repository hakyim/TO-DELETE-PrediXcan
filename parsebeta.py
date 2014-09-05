import argparse,gzip,os,sys,datetime,StringIO,math
import numpy as np
import MySQLdb as db 
import pandas as pd

#ABC11-48400900C8.1-WB.GTEx.txt
def parse_title(filename):
	step1 = filename.rsplit("-",1)
	#print step1
	gene = step1[0]
	step2 = step1[1].split(".")
	#print step2
	tissue = step2[0]
	study = step2[1]
	return (gene, tissue, study)

## MAIN ## 


parser = argparse.ArgumentParser(description="Parse input/output files.")
parser.add_argument("--betapath", help="path for betafiles")
parser.add_argument("--filelist", help="list of beta-files to process", default="trfilelist.txt")
args = parser.parse_args()

fl = args.filelist 
pathname = args.betapath

filelist = open(fl)

database = db.connect(host="localhost", # your host 
                     user="root", # your username
                      passwd="password", # your password
                      db="mysql") # name of the data base
cur = database.cursor()


for fname in filelist.readlines():
	fname = fname.strip('\n')
	fullpath = pathname + fname
	gene,tissue,study = parse_title(fname)
	try:
		snpframe = pd.read_table(fullpath)
	except:
		print "error in read_table", sys.exc_info()[0]
		print fullpath
		continue 
        
        #fuck. need to account for missing N. DAMN IT.
        if snpframe.cis.isnull()[0]:
                snpframe.cis = snpframe.N
                snpframe.N = snpframe.N.map(lambda x: np.nan if x == True else x)
                snpframe = snpframe.where(pd.notnull(snpframe), None)

	for row in snpframe.iterrows():
		statement = """INSERT INTO SNPs (rsnum, eff_allele, beta, p_value, N, cis, genename, tissue, study_name) VALUES ("%s","%s",%f,%f,%s,%r,"%s","%s","%s");""" % (row[1]["SNP"],row[1]["eff.allele"],row[1]["beta"],row[1]["p.value"],"Null",row[1]["cis"],gene,tissue,study)
		cur.execute(statement)

	database.commit()
	"""
	cur.execute("SELECT * FROM SNPs where study_name = 'GD' ;")
	print cur.fetchall()
	"""
	
database.close()
