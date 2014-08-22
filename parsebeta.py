import argparse,gzip,os,sys,datetime,StringIO,math
import numpy as np
import MySQLdb as db 
import pandas as pd


def parse_title(filename):
	step1 = filename.split(".")
	study = step1[1]
	step2 = step1[0].split("-")
	if len(step2) == 2:
		gene = step2[0]
		tissue = step2[1]
	elif len(step2) == 3: #to account for stupid gene names with '-' in them
		print step2 
		gene = step2[0] + '-' + step2[1]
		tissue = step2[2]
		print gene, tissue

	return (gene, tissue, study)

## MAIN ## 

parse_title("A2ML1-AS1-WB.GTEx.txt")
exit(1)

parser = argparse.ArgumentParser(description="Parse input/output files.")
#parser.add_argument("--betafile", help="betafile for processing ")
parser.add_argument("--filelist", help="list of beta-files to process", default="trfilelist.txt")
args = parser.parse_args()

fl = args.filelist 
filelist = open(filelist)

database = db.connect(host="localhost", # your host 
                     user="root", # your username
                      passwd="password", # your password
                      db="mysql") # name of the data base
cur = database.cursor()


for fname in filelist.readlines():
	fname = fname.strip('\n')
	gene,tissue,study = parse_title(fname)
	snpframe = pd.read_table(fname)

	for row in snpframe.iterrows():
		statement = """INSERT INTO SNPs (rsnum, eff_allele, beta, p_value, N, cis, genename, tissue, study_name) VALUES ("%s","%s",%f,%f,%d,%r,"%s","%s","%s");""" % (row[1]["SNP"],row[1]["eff.allele"],row[1]["beta"],row[1]["p.value"],row[1]["N"],row[1]["cis"],gene,tissue,study)
		cur.execute(statement)

	database.commit()
	"""
	cur.execute("SELECT * FROM SNPs;")
	print cur.fetchall()
	"""

database.close()