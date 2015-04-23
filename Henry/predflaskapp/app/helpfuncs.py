import os
import sys
import MySQLdb as db 

def _getTissueTypes(database):
	"""
	database = db.connect(host="192.170.232.66", # your host 
                     user='public', # your username
                      passwd='foobar', # your password
                      db="mysql",port=3306) # name of the data base
	"""
	cur = database.cursor()
	statement = "SELECT distinct tissue from SNPs;"
	cur.execute(statement)
	tissues = cur.fetchall()
	tissueList = [(t[0],t[0]) for t in tissues]
	return tissueList


def _getStudyNames(database):
	"""
	database = db.connect(host="192.170.232.66", # your host 
                     user='public', # your username
                      passwd='foobar', # your password
                      db="mysql",port=3306) # name of the data base
	"""
	cur = database.cursor()
	statement = "SELECT distinct study_name from SNPs;"
	cur.execute(statement)
	studies = cur.fetchall()
	studyList = [(t[0],t[0]) for t in studies] #(key,value)
	return studyList

def _generateCommand(phenopath,genedatapath,genotypeheader,genotypetail,tissue,study):
	prologue = "python compute.scores.transcriptome.py "
	pheno = "-pfn " + phenopath + " "
	gened = "-gdfn %s " % (genedatapath) 
	genohead = "--genoheader %s " % (genotypeheader)
	genotail = "--genotail %s " % (genotypetail)
	betahead = "--betaheader default/ "
	betatail = "--betatail -%s.%s.txt " % (tissue, study)
	db = "--db db --dbUser public --dbPass foobar "
	outfile = "--outfile computeOut.txt "
	return prologue + pheno + gened + genohead + genotail + betahead + betatail + outfile + db