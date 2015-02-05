import MySQLdb as db 

def _getTissueTypes():
	database = db.connect(host="192.170.232.66", # your host 
                     user='public', # your username
                      passwd='foobar', # your password
                      db="mysql",port=3306) # name of the data base
	cur = database.cursor()
	statement = "SELECT distinct tissue from SNPs;"
	cur.execute(statement)
	tissues = cur.fetchall()
	tissueList = [(t[0],t[0]) for t in tissues]
	return tissueList


def _getStudyNames():
	database = db.connect(host="192.170.232.66", # your host 
                     user='public', # your username
                      passwd='foobar', # your password
                      db="mysql",port=3306) # name of the data base
	cur = database.cursor()
	statement = "SELECT distinct study_name from SNPs;"
	cur.execute(statement)
	studies = cur.fetchall()
	studyList = [(t[0],t[0]) for t in studies] #(key,value)
	return studyList

def _generateCommand(phenopath,genedatapath,genotypepath,tissue,study):
	pass 