import predictX as px 
import MySQLdb as db 
from app import app
import time
import os


"""class wrapper to run the core prediction code"""
class prediction_maker:
	def __init__(self,gene_list,dosage_dir,dosage_prefix,dosage_buffer=None):
		self.GENE_LIST = gene_list
		self.DOSAGE_DIR = dosage_dir
		self.DOSAGE_PREFIX = dosage_prefix
		self.DOSAGE_BUFFER = int(dosage_buffer) if dosage_buffer else None
		self.database = px.WeightsDB()

	def do_predictions(self):
		get_applications_of = px.GetApplicationsOf()
		transcription_matrix = px.TranscriptionMatrix(self.GENE_LIST)
		
		for rsid, allele, dosage_row in px.get_all_dosages(self.DOSAGE_DIR,self.DOSAGE_PREFIX): 
			for gene, weight, ref_allele in get_applications_of(rsid):
				print "doing a rsid"
				transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
		fname = self.DOSAGE_PREFIX  + time.asctime().replace(" ", "") + ".txt"
		transcription_matrix.save(os.path.join(app.static_folder,fname)) #hard coded for the moment        
		return fname