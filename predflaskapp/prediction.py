import predictX as px 
import MySQLdb as db 

class prediction_maker:
	def __init__(self,gene_list,dosage_dir,dosage_prefix,dosage_buffer=None):
		self.GENE_LIST = gene_list
		self.DOSAGE_DIR = dosage_dir
		self.DOSAGE_PREFIX = dosage_prefix
		self.DOSAGE_BUFFER = int(dosage_buffer) if dosage_buffer else None
		self.database = px.WeightsDB()

	def predictions(self):
		get_applications_of = px.GetApplicationsOf()
		transcription_matrix = px.TranscriptionMatrix()
		
		for rsid, allele, dosage_row in pxget_all_dosages():
			for gene, weight, ref_allele in get_applications_of(rsid):
				transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
		transcription_matrix.save("PredXResults.txt")#hard coded for the moment`        
