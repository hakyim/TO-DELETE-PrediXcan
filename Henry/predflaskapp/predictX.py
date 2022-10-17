#!/usr/bin/env python

import argparse
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import MySQLdb as db
import sys

"""
parser = argparse.ArgumentParser()
parser.add_argument('--genelist', action="store", dest="genelist", default=None, help="Text file with chromosome, gene pairs.")
parser.add_argument('--dosages', action="store", dest="dosages", default="data/dosages", help="Path to a directory of gzipped dosage files.")
parser.add_argument('--dosages_prefix', dest="dosages_prefix", default="chr", action="store", help="Prefix of filenames of gzipped dosage files.")
parser.add_argument('--dosages_buffer', dest="dosages_buffer", default=None, action="store", help="Buffer size in GB for each dosage file (default: read line by line)")
parser.add_argument('--weights', action="store", dest="weights",default="data/weights.db",  help="SQLite database with rsid weights.")
parser.add_argument('--weights_on_disk', action="store_true", dest="weights_on_disk",  help="Don't load weights db to memory.")
parser.add_argument('--output', action="store", dest="output", default="output.txt", help="Path to the output file.")
args = parser.parse_args()
"""
 
#GENE_LIST = args.genelist
#DOSAGE_DIR = args.dosages
#DOSAGE_PREFIX = args.dosages_prefix
#DOSAGE_BUFFER = int(args.dosages_buffer) if args.dosages_buffer else None
BETA_FILE = "DGN-WB_0.5.db"
#PRELOAD_WEIGHTS = not args.weights_on_disk
UPLOAD_FOLDER = "./puploads."

def buffered_file(file,DOSAGE_BUFFER=None):
    if not DOSAGE_BUFFER:
        for line in file:
            yield line
    else:
        buf = ''
        while True:
            buf = buf + file.read(DOSAGE_BUFFER*(1024**3))
            if not buf:
                raise StopIteration
            last_eol = 0
            while True:
                next_eol = buf.find('\n', last_eol)
                if next_eol == -1: # No end of line here.
                    buf = buf[last_eol:]
                    break # keep the last fragment, read the next chunk
                else:
                    yield buf[last_eol:next_eol+1]
                    last_eol = next_eol + 1
                    if last_eol >= len(buf):
                        buf = ''
                        break
                            

def get_all_dosages(dosage_dir, dosage_prefix):
    for chrfile in [x for x in sorted(os.listdir(dosage_dir)) if x.startswith(dosage_prefix)]:
        print datetime.datetime.now(), "Processing %s"%chrfile
        for line in buffered_file(gzip.open(os.path.join(dosage_dir, chrfile))):
            arr = line.strip().split()
            rsid = arr[1]
            refallele = arr[4]
            dosage_row = np.array(map(float, arr[6:]))
            yield rsid, refallele, dosage_row


class WeightsDB:
    def __init__(self):
        self.conn = sqlite3.connect(BETA_FILE)

    def query(self, sql, args=None):
        c = self.conn.cursor()
        if args:
            for ret in c.execute(sql, args):
                yield ret
        else:
            for ret in c.execute(sql):
                yield ret

class GetApplicationsOf:
    def __init__(self,PRELOAD_WEIGHTS=False):
        self.db = WeightsDB()
        if PRELOAD_WEIGHTS:
            print datetime.datetime.now(), "Preloading weights..."
            self.tuples = defaultdict(list)
            for tup in self.db.query("SELECT rsid, gene, weight, eff_allele FROM weights"):
                self.tuples[tup[0]].append(tup[1:])
        else:
            self.tuples = None

    def __call__(self, rsid):
        if self.tuples:
            for tup in self.tuples[rsid]:
                yield tup
        else:

            for tup in self.db.query("SELECT gene, weight, eff_allele FROM weights WHERE rsid=?", args=[rsid]):
                yield tup

get_applications_of = GetApplicationsOf()

class TranscriptionMatrix:
    def __init__(self,gene_list=None):
        self.D = None
        self.gene_list = gene_list #a filename

    #this
    def get_gene_list(self):
        if self.gene_list:
            return list(sorted([line.strip().split()[-1] for line in open("puploads/"+self.gene_list)]))
        else:
            return [tup[0] for tup in WeightsDB().query("SELECT DISTINCT gene FROM weights ORDER BY gene")]

    def update(self, gene, weight, ref_allele, allele, dosage_row):
        if self.D is None:
            self.gene_list = self.get_gene_list()  
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            self.D = np.zeros((len(self.gene_list), len(dosage_row))) # Genes x Cases
            print "I GOT HERE"
            print self.D
        if gene in self.gene_index:            
            multiplier = 1 if ref_allele == allele else -1
            self.D[self.gene_index[gene],] += dosage_row * weight * multiplier # Update all cases for that gene


    def save(self,output_file):
        with open(output_file, 'w+') as outfile:
            outfile.write('\t'.join(self.gene_list) + '\n') # Nb. this lists the names of rows, not of columns
            for col in range(0, self.D.shape[1]):
                outfile.write('\t'.join(map(str, self.D[:,col]))+'\n')
"""
transcription_matrix = TranscriptionMatrix()
for rsid, allele, dosage_row in get_all_dosages():
    for gene, weight, ref_allele in get_applications_of(rsid):
        transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
transcription_matrix.save()        
"""
