#!/usr/bin/env python

import argparse
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--predict', action="store_false", dest="predict", help="Include to predict gene expression")
parser.add_argument('--assoc', action="store_false", dest="assoc", help="Include to perform association test")
parser.add_argument('--genelist', action="store", dest="genelist", default=None, help="Text file with chromosome, gene pairs.")
parser.add_argument('--dosages', action="store", dest="dosages", default="data/dosages", help="Path to a directory of gzipped dosage files.")
parser.add_argument('--dosages_prefix', dest="dosages_prefix", default="chr", action="store", help="Prefix of filenames of gzipped dosage files.")
parser.add_argument('--dosages_buffer', dest="dosages_buffer", default=None, action="store", help="Buffer size in GB for each dosage file (default: read line by line)")
parser.add_argument('--weights', action="store", dest="weights",default="data/weights.db",  help="SQLite database with rsid weights.")
parser.add_argument('--weights_on_disk', action="store_true", dest="weights_on_disk",  help="Don't load weights db to memory.")
parser.add_argument('--pheno', action="store", dest="pheno", default=None, help="Phenotype file")
parser.add_argument('--pheno-name', action="store", dest="pheno_name", default=None, help="Column name of the phenotype to perform association on.")
parser.add_argument('--filter', nargs=2, action="store", dest="fil", default=None, help="Takes two arguments. First is the name of the filter file, the second is a value to filter on.")
parser.add_argument('--mfilter', action="store", dest="mfil", default=None, help="Column number of filter file to filter on.  '1' specifies the first filter column")
parser.add_argument('--output', action="store", dest="output", default="output", help="Path to output directory")

args = parser.parse_args()

GENE_LIST = args.genelist
DOSAGE_DIR = args.dosages
DOSAGE_PREFIX = args.dosages_prefix
DOSAGE_BUFFER = int(args.dosages_buffer) if args.dosages_buffer else None
BETA_FILE = args.weights
PRELOAD_WEIGHTS = not args.weights_on_disk
PHENO_FILE = args.pheno
PHENO_NAME = args.pheno_name
FILTER = args.fil
MFILTER = args.mfil
OUTPUT_DIR = args.output

PRED_EXP_FILE = os.path.join(OUTPUT_DIR, "predicted_expression.txt")
ASSOC_FILE = os.path.join(OUTPUT_DIR, "association.txt")                           

def buffered_file(file):
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

def get_all_dosages():
    for chrfile in [x for x in sorted(os.listdir(DOSAGE_DIR)) if x.startswith(DOSAGE_PREFIX)]:
        print datetime.datetime.now(), "Processing %s" % chrfile
        for line in buffered_file(gzip.open(os.path.join(DOSAGE_DIR, chrfile))):
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

    def __init__(self):
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
            for tup in self.db.query("SELECT gene, weight, eff_allele FROM weights WHERE rsid=?", (rsid,)):
                yield tup


class TranscriptionMatrix:

    def __init__(self):
        self.D = None

    def get_gene_list(self):
        if GENE_LIST:
            return list(sorted([line.strip().split()[-1] for line in open(GENE_LIST)]))
        else:
            return [tup[0] for tup in WeightsDB().query("SELECT DISTINCT gene FROM weights ORDER BY gene")]

    def update(self, gene, weight, ref_allele, allele, dosage_row):
        if self.D is None:
            self.gene_list = self.get_gene_list()  
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            self.D = np.zeros((len(self.gene_list), len(dosage_row))) # Genes x Cases
        if gene in self.gene_index:            
            multiplier = 1 if ref_allele == allele else -1
            self.D[self.gene_index[gene],] += dosage_row * weight * multiplier # Update all cases for that gene

    def save(self):
        with open(PRED_EXP_FILE, 'w+') as outfile:
            outfile.write('\t'.join(self.gene_list) + '\n') # Nb. this lists the names of rows, not of columns
            for col in range(0, self.D.shape[1]):
                outfile.write('\t'.join(map(str, self.D[:,col]))+'\n')

    def read(self):
        # Creates the transcription matrix from a file.
        with open(PRED_EXP_FILE, 'r') as infile:
            self.gene_list = infile.readline().rstrip(['\n']).split('\t')
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            levels = []
            for line in infile:
                parsed_line = line.rstrip(['\n']).split('\t')
                row = [float(level) for level in parsed_line]
                levels.append(row)
        self.D = np.array(levels)


class PhenotypeArray:

    def __init__(self, pheno_name=None, mpheno=None):
        
        self.phenodict = {}
        self.filterdict = {}
        self.pheno_name = PHENO_NAME

        if FILTER:
            # Specify which column to be selecting from.
            filter_column = MFILTER + 1 if MFILTER else -1
            with open(FILTER[0], 'r') as infile:
                for line in infile:
                    parsed_line = line.rstrip(['\n']).split()
                    self.filterdict[(parsed_line[0], parsed_line[1])] = parsed_line[filter_column]

        column_index = mpheno + 1 if mpheno else -1

        with open(PHENO_FILE, 'r') as infile:
            for line in infile:
                parsed_line = line.rstrip(['\n']).split()
                # Check if there's a header in the file.
                if parsed_line[0] == 'FID' and parsed_line[1] == 'IID':
                    if self.pheno_name:
                        column_index = parsed_line.index(self.pheno_name)
                    continue

                fid = parsed_line[0]
                iid = parsed_line[1]
                try:
                    # Only include cases matching the filter value.
                    if not FILTER:
                        self.phenodict[(fid, iid)] = parsed_line[column_index]
                    else:
                        if self.filterdict[(fid, iid)] == FILTER[1]:
                            self.phenodict[(fid, iid)] = parsed_line[column_index]
                except KeyError as e:
                    # Person in pheno file not found in filter file.
                    # Do not include this person in association.
                        print "Case %s %s from pheno file not found in filter.  Not including in association." % fid, iid
                        pass

class Association:

    def __init__(self):
        self.D = np.zeros(len(self.gene_list), 5)

    def association_test(self, transcription_matrix, phenotype, model_type='logistic'):
        pass

    def get_significant(self, alpha=0.01):
        '''
        Returns all genes with p-values under the specifed alpha value.
        '''
        pass

    def save(self):
        pass

if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)    
if not os.path.exists(PRED_EXP_FILE):
    get_applications_of = GetApplicationsOf()
    transcription_matrix = TranscriptionMatrix()
    for rsid, allele, dosage_row in get_all_dosages():
        for gene, weight, ref_allele in get_applications_of(rsid):
            transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
    transcription_matrix.save()
