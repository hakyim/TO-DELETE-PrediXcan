#!/usr/bin/env python

import argparse
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import sys
import subprocess


def buffered_file(file, dosage_buffer=None):
    if not dosage_buffer:
        for line in file:
            yield line
    else:
        buf = ''
        while True:
            buf = buf + file.read(dosage_buffer*(1024**3))
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


def get_all_dosages(dosage_dir, dosage_prefix, dbuffer=None):
    for chrfile in [x for x in sorted(os.listdir(dosage_dir)) if x.startswith(dosage_prefix) and x.endswith(".gz")]:
        print("{} Processing {}".format(datetime.datetime.now(), chrfile))
        for line in buffered_file(gzip.open(os.path.join(dosage_dir, chrfile)), dosage_buffer=dbuffer):
            arr = line.decode('utf-8').strip().split()
            rsid = arr[1]
            eff_allele = arr[4]
            dosage_row = np.array(arr[6:], dtype=np.float64)
            yield rsid, eff_allele, dosage_row


class WeightsDB:
    def __init__(self, beta_file):
        self.conn = sqlite3.connect(beta_file)

    def query(self, sql, args=None):
        c = self.conn.cursor()
        if args:
            for ret in c.execute(sql, args):
                yield ret
        else:
            for ret in c.execute(sql):
                yield ret


class GetApplicationsOf:
    def __init__(self, beta_file, preload_weights=True):
        self.db = WeightsDB(beta_file)
        if preload_weights:
            print("{} Preloading weights...".format(datetime.datetime.now()))
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
    def __init__(self, beta_file, sample_file, gene_file=None):
        self.D = None
        self.beta_file = beta_file
        self.gene_file = gene_file
        self.sample_file = sample_file
        self.complements = {"A":"T","C":"G","G":"C","T":"A"}

    def get_gene_list(self):
        if self.gene_file:
            return list(sorted([line.strip().split()[-1] for line in open(self.gene_file)]))
        else:
            return [tup[0] for tup in WeightsDB(self.beta_file).query("SELECT DISTINCT gene FROM weights ORDER BY gene")]

    def update(self, gene, weight, ref_allele, allele, dosage_row):
        if self.D is None:
            self.gene_list = self.get_gene_list()
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            self.D = np.zeros((len(self.gene_list), len(dosage_row))) # Genes x Cases
        if gene in self.gene_index: #assumes dosage coding 0 to 2           
            if ref_allele == allele or self.complements[ref_allele] == allele: # assumes non-ambiguous SNPs to resolve strand issues: 
                self.D[self.gene_index[gene],] += dosage_row * weight
            else:
                self.D[self.gene_index[gene],] += (2-dosage_row) * weight # Update all cases for that gene 


    def get_samples(self):
        with open(self.sample_file, 'r') as samples:
            for line in samples:
                yield [line.split()[0], line.split()[1]]
                    
    def save(self, pred_exp_file):
        sample_generator = self.get_samples()
        with open(pred_exp_file, 'w+') as outfile:
            outfile.write('FID\t' + 'IID\t' + '\t'.join(self.gene_list) + '\n') # Nb. this lists the names of rows, not of columns
            for col in range(0, self.D.shape[1]):
                try:
                    outfile.write('\t'.join(next(sample_generator)) + '\t' + '\t'.join(map(str, self.D[:,col]))+'\n')
                except StopIteration:
                    print("ERROR: There are not enough rows in your sample file!")
                    print("Make sure dosage files and sample files have the same number of individuals in the same order.")
                    os.remove(pred_exp_file)
                    sys.exit(1)
            try:
                next(sample_generator)
            except StopIteration:
                print("Predicted expression file complete!")
            else:
                print("ERROR: There are too many rows in your sample file!")
                print("Make sure dosage files and sample files have the same number of individuals in the ame order.")
                os.remove(pred_exp_file)
                sys.exit(1)

def check_out_file(out_file):
    try:
        test_fo = open(out_file, 'w')
        test_fo.close()
    except IOError:
        print("ERROR: Cannot open {} for writing. ".format(out_file) +
              "Make sure path to file exists.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--predict', action="store_true", dest="predict", default=False, help="Include to predict gene expression")
    parser.add_argument('--assoc', action="store_true", dest="assoc", default=False, help="Include to perform association test")
    parser.add_argument('--genelist', action="store", dest="genelist", default=None, help="Text file with chromosome, gene pairs.")
    parser.add_argument('--dosages', action="store", dest="dosages", default="data/dosages", help="Path to a directory of gzipped dosage files.")
    parser.add_argument('--dosages_prefix', dest="dosages_prefix", default="chr", action="store", help="Prefix of filenames of gzipped dosage files.")
    parser.add_argument('--dosages_buffer', dest="dosages_buffer", default=None, action="store", help="Buffer size in GB for each dosage file (default: read line by line)")
    parser.add_argument('--samples', dest='sample_file', default="samples.txt", action="store", help="File in dosages directory with individual ids.  Must be in same order as columns for dosages")
    parser.add_argument('--weights', action="store", dest="weights",default="data/weights.db", help="SQLite database with rsid weights.")
    parser.add_argument('--weights_on_disk', action="store_true", dest="weights_on_disk", help="Don't load weights db to memory.")
    parser.add_argument('--pheno', action="store", dest="pheno", default=None, help="Phenotype file")
    parser.add_argument('--mpheno', action="store", dest="mpheno", default=None, help="Specify which phenotype column if > 1")
    parser.add_argument('--pheno_name', action="store", dest="pheno_name", default=None, help="Column name of the phenotype to perform association on.")
    parser.add_argument('--missing-phenotype', action="store", dest="missing_phenotype",  default='NA', help="Specify code for missing phenotype information.  Default is NA")
    parser.add_argument('--filter', nargs=2, action="store", dest="fil", default=None, help="Takes two arguments. First is the name of the filter file, the second is a value to filter on.")
    parser.add_argument('--mfilter', action="store", dest="mfil", default=None, help="Column number of filter file to filter on.  '1' specifies the first filter column")
    parser.add_argument('--output_dir', action="store", dest="output", default=None, help="This option is deprecated. Use --output_prefix instead.")
    parser.add_argument('--output_prefix', action="store", dest="output_prefix", default=None, help="Optional prefix for output files. Will concatenate specified prefix with underscore.")
    parser.add_argument('--pred_exp', action="store", dest="pred_exp", default=None, help="Predicted expression file from earlier run of PrediXcan")
    parser.add_argument('--logistic', action="store_true", dest="logistic", default=False, help="Include to perform a logistic regression")
    parser.add_argument('--linear', action="store_true", dest="linear", default=False, help="Include to perform a linear regression")
    parser.add_argument('--survival', action="store_true", dest="survival", default=False, help="Include to perform survival analysis")

    args = parser.parse_args()

    PREDICT = args.predict
    GENE_LIST = args.genelist
    DOSAGE_DIR = args.dosages
    DOSAGE_PREFIX = args.dosages_prefix
    DOSAGE_BUFFER = int(args.dosages_buffer) if args.dosages_buffer else None
    SAMPLE_FILE = os.path.join(DOSAGE_DIR, args.sample_file)
    BETA_FILE = args.weights
    PRELOAD_WEIGHTS = not args.weights_on_disk
    ASSOC = args.assoc
    PHENO_FILE = args.pheno
    MPHENO = str(int(args.mpheno) + 2) if args.mpheno else 'None'
    PHENO_NAME = args.pheno_name if args.pheno_name else 'None'
    MISSING_PHENOTYPE = args.missing_phenotype
    FILTER_FILE, FILTER_VAL = args.fil if args.fil else ('None', '1')
    MFILTER = args.mfil if args.mfil else 'None'
    OUT_EXP_NAME = args.output_prefix + "_predicted_expression.txt"  if args.output_prefix else "predicted_expression.txt"
    PRED_EXP_FILE = args.pred_exp if args.pred_exp else OUT_EXP_NAME
    ASSOC_FILE= args.output_prefix + "_association.txt" if args.output_prefix else "association.txt"
    if args.logistic:
        TEST_TYPE = "logistic"
    elif args.survival:
        TEST_TYPE = "survival"
    else:
        TEST_TYPE = "linear"

    if not PREDICT and not ASSOC:
        print("Error: User did not specify --predict or --assoc. Please specify one or both options.")
        sys.exit(1)

    if not args.output is None:
        print("Error: --output_dir deprecated. Use --output_prefix instead.")
        sys.exit(1)

    if PREDICT:
        check_out_file(PRED_EXP_FILE)
        get_applications_of = GetApplicationsOf(BETA_FILE, PRELOAD_WEIGHTS)
        transcription_matrix = TranscriptionMatrix(BETA_FILE, SAMPLE_FILE, GENE_LIST)
        for rsid, allele, dosage_row in get_all_dosages(DOSAGE_DIR, DOSAGE_PREFIX, DOSAGE_BUFFER):
            for gene, weight, ref_allele in get_applications_of(rsid):
                transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
        transcription_matrix.save(PRED_EXP_FILE)
    if ASSOC:
        check_out_file(ASSOC_FILE)
        subprocess.call(
            ["./PrediXcanAssociation.R",
            "PRED_EXP_FILE", PRED_EXP_FILE,
            "PHENO_FILE", PHENO_FILE,
            "PHENO_COLUMN", MPHENO,
            "PHENO_NAME", PHENO_NAME,
            "MISSING_PHENOTYPE", MISSING_PHENOTYPE,
            "FILTER_FILE", FILTER_FILE,
            "FILTER_VAL", FILTER_VAL,
            "FILTER_COLUMN", MFILTER,
            "TEST_TYPE", TEST_TYPE,
            "OUT", ASSOC_FILE])


if __name__ == '__main__':
    main()
