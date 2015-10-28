#!/usr/bin/env python

import sys
import os
import subprocess
import gzip
from argparse import ArgumentParser

if __name__ == "__main__":
  parser = ArgumentParser(
    description="Converts plink file format to dosage format for prediXcan"
  )
  parser.add_argument(
    "-b", "--bfile", dest="bfile", required=True,
    help="prefix for binary ped, bim, and bam files."
  )
  parser.add_argument(
    "-o", "--out", dest="out", required=False,
    help="prefix for output files", default="chr"
  )

  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

  args = parser.parse_args()

  # First we get the minor allele dosages for *all* chromosomes:
  subprocess.call([
    'plink2', '--bfile', args.bfile,
    '--recode', 'A-transpose', '--out', args.out
  ])

  # Now calculate the minor allele frequency:
  subprocess.call([
    'plink2', '--bfile', args.bfile, '--freq', '--out', args.out
  ])

  # now we need to process the files line by line to create the dosage files
  with open(args.out + ".traw") as dfile, open(args.out + ".frq") as ffile, open(args.bfile + ".bim") as bfile:
    i = 1
    curCHR = -1
    while True:
      # skip the headers
      if i == 1:
        dfile.readline()
        ffile.readline()
        i += 1
      else:
        # First check that we're not at the end of the files
        dline = dfile.readline()
        fline = ffile.readline()
        bline = bfile.readline()
        if not dline: break
        if not fline: break
        if not bline: break
        # Split into columns by whitespace
        dcols = dline.split()
        fcols = fline.split()
        bcols = bline.split()
        # Make sure rsIDs match in maf frequency file -- note this cannot
        # be filtered unlike the --recode 'd file
        while True:
          if fcols[1] != dcols[1]:
            fcols = ffile.readline().split()
          else:
            break
        # Combine columns as per 'dosage' format. Impute missing data as 2*MAF.
        nline = fcols[0:2] + [bcols[3]]+ fcols[2:5] + [fcols[4]*2 if e == "NA" else e for e in dcols[6:]]
        # Write out to the appropriate file for that chromosome
        if fcols[0] != curCHR:
          if curCHR != -1:
            ofile.close()
          ofile = open(args.out + fcols[0] + ".txt", "a")
        ofile.write(" ".join(nline) + "\n")
  ofile.close()

  # now we need to gzip the files
  out_dir = os.path.dirname(args.out)
  prefix = os.path.basename(args.out)
  for chrfile in [x for x in os.listdir(out_dir) if (x.startswith(prefix) and x.endswith(".txt"))]:
    f = open(os.path.join(out_dir, chrfile), 'rb')
    ofile = gzip.open(os.path.join(out_dir, chrfile + ".gz"), "wb")
    ofile.writelines(f)
    f.close()
    ofile.close()
    # remove non-zipped file
    os.remove(os.path.join(out_dir, chrfile))

  # Remove left over files
  os.remove(os.path.join(out_dir, prefix + ".frq"))
  os.remove(os.path.join(out_dir, prefix + ".traw"))
  os.remove(os.path.join(out_dir, prefix + ".log"))

