Frequently Asked Questions

- What does TW mean in the model names?
  - We use TW to refer to the models trained with the original tissue
  expression. This makes sense in relation to the tissue specific
  expression models (TS*) we have derived using Orthogonal Tissue
  Decomposition [See Details Here](http://biorxiv.org/content/early/2016/03/15/043653.1)

- Are performance measures of prediction available?
  - Yes, they can be found in the `extra` table in each sqlite db.

- How many samples were used to train each model?
  - This can be found in the `sample_info` table in each sqlite db.

- How do I query the prediction model db?
  - Find example queries in
  [here](https://github.com/hakyimlab/PrediXcan/blob/master/Software/query-db.Rmd)
  
- PrediXcan says my samples file has too few/too many rows.
  - In order for the predicted expression file to be correct, the dosage
  files and samples file must correspond to the same individual.  The
  columns for the dosage files are columns are snpid rsid position
  allele1 allele2 MAF id1 ..... idn and the samples file lists out the
  information for id1 ... idn row by row.  As the id numbers are not
  included in the dosage file, it is critically important that the
  samples file has the same number of individuals in the same order as
  the dosages files.  Otherwise later association tests will be invalid.
