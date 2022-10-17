#!/usr/bin/env python

'''make a run script for each chr and output a qsub file'''

qsubfile = open('qsub.txt','w')
prescript = '03_imputedDGN-WB_calc_topSNP_R2'

for i in range(1,23):
    outfilename = 'runscripts/run_' + prescript + '_chr' + str(i) + '.sh'
    outfile = open(outfilename,'w')
    output = '''#!/bin/bash
#PBS -N R.ps.''' + str(i) +'''\n#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R

time R --vanilla < ''' + prescript + '.r --args ' + str(i) + '\n'
    outfile.write(output)
    qsubfile.write('qsub ' + outfilename + '\nsleep 3\n')


