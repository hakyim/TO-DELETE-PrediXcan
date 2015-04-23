use warnings;
use strict;

my $prescript = "01_imputedDGN-WB_CV_elasticNet";
open(QSUB, ">qsub.txt");

for(my $i = 1; $i<=22; $i++){
#    for(my $j = 1; $j <= 2; $j++){ ##alpha=0.5 & 1
    for(my $j = 1; $j <= 1; $j++){ ##alpha = 0.5 only
	my $k = $j/2;
	my $outfile = "runscripts/run_" . $prescript . "_chr" . $i . "_alpha" . $k . ".sh";
	print QSUB "qsub $outfile\nsleep 4\n";
	open(OUT, ">$outfile");
	print OUT "#!/bin/bash\n#PBS -N R.glmnet\n#PBS -S /bin/bash\n#PBS -l walltime=100:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=10gb\n#PBS -o logs/\$\{PBS_JOBNAME\}.o\$\{PBS_JOBID\}.log\n#PBS -e logs/\$\{PBS_JOBNAME\}.e\$\{PBS_JOBID\}.err\n\ncd \$PBS_O_WORKDIR\nmodule load R\n\ntime R --no-save < $prescript.r --args $i $k\n";
    }
}
