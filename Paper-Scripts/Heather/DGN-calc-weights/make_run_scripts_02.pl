use warnings;
use strict;

my $prescript = "02_imputedDGN-WB_CV_polyscore";
open(QSUB, ">qsub.txt");

for(my $i = 1; $i<=22; $i++){
    my @pvals = (0.0001,0.001,0.01,0.05,0.5,1);
    for(my $j = 0; $j <= 5; $j++){
	my $k = $pvals[$j];
	my $outfile = "runscripts/run_" . $prescript . "_chr" . $i . "_Pthresh" . $k . ".sh";
	print QSUB "qsub $outfile\nsleep 4\n";
	open(OUT, ">$outfile");
	print OUT "#!/bin/bash\n#PBS -N R.pscore\n#PBS -S /bin/bash\n#PBS -l walltime=72:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=10gb\n#PBS -o logs/\$\{PBS_JOBNAME\}.o\$\{PBS_JOBID\}.log\n#PBS -e logs/\$\{PBS_JOBNAME\}.e\$\{PBS_JOBID\}.err\n\ncd \$PBS_O_WORKDIR\nmodule load R\n\ntime R --no-save < $prescript.r --args $i $k\nchmod 750 *\n";
    }
}
