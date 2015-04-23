use warnings;
use strict;

open(A, "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB.exp.ID.list");

my %ids;
while(<A>){
    chomp;
    my ($id) = split(/\n/);
    $ids{$id} = 1;
}

open(B, "DGN.hapmap2.fam");
open(OUT, ">DGN-WB.exp.FID.IID");

while(<B>){
    chomp;
    my ($fid, $iid, $d, $m, $s, $p) = split(/\s+/);
    my ($a, $b, $c) = split(/_/, $fid);
    if(defined($ids{$c})){
	print OUT "$fid $iid\n";
    }
}
