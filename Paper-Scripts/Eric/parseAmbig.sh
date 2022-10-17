#!/bin/bash
# Generate ambiguous SNPs for imputation 

echo "$1"

perl -ane 'print "$F[1]\t$F[4]\t$F[5]\n" if($F[4] eq "A" and $F[5] eq "T")' "$1"
perl -ane 'print "$F[1]\t$F[4]\t$F[5]\n" if($F[4] eq "T" and $F[5] eq "A")' "$1"
perl -ane 'print "$F[1]\t$F[4]\t$F[5]\n" if($F[4] eq "C" and $F[5] eq "G")' "$1"
perl -ane 'print "$F[1]\t$F[4]\t$F[5]\n" if($F[4] eq "G" and $F[5] eq "C")' "$1"

exit 0
