# Eric Gamazon
# eQTL mapping per chromosome on DGN expression data

# cat /nas40t0/egamazon/VANDY/PREDIXCAN/snp_location_info.MOD | perl -ane '$chr = $F[1]; $chr =~ s/chr//; print $F[0] . "\t" . $F[1] . "\t" . $F[2] . "\n" if ($chr >= 1 and $chr <= 22) ' > /nas40t0/egamazon/VANDY/PREDIXCAN/snp_location_info.MOD.v1 &

R --no-save --args  22 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r

nohup R --no-save --args  21 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r &

nohup   R --no-save --args  20 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup20.out &
nohup   R --no-save --args  19 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup19.out &
nohup   R --no-save --args  18 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup18.out &


nohup   R --no-save --args  17 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup17.out &
nohup   R --no-save --args  16 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup16.out &
nohup   R --no-save --args  15 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup15.out &
nohup   R --no-save --args  14 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup14.out &


nohup   R --no-save --args  13 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup13.out &
nohup   R --no-save --args  12 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup12.out &
nohup   R --no-save --args  11 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup11.out &
nohup   R --no-save --args  10 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r    > nohup10.out &


nohup   R --no-save --args  9 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup9.out &
nohup   R --no-save --args  8 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup8.out &
nohup   R --no-save --args  7 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup7.out &
nohup   R --no-save --args  6 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup6.out &


nohup   R --no-save --args  5 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup5.out &
nohup   R --no-save --args  4 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup4.out &
nohup   R --no-save --args  3 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup3.out &
nohup   R --no-save --args  2 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup2.out &
nohup   R --no-save --args  1 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.r     > nohup1.out &

#nohup   R --no-save --args  1 < /nas40t0/egamazon/VANDY/PREDIXCAN/mapeQTLs.TMP.r         > nohup1.tmp.out &


