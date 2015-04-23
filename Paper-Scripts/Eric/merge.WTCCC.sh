# Eric Gamazon
PERL5LIB=/nas40t0/egamazon/VANDY/VCF/vcftools_0.1.12b/perl/; export PERL5LIB
# RA
      for i in `seq 22 22`;
      do
             tabix AffyV4.Plinkchecked.58C+NBS.excluded.samples.snps.$i.vcf.gz
             tabix AffyV4.Plinkchecked.RA.excluded.samples.snps.$i.vcf.gz
             /nas40t0/egamazon/VANDY/VCF/vcftools_0.1.12b/bin/vcf-merge AffyV4.Plinkchecked.58C+NBS.excluded.samples.snps.$i.vcf.gz AffyV4.Plinkchecked.RA.excluded.samples.snps.$i.vcf.gz | bgzip -c > OUT/Combined.RAControls.$i.vcf.gz
             echo "/nas40t0/egamazon/VANDY/VCF/vcftools_0.1.12b/bin/vcf-merge AffyV4.Plinkchecked.58C+NBS.excluded.samples.snps.$i.vcf.gz AffyV4.Plinkchecked.RA.excluded.samples.snps.$i.vcf.gz | bgzip -c > OUT/Combined.RAControls.$i.vcf.gz"
      done

