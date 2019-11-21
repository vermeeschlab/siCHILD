#!/bin/bash

# run by: ~/git_gcap/sichild/comb_GBS_gvcf.sh GBS_PGD440_C1

wd='/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/'

a='gbs_run_'
sh="$a$1"


o1='GBS_PGD_output/'
o2="$wd$o1$1"
o="$o2/merge/"
mkdir -p "$o"

#for j in {1..22}; do echo -n "/uz/data/hydra/shared_app/apps/jdk/jdk1.8.0_131/bin/java -Xmx16G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=tmp -jar /uz/data/hydra/shared_app/apps/gatk/3.8/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /uz/data/shortcuts/genomicscore/folders/raw/bio/resources/species/homo_sapiens/hg19/picard/2.17.10/genome.fa --out $o/$1.$j.merge.vcf " >> $wd$sh$j;  for i in $(ls /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/$1*/gvcf_variants_gatk/current/result/*.$j.HaplotypeCaller.raw.g.vcf_filter.vcf ); do echo -n "-V $i " >> $wd$sh$j ; done; done;
#for j in {"X","Y"}; do echo -n "/uz/data/hydra/shared_app/apps/jdk/jdk1.8.0_131/bin/java -Xmx16G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=tmp -jar /uz/data/hydra/shared_app/apps/gatk/3.8/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /uz/data/shortcuts/genomicscore/folders/raw/bio/resources/species/homo_sapiens/hg19/picard/2.17.10/genome.fa --out $o/$1.$j.merge.vcf " >> $wd$sh$j;  for i in $(ls /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/$1*/gvcf_variants_gatk/current/result/*.$j.HaplotypeCaller.raw.g.vcf_filter.vcf ); do echo -n "-V $i " >> $wd$sh$j ; done; done;
#chmod 755 $wd$sh*

echo $wd$sh$j
#for j in {1..22}; do echo -n "/uz/data/hydra/shared_app/apps/freebayes/freebayes/bin/freebayes -f /uz/data/shortcuts/genomicscore/folders/raw/bio/resources/species/homo_sapiens/hg19/picard/2.17.10/genome.fa --ploidy 2 -r $j" >> $wd$sh$j;  for i in $(ls /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/$1*/align_bwa/current/result/*REF_$j.bam); do echo -n " $i " >> $wd$sh$j; done; echo " > $o/$1.$j.merge.freebayers.vcf" >> $wd$sh$j; done;
chmod 755 $wd$sh*
for j in {"X","Y"}; do echo -n "/uz/data/hydra/shared_app/apps/freebayes/freebayes/bin/freebayes -f /uz/data/shortcuts/genomicscore/folders/raw/bio/resources/species/homo_sapiens/hg19/picard/2.17.10/genome.fa --ploidy 2 -r $j" >> $wd$sh$j;  for i in $(ls /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/$1*/align_bwa/current/result/*REF_$j.bam); do echo -n " $i " >> $wd$sh$j; done; echo " > $o/$1.$j.merge.freebayers.vcf" >> $wd$sh$j; done;
chmod 755 $wd$sh*


# ls -1  /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/gbs_run_GBS_PGD440_C1* > /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/gbs_run_GBS_PGD440_C1.sh

