#!/bin/bash
#SBATCH -J knn
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o %j.out
#SBATCH -e %j.err
samples=tp_chr15_within
chr=15
vcf=14pops_rm_allhet_chr15.recode.vcf.gz
prefix=tp


bcftools view -S $samples -r $chr $vcf | vcftools --vcf -  --maf 0.05 --thin 100 -c --geno-r2 | awk '!/-na/' | perl /work/share/kdy_caolijun/software/reformat/emerald2windowldcounts.pl | gzip >${prefix}.thin100.maf5.Chr${chr}.${samples}.windows.ld.gz
