#!/bin/bash
#SBATCH -J MLJsplitPop
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o %j.out
#SBATCH -e %j.err


popName=popList
vcfFile=/public/home/malijun/malijun/01_data/1vcffiles/input.recode.vcf
prefix=tp
for indList in $(cat ${popName})
do vcftools --vcf $vcfFile --keep $indList --out ${prefix}_$indList --recode
vcftools --vcf ${prefix}_${indList}.recode.vcf --plink --out ${prefix}_${indList}
awk -F" " '{print $2,$4}' ${prefix}_${indList}.map > ${prefix}_${indList}.info
done

