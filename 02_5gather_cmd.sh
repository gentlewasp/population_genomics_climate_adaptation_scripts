#!/bin/bash
#SBATCH -J gather_vcf
#SBATCH -p i01203share
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=84g
#SBATCH -o %j.out
#SBATCH -e %j.err

###注意染色体的编号需要按照基因组的排序顺序一样
gatk GatherVcfs  --java-options "-Xmx84g -Djava.io.tmpdir=./tmp" \
-I vcf_by_chr/final_chr01.vcf \
-I vcf_by_chr/final_chr02.vcf \
-I vcf_by_chr/final_chr03.vcf \
-I vcf_by_chr/final_chr04.vcf \
-I vcf_by_chr/final_chr05.vcf \
-I vcf_by_chr/final_chr06.vcf \
-I vcf_by_chr/final_chr07.vcf \
-I vcf_by_chr/final_chr08.vcf \
-I vcf_by_chr/final_chr09.vcf \
-I vcf_by_chr/final_chr10.vcf \
-I vcf_by_chr/final_chr11.vcf \
-I vcf_by_chr/final_chr12.vcf \
-I vcf_by_chr/final_chr13.vcf \
-I vcf_by_chr/final_chr14.vcf \
-I vcf_by_chr/final_chr15.vcf \
-I vcf_by_chr/final_chr16.vcf \
-O vcf_by_chr/19pop_285samples_gather.vcf
