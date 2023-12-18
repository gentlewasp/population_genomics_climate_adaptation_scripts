#!/bin/bash 
#SBATCH -J mljadmixture 
#SBATCH -p a01216amd
#SBATCH -N 1  
#SBATCH -n 10 
#SBATCH -o %j.out 
#SBATCH -e %j.err 

bed_prefix=tp.rad
popmap=tp.popmap
popOrder=tp.poporder
pdf_prefix=tp
for K in  $(seq 10); \
do admixture --cv=10 $bed_prefix.bed $K | tee log${K}.out; done

grep -h CV log*.out > cv.results

for i in  $(seq 10); \
do python distruct2.2_ocwa.py -K $i --input=$bed_prefix --output=${pdf_prefix}${i}.pdf --title=K=$i --popfile=$popmap --poporder=$popOrder; done
