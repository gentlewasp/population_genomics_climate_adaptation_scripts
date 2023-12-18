#!/bin/bash 
#SBATCH -J smcYNXS 
#SBATCH -p xhacnormalb
#SBATCH -N 1  
#SBATCH -n 15 
#SBATCH -o %j.out 
#SBATCH -e %j.err 

vcf=/work/home/malijun/malijun/01_data/1vcffiles/INPUT.vcf.gz
SMC=$(which smc++)
pop=YNXS
mkdir ${pop}_out
#bcftools sort $vcf -O z -o ${vcf}.gz
#bcftools index ${vcf}.gz
#for i in {YNXS01,YNXS02,YNXS03,YNXS04,YNXS05,YNXS06,YNXS07,YNXS08,YNXS09,YNXS10,YNXS11,YNXS12,YNXS13,YNXS14,YNXS15}
#  do
#  for j in {chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12,chr13,chr14,chr15,chr16}
#    do $SMC vcf2smc --cores 50 -m genome_mask_output.gz -d $i $i ${vcf}.gz ${pop}_out/${pop}_${j}_${i}.smc.gz $j YNXS:YNXS01,YNXS02,YNXS03,YNXS04,YNXS05,YNXS06,YNXS07,YNXS08,YNXS09,YNXS10,YNXS11,YNXS12,YNXS13,YNXS14,YNXS15
#  done
#done
smc++ estimate --thinning 2000 --em-iterations 50 --base ${pop}_1000_8.4 -o ${pop}_out/ 8.4e-9 ${pop}_out/${pop}_*_*.smc.gz --cores 50 --spline cubic --timepoints 1 1.5e6

#smc++ plot -g 0.25 fww_model3_plot_years.pdf BJPG_model3.final.json YNHH_model3.final.json SCCD_model3.final.json  -x 1 1e6
