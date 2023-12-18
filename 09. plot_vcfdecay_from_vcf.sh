#!/bin/bash
#SBATCH -J MLJ_LDdecay
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o %j.out
#SBATCH -e %j.err



popName=popList
prefix=tp

for pop in $(cat $popName)
do
/public/share/kdy_caolijun/shared/software/PopLDdecay/bin/PopLDdecay -InVCF ${prefix}_${pop}.recode.vcf -MAF 0.05  -OutStat ${prefix}_${pop}.LDdecay
done
ls *.LDdecay.stat.gz > ld.results.ls
paste ld.results.ls popList > Pop.ReslutPath.list
perl /public/share/kdy_caolijun/shared/software/PopLDdecay/bin/Plot_MultiPop.pl  -inList  Pop.ReslutPath.list -output  tp.19pop.LDdecay

