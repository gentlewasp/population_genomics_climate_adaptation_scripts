#!/bin/bash
#SBATCH -J mlj_populations
#SBATCH -p a01216amd
#SBATCH -N 1 
#SBATCH -n 32
#SBATCH -o %j.out
#SBATCH -e %j.err

###############
###############
##populations##
###############
###############


#01_BJDX_15sample
populations -V BJDX_15sample_hardfilter_snp_mis_1.recode.vcf -M BJDX_pop.list -t 32 -p 1 -O ./ --fstatss

#02_BJCY_15sample
populations -V BJCY_15sample_hardfilter_snp_mis_1.recode.vcf -M BJCY_pop.list -t 32 -p 1 -O ./ --fstats

#03_GDSZ_15sample
populations -V GDSZ_15sample_hardfilter_snp_mis_1.recode.vcf -M GDSZ_pop.list -t 32 -p 1 -O ./ --fstats

#04_GDZQ_15sample
populations -V GDZQ_15sample_hardfilter_snp_mis_1.recode.vcf -M GDZQ_pop.list -t 32 -p 1 -O ./ --fstats

#05_HNCS_15sample
populations -V HNCS_15sample_hardfilter_snp_mis_1.recode.vcf -M HNCS_pop.list -t 32 -p 1 -O ./ --fstats

#06_HNS1_15sample
populations -V HNS1_15sample_hardfilter_snp_mis_1.recode.vcf -M HNS1_pop.list -t 32 -p 1 -O ./ --fstats

#07_HNSY_15sample
populations -V HNSY_15sample_hardfilter_snp_mis_1.recode.vcf -M HNSY_pop.list -t 32 -p 1 -O ./ --fstats

#08_JANP_15sample
populations -V JANP_15sample_hardfilter_snp_mis_1.recode.vcf -M JANP_pop.list -t 32 -p 1 -O ./ --fstats

#09_LNAS_15sample
populations -V LNAS_15sample_hardfilter_snp_mis_1.recode.vcf -M LNAS_pop.list -t 32 -p 1 -O ./ --fstats

#10_NMHS_15sample
populations -V NMHS_15sample_hardfilter_snp_mis_1.recode.vcf -M NMHS_pop.list -t 32 -p 1 -O ./ --fstats

#11_SCCD_15sample
populations -V SCCD_15sample_hardfilter_snp_mis_1.recode.vcf -M SCCD_pop.list -t 32 -p 1 -O ./ --fstats

#12_SCDY_15sample
populations -V SCDY_15sample_hardfilter_snp_mis_1.recode.vcf -M SCDY_pop.list -t 32 -p 1 -O ./ --fstats

#13_SDS1_14sample
populations -V SDS1_14sample_hardfilter_snp_mis_1.recode.vcf -M SDS1_pop.list -t 32 -p 1 -O ./ --fstats

#14_SDS3_15sample
populations -V SDS3_15sample_hardfilter_snp_mis_1.recode.vcf -M SDS3_pop.list -t 32 -p 1 -O ./ --fstats

#15_SDSG_15sample
populations -V SDSG_15sample_hardfilter_snp_mis_1.recode.vcf -M SDSG_pop.list -t 32 -p 1 -O ./ --fstats

#16_SJY1_15sample
populations -V SJY1_15sample_hardfilter_snp_mis_1.recode.vcf -M SJY1_pop.list -t 32 -p 1 -O ./ --fstats

#17_SJY4_11sample
populations -V SJY4_15sample_hardfilter_snp_mis_1.recode.vcf -M SJY4_pop.list -t 32 -p 1 -O ./ --fstats

#18_YNBN_15sample
populations -V YNBN_15sample_hardfilter_snp_mis_1.recode.vcf -M YNBN_pop.list -t 32 -p 1 -O ./ --fstats

#19_YNXS_15sample
populations -V YNXS_15sample_hardfilter_snp_mis_1.recode.vcf -M YNXS_pop.list -t 32 -p 1 -O ./ --fstats
