#!/bin/bash
#SBATCH -J mlj
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -o %j.out
#SBATCH -e %j.err


pixy --stats dxy\
--vcf input.vcf.gz \
--populations pop \
--window_size 2000 \
--n_cores 30 \
--bypass_invariant_check 'no'
~                                                                                                                                                                                         
~ 