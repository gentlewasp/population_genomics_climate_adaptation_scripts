#!/bin/bash
#SBATCH -J mlj_03rda
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -o %j.out
#SBATCH -e %j.err


#Rscript 01_PCs_4fds.R
Rscript 02rda_tp_16pops.R 
