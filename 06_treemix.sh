#!/bin/bash
#SBATCH -J mlj_3mig
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -o %j.out
#SBATCH -e %j.err
module load compiler/gcc/10.2.0

treemix=/work/share/kdy_caolijun/software/treemix-1.13/bin/treemix
# Running Treemix
for j in {1..100} 
    do
   $treemix -i ../16.treemix.txt.gz -bootstrap -k 100 -m 3 -root YNBN -o thrips_tree_3mig_bootstrap_${j} # -m needs change from 1 to 10
done
