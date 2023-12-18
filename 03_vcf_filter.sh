rawvcf=mlj_19pop_277sample_hardfilter_snp.recode.vcf1
outprefix=19pops
##### 1.Filter on minimum read depth (DP) and genotype quality (GQ)
bcftools filter $rawvcf -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o 19pops_DP3.GQ20.vcf.gz 
## check
bcftools query -i 'FMT/DP<3' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' 19pops_DP3.GQ20.vcf.gz | less -S
bcftools query -i 'FMT/DP=3' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' 19pops_DP3.GQ20.vcf.gz | less -S
##### 2.Remove multiallelic SNPs and indels, monomorphic SNPs, and SNPs in the close proximity of indels
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 19pops_DP3.GQ20.vcf.gz | \
bcftools view -m2 -M2 -v snps -O z -o 19pops_DP3.GQ20.allele.vcf.gz 
## check
bcftools view -H 19pops_DP3.GQ20.vcf.gz | less -S
bcftools view -H 19pops_DP3.GQ20.allele.vcf.gz | less -S
zless -S 19pops_DP3.GQ20.vcf.gz
zless -S 19pops_DP3.GQ20.allele.vcf.gz
## How many variants are left in our data file after all the above filtering steps have been applied?
bcftools view -H 19pops_DP3.GQ20.allele.vcf.gz | wc -l
##### 3. Remove individuals with a high amount of missing data
bcftools stats -s - 19pops_DP3.GQ20.allele.vcf.gz | grep -E ^PSC | cut -f3,14 > 19pops_DP3.GQ20.allele.imiss
###plot
Rscript missing_persample.R
# remove individuals
bcftools view 19pops_DP3.GQ20.allele.vcf.gz -s ^SDS313,NMHS15,SDS306,GDZQ06,GDZQ04,GDZQ03,SDS305,YNBN03,YNBN01,SDS105,SDS304,SDS301 \
-O z -o 19pops_DP3.GQ20.allele.missi.vcf.gz 
# Count the number of samples
bcftools query -l 19pops_DP3.GQ20.allele.missi.vcf.gz | wc -l
##### 4. Remove variants with a high amount of missing genotypes and filter on minor allele frequency
bcftools filter 19pops_DP3.GQ20.allele.missi.vcf.gz -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o 19pops_DP3.GQ20.allele.missi.miss20.maf0.05.vcf.gz 
# check
bcftools query 19pops_DP3.GQ20.allele.missi.miss20.maf0.05.vcf.gz  -f'%AC\n' | sort -g | head
bcftools view -H 19pops_DP3.GQ20.allele.missi.miss20.maf0.05.vcf.gz | wc -l

# check depth of each sample
bcftools stats -s - 19pops_DP3.GQ20.allele.missi.miss20.maf0.05.vcf.gz | grep -E ^PSC | cut -f3,14 > 19pops_DP3.GQ20.allele.missi.miss20.maf0.05.imiss
Rscript missing_persample.R