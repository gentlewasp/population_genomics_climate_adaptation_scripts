#! /usr/bin/perl -w
use strict;

system('mkdir -p gvcf');
system('mkdir -p callsnp_log');
#need change
my @list1s=`find /public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling/0_mapping/tmp_pipe_data -name "*.sorted.rmdup.bam"`;

my(@lists);
push @lists,@list1s;
print "$#lists\n";

my($core,$genome,$sample,$IN,$OUT);
$OUT='/public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling';
$genome="/public1/home/caolijun/yuelei/malijun/2_genome/genome.fasta";
$core =6;
my $tmp_pipe_data="/public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling/0_mapping/tmp_pipe_data";
my $q="i01203share";
&createJob(@lists);
&submitJob(@lists);
sub submitJob{
        my @arr=@_;
        foreach (@arr){
           chomp;
           my @arrays = split/\//;
           my $sample = $arrays[-1];
           print $sample,"\n";
           $sample =~s/.sorted.rmdup.bam//; 
           my $Jobs=`squeue | grep "callsnp" | wc -l`;
           while( $Jobs > 98 ){
              print "JOB Remain $Jobs\n";
              sleep 98;
              $Jobs=`squeue | grep "callsnp" | wc -l`;
           }
           system("sbatch $OUT/gvcf/$sample.pbs >$OUT/callsnp_log/$sample.log");
        }
}
sub createJob{
        my @arr=@_;
        foreach (@arr){
                chomp;
                my $sra = $_;
                my @arrays = split/\//,$sra;
                my $sample = $arrays[-1];
                $sample =~s/.sorted.rmdup.bam//; #need change
                open OUT,">$OUT/gvcf/$sample.pbs";
                print OUT <<SET;
#!/bin/bash
#SBATCH -J TPcallsnp_$sample
#SBATCH -p $q
#SBATCH -n $core
#SBATCH -N 1
#SBATCH --mem=18g
#SBATCH -o out$sample.%j
#SBATCH -e err$sample.%j
date
#callsnp
/public1/home/caolijun/yuelei/miniconda3/bin/gatk --java-options "-Xmx18g -Djava.io.tmdir=./tmp" HaplotypeCaller -R $genome -I $tmp_pipe_data/$sample.sorted.rmdup.bam -ERC GVCF -O $OUT/gvcf/$sample.g.vcf 1>$sample.callsnp.log 2>&1

date
SET
        close OUT;
        }
}
