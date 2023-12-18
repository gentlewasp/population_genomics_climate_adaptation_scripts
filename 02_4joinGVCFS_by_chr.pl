#! /usr/bin/perl -w
use strict;

system('mkdir -p joinGVCF_log');
system('mkdir -p vcf_by_chr');
##system('mkdir -p sam');
##system('mkdir -p un_conc');
open(DATA, "</public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling/tp_chr_list");
my @list1s = <DATA>;
close DATA;
#my @list1s=`find /public1/home/caolijun/caolijun/dbm_chenmingzhu -name "*.rem.1.fq.gz"`;

my(@lists);
push @lists,@list1s;
#push @lists,@list2s;
print "@lists\n";

my($core,$IN,$OUT,$GVCFS);
$OUT='/public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling';   #need change
$GVCFS='/public1/home/caolijun/yuelei/malijun/tp_rawdata2/3_snp_calling/gvcf.list';
my $REF='/public1/home/caolijun/yuelei/malijun/2_genome/genome.fasta';  #need change
$core =28;
my $q="i01203share";
my $mem = $core * 3;
&createJob(@lists);
&submitJob(@lists); ##if sub multipbs,and need remove #

sub submitJob{
	my @arr=@_;
	foreach (@arr){
	  chomp;
          ##my @arrays = split/\//;
           my $chr = $_;
##           print $arr,"\n";
  ##         $chr = $_; #need change
           #$chr =~s/_1.clean.fq.gz//; #need change
          # my $Jobs=`squeue | grep "bowtie" | wc -l`;
           #while( $Jobs > 60 ){
             ##print "JOB Remain $Jobs\n";
             ## sleep 60;
              ##$Jobs=`squeue | grep "bowtie" | wc -l`;
          ## }
	  system("sbatch joinGVCF_log/$chr".".pbs >$OUT/joinGVCF_log/$chr".".log");
          system('sleep 5');
	}
}


sub createJob{
	my @arr=@_;
	foreach (@arr){
		chomp;
                my $chr = $_;
                my @job = split /\_/,$chr;
                #my $jobname = $job[-2];
               ## open (VCF, "<vcf_3.list");
               ## my @VCFS = <VCF>;
                ##close VCFS;
                ##my @arrays = split/\//,$sra;
               ###my $chr = $arrays[-1];
		###$chr =~s/\.rem.1.fq.gz//; #need change
	        #$chr =~s/_1.clean.fq.gz//; #need change
		###my $left = $sra;
		##my $right = $sra;
               ## my $left_raw = $sra;
		###my $right_raw = $sra;	
		##$left_raw =~s/\.rem.1.fq.gz/\.1.fq.gz/;#need change
		###$right =~s/\.rem.1.fq.gz/\.rem.2.fq.gz/;#need change
                ##$right_raw =~s/\.rem.1.fq.gz/\.2.fq.gz/;#need change
                open OUT,">$OUT/joinGVCF_log/$chr".".pbs";
		print OUT <<SET;
#!/bin/bash
#SBATCH -J TPjoin_$chr
#SBATCH -p $q
#SBATCH -n $core
#SBATCH -N 1
#SBATCH --mem=$mem\g
#SBATCH -o out$chr.%j
#SBATCH -e err$chr.%j

export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.131-11.b12.el7.x86_64/
export CLASSPATH=.:\$JAVA_HOME/lib:\$JAVA_HOME/jre/lib
export PATH=\$JAVA_HOME/bin:\$JAVA_HOME/jre/bin:\$PATH

date
#GenomicsDBImport
gatk --java-options "-Xmx$mem\g" GenomicsDBImport --sample-name-map $GVCFS --genomicsdb-workspace-path joinGVCF_log/my_database_$chr --intervals $chr 

##joint genotyping
gatk GenotypeGVCFs -R $REF -V gendb://joinGVCF_log/my_database_$chr -O vcf_by_chr/final_$chr.vcf

date
SET
	close OUT;
	}
}

