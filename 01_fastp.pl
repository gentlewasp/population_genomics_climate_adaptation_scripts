#! /usr/bin/perl -w
use strict;
use Getopt::Long;
my($input_dir,$outdir,$suffix,$first_position,$last_position,$suffix_out);
GetOptions(
        "input_dir:s"           =>\$input_dir,
        "outdir:s"              =>\$outdir,
        "suffix:s"                  =>\$suffix,
	"first_position:s"	=>\$first_position,
	"last_position:s"	=>\$last_position,
	"suffix_out:s"		=>\$suffix_out
        
)||&help;
&help unless ($input_dir && $outdir && $suffix && $first_position && $last_position && $suffix_out);
sub help{
print<<"Usage.end";

        Description:
        Function: Pick up mapping ratio based on samtools flagstet file.
        Usage:
                -input_dir      the mapinfo dir         must be given
                -outdir         output dir              must be given
                -suffix         the suffix of the name     muse be given
		-first_position  	the first posiion	must be given
		-last_position		the last position	must be given
                -suffix_out		the suffix of the outfiles must be given
Usage.end
        exit;
}

system('mkdir -p pbs_folder');
system('mkdir -p fastp_log');
#need change
my @list1s=`find $input_dir -name "*$suffix"`;

my(@lists);
push @lists,@list1s;
#push @lists,@list2s;
print "$#lists\n";

my($core,$threads,$left_raw,$right_raw,$left_clean,$right_clean,$IN,$OUT);
$OUT=$outdir;
$core =12;
$threads =24;
my $q="xhacnormalb";
&createJob(@lists);
&submitJob(@lists);
sub submitJob{
        my @arr=@_;
        foreach (@arr){
           chomp;
           my @arrays = split/\//;
           my $sample = $arrays[-1];
           print $sample,"\n";
           $sample =~s/$suffix//; 
           system("sbatch $OUT/pbs_folder/$sample".".pbs >$OUT/fastp_log/$sample".".fastp.log");
        }
}

sub createJob{
        my @arr=@_;
        foreach (@arr){
                chomp;
                my $sra = $_;
                my @arrays = split/\//,$sra1;
                my $sample = $arrays[-1];
                $sample =~s/$suffix//;
                my $left_raw = $sra;
                my $right_raw = $sra;
                $left_raw =~s/.R1/.R1/;
                $right_raw =~s/.R1/.R2/;
		my $sample_out= substr($sample,$first_position - 1,$last_position - $first_position+1);
		my $left_paired = $sample_out .".R1".".".$suffix_out;
                my $right_paired = $sample_out .".R2".".".$suffix_out;
                open OUT,">$OUT/pbs_folder/$sample".".pbs";
                print OUT <<SET;
#!/bin/bash
#SBATCH -J trim_$sample
#SBATCH -p $q
#SBATCH -n $core
#SBATCH --mem=36g
#SBATCH -o out$sample.%j
#SBATCH -e err$sample.%j
date
#fastp
fastp -i $left_raw -I $right_raw -o $left_paired -O $right_paired
mv $left_paired $outdir
mv $right_paired $outdir
date
SET
        close OUT;
        }
}
