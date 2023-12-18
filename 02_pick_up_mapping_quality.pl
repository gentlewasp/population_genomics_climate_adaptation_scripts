#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $ver	=	"1.0";
my $Writer =	"wangzhiwei <wangzhiwei1\@genomics.cn>";
my $Data = 	"2019/7/9";

my($input_dir,$outdir,$PE,$ref_genome_size);
GetOptions(
	"input_dir:s"		=>\$input_dir,
	"outdir:s"		=>\$outdir,
	"pe:s"			=>\$PE,
	"size:s"		=>\$ref_genome_size
)||&help;
&help unless ($input_dir && $outdir && $PE && $ref_genome_size);
sub help{
print<<"Usage.end";

	Description:
	Writer	: $Writer
        Data    : $Data
        Version : $ver
	Function: Pick up mapping ratio based on samtools flagstet file.
	Usage:
		-input_dir	the mapinfo dir		must be given
		-outdir		output dir 		must be given
		-pe		base number of read     muse be given
		-size 		Ref genome size		must be given[bp]
Usage.end
	exit;
}
my @mapinfo_list = `find $input_dir -name "*.source.mapinfo"`;

open (OUT,">$outdir/mapping_ratio.quality");
print OUT "Sample\tMapping read\tMapping ratio\tClean read\tClean base\tDepth\tProperly paired read\tProperly paired ratio\n";
my($samp,$number,$free_number,$type);
for (my $i=0;$i<scalar (@mapinfo_list);$i++){
	open (IN,"$mapinfo_list[$i]");
	while (<IN>){
	chomp;
	next if (/^\s+$/);
	($number,undef,$free_number,$type)=split(/\s/,$_,4);
	if ($type =~ /mapped/){
		if ($type !~ /mate/){
	$samp=$mapinfo_list[$i];
        $samp =~s/\.\///g;
        $samp =~s/\.source\.mapinfo//g;
        chomp($samp);
	$type=~s/mapped\s\(//g;
	$type=~s/\s\:\sN\/A\)//g;
	print OUT "$samp\t$number\t$type\t";
}
}
	if ($type =~ /paired in sequencing/){
        print OUT "$number\t";
        my$basenum=$number*$PE;
        print OUT "$basenum\t";
	my$depth=sprintf "%.3f",$basenum/$ref_genome_size;
	my$depth_out=$depth."X";
	print OUT "$depth_out\t";
        }
	if ($type=~ /properly paired/){
	$type =~ s/properly\spaired \(//g;
	$type =~ s/\s\:\sN\/A\)//g;
	print OUT "$number\t$type\n";
	}
}
}
close IN;
close OUT;
print "Pick up mapping info is ok.\n";
#################End####################
