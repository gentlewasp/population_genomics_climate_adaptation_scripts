#########################################################################
#	File Name: getNM.pl
#	> Author: qgao
#	> Mail: qgao@genetics.ac.cn 
#	Created Time: Sun 12 Jun 2016 01:12:58 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my ($file)=@ARGV;
if(!$file){
	print "USAGE:\nperl $0 sample_1.bam $file\n";
	exit;
}
my ($out)=$file;
open(OUT,">$out.NM");
print OUT "$file\tNumber\tRatio\n";
if($file=~/\.bam$/){
	open(IN,"samtools view $file|");
}else{
	open(IN,"$file");
}
my %hash;
my $all;
while(<IN>){
	next if($_=~/^\@/);
	next if($_=~/\s\*\s/);
	my ($num)=$_=~/NM\:i\:(\d+)/;
	#print "$num\n";
	$all+=1;
	$hash{$num}+=1;
}
close IN;
for(my $i=0;$i<100;$i++){
	my $now=0;
	if(exists $hash{$i}){
		$now=$hash{$i};
	}
	my $ra=int($now/$all*10000)/100;
	print OUT "$i\t$now\t$ra\%\n";
}
close OUT;
