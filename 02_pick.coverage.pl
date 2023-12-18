#!/usr/bin/perl

my$input_dir=$ARGV[0];
die "perl $0 <input_dir> "if (scalar(@ARGV) !=1);


my%all_1X_hash;
my%all_2X_hash;
my%all_4X_hash;
my%all_6X_hash;
my%all_8X_hash;
my%all_10X_hash;
my%all_20X_hash;
my%all_30X_hash;
my%all_40X_hash;
my%all_50X_hash;
my$file_all=`find $input_dir -name "*sorted.bam"| grep -v 'log' |grep -v 'err'|grep -v 'pbs'| grep -v 'slurm'`;
my@file_single=split("\n",$file_all);
for(my$a=0;$a<scalar(@file_single);$a++){
	my($name)=$file_single[$a]=~/.*\/(.*?)\.sorted\.bam/;
	$name=~s/Clean_//g;
	$name=~s/\.sam//g;
	open (A,"<$file_single[$a]/genome_results.txt") or die "cannot open raw result files\n";
	while(<A>){
		chomp;
		next if (/^\s+$/);
		if (/\s+There .* >= 1X/){
			my(undef,undef,undef,undef,$coverage_1X,undef)=split(/\s+/,$_,6);
			$all_1X_hash{$name}=$coverage_1X;
		}
		if (/\s+There .* >= 2X/){
			my(undef,undef,undef,undef,$coverage_2X,undef)=split(/\s+/,$_,6);
			$all_2X_hash{$name}=$coverage_2X;
		}
		if (/\s+There .* >= 4X/){
			my(undef,undef,undef,undef,$coverage_4X,undef)=split(/\s+/,$_,6);
			$all_4X_hash{$name}=$coverage_4X;
			}
		if (/\s+There .* >= 6X/){
			my(undef,undef,undef,undef,$coverage_6X,undef)=split(/\s+/,$_,6);
			$all_6X_hash{$name}=$coverage_6X;
		}
		if (/\s+There .* >= 8X/){
			my(undef,undef,undef,undef,$coverage_8X,undef)=split(/\s+/,$_,6);
			$all_8X_hash{$name}=$coverage_8X;
		}
		if (/\s+There .* >= 10X/){
			my(undef,undef,undef,undef,$coverage_10X,undef)=split(/\s+/,$_,6);
			$all_10X_hash{$name}=$coverage_10X;
			}
		if (/\s+There .* >= 20X/){
                        my(undef,undef,undef,undef,$coverage_20X,undef)=split(/\s+/,$_,6);
                        $all_20X_hash{$name}=$coverage_20X;
                        }
		if (/\s+There .* >= 30X/){
                        my(undef,undef,undef,undef,$coverage_30X,undef)=split(/\s+/,$_,6);
                        $all_30X_hash{$name}=$coverage_30X;
                        }
		if (/\s+There .* >= 40X/){
                        my(undef,undef,undef,undef,$coverage_40X,undef)=split(/\s+/,$_,6);
                        $all_40X_hash{$name}=$coverage_40X;
                        }
		if (/\s+There .* >= 50X/){
                        my(undef,undef,undef,undef,$coverage_50X,undef)=split(/\s+/,$_,6);
                        $all_50X_hash{$name}=$coverage_50X;
                        }	
	}
	close A;
}

open (OUT,">$input_dir/Mapping.coverage") or die "cannot open Mapping.coverage file .\n";
print OUT "Sample\t1X_depth_coverage\t2X_depth_coverage\t4X_depth_coverage\t6X_depth_coverage\t8X_depth_coverage\t10X_depth_coverage\t20X_depth_coverage\t30X_depth_coverage\t40X_depth_coverage\t50X_depth_coverage\n";
foreach my $key( keys %all_2X_hash){
	print OUT "$key\t$all_1X_hash{$key}\t$all_2X_hash{$key}\t$all_4X_hash{$key}\t$all_6X_hash{$key}\t$all_8X_hash{$key}\t$all_10X_hash{$key}\t$all_20X_hash{$key}\t$all_30X_hash{$key}\t$all_40X_hash{$key}\t$all_50X_hash{$key}\n";
}
close OUT;
chdir $input_dir;
`rm -r */css`;
`rm -r */images_qualimapReport`;
`rm -r */raw_data_qualimapReport`;
`rm -r */qualimapReport.html`;
print "Thanks.\n";
