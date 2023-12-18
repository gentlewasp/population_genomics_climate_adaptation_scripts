########################################################################
#	File Name: globBWAMEM.pl
#	> Author: QiangGao
#	> Mail: qgao@genetics.ac.cn 
#	Created Time: Fri 11 Nov 2016 11:22:41 AM CST
#	Modified time: Mon 30 Dec 2019 16:43 PM BGI
#	Modified Author: wangzhiwei
#########################################################################

############nohup to submit the mission


#!/usr/bin/perl
######data setting
use strict;

my $queue="xhacnormalb";
my $CORE=14;
my $pe=150;
my $REF="/work/home/malijun/renyajing/00_genome/genome.fasta";
my $input_raw="/work/home/malijun/renyajing/02_clean_fastp_data/0_ML"; ###including all samples in the path, and withou sub-folds
my $prefix="TP";
my $suffix=".R1.fastp.fq.gz";
&help unless defined ($prefix && $REF && $input_raw && $suffix);
#######default software diray
#######default software diray
my$script_dir="/work/home/malijun/renyajing/03_mapping"; #done
my$samtool_dir="/work/share/kdy_caolijun/software/samtools/samtools-1.14"; #done
my$bwa_dir="/work/share/kdy_caolijun/software/bwa"; #done
my$picard_dir="/work/share/kdy_caolijun/software/picard.jar"; #done
my$seqkit_dir="/work/share/kdy_caolijun/software/seqkit";  #done
my$Qualimap_dir="/work/share/kdy_caolijun/software/qualimap_v2.2.1/qualimap1"; #done
if ($ARGV[0] ne "yes" && $ARGV[0] ne "YES"){ &help;die;}
sub help{
print "\n";
print "Mission Prefix		= $prefix\n";
print "suffix			= $suffix\n";
print "queue 			= $queue\n";
print "CPU   			= $CORE\n";
print "length			= PE150\n";
print "REF  			= $REF\n";
print "input fastq dir		= $input_raw\n";
print "script			= $script_dir\n";
print "samtools 		= $samtool_dir\n";
print "bwa			= $bwa_dir\n";
print "seqkit			= $seqkit_dir\n";
print "qualimap		= $Qualimap_dir\n\n";

print "if all parameter is OK, please type \n\n";
print "\t nohup perl $0 yes & \n\n";
}

#################Start
my @file=`find $input_raw -name "*$suffix" `;


if(-e "tmp_pipe_data"){
	
}else{
	my $cc=`mkdir -p tmp_pipe_data`;
}

if(-e "tmp_pipe_log"){
	
}else{
	my $cc=`mkdir -p tmp_pipe_log`;
}
my$bwa_mem=$CORE*3;
$bwa_mem.="g";
my$name;
my$total_number=scalar(@file);
open(USED,">used.data");
foreach(@file){
	chomp $_;
	($name)=$_=~/.*\/(.*?)$suffix/;
        print $name;
	my $R2=$_;
	$R2=~s/R1/R2/;
	open(OUT,">tmp_pipe_log/$name.map.slurm");
print OUT <<EOF;
#!/bin/bash
#SBATCH -J $prefix\_bwa$name
#SBATCH -p $queue
#SBATCH -n $CORE
#SBATCH -N 1
#SBATCH --mem=$bwa_mem
#SBATCH -o tmp_pipe_log/$name.slurm.log
#SBATCH -e tmp_pipe_log/$name.slurm.err

export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.131-11.b12.el7.x86_64/
                              java-1.8.0-openjdk-1.8.0.181-7.b13.el7.x86_64
export CLASSPATH=.:\$JAVA_HOME/lib:\$JAVA_HOME/jre/lib
export PATH=\$JAVA_HOME/bin:\$JAVA_HOME/jre/bin:\$PATH

$bwa_dir/bwa mem -M -t $CORE -R \"\@RG\\tID:$name\\tLB:$name\\tSM:$name\\tPL:illumina\\tPU:$name\" $REF $_ $R2| $samtool_dir/samtools view -@ $CORE -bS - > tmp_pipe_data/$name.source.bam
$samtool_dir/samtools sort -m 3g -@ $CORE tmp_pipe_data/$name.source.bam -o tmp_pipe_data/$name.noq30.sorted.bam
$samtool_dir/samtools stats tmp_pipe_data/$name.noq30.sorted.bam > tmp_pipe_data/$name.noq30.sorted.bam.stat
perl $script_dir/02_getNM.pl tmp_pipe_data/$name.noq30.sorted.bam > tmp_pipe_data/$name.noq30.sorted.bam.NM
$samtool_dir/samtools flagstat tmp_pipe_data/$name.source.bam  >tmp_pipe_data/$name.source.mapinfo
rm tmp_pipe_data/$name.source.bam
$samtool_dir/samtools view -@ $CORE -h -b -q30 tmp_pipe_data/$name.noq30.sorted.bam > tmp_pipe_data/$name.sorted.bam
$picard_dir/picard -Xmx6g MarkDuplicates I=tmp_pipe_data/$name.sorted.bam O=tmp_pipe_data/$name.sorted.rmdup.bam CREATE_INDEX=true REMOVE_DUPLICATES=true M=tmp_pipe_data/$name.marked_dup_metrics.txt
###########
######
rm tmp_pipe_data/$name.sorted.bam
EOF
  close OUT;
  next if(-e "tmp_pipe_data/$name.source.mapinfo");
#  my $cc=`sbatch < tmp_pipe_log/$name.map.slurm >> bwa.monitor.list`;
  print "run $name\n";
}
close USED;
#####Genome length
my$ref_length=`$seqkit_dir/seqkit stat ${REF}.fasta |tail -1 |awk '{OPS="\s+";print \$5}'`;
$ref_length =~ s/,//g;
my$now_path=`pwd`;
###############monitor

#sleep 60;
#my$bwa_mission_all = `cat bwa.monitor.list|grep -P '\d\+'|cut -d ' ' -f 4 `;
#print "$bwa_mission_all\n";
#my@bwa_mission=split(/\n/,$bwa_mission_all);
#my$bwa_run_number=0;
#for (my$b=0;$b<scalar(@bwa_mission);$b++){
#        my$count_number=`squeue |awk '{print \$1}'|grep $bwa_mission[$b]|wc -l`;
#        $bwa_run_number+=$count_number;
#}
#until($bwa_run_number == 0){
#        sleep 300;
#        $bwa_run_number=0;
#        for (my$b=0;$b<scalar(@bwa_mission);$b++){
#        my$count_number=`squeue |awk '{print \$1}'|grep $bwa_mission[$b]|wc -l`;
#        $bwa_run_number+=$count_number;
#        }
#}


sleep 60;
my$mapinfo_file=`ls tmp_pipe_data/*mapinfo |wc -l `;
chomp($mapinfo_file);
until($mapinfo_file == $total_number){
        sleep 300;
        print "circle is on .\n";
        $mapinfo_file=`ls tmp_pipe_data/*mapinfo |wc -l`;
        chomp ($mapinfo_file);
        print "$total_number\t $mapinfo_file\t .\n";
}
my$mapinfo_all=`find tmp_pipe_data -name "*mapinfo"`;
my@mapinfo_single=split(/\n/,$mapinfo_all);
my$tem_count=0;
for(my$e=0;$e<scalar(@mapinfo_single);$e++){
        if (-z $mapinfo_single[$e]){
        $tem_count++;
        }
}
until($tem_count ==0){
        sleep 60;
        $tem_count=0;
        for(my$e=0;$e<scalar(@mapinfo_single);$e++){
        if (-z $mapinfo_single[$e]){
        $tem_count++;
                }
        }
}

if ($tem_count != 0 ){die "$tem_count\tmapinfo is zero file,error.\n";}
`rm tmp_pipe_data/*source.bam`;
#####check sorted.bam be killed or not
#my$kill=`grep -E "kill|cound't|out-of-memory" tmp_pipe_log/*err `;
#my$kill_mission=`grep -E "kill|cound't|out-of-memory" tmp_pipe_log/*err |wc -l`;
#if ($kill_mission != 0 ){
#print "$kill\n";
#die "some bwa mission was killed, there are $kill_mission mission was killer.\n";}

######pick up
`perl $script_dir/02_pick_up_mapping_quality.pl -input tmp_pipe_data -outdir tmp_pipe_data -pe $pe -size $ref_length`;
print "pick up is complete.\n";
####modified by wangzw
if(-e "tmp_pipe_data/Qualimap_data"){

}else{
        my $dd=`mkdir -p tmp_pipe_data/Qualimap_data`;
}
my@bam_list = `find tmp_pipe_data/ -name "*.noq30.sorted.bam"|grep -v "Qualimap"`;
for(my $i=0;$i<scalar(@bam_list);$i++){
        my $single_slurm=$bam_list[$i];
        chomp($single_slurm);
        $single_slurm =~ s/.*\///g;
	open(SLURM,">tmp_pipe_data/Qualimap_data/$single_slurm\.slurm");
	print SLURM <<"SLURM";
#!/bin/bash
#SBATCH -J $prefix\_qualimap_$single_slurm
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=36g
#SBATCH -p $queue
#SBATCH -o tmp_pipe_data/Qualimap_data/$single_slurm\.log
#SBATCH -e tmp_pipe_data/Qualimap_data/$single_slurm\.err

$Qualimap_dir/qualimap bamqc -bam tmp_pipe_data/$single_slurm -nt 12 -c --java-mem-size=4G -outdir tmp_pipe_data/Qualimap_data/$single_slurm -outformat HTML
SLURM
`sbatch < tmp_pipe_data/Qualimap_data/$single_slurm\.slurm >>qualimap.monitor.list`;
}
sleep 30;
#my$qualimap_num=`grep "Finish" tmp_pipe_data/Qualimap_data/*log |wc -l`;
#chomp($qualimap_num);
#until($qualimap_num == $total_number){
#	sleep 100;
#	print "waiting qualimap result.\n";
#	$qualimap_num=`grep "Finished" tmp_pipe_data/Qualimap_data/*log |wc -l`;
#	chomp($qualimap_num);
#	}

###my$quali_mission_all = `cat qualimap.monitor.list|grep -P '\d\+'|cut -d ' ' -f 4 `;
###print "$quali_mission_all";
###my@quali_mission=split(/\n/,$quali_mission_all);
###my$quali_run_number=0;
###for (my$b=0;$b<scalar(@quali_mission);$b++){
###        my$count_number=`squeue |awk '{print \$1}'|grep $quali_mission[$b]|wc -l`;
###        $quali_run_number+=$count_number;
###}
###until($quali_run_number == 0){
###        sleep 300;
###        $quali_run_number=0;
###        for (my$b=0;$b<scalar(@quali_mission);$b++){
###        my$count_number=`squeue |awk '{print \$1}'|grep $quali_mission[$b]|wc -l`;
###        $quali_run_number+=$count_number;
###        }
###}


chdir "tmp_pipe_data" or die "cannot cd dir tmp_pipe_data";
`rm *.noq30.sorted.bam*`;

#`perl $script_dir/02_pick.coverage.pl Qualimap_data`;
#`mv Qualimap_data/Mapping.coverage .`;
#`rm -rf Qualimap_data`;
print "Thanks.\n";
