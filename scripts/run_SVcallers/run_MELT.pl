#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $input_bam = '';
my $ref = '';
my $target_chr = 'ALL';
my $out_prefix = '';
my $melt_jar = '';
my $melt_lib = '';
my $melt_gene = '';
my $melt_targets = '';
my $java_path = '';
my $bowtie_path = '';
my $excl_chr_list = '';
my $non_human = 0;
my $MOP_dir = '';

my $cores = 1;
my $read_length = 150;
my $read_coverage = 30;

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'melt_jar|mj=s' => \$melt_jar,
    'melt_lib|ml=s' => \$melt_lib,
    'melt_gene|mg=s' => \$melt_gene,
    'melt_target|mt=s' => \$melt_targets,
    'java_path|jp=s' => \$java_path,
    'bowtie_path|bp=s' => \$bowtie_path,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'excl_chr|ex=s' => \$excl_chr_list,
    'threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_MELT.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -mj <MELT.jar path> -ml <MELT reference library directory> -mg <MELT gene annotation bed> -n <threads> (-mt <MELT target list>)

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory and contain RG tag [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --melt_jar or -mj <STR>    MELT jar executable path (e.g., ..../MELTv2.2.2/MELT.jar) [mandatory]
   --melt_lib or -ml <STR>    MELT mobile element reference data directory, containing *_MELT.zip files (* corresponds to target element specified with -mt option) 
                              (e.g., ..../MELTv2.2.2/me_refs/Hg38) [mandatory]
   --melt_gene or -mg <STR>   MELT gene annotation file (e.g., ..../MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed) [mandatory]
   --prefix or -p <STR>       outpout prefix [mandatory]
   --melt_target or -mt <STR> MELT mobile element or other sequence target comma-separated list [default: ALU,LINE1,SVA,HERVK]
   --java_path or -jp <STR>   directory containing java (v1.8 or later) executable, unless the corrsponding version of java is set in $PATH
   --bowtie_path or -bp <STR> directory containing bowtie2 executable, unless bowtie2 is set in $PATH
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --excl_chr or -ex <STR>    exclude chromosome list separated by '/' (MELT excludes < 1 Mb reference chromosomes by default) [default for human: decoy contigs]
   --threads or -n <INT>      number of threads to be used [default: 1]
   --help or -h               output help message
   
=cut


$ENV{PATH} = "$bowtie_path:$ENV{PATH}" if ($bowtie_path ne '');
$ENV{PATH} = "$java_path:$ENV{PATH}" if ($java_path ne '');

my $java_path2 = `which java`;
print STDERR "$java_path2\n";
my $bowtie_path2 = `which bowtie2`;
print STDERR "$bowtie_path2\n";

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

my $ref_index = "$ref.fai";
my $min_chr_size = 1000000;

if (($non_human == 0) and ($excl_chr_list eq '')){
    my @excl_list;
    open (FILE, $ref_index) or die "$ref_index is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $size) = split (/\t/, $line);
        if ($target_chr eq 'ALL'){
            push @excl_list, $chr if ($chr !~ /^c*h*r*[\dXY]+$/i) and ($size >= 1000000);
        }
        else{
            push @excl_list, $chr if (!exists $target_chr{$chr});
        }
        
    }
    close (FILE);
    $excl_chr_list = join ('/', @excl_list) if (@excl_list > 0);
}
else{
    my %excl_chr;
    if ($excl_chr_list ne ''){
        my @excl_list = split (/\//, $excl_chr_list);
        map{$excl_chr{$_} = 1} @excl_list;
    }
    open (FILE, $ref_index) or die "$ref_index is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $size) = split (/\t/, $line);
        next if (exists $excl_chr{$chr});
        $min_chr_size = $size if ($size < $min_chr_size);
        if ($target_chr ne 'ALL'){
            $excl_chr_list .= "$chr/" if (!exists $target_chr{$chr});
        }
    }
    close (FILE);
    $min_chr_size -= 1000;
    $excl_chr_list =~ s/\/$// if ($excl_chr_list =~ /\/$/);
}

my %chr_len;
my @read_dist;
my @read_len;
my $read_num = 0;
my $ave_dist = 480;
my $coverage = 30;
my $first_chr = '';

open (FILE, $input_bam) or die "$input_bam: $!" if ($input_bam =~ /\.sam$/);
open (FILE, "samtools view -h $input_bam |") or die "$input_bam: $!" if ($input_bam =~ /\.bam$/);
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^\@/){
        if ($line =~ /^\@SQ\tSN:(\S+)\tLN:(\d+)/){
            $chr_len{$1} = $2;
        }
        next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[2];
    $first_chr = $chr if ($first_chr eq '');
    next if ($chr eq '*');
    if ($chr eq $first_chr){
        my $tag = $line[1];
        my $mapq = $line[4];
        my $distance = $line[8];
        my $readlen = length $line[9];
	push @read_len, $readlen;
        $read_num ++;
        if (($tag == 99) or ($tag == 163)){
            if ($distance > 0){
                if (($mapq >= 20) and ($distance < 1000)){
                    push @read_dist, $distance;
                }
            }
        }
    }
    else{
        last;
    }
}
close (FILE);


my $first_chr_len = $chr_len{$first_chr};
$coverage = int ($read_num * $read_length / $first_chr_len);
$coverage = 30 if ($coverage == 0);

print STDERR "$first_chr: LEN-$first_chr_len\treads: $read_num\tread_len: $read_length\tcoverage: $coverage\n";

if (@read_dist > 1){
    my $sum_insert = 0;
    foreach (@read_dist){
	$sum_insert += $_;
    }
    $ave_dist = int ($sum_insert / @read_dist);
}
my $sum_read_len = 0;
foreach (@read_len){
    $sum_read_len += $_;
}
$read_length = int ($sum_read_len / @read_len + 0.5);

my $bam_base = basename ($input_bam);
my $input_bai = "$input_bam.bai";
if (!-f $input_bai){
    my $bai2 = $input_bai;
    my $bam_base2 = $1 if ($input_bam =~ /(.+)\.bam/);
    $input_bai = "$bam_base2.bai";
    if (!-f $input_bai){
        die "bam index file ($bai2 or $input_bai) is not found: \n";
    }
}

system ("ln -s $input_bam") if (!-f $bam_base);
system ("ln -s $input_bai") if (!-f "$bam_base.bai");

my $command_preprocess = "java -jar $melt_jar Preprocess -bamfile $bam_base -h $ref 2>$out_prefix.preprocess.log";

open (OUT, ">> $out_prefix.command.log");
print OUT "MELT preprocess: $command_preprocess\n";
close (OUT);

system ("rm $bam_base.disc") if (-f "$bam_base.disc");
system ("rm $bam_base.disc.bai") if (-f "$bam_base.disc.bai");
system ("rm $bam_base.fq") if (-f "$bam_base.fq");

system ("$command_preprocess");

my ($cur_dir) = `pwd`;
chomp $cur_dir;

my @MELT_target = split (/,/, $melt_targets);

my @out_files;

foreach my $target (@MELT_target){
    my $target_zip = $target . '_MELT.zip';
    my $target_ref = "$melt_lib/$target_zip";
    if (!-f $target_ref){
        die "$target_ref is not found:\n";
    }
    my $command_melt = "java -Xmx36G -jar $melt_jar Single -bamfile $cur_dir/$bam_base -h $ref -n $melt_gene -t $target_ref -w . -r $read_length -e $ave_dist -d $min_chr_size -c $coverage -b $excl_chr_list 2>$out_prefix.MELT-$target.log";
    $command_melt = "java -Xmx36G -jar $melt_jar Single -bamfile $cur_dir/$bam_base -h $ref -n $melt_gene -t $target_ref -w . -r $read_length -e $ave_dist -d $min_chr_size -c $coverage 2>$out_prefix.MELT-$target.log" if ($excl_chr_list eq '');
    open (OUT, ">> $out_prefix.command.log");
    print OUT "MELT $target: $command_melt\n";
    close (OUT);
    my $out_vcf = "$target.final_comp.vcf";
    system ("$command_melt") if (!-f $out_vcf) or (-z $out_vcf);
    push @out_files, $out_vcf;
}

my $str = join (' ', @out_files);
if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_MELT_vcf.pl $str > MELT.$out_prefix.vcf");
}
else{
    system ("$Bin/convert_MELT_vcf.pl $str > MELT.$out_prefix.vcf");
}

print STDERR "MELT run was completed\n";
