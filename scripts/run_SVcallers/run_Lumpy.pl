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
my $lumpy_path = '';
my $lumpy_scripts = '';
my $samtools_path = '';
my $non_human = 0;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'command_path|cp=s' => \$lumpy_path,
    'lumpy_scripts|ls=s' => \$lumpy_scripts,
    'samtools_path|sp=s' => \$samtools_path,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);


=head1 SYNOPSIS

  run_Lumpy.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -cp <lumpy command path> -ls <lumpy scripts directory>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory and contain RG tag [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR>  lumpy executable path (e.g., ..../lumpy-sv/bin/lumpy) [mandatory]
   --lumpy_scripts or -ls <STR> lumpy scripts directory (e.g., ..../lumpy-sv/scripts) [mandatory]
   --samtools_path or -bp <STR> samtools path unless samtools is set in $PATH
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --help or -h               output help message
   
=cut

$ENV{PATH} = "$samtools_path:$ENV{PATH}" if ($samtools_path ne '');

my $input_base = basename ($input_bam);

my %chr_len;
my @read_dist;
my @read_len;
my $read_num = 0;
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

my $sum_readlen = 0;
my $sum_inssize = 0;
my $ave_readlen = 150;
my $ave_inssize = 500;

map{$sum_readlen += $_} @read_len;
map{$sum_inssize += $_} @read_dist;

$ave_readlen = int ($sum_readlen / $read_num + 0.5);
$ave_inssize = int ($sum_inssize / @read_dist + 0.5);

my $sum_diff = 0;
map{$sum_diff += ($_ - $ave_inssize) ** 2} @read_dist;
my $ins_sd = int (($sum_diff / @read_dist) ** 0.5 + 0.5);

system ("samtools view -uF 0x0002 $input_bam | samtools view -uF 0x100 - | samtools view -uF 0x0004 - | samtools view -uF 0x0008 - | samtools view -bF 0x0400 - | samtools sort - -T $out_prefix.discordant.sort -o $out_prefix.discordant.sort.bam");

system ("samtools view -h $input_bam | $lumpy_scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - | samtools sort - -T $out_prefix.sr.sort -o $out_prefix.sr.sort.bam");

system ("samtools view $input_bam | tail -n+100000 | python $lumpy_scripts/pairend_distro.py -r $ave_readlen -X $ins_sd -N 10000 -o $out_prefix.histo");

my $command = "$lumpy_path -mw 4 -tt 0.0 -pe bam_file:$out_prefix.discordant.sort.bam,histo_file:$out_prefix.histo,mean:$ave_inssize,stdev:$ins_sd,read_length:$ave_readlen,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:20 -sr bam_file:$out_prefix.sr.sort.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:20 > $out_prefix.vcf 2>$out_prefix.Lumpy.log";

open (OUT, ">> $out_prefix.command.log");
print OUT "samtools view $input_bam | tail -n+100000 | python $lumpy_scripts/pairend_distro.py -r $ave_readlen -X $ins_sd -N 10000 -o $out_prefix.histo\n";
print OUT "Lumpy command: $command\n";
close (OUT);

system ("$command");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_Lumpy_vcf.pl $out_prefix.vcf $non_human $target_chr > Lumpy.$out_prefix.vcf");
}
else{
    system ("$Bin/convert_Lumpy_vcf.pl $out_prefix.vcf $non_human $target_chr > Lumpy.$out_prefix.vcf");
}

print STDERR "Lumpy run was completed\n";
