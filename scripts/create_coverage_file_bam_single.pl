#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $bam = '';
my $sample_name = '';
my $sample_dir = '.';
my $ref = '';
my $non_human = 0;
my $out_prefix = '';
my $read_length = 0;
my $exclude_discread = 0;
my $min_MAPQ = 1;
my $bin_size = 50;

my $exclude_mit = 1;

my $help;
 
GetOptions(
    'bam|b=s' => \$bam,
    'sample_name|sn=s' => \$sample_name,
    'sample_dir|sd=s' => \$sample_dir,
    'ref|r=s' => \$ref,
    'non_human|nh=i' => \$non_human,
    'bin_size|bs=i' => \$bin_size,
    'mapq|q=i' => \$min_MAPQ,
    'ex_disc|xd' => \$exclude_discread,
    'read_len|rl=i' => \$read_length,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  create_coverage_file_bam_single.pl -b <bam file path> -sn <sample name> -sd <sample dir> -r <ref fasta> -rl <read length> 
  (-nh 1 if sample is a non-human species)
  (create files with coverage and split read information of read alignment at 50-bp window for a single bam)
  outout Cov directory will be created in the working directory or under ${sample_name} directory if sample-bam/cram table or sample list file is provided

  Options:
   --bam or -b <STR>        file path of input bam. [mandatory]
   --sample_name or -sn <STR> sample name (ID) [mandatory]
   --sample_dir or -sd <STR>  sample directory path [default: ./]
   --ref or -r <STR>        reference fasta file used for creating cram (only needed for cram files as input) [mandatory]
   --read_len or -rl <INT>  read length [mandatory]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --bin_size or -bs <INT>  bin size (bp) for calculating coverage [default: 50]
   --mapq or -q <INT>       minimum mapping quality [default: 1]
   --ex_disc or -xd         exclude read pairs with discordant insert size [default: false]
   --prefix or -p <STR>     output prefix
   --help or -h             output help message
   
=cut

my @chr_list;

die "bam file is not specified:" if ($bam eq '');
die "Reference index file: $ref.fai is not found:$!\n" if (!-f "$ref.fai");
die "sample name is not specified: \n" if ($sample_name eq '');
die "read length is not specified or too short:\n" if ($read_length < 10);

open (FILE, "$ref.fai") or die "$ref.fai is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $size) = split (/\t/, $line);
    next if ($chr =~ /^c*h*r*Mit$|^c*h*r*Mt$|^chrM$/i) and ($exclude_mit == 1);
    if ($non_human == 0){
        if ($chr =~ /^c*h*r*[\dXY]+$/){
            push @chr_list, $chr;
        }
    }
    else{
        push @chr_list, $chr;
    }
}
close (FILE);

foreach my $chr (@chr_list){
    my $command = "$Bin/create_coverage_file_bam_chr.pl -b $bam -c $chr -rl $read_length -r $ref -sn $sample_name -sd $sample_dir -bs $bin_size -rl $read_length -q $min_MAPQ";
    $command = "$Bin/create_coverage_file_bam_chr.pl -b $bam -c $chr -rl $read_length -r $ref -sn $sample_name -sd $sample_dir -bs $bin_size -rl $read_length -q $min_MAPQ -xd" if ($exclude_discread == 1);
    
    open (OUT, ">> $out_prefix.command.log");
    print OUT "$command\n";
    close (OUT);
    
    system ("$command");
}

print STDERR "Coverage calculation completed:\n";
