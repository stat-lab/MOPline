#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

# this script assumes that the current working directory contains sample directories.

my $mopline_dir = '';
my $sample_list = '';
my $ref = '';
my $non_human = 0;
my $read_length = 0;
my $help;

GetOptions(
    'mopline_dir|md=s' => \$mopline_dir,
    'sample|s=s' => \$sample_list,
    'ref|r=s' => \$ref,
    'read_len|rl=i' => \$read_length,
    'non_human|nh=i' => \$non_human,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_create_cov_batch.pl -s <sample list file> -r <ref fasta> -rl <read length> 
  (-nh 1 if sample is a non-human species)
  This script assumes that sample directories (${sample_name})are under the working directory and each bam file ($sample_name}.bam) is present in each sample directory.)

  Options:
   --sample or -s <STR>     sample list file [mandatory]
   --ref or -r <STR>        absolute path of reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --read_len or -rl <INT>  read length [mandatory for create_cov]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --help or -h             output help message
   
=cut

my $cur_dir = `pwd`;
chomp $cur_dir;

if (!-f "$Bin/create_coverage_file_bam_single.pl"){
    die "$Bin/create_coverage_file_bam_single.pl is not found:\n";
}

open (FILE, $sample_list) or die "$sample_list is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^$|^#/);
    my ($sample_name) = split (/\s+/, $line);
    chdir $sample_name;
    my $bam = "$cur_dir/$sample_name/$sample_name.bam";
    die "$bam is not found:\n" if (!-f $bam);
    my $run_command = '';
    my $run_script = "$Bin/create_coverage_file_bam_single.pl";
    die "reference fasta not specified:\n" if ($ref eq '');
    die "$ref is not found:\n" if (!-f $ref);
    die "read length not specified:\n" if ($read_length <= 0);
    $run_command = "$run_script -b $bam -r $ref -sn $sample_name -rl $read_length";
    $run_command .= " -nh 1" if ($non_human == 1);
    $run_command .= " 2>$sample_name.create_cov.log";
    system ("$run_command");
    chdir '..';
}
close (FILE);
