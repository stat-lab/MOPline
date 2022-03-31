#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $bam_list = '';
my $bam_dir = '';

my $help;

GetOptions(
    'bam_list|b=s' => \$bam_list,
    'bam_dir|bd=s' => \$bam_dir,
    'help' => \$help
) or pod2usage(-verbose => 0);


=head1 SYNOPSIS

  run_Lumpy.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -cp <lumpy command path> -ls <lumpy scripts directory>

  Options:
   --bam_list or -b <STR>     list file of bam file names (not bam file path) [mandatory]
   --bam_dir or -bd <STR>     directory path containing bam files [mandatory]
   --help or -h               output help message
   
=cut

die "-b or -bd not specified:\n" if ($bam_list eq '') or ($bam_dir eq '');

open (FILE, $bam_list) or die "$bam_list is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my $bam_name = $line;
    my $sample_name = '';
    $sample_name = $1 if ($bam_name =~ /(.+)\.bam/);
    my $bam = "$bam_dir/$bam_name";
    mkdir $sample_name if (!-d $sample_name);
    chdir $sample_name;
    system ("ln -s $bam");
    system ("ln -s $bam.bai");
    chdir '..';
}
close (FILE);
