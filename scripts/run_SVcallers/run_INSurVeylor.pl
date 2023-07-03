#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

# run INSurVeylor using singularity image file

my $input_bam = '';
my $ref = '';
my $out_prefix = '';
my $sif_path = '';
my $temp_dir = '';
my $bind_dir = '';
my $no_home = 0;
my $cores = 1;
my $non_human = 0;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'command_path|cp=s' => \$sif_path,
    'prefix|p=s' => \$out_prefix,
    'temp_dir|td=s' => \$temp_dir,
    'bind_dir|bd=s' => \$bind_dir,
    'no_home|noh' => \$no_home,
    'mopdir|md=s' => \$MOP_dir,
    'threads|n=i' => \$cores,
    'non_human|nh=i' => \$non_human,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_INSurVeylor.pl -cp <singularity sif file path> -b <full path of input bam file> -r <full path of reference fasta> -p <output prefix> -td <full path of temp directory> -n <threads>

  Options:
   --bam or -b <STR>          absolute path of input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>          absolute path of reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --command_path or -cp <STR>  INSurVeylor singularity sif/container file path [mandatory]
   --temp_dir or -td <STR>  absolute path of tmp directory on the host [mandatory]
   --bind_dir or -bd <STR>  comma-separated list of absolute path (except for the working directory and directories containing input bam and reference fasta) on the host to be added to singularity container (optional)
   --no_home or -noh <BOOLEAN>  do not add $HOME on the host to singularity container [default: false]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --threads or -n <INT>      number of threads to be used [default: 1]
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --help or -h               output help message
   
=cut

my $work_dir = `pwd`;
chomp $work_dir;

my $bam_dir = '';
my $ref_dir = '';

$bam_dir = $1 if ($input_bam =~ /^(.+)\/.+\.bam$/);
$bam_dir = $1 if ($bam_dir eq '') and ($input_bam =~ /^(.+)\/.+\.cram$/);

$ref_dir = $1 if ($ref =~ /^(.+)\//);

$bind_dir = $work_dir;
$bind_dir .= ",$temp_dir" if ($temp_dir ne '') and ($bind_dir ne $work_dir);

$bind_dir .= ",$bam_dir" if ($bam_dir ne '') and ($bam_dir ne $work_dir);
$bind_dir .= ",$ref_dir" if ($ref_dir ne '') and ($ref_dir ne $work_dir) and ($ref_dir ne $bam_dir);


my $command = "singularity run --bind $bind_dir $sif_path --threads $cores $input_bam . $ref 2>$out_prefix.insurveyor.log";
$command = "singularity run --bind $bind_dir $sif_path --no-home --threads $cores $input_bam . $ref 2>$out_prefix.insurveyor.log" if ($no_home == 1);
    
open (OUT, ">> $out_prefix.command.log");
print OUT "INSurVeylor command: $command\n";
close (OUT);

system ("$command");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_INSurVeylor_vcf.pl out.pass.vcf.gz $non_human > INSurVeylor.$out_prefix.vcf");
}
else{
  system ("$Bin/convert_INSurVeylor_vcf.pl out.pass.vcf.gz $non_human > INSurVeylor.$out_prefix.vcf");
}

print STDERR "INSurVeylor run was completed\n";
