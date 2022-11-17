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
my $matchclip_path = '';
my $cores = 1;
my $non_human = 0;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c=s' => \$target_chr,
    'command_path|cp=s' => \$matchclip_path,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_MATCHCLIP.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -ms <matchclips command path> -n <threads>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR>  MATCHCLIP executable path, matchclips (e.g., ..../matchclips2/matchclips) [mandatory]
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --threads or -n <INT>      number of threads to be used [default: 1]
   --help or -h               output help message
   
=cut

my $command = "$matchclip_path -f $ref -b $input_bam -t $cores -o $out_prefix.matchclip.out 2>$out_prefix.run.log";
    
open (OUT, ">> $out_prefix.command.log");
print OUT "Matchclip command: $command\n";
close (OUT);

system ("$command");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_MATCHCLIP_vcf.pl $out_prefix.matchclip.out $non_human $target_chr > MATCHCLIP.$out_prefix.vcf");
}
else{
  system ("$Bin/convert_MATCHCLIP_vcf.pl $out_prefix.matchclip.out $non_human $target_chr > MATCHCLIP.$out_prefix.vcf");
}

print STDERR "MATCHCLIP run was completed\n";
