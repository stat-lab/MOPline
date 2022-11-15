#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

# running bsub command passed from job_sv_call_lsf_1.pl script

my $input_bam = '';
my $ref = '';
my $target_chr = 'ALL';
my $out_prefix = '';
my $manta_script = '';
my $cores = 1;
my $non_human = 0;
my $build = 37;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'command_path|cp=s' => \$manta_script,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'build=s' => \$build,
    'mopdir|md=s' => \$MOP_dir,
    'threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_Manta.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -ms <configManta.py command path> -n <threads>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR>   Manata configManta.py command path (e.g., ..../manta-1.6.0.centos6_x86_64/bin/configManta.py) [mandatory]
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --build <STR>              reference build (GRCh37, GRCh38) number (37, 38, or T2T) when sample is human [default: 37]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --threads or -n <INT>      number of threads to be used [default: 1]
   --help or -h               output help message
   
=cut

my $chr_bed = '';
if ($MOP_dir ne ''){
  $chr_bed = "$MOP_dir/Data/hs38.bed.gz";
}
else{
  $chr_bed = "$Bin/../../Data/hs38.bed.gz";
}

my $command_config = "$manta_script --bam $input_bam --referenceFasta $ref --runDir ./ 2>$out_prefix.config.log";
$command_config = "$manta_script --bam $input_bam --referenceFasta $ref --runDir ./ --callRegions $chr_bed 2>$out_prefix.config.log" if ($non_human == 0) and ($build eq '38');

my $command_run = "./runWorkflow.py -m local -j $cores -g 100 2>$out_prefix.find.log";

system ("rm -r results") if (-d 'results');
system ("rm -r workspace") if (-d 'workspace');

open (OUT, ">> $out_prefix.command.log");
print OUT "Manta config command: $command_config\n";
close (OUT);

system ("$command_config");

open (OUT, ">> $out_prefix.command.log");
print OUT "Manta run command: $command_run\n";
close (OUT);

system ("$command_run");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_Manta_vcf.pl results/variants/diploidSV.vcf.gz $non_human $target_chr > Manta.$out_prefix.vcf");
}
else{
  system ("$Bin/convert_Manta_vcf.pl results/variants/diploidSV.vcf.gz $non_human $target_chr > Manta.$out_prefix.vcf");
}

print STDERR "Manta run was completed\n";

