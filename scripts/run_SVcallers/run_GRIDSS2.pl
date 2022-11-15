#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);


my $input_bam = '';
my $ref = '';
my $target_chr = 'ALL';
my $gridss_jar = '';
my $gridss_sh = '';
my $black_list = '';
my $java_dir = '';
my $R_dir = '';
my $out_prefix = '';
my $cores = 1;
my $non_human = 0;
my $MOP_dir = '';
my $heap_size = '15g';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'gridss_jar|gj=s' => \$gridss_jar,
    'gridss_sh|gs=s' => \$gridss_sh,
    'black_list|bl=s' => \$black_list,
    'java_path|jp=s' => \$java_dir,
    'r_path|rp=s' => \$R_dir,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_GRIDSS2.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -gj <gridss jar command path> -gs <gridss.sh path> -bl <black list> -n <threads>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>          reference fasta file [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --gridss_jar or -gj <STR>  GRIDSS jar file (e.g., ..../gridss-2.13.1/gridss-2.13.1-gridss-jar-with-dependencies.jar) [mandatory]
   --gridss_sh or -gs <STR>   GRIDSS sh script (e.g., ..../gridss-2.13.1/scripts/gridss.sh) [mandatory]
   --black_list or -bl <STR>  black list bed file (e.g., wgEncodeDacMapabilityConsensusExcludable.b37.bed for human build 37 in GRIDSS package) [optional]
   --java_path or -jp <STR>   directory containing java (1.8 or later) executable, unless the corrsponding version of java is set in $PATH
   --r_path or -rp <STR>      directory containing R (4.0 or later) executable, unless the corrsponding version of R is set in $PATH
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --num_threads or -n <INT>  number of threads to be used [default: 1]
   --help or -h               output help message
   
=cut

$ENV{PATH} = "$java_dir:" . $ENV{PATH} if ($java_dir ne '');
$ENV{PATH} = "$R_dir:" . $ENV{PATH} if ($R_dir ne '');
my $java_opt = '-XX:+UseSerialGC -Xmx' . $heap_size . ' -Xms' . $heap_size;
$ENV{JAVA_TOOL_OPTIONS} = $java_opt;

my $java_pass = `which java`;
chomp $java_pass;
print STDERR "java_pass: $java_pass\n";
my $R_pass = `which R`;
chomp $R_pass;
print STDERR "R_pass: $R_pass\n";

my @working_dir = <./*.working>;
my @log_file = <./gridss.*>;

foreach (@working_dir){
	system ("rm -r $_");
}
foreach (@log_file){
	system ("rm $_");
}

my $ref_base = basename ($ref);
system ("ln -s $ref") if ($ref ne $ref_base);

my $command = "$gridss_sh --steps All --jvmheap $heap_size --otherjvmheap $heap_size --reference $ref_base --output $out_prefix.vcf.gz --assembly $out_prefix.gridss.assembly.bam --jar $gridss_jar --blacklist $black_list --threads $cores $input_bam 2>$out_prefix.gridss.log" if ($gridss_sh !~ /\.sh$/);
$command = "$gridss_sh --steps All --jvmheap $heap_size --otherjvmheap $heap_size --reference $ref_base --output $out_prefix.vcf.gz --assembly $out_prefix.gridss.assembly.bam --jar $gridss_jar --threads $cores $input_bam 2>$out_prefix.gridss.log" if ($black_list eq '') and ($gridss_sh !~ /\.sh$/);
$command = "$gridss_sh --steps All --jvmheap $heap_size --reference $ref_base --output $out_prefix.vcf.gz --assembly $out_prefix.gridss.assembly.bam --jar $gridss_jar --blacklist $black_list --threads $cores $input_bam 2>$out_prefix.gridss.log" if ($gridss_sh =~ /\.sh$/);
$command = "$gridss_sh --steps All --jvmheap $heap_size --reference $ref_base --output $out_prefix.vcf.gz --assembly $out_prefix.gridss.assembly.bam --jar $gridss_jar --threads $cores $input_bam 2>$out_prefix.gridss.log" if ($black_list eq '') and ($gridss_sh =~ /\.sh$/);


open (OUT, ">> $out_prefix.command.log");
print OUT "GRIDSS command: $command\n";
close (OUT);

system ($command);

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_GRIDSS_vcf.pl $out_prefix.vcf.gz $non_human $target_chr > GRIDSS.$out_prefix.vcf");
}
else{
  system ("$Bin/convert_GRIDSS_vcf.pl $out_prefix.vcf.gz $non_human $target_chr > GRIDSS.$out_prefix.vcf");
}

print STDERR "GRIDSS was completed\n";
