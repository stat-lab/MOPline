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
my $wham_path = '';
my $picard_jar = '';
my $samtool_path = '';
my $non_human = 0;
my $cores = 1;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'command_path|cp=s' => \$wham_path,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'samtool_path|sp=s' => \$samtool_path,
    'picard_jar|pj=s' => \$picard_jar,
    'mopdir|md=s' => \$MOP_dir,
    'threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);


=head1 SYNOPSIS

  run_Wham.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -wp <wham command path> -n <threads>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory and contain RG tag [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR> Wham executable path (e.g., ..../wham/bin/whamg) [mandatory]
   --prefix or -p <STR>       outpout prefix [mandatory]
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --samtool_path or sp <STR> samtools path unless samtools is not set in $PATH
   --picard_jar or -pp <STR>  Picard jar path (necessary if input bam file has no RG tag)
   --threads or -n <INT>      number of threads to be used [default: 1]
   --help or -h               output help message
   
=cut

$ENV{PATH} = "$samtool_path:$ENV{PATH}" if ($samtool_path ne '');

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

my $ref_index = "$ref.fai";
my @chr_list;

open (FILE, $ref_index) or die "$ref_index is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr) = split (/\t/, $line);
    if ($target_chr eq 'ALL'){
        if ($non_human == 0){
            push @chr_list, $chr if ($chr =~ /^c*h*r*[\dXY]+$/);
        }
        else{
            push @chr_list, $chr;
        }
    }
    else{
        push @chr_list, $chr if (exists $target_chr{$chr});
    }
}
close (FILE);

my $chr_list = join (',', @chr_list);

my $RG_flag = 0;

open (FILE, "samtools view -H $input_bam |") or die "$input_bam is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^\@RG/){
        $RG_flag = 1;
    }
}
close (FILE);

if ($RG_flag == 0){
    my $out_bam = "$out_prefix.RG.bam";
    my $command_picard = "java -jar -Xmx8g $picard_jar AddOrReplaceReadGroups I=$input_bam O=$out_bam ID=1 SM=$out_prefix LB=PE PL=illumina PU=barcode 2>$out_prefix.picard.log";
    open (OUT, ">> $out_prefix.command.log");
    print OUT "picard add RG command: $command_picard\n";
    close (OUT);
    system ("$command_picard");
    system ("samtools index $out_bam");
    $input_bam = $out_bam;
}

my $command = "$wham_path -f $input_bam -a $ref -x $cores -c $chr_list > $out_prefix.vcf 2> $out_prefix.run.log";

open (OUT, ">> $out_prefix.command.log");
print OUT "whamg command: $command\n";
close (OUT);

system ("$command");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_Wham_vcf.pl $out_prefix.vcf $non_human > Wham.$out_prefix.vcf");
}
else{
    system ("$Bin/convert_Wham_vcf.pl $out_prefix.vcf $non_human > Wham.$out_prefix.vcf");
}

print STDERR "Wham run was completed\n";
