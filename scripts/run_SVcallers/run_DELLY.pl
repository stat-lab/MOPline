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
my $delly_path = '';
my $exclude_list_file = '';
my $bcftools_path = '';
my $non_human = 0;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'command_path|cp=s' => \$delly_path,
    'exclude|ex=s' => \$exclude_list_file,
    'bcftools_path|bp=s' => \$bcftools_path,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);


=head1 SYNOPSIS

  run_DELLY.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -cp <DELLY command path> -n <threads>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory and contain RG tag [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR>  delly executable path (e.g., ..../delly-0.8.7/delly) [mandatory]
   --exclude or -ex <STR>     tsv file indicating exclude regions (e.g., ..../delly-0.8.7/excludeTemplates/human.hg19.excl.tsv)
   --bcftools_path or -bp <STR> bcftools path unless bcftools is set in $PATH
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --help or -h               output help message
   
=cut

$ENV{PATH} = "$bcftools_path:$ENV{PATH}" if ($bcftools_path ne '');

my $command_del = "$delly_path call -g $ref -o $out_prefix.del.bcf -t DEL $input_bam -x $exclude_list_file 2>$out_prefix.Delly.DEL.log";
$command_del = "$delly_path call -g $ref -o $out_prefix.del.bcf -t DEL $input_bam 2>$out_prefix.Delly.DEL.log" if ($exclude_list_file eq '');
my $command_dup = "$delly_path call -g $ref -o $out_prefix.dup.bcf -t DUP $input_bam -x $exclude_list_file 2>$out_prefix.Delly.DUP.log";
$command_dup = "$delly_path call -g $ref -o $out_prefix.dup.bcf -t DUP $input_bam 2>$out_prefix.Delly.DUP.log" if ($exclude_list_file eq '');
my $command_inv = "$delly_path call -g $ref -o $out_prefix.inv.bcf -t INV $input_bam -x $exclude_list_file 2>$out_prefix.Delly.INV.log";
$command_inv = "$delly_path call -g $ref -o $out_prefix.inv.bcf -t INV $input_bam 2>$out_prefix.Delly.INV.log" if ($exclude_list_file eq '');
#my $command_tra = "~/tools/SV_callers/delly call -g $ref -o $out_prefix.tra.vcf -t TRA $input_bam -x $exclude_list 2>$out_prefix.Delly.TRA.log";

my @out_files;

open (OUT, ">> $out_prefix.command.log");
print OUT "DEL command: $command_del\n";
close (OUT);

system ("$command_del");
push @out_files, $out_prefix . '.del.bcf';

sleep 300;

open (OUT, ">> $out_prefix.command.log");
print OUT "DUP command: $command_dup\n";
close (OUT);

system ("$command_dup");
push @out_files, $out_prefix . '.dup.bcf';

sleep 300;

open (OUT, ">> $out_prefix.command.log");
print OUT "INV command: $command_inv\n";
close (OUT);

system ("$command_inv");
push @out_files, $out_prefix . '.inv.bcf';

sleep 300;

my $out_str = join (',', @out_files);

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_DELLY_vcf.pl $out_str $non_human $target_chr > DELLY.$out_prefix.vcf") if ($bcftools_path eq '');
    system ("$MOP_dir/scripts/run_SVcallers/convert_DELLY_vcf.pl $out_str $non_human $target_chr $bcftools_path > DELLY.$out_prefix.vcf") if ($bcftools_path ne '');
}
else{
  system ("$Bin/convert_DELLY_vcf.pl $out_str $non_human $target_chr > DELLY.$out_prefix.vcf") if ($bcftools_path eq '');
  system ("$Bin/convert_DELLY_vcf.pl $out_str $non_human $target_chr $bcftools_path > DELLY.$out_prefix.vcf") if ($bcftools_path ne '');
}

print STDERR "DELLY run was completed\n";
