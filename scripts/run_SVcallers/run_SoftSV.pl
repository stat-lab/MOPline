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
my $softsv_path = '';
my $bamtools_lib = '';
my $non_human = 0;
my $MOP_dir = '';

my $cores = 1;

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c' => \$target_chr,
    'command_path|cp=s' => \$softsv_path,
    'bamtools_lib|bl=s' => \$bamtools_lib,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);


=head1 SYNOPSIS

  run_SoftSV.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -cp <SoftSV command path> -bl <bamtools lib directory>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory and contain RG tag [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR> SoftSV executable path (e.g., ..../SoftSV_1.4.2/SoftSV) [mandatory]
   --bamtools_lib or -bl <STR> bamtools lib path unless bamtools lib directory is set in $LD_LIBRARY_PATH
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --help or -h               output help message
   
=cut

my $input_base = basename ($input_bam);
my $input_base2 = $1 if ($input_base =~ /(.+)\.bam$/);

$ENV{LD_LIBRARY_PATH} = "$bamtools_lib:$ENV{LD_LIBRARY_PATH}";

#$ENV{LD_LIBRARY_PATH} = '/opt/local/python/2.7/lib:/home/gluster/sg/kosugi/tools/boost_1_49_0/stage/lib:/home/gluster/sg/kosugi/tools/lib/openmpi:/home/gluster/sg/kosugi/tools/gmp-6.0.0/lib:/home/gluster/sg/kosugi/tools/mpfr-3.1.2/lib:/home/gluster/sg/kosugi/tools/isl-0111/lib:/home/gluster/sg/kosugi/tools/mpc-1.0.3/lib:/home/gluster/sg/kosugi/tools/ppl-0.11/lib:/opt/local/samtools/0.1.19/lib:/opt/local/py-tcdb/0.4/lib:/opt/local/pyfasta/0.5.2/lib:/opt/local/pash/3.0.6.2/lib:/opt/local/methylcoder/59ac293/lib:/opt/local/mcr/v79/runtime/glnxa64:/opt/local/bsmap/2.74/lib:/opt/local/brat_bw/2.0.1/lib:/opt/local/brat/1.2.4/lib:/opt/local/BS_SEEKER/1.0/lib:/opt/local/BMap/1.0/lib:/opt/local/zsh/5.0.5/lib:/opt/local/xz/5.0.5/lib:/opt/local/sqlite/3.8.2/lib:/opt/local/pgplot/5.2:/opt/local/libevent/2.0.21/lib:/opt/local/lapack-gcc/3.4.2/lib:/opt/local/jdk/6/lib:/opt/local/gsl-gcc/1.16/lib:/opt/local/cmake/2.8.12.1/lib:/opt/local/boost-gcc/1.54.0/lib:/opt/local/blas-gcc/20110419/lib:/opt/local/r/3.0.2/lib64/R/lib:/opt/local/atlas-gcc/3.10.1/lib:/usr/share/lava/1.0/linux2.6-glibc2.12-x86_64/lib:/opt/intel/impi/4.1.3.048/intel64/lib:/opt/intel/composer_xe_2013.5.192/compiler/lib/intel64:/opt/intel/composer_xe_2013.5.192/mpirt/lib/intel64:/opt/intel/composer_xe_2013.5.192/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2013.5.192/ipp/lib/intel64:/opt/intel/composer_xe_2013.5.192/compiler/lib/intel64:/opt/intel/composer_xe_2013.5.192/mkl/lib/intel64:/opt/intel/composer_xe_2013.5.192/tbb/lib/intel64/gcc4.4:/home/gluster/sg/kosugi/tools/MCR_R2013b_glnxa64_installer/v82/runtime/glnxa64:/home/gluster/sg/kosugi/tools/MCR_R2013b_glnxa64_installer/v82/bin/glnxa64:/home/gluster/sg/kosugi/tools/MCR_R2013b_glnxa64_installer/v82/sys/os/glnxa64:/home/gluster/sg/kosugi/tools/re2-master/usr/local/lib:/home/gluster/sg/kosugi/tools/samtools-0.1.19:/home/gluster/sg/kosugi/tools/bamtools:/home/gluster/sg/kosugi/tools/lib:/home/gluster/sg/kosugi/tools/hdf5-1.8.15/lib:/home/gluster/sg/kosugi/tools/bamtools-master/lib';


my $command = "$softsv_path -i $input_bam 2>$out_prefix.SoftSV.log";

open (OUT, ">> $out_prefix.command.log");
print OUT "SoftSV command: $command\n";
close (OUT);

system ("$command");

my @out;

push @out, 'deletions_small.txt' if (-f 'deletions_small.txt');
push @out, 'deletions.txt' if (-f 'deletions.txt');
push @out, 'insertion_small.txt' if (-f 'insersion_small.txt');
push @out, 'insersion.txt' if (-f 'insersion.txt');
push @out, 'insersions.txt' if (-f 'insersions.txt');
push @out, 'inversions_small.txt' if (-f 'inversions_small.txt');
push @out, 'inversions.txt' if (-f 'inversions.txt');
push @out, 'tandems_small.txt' if (-f 'tandems_small.txt');
push @out, 'tandems.txt' if (-f 'tandems.txt');
#push @out, 'translocations_inverted.txt' if (-f 'translocations_inverted.txt');
#push @out, 'translocations.txt' if (-f 'translocations.txt');

my $file_str = join (',', @out);

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_SoftSV_vcf.pl $file_str $non_human $target_chr > SoftSV.$out_prefix.vcf");
}
else{
  system ("$Bin/convert_SoftSV_vcf.pl $file_str $non_human $target_chr > SoftSV.$out_prefix.vcf");
}


print STDERR "SoftSV run was completed\n";
