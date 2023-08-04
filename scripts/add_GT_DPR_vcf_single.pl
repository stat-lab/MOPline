#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

# running bsub command passed from job_sv_call_lsf_1.pl script

my $vcf = '';
my $sample_name = '';
my $sample_dir = '';
my $out_prefix = '';
my $tool_set = '7tools';
my $min_score = 85;
my $disable_genotype = 0;

my $cov_dir = 'Cov';
my $flank_len = 1000;
my $flank_rate = 0.1;
my $bin_size = 50;
my $non_human = 0;
my $build = 37;
my $ref_index = '';
my $gap_bed = '';

my $data_dir = "$Bin/../Data";
$data_dir = "$ENV{MOPLINE_DIR}/Data" if (exists $ENV{MOPLINE_DIR});

my $help;

GetOptions(
    'vcf|s=s' => \$vcf,
    'sample_name|sn=s' => \$sample_name,
    'sample_dir|sd=s' => \$sample_dir,
    'non_human|nh=i' => \$non_human,
    'build|b=s' => \$build,
    'ref_index|r=s' => \$ref_index,
    'gap_bed|gap=s' => \$gap_bed,
    'toolset|ts=s' => \$tool_set,
    'min_score|ms=i' => \$min_score,
    'disable_genotype|dg' => \$disable_genotype,
    'cov_dname|cd=s' => \$cov_dir,
    'flank|f=i' => \$flank_len,
    'bin=i' => \$bin_size,
    'help|h' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  add_GT_DPR_vcf_single.pl -v <vcf file path> -sn <sample_name> -sd <sample_directory> -ts <tool_set> 
  (-ri <reference_index> -gap <gap_bed> -nh 1 if sample is a non-human species)
  (Add genotype and alignment information (GT/DR/DS/SR tags) to the genotype field of a single 
   MOP-merged vcf file with cov files that were generated with add_GT_DPR_vcf_thread.pl)

  Options:
   --vcf or -v <STR>        vcf file path generated with MOP-based merging strategy [mandatory]
   --sample_name or -sn <STR> sample name/ID [mandatory]
   --sample_dir or -sd <STR>  sample directory path [mandatory]
   --non_human or -nh <INT> sample is non-human species (human: 0. non-human: 1) [default: 0]
   --build or -b <STR>      human reference build (GRCh37, GRCh38, T2T-CHM13) (37, 38, or T2T) [default: 37]
   --ref_index or -r <STR>  reference fasta index file [mandatory if non-human species]
   --gap_bed or -gb <STR>   gap bed file indicating reference gap regions 
                            (1st column: chr-name, 2nd and 3rd columns: start and end positions of a gap) [optional, automatically selected for human]
   --help or -h             output help message

   ## Options for adding SV genotype from SV callers to vcf ##
   --toolset or -ts <STR>    a preset toolset of SV detection tools used (6tools_1, 6tools_2, 7tools, 7tools_2, 8tools, 9tools, or 11tools) 
                             or a list file describing tool names in each line [default: 7tools]
   --min_score or -ms <INT> minimum score to accept SV genotype 
                            (if multiple tools call the same genotype at a site, a summed score from these tools is considered) [default: 85]
   --disable_genotype or -dg <BOOLEAN> disable SV geneotype addition [default: false]

   ## Options for adding DPR/SR alignment information to vcf ##
   --cov_dname or -cd <STR> directory name containing coverage data [default: Cov]
   --flank or -f <INT>      flanking size of SV (DEL and DUP) breakpoints to be analyzed for DPR/SR [default: 1000]
   --bin <INT>              bin size of coverage data (Cov/*.cov file) [default: 50]
   
=cut

die "vcf file is not specified: \n" if ($vcf eq '');
die "sample name is not specified: \n" if ($sample_name eq '');
die "sample directory is not specified: \n" if ($sample_dir eq '');

if (($non_human == 0) and ($ref_index eq '')){
    $ref_index = "$data_dir/hs37.fa.fai" if ($build eq '37');
    $ref_index = "$data_dir/hs38.fa.fai" if ($build eq '38');
    $ref_index = "$data_dir/chm13v2.0.fa.fai" if ($build eq 'T2T');
}
if (($non_human == 0) and ($gap_bed eq '')){
    $gap_bed = "$data_dir/gap.bed" if ($build eq '37');
    $gap_bed = "$data_dir/gap.b38.bed" if ($build eq '38');
}

my $work_dir = `pwd`;
chomp $work_dir;

my $vcf_base = $1 if ($vcf =~ /(.+)\.vcf$/);
$vcf_base = $1 if ($vcf =~ /(.+)\.vcf\.gz$/);
my $vcf2 = "$vcf_base.noAdd.vcf";
if ((-f $vcf2) and (!-z $vcf2)){
    system ("mv -f $vcf2 $vcf");
}
die "$vcf is not found:\n" if (!-f $vcf);

my $cov_dir_path = '';
$cov_dir_path = "$sample_name/$cov_dir" if ($sample_name ne '');
$cov_dir_path = "$sample_dir/$cov_dir" if ($sample_dir ne '');
die "$cov_dir_path directory is not found:\n" if (!-d $cov_dir_path);
my @cov_files = <$cov_dir_path/$sample_name.chr*>;
die "$cov_dir_path/$sample_name.chr\*\.cov.gz files not found:\n" if (@cov_files == 0);

system ("cp $vcf $vcf2");
if ($disable_genotype == 0){
    system ("$Bin/add_genotype_SV_callers_vcf.pl -v $vcf -sn $sample_name -sd $sample_dir -ms $min_score -ts $tool_set > $vcf.gt");
    if ((-f "$vcf.gt") and (!-z "$vcf.gt")){
        system ("mv -f $vcf.gt $vcf");
    }
    else{
        print STDERR "$vcf genotype addition failed: \n";
    }
}

my $dp_command = "$Bin/add_coverage_rate_vcf.pl -v $vcf -cd $cov_dir -sn $sample_name -sd $sample_dir -nh $non_human -b $bin_size -f $flank_len -ri $ref_index";
$dp_command = "$Bin/add_coverage_rate_vcf.pl -v $vcf -cd $cov_dir -sn $sample_name -sd $sample_dir -nh $non_human -b $bin_size -f $flank_len -ri $ref_index -gap $gap_bed" if ($gap_bed ne '');
system ("$dp_command");

if ((-f "$vcf.cov") and (!-z "$vcf.cov")){
    system ("mv -f $vcf.cov $vcf");
}
else{
    print STDERR "$vcf coverage addition failed: \n";
}

print STDERR "$vcf genotype/covergae addition finished:\n";

