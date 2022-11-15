#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use threads;

my $sample_list = '';
my $vcf_dir = 'Merge_7tools';
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

my $cores = 1;
my $help;

GetOptions(
    'sample_list|s=s' => \$sample_list,
    'vcf_dname|vd=s' => \$vcf_dir,
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
    'threads|n=i' => \$cores,
    'help|h' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  add_GT_DPR_vcf.pl -s <sample_list> -vd <vcf_directory> -ts <tool_set> -n <num_threads> 
  (-ri <reference_index> -gap <gap_bed> -nh 1 if sample is a non-human species)
  (Add genotype and alignment information (GT/DR/DS/SR tags) to the genotype field of MOP-merged 
   vcf files with cov files that were generated with add_GT_DPR_vcf_thread.pl, using multiple threads)

  Options:
   --sample_list or -s <STR> sample list file or bam list file. [mandatory]
                            The working directory should contain (a) sample directories (i.e., ${sample_name}), which contain a bam file (.i.e., ${sample_name}.bam)
   --vcf_dname or -vd <STR> vcf directory name containing a MOP-based merged vcf file from multiple tools, 
                            which should be under the sample directory [default: Merge_7tools]
   --non_human or -nh <INT> sample is non-human species (human: 0. non-human: 1) [default: 0]
   --build or -b <STR>      human reference build (GRCh37, GRCh38, T2T-CHM13) (37, 38, or T2T) [default: 37]
   --ref_index or -r <STR>  reference fasta index file [mandatory if non-human species]
   --gap_bed or -gb <STR>   gap bed file indicating reference gap regions 
                            (1st column: chr-name, 2nd and 3rd columns: start and end positions of a gap) [optional, automatically selected for human]
   --threads or -n <INT>    number of threads [default: 1]
   --help or -h             output help message

   ## Options for adding SV genotype from SV callers to vcf ##
   --toolset or -ts <STR>    a preset toolset of SV detection tools used (6tools_1, 6tools_2, 7tools, 9tools, or 11tools) 
                             or a list file describing tool names in each line [default: 7tools]
   --min_score or -ms <INT> minimum score to accept SV genotype (if multiple tools call the same genotype at a site, 
                            a summed score from these tools is considered) [default: 85]
   --disable_genotype or -dg <BOOLEAN> disable SV geneotype addition [default: false]

   ## Options for adding DPR/SR alignment information to vcf ##
   --cov_dname or -cd <STR> directory name containing coverage data [default: Cov]
   --flank or -f <INT>      flanking size of SV (DEL and DUP) breakpoints to be analyzed for DPR/SR [default: 1000]
   --bin <INT>              bin size of coverage data (Cov/*.cov file) [default: 50]
   
=cut

die "sample list file is not specified: \n" if ($sample_list eq '');
die "ref index *.fa.fai is not specified: \n" if ($non_human == 1) and ($ref_index eq '');

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

my $count = 0;
my @jobs;

open (FILE, $sample_list) or die "$sample_list is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\s+/, $line);
    my $ID = basename ($line[0]);
    $ID = $1 if ($ID =~ /(.+?)\./);
    my $merged_vcf = "$ID/$vcf_dir/$ID.Merge.ALL.vcf";
    my $merged_vcf2 = "$ID/$vcf_dir/$ID.Merge.ALL.noAdd.vcf";
    if ((-f $merged_vcf2) and (!-z $merged_vcf2)){
        system ("rm $merged_vcf") if (-f $merged_vcf);
        system ("mv $merged_vcf2 $merged_vcf");
    }
    die "$merged_vcf is not found:\n" if (!-f $merged_vcf);
    my $cov_dir_path = "$ID/$cov_dir";
    die "$cov_dir_path directory is not found:\n" if (!-d $cov_dir_path);
    my @cov_files = <$cov_dir_path/$ID.chr*>;
    die "$cov_dir_path/$ID.chr\*\.cov.gz files not found:\n" if (@cov_files == 0);
    system ("cp $merged_vcf $merged_vcf2");
    $count ++;
    my ($thread_t) = threads->new(\&addGT_cov, $ID, $merged_vcf);
    push @jobs, $thread_t;
    if ($count == $cores){
        foreach (@jobs){
            my ($id2) = $_->join;
            print STDERR "$id2 SVs genotype and coverage info adding completed:\n";
        }
        $count = 0;
        undef @jobs;
    }
}
if (@jobs > 0){
    foreach (@jobs){
        my ($id2) = $_->join;
        print STDERR "$id2 SVs genotype and coverage info adding completed:\n";
    }
    $count = 0;
    undef @jobs;
}
close (FILE);

print STDERR "SV Genotype and Coverage info addition completed for all samples:\n";


sub addGT_cov {
    my ($ID, $merged_vcf) = @_;
    
    my $merged_vcf_base = $1 if ($merged_vcf =~ /(.+)\.vcf$/);
    my $merged_vcf_dir = $1 if ($merged_vcf =~ /(.+)\//);
    if ($disable_genotype == 0){
        system ("$Bin/add_genotype_SV_callers_vcf.pl -v $merged_vcf -sn $ID -ms $min_score -ts $tool_set > $merged_vcf.gt");
        if ((-f "$merged_vcf.gt") and (!-z "$merged_vcf.gt")){
            system ("mv -f $merged_vcf.gt $merged_vcf");
        }
        else{
            print STDERR "$merged_vcf genotype addition failed: \n";
            threads->yield();
            return ("$ID-01");
        }
    }
    
    my $dp_command = "$Bin/add_coverage_rate_vcf.pl -v $merged_vcf -cd $cov_dir -sn $ID -nh $non_human -b $bin_size -f $flank_len -ri $ref_index";
    $dp_command = "$Bin/add_coverage_rate_vcf.pl -v $merged_vcf -cd $cov_dir -sn $ID -nh $non_human -b $bin_size -f $flank_len -ri $ref_index -gap $gap_bed" if ($gap_bed ne '');
    system ("$dp_command");

    if ((-f "$merged_vcf.cov") and (!-z "$merged_vcf.cov")){
        system ("mv -f $merged_vcf.cov $merged_vcf");
    }
    else{
        print STDERR "$merged_vcf coverage addition failed: \n";
        threads->yield();
        return ("$ID-02");
    }
    
    threads->yield();
    sleep 1;
    
    return ("$ID");
}

