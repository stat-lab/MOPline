#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use threads;

# running bsub command passed from job_sv_call_lsf_1.pl script

my $sample_list = '';
my $GroupID = '';
my $out_prefix = '';
my $non_human = 0;
my $build = 37;
my $tool_set = '7tools';
my $merge_dir = '';
my $cores = 1;

my $help;

GetOptions(
    'sample|s=s' => \$sample_list,
    'group|g=s' => \$GroupID,
    'toolset|ts=s' => \$tool_set,
    'input_dir|d=s' => \$merge_dir,
    'nonhuman|nh=i' => \$non_human,
    'build=i' => \$build,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  add_coverage_rate_2vcf.pl

  Options:
   --sample or -s <STR>     sample list file or a single sample name
   --group or -g <STR>      group ID if group ID is contained in the input vcf file name (i.e., ${groupID}.${sample_name}.Merge.ALL.vcf)
   --toolset or -ts <STR>   preset tool set used for merging/selecting of SV calls from multiple SV call sets 
                            (7tools|6tools_1|6tools_2|9tools|11tools) [default: 7tools]
   --input_dir or -d <STR>  input directory name containing input merged/slected vcf files [default: Merge_ + $toolset, e.g., Merge_7tools]
   --non_human or -nh <INT> samples are non-human species (human: 0, non-human: 1) [default: 0]
   --build <INT>            human reference build build 37 or build 38 (37|38) [default: 37] (adapted for only human)
   --prefix or -p <STR>     output prefix
   --help or -h             output help message
   
=cut

$merge_dir = 'Merge_' . $tool_set if ($merge_dir eq '');

my $work_dir = `pwd`;
chomp $work_dir;

my @sample_list;

if (-f $sample_list){
    open (FILE, $sample_list) or die "$sample_list is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my $sample = $1 if ($line =~ /(\S+)/);
        push @sample_list, $sample;
    }
    close (FILE);
}
else{
    push @sample_list, $sample_list;
}

my $count = 0;
my @jobs;

foreach my $ID (@sample_list){
    my $vcf = "$ID/$merge_dir/$ID.Merge.ALL.vcf";
    $vcf = "$ID/$merge_dir/$GroupID.$ID.Merge.ALL.vcf" if ($GroupID ne '');
    if (!-f $vcf){
        print STDERR "$vcf merged vcf file is not found: \n";
        print STDERR "MOPline expects the above input vcf file in sample_name/input_dir/\n";
        next
    }
    $count ++;
    my ($thread_t) = threads->new(\&addGT_cov, $ID, $vcf);
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

print STDERR "$GroupID $out_prefix SV Genotype and Coverage info addition completed:\n";


sub addGT_cov {
    my ($ID, $vcf) = @_;
    my $vcf_base = $1 if ($vcf =~ /(.+)\.vcf$/);
    my $vcf_dir = $1 if ($vcf =~ /(.+)\//);
    my $flag = 0;
    if ((-f "$vcf_base.noGT.vcf") and (!-z "$vcf_base.noGT.vcf")){
        system ("rm $vcf") if (-f $vcf);
        system ("mv -f $vcf_base.noGT.vcf $vcf");
        $flag ++;
    }
    if ($tool_set eq ''){
        system ("add_genotype_SV_callers_vcf.pl -v $vcf -s $ID > $vcf.gt");
    }
    else{
        system ("add_genotype_SV_callers_vcf.pl -v $vcf -s $ID -ts $tool_set > $vcf.gt");
    }
    my $vcf_base2 = basename ($vcf_base);
    my $vcf2 = "$vcf_base.noGT.vcf";
    if ((-f "$vcf.gt") and (!-z "$vcf.gt")){
        system ("mv -f $vcf $vcf2");
        system ("mv $vcf.gt $vcf");
    }
    else{
        print STDERR "$vcf.gt file is not found: \n";
        threads->yield();
        return ("$ID-01");
    }
    
    my $cov_file = "$ID/Cov/$ID";
    $cov_file = "Cov/$ID" if (!-f "$cov_file.chr1.cov.gz");
    system ("add_coverage_rate_vcf.pl -v $vcf -c $cov_file --build $build");
    my $vcf_out = $vcf;
    if ((-f "$vcf.cov") and (!-z "$vcf.cov")){
        system ("mv -f $vcf.cov $vcf_out");
    }
    else{
        print STDERR "$vcf.cov file is not found: \n";
        threads->yield();
        return ("$ID-02");
    }
    
    threads->yield();
    sleep 1;
    
    return ("$ID");
}

