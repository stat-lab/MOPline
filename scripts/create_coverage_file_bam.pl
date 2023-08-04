#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use threads;

my $bam_list = '';
my $ref = '';
my $non_human = 0;
my $out_prefix = '';
my $read_length = 0;
my $exclude_discread = 0;
my $min_MAPQ = 1;
my $bin_size = 50;
my $cov_dir = 'Cov';
my $cores = 1;

my $exclude_mit = 1;

my $help;
 
GetOptions(
    'bam_list|b=s' => \$bam_list,
    'ref|r=s' => \$ref,
    'non_human|nh=i' => \$non_human,
    'bin_size|bs=i' => \$bin_size,
    'mapq|q=i' => \$min_MAPQ,
    'ex_disc|xd' => \$exclude_discread,
    'read_len|rl=i' => \$read_length,
    'prefix|p=s' => \$out_prefix,
    'num_threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  create_coverage_file_bam.pl -b <bam_list or a bam file path> -r <ref fasta> -rl <read length> -n <threads> 
  (-nh 1 if sample is a non-human species)
  (create files with coverage and split read information of read alignment at 50-bp window using mutiple threads)
  outout Cov directory will be created in the working directory or under ${sample_name} directory if sample-bam/cram table or sample list file is provided

  Options:
   --bam_list or -b <STR>   sample-bam/cram table file. Table file contains sample directory and bam file name in the 1st and 2nd columns, respectively.
                            or a sample (bam) list file if the working directory contains ${sample_name}/${sample_name}.bam [mandatory]
   --ref or -r <STR>        reference fasta file used for creating bam/cram [mandatory]
   --read_len or -rl <INT>  read length [mandatory]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --bin_size or -bs <INT>  bin size (bp) for calculating coverage [default: 50]
   --mapq or -q <INT>       minimum mapping quality [default: 1]
   --ex_disc or -xd         exclude read pairs with discordant insert size [default: false]
   --prefix or -p <STR>     output prefix
   --num_threads or -n <INT> number of threads to be used [default: 1]
   --help or -h             output help message
   
=cut

my @chr_list;

die "Reference index file: $ref.fai is not found:$!\n" if (!-f "$ref.fai");
die "read length is not specified or too short:\n" if ($read_length < 10);

open (FILE, "$ref.fai") or die "$ref.fai is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $size) = split (/\t/, $line);
    next if ($chr =~ /^c*h*r*Mit$|^c*h*r*Mt$|^chrM$/i) and ($exclude_mit == 1);
    if ($non_human == 0){
        if ($chr =~ /^c*h*r*[\dXY]+$/){
            push @chr_list, $chr;
        }
    }
    else{
        push @chr_list, $chr;
    }
}
close (FILE);

my @bam_files;

if ((-f $bam_list) and ($bam_list =~ /\.bam$|\.cram$/)){
    push @bam_files, $bam_list;
}
else{
    open (FILE, $bam_list) or die "$bam_list is not found:$!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my @line = split (/\s+/, $line);
        my $sample = basename ($line[0]);
        $sample = $1 if ($sample =~ /(.+?)\./);
        my $bam_base = '';
        if ($line[0] =~ /\.bam$|\.cram$/){
            $bam_base = basename ($line[0]);
        }
        else{
            $bam_base = "$sample.bam";
        }
        my $bam = '';
        if (@line == 1){
            if (-d $sample){
                $bam = "$sample/$bam_base";
                $bam = "$sample/$sample.cram" if (!-f $bam);
            }
            else{
                $bam = $bam_base;
                $bam = "$sample.cram" if (!-f $bam);
            }

        }
        elsif (@line >= 2){
            $bam = "$line[0]/$line[1]";
        }
        if (!-f $bam){
            die "BAM/CRAM file for $sample: $bam is not found:\n";
        }
        my $bam_base2 = $1 if ($bam =~ /(.+)\.bam$/);
        if ((!-f "$bam.bai") and (!-f "$bam_base2.bai") and (!-f "$bam.crai") and (!-f "$bam_base2.crai")){
            die "BAM index file: $bam.bai is not found:\n" if ($bam =~ /\.bam$/);
            die "BAM index file: $bam.crai is not found:\n" if ($bam =~ /\.cram$/);
        }
        push @bam_files, $bam;
    }
    close (FILE);
}

foreach my $bam (@bam_files){
    my @jobs;
    my $count = 0;
    my $sample_dir = '.';
    my $sample_name = '';
    $sample_dir = $1 if ($bam =~ /(.+)\//);
    $sample_name = basename ($bam);
    $sample_name = $1 if ($sample_name =~ /(.+?)\./);
    $sample_name = $out_prefix if ($sample_name eq '');
#    system ("rm -r $cov_dir") if (-d $cov_dir) and ($sample_dir eq '');
#    system ("rm -r $sample_dir/$cov_dir") if (-d "$sample_dir/$cov_dir") and ($sample_dir ne '');
    system ("mkdir $cov_dir") if (!-d $cov_dir) and ($sample_dir eq '');
    system ("mkdir $sample_dir/$cov_dir") if (!-d "$sample_dir/$cov_dir") and ($sample_dir ne '');
    foreach my $chr (@chr_list){
        $count ++;
        my ($thread_t) = threads->new(\&cal_cov, $sample_dir, $sample_name, $bam, $chr);
        push @jobs, $thread_t;
        if ($count == $cores){
            foreach (@jobs){
                my ($sample, $chr2) = $_->join;
                print STDERR "$sample: $chr2 finished:\n";
            }
            $count = 0;
            undef @jobs;
        }
    }
    if (@jobs > 0){
        foreach (@jobs){
            my ($sample, $chr2) = $_->join;
            print STDERR "$sample: $chr2 finished:\n";
        }
        $count = 0;
        undef @jobs;
    }
}

print STDERR "Coverage calculation completed:\n";

# record coverage and the number of split reads in a 50-bp window size for a single chromosome
sub cal_cov{
    my ($sample_dir, $sample_name, $bam, $chr) = @_;
    my $command = "$Bin/create_coverage_file_bam_chr.pl -b $bam -c $chr -rl $read_length -r $ref -sn $sample_name -sd $sample_dir -bs $bin_size -rl $read_length -q $min_MAPQ";
    $command = "$Bin/create_coverage_file_bam_chr.pl -b $bam -c $chr -rl $read_length -r $ref -sn $sample_name -sd $sample_dir -bs $bin_size -rl $read_length -q $min_MAPQ -xd" if ($exclude_discread == 1);
    
    open (OUT, ">> $out_prefix.command.log");
    print OUT "$command\n";
    close (OUT);
    
    system ("$command");
    
    threads->yield();
    sleep 1;
    
    return ($sample_name, $chr);
}
