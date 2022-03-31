#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $bam = '';

my $sample_name = '';
my $sample_dir = '';

my $target_chr = '';

my $ref = '';

my $bin_size = 50;
my $Mbin_size = 1000000;

my $read_length = 0;

my $exclude_discread = 0;

my $min_MAPQ = 1;

my $cov_dir = 'Cov';

my $help;
 
GetOptions(
    'bam|b=s' => \$bam,
    'sample_name|sn=s' => \$sample_name,
    'sample_dir|sd=s' => \$sample_dir,
    'chr|c=s' => \$target_chr,
    'ref|r=s' => \$ref,
    'bin_size|bs=i' => \$bin_size,
    'mapq|q=i' => \$min_MAPQ,
    'ex_disc|xd' => \$exclude_discread,
    'read_len|rl=i' => \$read_length,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  create_coverage_file_bam_chr.pl -b <bam/cram file> -c <chr name> -sn <sample name> -sd <sample dir>
  (create a file with coverage and split read information of read alignment at 50-bp window for a single chromosome of a single bam)

  Options:
   --bam or -b <STR>        bam/cram file path [mandatory]
   --sample_name or -sn <STR>  sample name [mandatory]
   --sample_dir or -sd <STR>   sample directory (optional, Cov directory is created in the working directory if not specified)
   --chr or -c              target chromosome [mandatory]
   --bin_size or -bs <INT>  bin size (bp) for calculating coverage [default: 50]
   --mapq or -q <INT>       minimum mapping quality [default: 1]
   --read_len or -rl <INT>  read length [default: 0]
   --ex_disc or -xd         exclude read pairs with discordant insert size [default: false]
   --help or -h             output help message
   
=cut

my $max_bin = $Mbin_size / $bin_size;

system ("mkdir $cov_dir") if (!-d $cov_dir) and ($sample_dir eq '');
system ("mkdir $sample_dir/$cov_dir") if (!-d "$sample_dir/$cov_dir") and ($sample_dir ne '');

my $out_chr = "$sample_name.$target_chr.cov";
$out_chr = "$sample_name.chr$target_chr.cov" if ($target_chr =~ /^\d+/);

if ($sample_dir eq ''){
    $out_chr = "$cov_dir/$out_chr";
}
else{
    $out_chr = "$sample_dir/$cov_dir/$out_chr";
}

my $Mbin = 0;
my %DP2;
my %DP2_2;
my %split5;
my %split3;
open (OUT, "> $out_chr");
open (FILE, "samtools view $bam |") or die "$bam is not found: $!\n" if ($bam =~ /\.bam$/) and ($bam =~ /\.chr/);
open (FILE, "samtools view $bam $target_chr |") or die "$bam is not found: $!\n" if ($bam =~ /\.bam$/) and ($bam !~ /\.chr/);
open (FILE, "samtools view -T $ref $bam |") or die "$bam is not found: $!\n" if ($bam =~ /\.cram$/) and ($bam =~ /\.chr/);
open (FILE, "samtools view -T $ref $bam $target_chr |") or die "$bam is not found: $!\n" if ($bam =~ /\.cram$/) and ($bam !~ /\.chr/);
open (FILE, "samtools view -S $bam $target_chr |") or die "$bam is not found: $!\n" if ($bam =~ /\.sam$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^\@/);
    my @line = split (/\s+/, $line);
    my $chr2 = $line[2];
    my $cigar = $line[5];
    next if ($cigar eq '*');
    next if ($cigar =~ /H/);
    my $mapq = $line[4];
    my $pos = $line[3];
    my $len = $read_length if ($read_length > 0);
    $len = length $line[9] if ($read_length == 0);
    if ($cigar =~ /(\d+)S$/){
        $len -= $1;
    }
    if ($cigar =~ /^(\d+)S/){
        $len -= $1;
    }
    my $Mpos_div = int ($pos / $Mbin_size);
    my $Mpos = $pos % $Mbin_size;
    my $Bpos_div = int ($Mpos / $bin_size);
    my $Bpos = $Mpos % $bin_size;
    my $count = 0;
#print STDERR "$chr:$pos\t$Mpos_div-$Mpos\t$Bpos_div-$Bpos\n";
    while (1){
        my $Bpos_e = $bin_size - $Bpos;
        $count ++;
        if ($Bpos_e < $len){
            ${$DP2{$Mpos_div}}{$Bpos_div} += $Bpos_e;
            ${$DP2_2{$Mpos_div}}{$Bpos_div} += $Bpos_e if ($mapq >= $min_MAPQ);
            ${$split5{$Mpos_div}}{$Bpos_div} ++ if ($count == 1) and ($cigar =~ /^\d+S/);
            $len -= $Bpos_e;
            $Bpos = 0;
            $Bpos_div ++;
            if ($Bpos_div >= $max_bin){
                $Mpos_div ++;
            }
        }
        else{
            ${$DP2{$Mpos_div}}{$Bpos_div} += $len;
            ${$DP2_2{$Mpos_div}}{$Bpos_div} += $len if ($mapq >= $min_MAPQ);
            ${$split3{$Mpos_div}}{$Bpos_div} ++ if ($cigar =~ /\d+S$/);
            last;
        }
    }
    if (int ($pos / $Mbin_size) > $Mbin){
        foreach my $bin (sort {$a <=> $b} keys %{$DP2{$Mbin}}){
            my $pos2 = $Mbin * $Mbin_size + $bin * $bin_size;
            my $cov = int (${$DP2{$Mbin}}{$bin} / $bin_size * 10) / 10;
            my $cov2 = 0;
            $cov2 = int (${$DP2_2{$Mbin}}{$bin} / $bin_size * 10) / 10 if (exists ${$DP2_2{$Mbin}}{$bin});
            my $split5 = 0;
            my $split3 = 0;
            $split5 = ${$split5{$Mbin}}{$bin} if (exists ${$split5{$Mbin}}{$bin});
            $split3 = ${$split3{$Mbin}}{$bin} if (exists ${$split3{$Mbin}}{$bin});
            if (int ($pos2 / $Mbin_size) == $Mbin){
                print OUT "$target_chr\t$pos2\t$cov\t$cov2\t$split5\t$split3\n";
            }
            else{
                my $Mpos_div2 = int ($pos2 / $Mbin_size);
                my $Mpos2 = $pos2 % $Mbin_size;
                my $Bpos_div2 = int ($Mpos2 / $bin_size);
                ${$DP2{$Mpos_div2}}{$Bpos_div2} += ${$DP2{$Mbin}}{$bin};
                ${$DP2_2{$Mpos_div2}}{$Bpos_div2} += ${$DP2_2{$Mbin}}{$bin} if (exists ${$DP2_2{$Mbin}}{$bin});
                ${$split5{$Mpos_div2}}{$Bpos_div2} += $split5;
                ${$split3{$Mpos_div2}}{$Bpos_div2} += $split3; 
            }
        }
        delete $DP2{$Mbin};
        delete $DP2_2{$Mbin};
        delete $split5{$Mbin};
        delete $split3{$Mbin};
        $Mbin ++;
    }
}
if (scalar keys %DP2 > 0){
    foreach my $mbin (sort {$a <=> $b} keys %DP2){
        foreach my $bin (sort {$a <=> $b} keys %{$DP2{$mbin}}){
            my $pos2 = $mbin * $Mbin_size + $bin * $bin_size;
            my $cov = int (${$DP2{$mbin}}{$bin} / $bin_size * 10) / 10;
            my $cov2 = 0;
            $cov2 = int (${$DP2_2{$mbin}}{$bin} / $bin_size * 10) / 10 if (exists ${$DP2_2{$mbin}}{$bin});
            my $split5 = 0;
            my $split3 = 0;
            $split5 = ${$split5{$mbin}}{$bin} if (exists ${$split5{$mbin}}{$bin});
            $split3 = ${$split3{$mbin}}{$bin} if (exists ${$split3{$mbin}}{$bin});
            print OUT "$target_chr\t$pos2\t$cov\t$cov2\t$split5\t$split3\n";
        }
    }
}
close (FILE);
close (OUT);

system ("bgzip -f $out_chr");
#system ("tabix -f -b 2 -e 2 $out_chr.gz");
