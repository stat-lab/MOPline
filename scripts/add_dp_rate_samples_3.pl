#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);

# script for add_dp_rate subroutine of merge_SV_calls_Missing_calls_mthreads_2.pl

my $sample_list = shift @ARGV;
my $out_file = shift @ARGV;
my $temp_dir = shift @ARGV;
my $min_split_diff = shift @ARGV;

my $cur_dir = `pwd`;
chomp $cur_dir;

system ("mkdir $temp_dir") if (!-d $temp_dir);

my @out_files;

open (FILE, $sample_list) or die "$sample_list is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my ($gr, $id, $pos_file, $cov_dir) = split (/\t/, $line);
    system ("mkdir $temp_dir/$gr") if ($gr ne '') and (!-d "$temp_dir/$gr");
    my $out_dprate = "$temp_dir/$gr/$id.dprate" if ($gr ne '');
    $out_dprate = "$temp_dir/$gr/$id.dprate" if ($gr eq '');
    my $command = "$Bin/add_dp_rate_4.pl $id $pos_file $cov_dir $min_split_diff > $out_dprate";
    system ($command);
    push @out_files, $out_dprate;
    print STDERR "$id dprate calc\n";
}
close (FILE);

open (OUT, "> $out_file");
foreach my $file (@out_files){
    open (FILE, $file) or die "$file is not found:$!\n";
    while (my $line = <FILE>){
        print OUT $line;
    }
    close (FILE);
}
close (OUT);
