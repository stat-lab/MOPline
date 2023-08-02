#!/usr/bin/perl -w
use strict;

# covert trf data output file to UCSC SimpleRepeat format file

my $data_file = shift @ARGV;

my $chr = '';

open (FILE, $data_file) or die "$data_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^Sequence:\s(\S+)/){
        $chr = $1;
        next;
    }
    next if ($line =~ /^#/);
    next if ($line !~ /^\d/);
    my ($start, $end) = split (/\s+/, $line);
    my @line = split (/\s+/, $line);
    my @info = splice (@line, 2, 12);
    my $info = join ("\t", @info);
    print "585\t$chr\t$start\t$end\ttrf\t$info\n";
}
close (FILE);
