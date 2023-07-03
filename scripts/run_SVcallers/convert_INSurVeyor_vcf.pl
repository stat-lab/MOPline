#!/usr/bin/perl -w
use strict;

# covert INSurVeylor vcf

my $vcf = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $min_len = 40;

open (FILE, $vcf) or die "$vcf is not found: $!\n" if ($vcf !~ /\.gz$/);
open (FILE, "gzip -dc $vcf |") or die "$vcf is not found: $!\n" if ($vcf =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    next if ($line[3] =~ /NB=Y/);
    my $chr = $line[0];
    next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXYMT]+$/);
    my $pos = $line[1];
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
    next if ($len > 0) and ($len < $min_len);
    my $read = 0;
    if ($line[7] =~ /SPANNING_READS=(\d*),(\d+)/){
        my $read1 = $1;
        my $read2 = $2;
        $read = int (($read1 + $read2) * 0.5 + 0.5);
    }
    elsif ($line[7] =~ /SPANNING_READS=(\d*)/){
        $read = $1;
    }
    $line[7] = "SVTYPE=INS;SVLEN=$len;READS=$read";
    splice (@line, 8, 2);
    print join ("\t", @line), "\n";
}
close (FILE);