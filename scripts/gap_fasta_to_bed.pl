#!/usr/bin/perl -w
use strict;

my $genome_fasta = shift (@ARGV);

my $min_gap = 1000;

$min_gap = shift @ARGV if (@ARGV > 0);

print STDERR "Minimum Gap length: $min_gap\n";

my %gap;
my $seq = '';
my $chr = '';

open (FILE, $genome_fasta) or die "$genome_fasta is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        if ($seq ne ''){
            my $pos = 0;
            while ($seq =~ /(.*?)([Nn]{$min_gap,})/g){
                my $len = length $1;
                my $glen = length $2;
                $pos += $len;
                my $end = $pos + $glen;
                ${$gap{$chr}}{$pos} = $end;
                $pos += $glen;
            }
            $seq = '';
        }
        $chr = $1;
    }
    else{
        $seq .= $line;
    }
}
if ($seq ne ''){
    my $pos = 0;
    while ($seq =~ /(.*?)([Nn]{$min_gap,})/g){
        my $len = length $1;
        my $glen = length $2;
        $pos += $len;
        my $end = $pos + $glen;
        ${$gap{$chr}}{$pos} = $end;
        $pos += $glen;
    }
    $seq = '';
}
close (FILE);

my $total_gap = 0;
my $total_gaplen = 0;

foreach my $chr (sort keys %gap){
    foreach my $pos (sort {$a <=> $b} keys %{$gap{$chr}}){
        my $end = ${$gap{$chr}}{$pos};
        my $len = $end - $pos;
        $total_gap ++;
        $total_gaplen += $len;
        print "$chr\t$pos\t$end\t$len\n";
    }
}

print STDERR "Total Gap sites (>= $min_gap bp): $total_gap\n";
print STDERR "Total Gap length (>= $min_gap bp): $total_gaplen\n";
