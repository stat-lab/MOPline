#!/usr/bin/perl -w
use strict;

# stat INS-read number in inGAP vcf

my $vcf_file = shift @ARGV;

my $num1 = shift @ARGV;

my $num2 = shift @ARGV;

my %ins;

my $total_ins = 0;

open (FILE, $vcf_file) or die "$vcf_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    if ($line =~ /SVTYPE=INS/){
        my $read = $1 if ($line =~ /READS=(\d+)/);
        $read = 100 if ($read > 100);
        for (my $i = 1; $i <= $read; $i ++){
            $ins{$i} ++;
        }
        $total_ins ++;
    }
}
close (FILE);

my $mread1 = 0;
my $mread2 = 0;

foreach my $read (sort {$a <=> $b} keys %ins){
    my $num = $ins{$read};
#print STDERR "$read $num\n";
    if (($num <= $num1) and ($mread1 == 0)){
        $mread1 = $read;
    }
    if (($num <= $num2) and ($mread2 == 0)){
        $mread2 = $read;
    }
}

print "$mread1,$mread2\n";
