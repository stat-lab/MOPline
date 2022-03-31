#!/usr/bin/perl -w
use strict;

# covert Wham output files to vcf

my $var_file = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $min_del_size = 40;
my $min_sv_size = 40;
my $max_sv_size = 2000000;

my $exclude_mit = 1;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
        print "$line\n";
        next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    next if ($chr =~ /^c*h*r*Mit$|^c*h*r*Mt$|^chrM$/i) and ($exclude_mit == 1);
    next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXYMT]+$/);
    my $pos = $line[1];
    my $type = '';
    $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type eq 'BND');
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($len > $max_sv_size);
    next if ($type eq 'DEL') and ($len < $min_del_size);
    next if ($type ne 'DEL') and ($len < $min_sv_size);
    my $read = 0;
    $read = $1 if ($line[7] =~ /A=(\d+);/);
    next if ($read < 3);
    my $cw = $1 if ($line[7] =~ /CW=(.+?);/);
    my @cw = split (/,/, $cw);
    next if ($cw[4] > 0.2);
    if ($type eq 'DEL'){
        next if ($cw[0] < 0.2);
    }
    elsif ($type eq 'DUP'){
        next if ($cw[1] < 0.2);
    }
    elsif ($type eq 'INV'){
        next if ($cw[2] < 0.2);
    }
    elsif ($type eq 'INS'){
        next if ($cw[3] < 0.2);
    }
    print "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read\n";
}
close (FILE);
