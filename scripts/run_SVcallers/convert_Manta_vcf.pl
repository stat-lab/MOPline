#!/usr/bin/perl -w
use strict;

# covert Manta output files to vcf

my $var_file = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $target_chr = 'ALL';

$target_chr = shift @ARGV if (@ARGV > 0);

#my $min_sv_len = 10;
my $min_sv_len = 40;
my $max_del_len = 20000000;
my $max_dup_len = 20000000;

my $exclude_mit = 1;

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

open (FILE, "gzip -dc $var_file |") or die "$var_file is not found: $!\n" if ($var_file =~ /\.gz$/);
open (FILE, $var_file) or die "$var_file is not found: $!\n" if ($var_file !~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
    	print "$line\n";
    	next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    if ($target_chr eq 'ALL'){
        next if ($chr =~ /^c*h*r*Mit$|^c*h*r*Mt$|^chrM$/i) and ($exclude_mit == 1);
        next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXYMT]+$/);
    }
    else{
        next if (!exists $target_chr{$chr});
    }
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($type ne 'DEL') and ($type ne 'DUP') and ($type ne 'INS') and ($type ne 'INV');
    my $len = 0;
    $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    next if ($len < $min_sv_len) and ($type ne 'INS');
    next if ($len > $max_del_len) and ($type eq 'DEL');
    next if ($len > $max_dup_len) and ($type eq 'DUP');
    my $reads = 0;
    my $gt = $1 if ($line[9] =~ /^(.+?):/);
    my @info = split (/:/, $line[9]);
    my ($ref_r1, $alt_r1) = split (/,/, $info[-1]);
    if ($alt_r1 == 0){
    	my ($ref_r2, $alt_r2) = split (/,/, $info[-2]);
    	$reads = $alt_r2;
    }
    else{
	   $reads = $alt_r1;
    }
    next if ($line[6] ne 'PASS');
    $line[2] = $type;
    $line[7] = "SVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt";
    splice (@line, 8, 2);
    $line = join ("\t", @line);
    print "$line\n";
}
close (FILE);
