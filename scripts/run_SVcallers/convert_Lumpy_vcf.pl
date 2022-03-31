#!/usr/bin/perl -w
use strict;

# covert Lumpy output files to vcf

my $var_file = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $target_chr = 'ALL';

$target_chr = shift @ARGV if (@ARGV > 0);

my $min_sv_len = 40;

my $min_reads = 2;

my $exclude_mit = 1;

my %vcf;

my $annot1 = 0;
my $annot2 = 0;

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

open (FILE, $var_file) or die "$var_file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
	    next;
	}	
	my @line = split (/\t/, $line);
    if (@line >= 11){
        if ($line[9] =~ /:0:0:0$/){
            $annot1 ++;
        }
        if ($line[10] =~ /:0:0:0$/){
            $annot2 ++;
        }
    }
}
close (FILE);

open (FILE, $var_file);
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
#	    print $line, "\n" if ($count == 1);
	    next;
	}	
	my @line = split (/\t/, $line);
    if (@line >= 11){
        if (($line[9] =~ /:0:0:0$/) and ($annot1 > $annot2)){
            next;
        }
        if (($line[10] =~ /:0:0:0$/) and ($annot1 < $annot2)){
            next;
        }
    }
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
	next if ($type eq 'BND');
	my $len = 0;
	$len = $1 if ($line[7] =~ /SVLEN=-*(\d+);/);
	next if ($len < $min_sv_len) and ($len > 0);
	my $end = $pos + $len - 1;
	my $reads = 0;
	$reads = $1 if ($line[7] =~ /SU=(\d+);/);
	my $gt = '';
	$gt = $1 if ($line[-1] =~ /^(.+?):/);
	my $chr_02d = $chr;
	$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt" if ($gt ne '');
	${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$reads" if ($gt eq '');
}
close (FILE);

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
		    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
		}
    }
}
