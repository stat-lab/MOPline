#!/usr/bin/perl -w
use strict;

# covert Delly output files to vcf

my $bcf_files = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $target_chr = 'ALL';

$target_chr = shift @ARGV if (@ARGV > 0);

my $bcftools = 'bcftools';

$bcftools = shift @ARGV if (@ARGV > 0);

my $min_sv_len = 40;
my $max_sv_len = 2000000;

my $qual_filter = 0;

my $exclude_mit = 1;

my @bcf_files = split (/,/, $bcf_files);

my %vcf;
my $count = 0;

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

foreach my $file (@bcf_files){
    if (!-e $file){
		print STDERR "$file is not found\n";
		next;
    }
    $count ++;
    open (FILE, $file) if ($file =~ /\.vcf$/);
    open (FILE, "$bcftools view $file |") if ($file =~ /\.bcf$/);
    while (my $line = <FILE>){
		chomp $line;
		if ($line =~ /^#/){
		    print $line, "\n" if ($count == 1);
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
		my $end = $1 if ($line[7] =~ /END=(\d+);/);
		my $reads = $1 if ($line[7] =~ /PE=(\d+);/);
		my $len = $end - $pos + 1;
		next if ($len < $min_sv_len);
		my $qual = $line[6];
		my $gt = $1 if ($line[9] =~ /^(.+?):/);
		if ($qual_filter == 1){
		    next if ($qual ne 'PASS');
		}
		my $chr2 = '';
		my $pos2 = 0;
		if ($type eq 'TRA'){
		    $chr2 =  $1 if ($line[7] =~ /CHR2=(.+?);/);
		    $pos2 = $end;
		    $end = 0;
		    $len = 0;
		}
		my $chr_02d = $chr;
		$chr_02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
		${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\t$qual\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt" if ($type ne 'TRA');
		${${$vcf{$chr_02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\t$qual\tSVTYPE=$type;SVLEN=$len;READS=$reads;GT=$gt\tCHR2=$chr2;POS2=$pos2" if ($type eq 'TRA');
    }
    close (FILE);
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
		    print ${${$vcf{$chr}}{$pos}}{$type}, "\n";
		}
    }
}
