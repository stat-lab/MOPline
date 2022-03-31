#!/usr/bin/perl -w
use strict;

# covert GRIDSS output vcf to vcf

my $var_file = shift @ARGV;

my $non_human = 0;

$non_human = shift @ARGV if (@ARGV > 0);

my $target_chr = 'ALL';

$target_chr = shift @ARGV if (@ARGV > 0);

my $min_ins_len = 30;
my $min_sv_len = 40;
my $max_sv_len = 20000000;

my $exclude_mit = 1;

my %used;

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

open (FILE, $var_file) or die "$var_file is not found: $!\n" if ($var_file !~ /\.gz$/);
open (FILE, "gzip -dc $var_file |") or die "$var_file is not found: $!\n" if ($var_file =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    if ($target_chr eq 'ALL'){
        next if ($chr =~ /^c*h*r*Mit$|^c*h*r*Mt$|^chrM$/i) and ($exclude_mit == 1);
        next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXYMT]+$/);
    }
    else{
        next if (!exists $target_chr{$chr});
    }
    my $pos = $line[1];
    my $qual = $line[5];
    my $alt = $line[4];
    my $ref = $line[3];
    my $chr_pos = "$chr:$pos";
    my $pass = $line[6];
    next if ($pass eq 'LOW_QUAL;NO_ASSEMBLY');
    next if (exists $used{$chr_pos});
    my $type = '';
    my $len = 0;
    if ($alt =~ /([\[\]])(.+)([\[\]])/){
        my $blnk1 = $1;
        my $blnk2 = $3;
        my $chrpos = $2;
        my ($chr2, $pos2) = split (/:/, $chrpos);
        next if ($chr ne $chr2);
        if (($blnk1 eq '[') and ($blnk2 eq '[')){
            if ($alt =~ /\[$ref$/){
                $type = 'INV';
                $len = $pos2 - $pos - 1;
            }
            elsif ($alt =~ /^([ACGTN)]+)/){
                if (length $1 == 1){
                    $type = 'DEL';
                    $len = $pos2 - $pos - 1;
                }
                else{
                    $len = length ($1) - 2;
                    my $distance = $pos2 - $pos;
                    if ($distance < 0){
                        $type = 'DEL';
                        $len += $distance;
                    }
                    elsif ($distance < 30){
                        $len -= $distance;
                        $type = 'INS';
                    }
                    else{
                        $type = 'DEL';
                        $len = $distance - $len;
                    }
                }
            }
            elsif ($alt =~ /([ACGTN)]{2,})$/){
                $type = 'INV';
                $len = $pos2 - $pos - 1;
                $len += length ($1) - 2;
            }
        }
        elsif (($blnk1 eq ']') and ($blnk2 eq ']')){
            if ($alt =~ /\]([ACGTN]+)$/){
                $len = length ($1) - 1;
                my $distance = $pos2 - $pos;
                if ($distance < 0){
                    $type = 'DEL';
                    $len += $distance;
                }
                elsif ($distance < 30){
                    $len -= $distance;
                    $type = 'INS';
                }
                else{
                    $type = 'DUP';
                    $len += $distance;
                }
            }
            elsif ($alt =~ /\]$/){
                $type = 'INV';
                $len = $pos2 - $pos - 1;
            }
        }
        $used{$chrpos} = 1;
        $used{$chr_pos} = 1;
        next if ($len < $min_sv_len) and ($type ne 'INS');
        next if ($len < $min_ins_len) and ($type eq 'INS');
        next if ($len > $max_sv_len);
        my $qual = 0;
        $qual = $line[5];
        $qual = int ($qual + 0.5);
        my $read = 2;
        if ($qual < 10){
            $read = 2;
        }
        elsif ($qual < 50){
            $read = 3;
        }
        elsif ($qual < 100){
            $read = 4;
        }
        elsif ($qual < 300){
            $read = 5;
        }
        elsif ($qual < 500){
            $read = 6;
        }
        elsif ($qual < 700){
            $read = 7;
        }
        elsif ($qual < 1000){
            $read = 8;
        }
        elsif ($qual < 1500){
            $read = 9;
        }
        elsif ($qual < 2000){
            $read = 10;
        }
        else{
            $read = 12;
        }
        print "$chr\t$pos\t.\t.\t.\t$qual\t$pass\tSVTYPE=$type;SVLEN=$len;READS=$read\n";
    }
}
close (FILE);
