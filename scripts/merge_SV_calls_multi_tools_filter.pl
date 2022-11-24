#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

# select high quality SVs (overlap calls and single calls) based on the specified list of algorithms and RSS

my $help;

my $sv_type = 'ALL';

my $vcf = '';

my @del;
my @dup;
my @ins;
my @inv;
my @del_ss;
my @del_s;
my @del_m;
my @del_l;
my @dup_s;
my @dup_m;
my @dup_l;
my @inv_s;
my @inv_m;
my @inv_l;

GetOptions(
    'vcf|v=s' => \$vcf,
    'del=s{,}' => \@del,
    'dup=s{,}' => \@dup,
    'ins=s{,}' => \@ins,
    'inv=s{,}' => \@inv,
    'del_ss=s{,}' => \@del_ss,
    'del_s=s{,}' => \@del_s,
    'del_m=s{,}' => \@del_m,
    'del_l=s{,}' => \@del_l,
    'dup_s=s{,}' => \@dup_s,
    'dup_m=s{,}' => \@dup_m,
    'dup_l=s{,}' => \@dup_l,
    'inv_s=s{,}' => \@inv_s,
    'inv_m=s{,}' => \@inv_m,
    'inv_l=s{,}' => \@inv_l,
    'sv_type|t=s' => \$sv_type,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;


=head1 SYNOPSIS
  
  # Example for DEL
  merge_SV_calls_multi_tools_filter.pl -v <input merged vcf file> -t DEL --del_s <list of tool name:RSS psir(s)> --del_m <list of tool name:RSS psir(s)> --del_l <list of tool name:RSS psir(s)>  > [output vcf]

  Options:
   --vcf or -v <STR>        input vcf file (generated with merge_SV_calls_multi_tools.pl or merge_SV_calls_multi_tools_eachSize.pl)
   --del <STR>              list (space-delimited) of tool name:RSS psir(s) to select overlap DEL calls (e.g., Pindel:3=Delly:5 CNVnator:3=Lumpy:7)
   --dup <STR>              list (space-delimited) of tool name:RSS psir(s) to select overlap DUP calls
   --ins <STR>              list (space-delimited) of tool name:RSS psir(s) to select overlap INS calls
   --inv <STR>              list (space-delimited) of tool name:RSS psir(s) to select overlap INV calls
   --del_ss <STR>           list (space-delimited) of tool name:RSS psir(s) to select tiny size (<= 100 bp) of overlap DEL calls
   --del_s <STR>            list (space-delimited) of tool name:RSS psir(s) to select small size (50-1000 bp) of overlap DEL calls
   --del_m <STR>            list (space-delimited) of tool name:RSS psir(s) to select midlde size (1-100 Kb) of overlap DEL calls
   --del_l <STR>            list (space-delimited) of tool name:RSS psir(s) to select large size (> 100 Kb) of overlap DEL calls
   --dup_s <STR>            list (space-delimited) of tool name:RSS psir(s) to select small size (50-1000 bp) of overlap DUP calls
   --dup_m <STR>            list (space-delimited) of tool name:RSS psir(s) to select midlde size (1-100 Kb) of overlap DUP calls
   --dup_l <STR>            list (space-delimited) of tool name:RSS psir(s) to select large size (> 100 Kb) of overlap DUP calls
   --inv_s <STR>            list (space-delimited) of tool name:RSS psir(s) to select small size (50-1000 bp) of overlap INV calls
   --inv_m <STR>            list (space-delimited) of tool name:RSS psir(s) to select midlde size (1-100 Kb) of overlap INV calls
   --inv_l <STR>            list (space-delimited) of tool name:RSS psir(s) to select large size (> 100 Kb) of overlap INV calls
   --sv_type or -t <STR>    SV type, DEL|DUP|INS|INV|TRA|MEI|VEI|NUMT|ALL [default: all]
   --help or -h             output help message
   
=cut

my %tool_pairs;

my %svtype = ('DEL-ALL' => \@del, 'DUP-ALL' => \@dup, 'INS-ALL' => \@ins, 'INV-ALL' => \@inv, 'DEL-SS' => \@del_ss, 'DEL-S' => \@del_s, 'DEL-M' => \@del_m, 'DEL-L' => \@del_l,
              'DUP-S' => \@dup_s, 'DUP-M' => \@dup_m, 'DUP-L' => \@dup_l, 'INV-S' => \@inv_s, 'INV-M' => \@inv_m, 'INV-L' => \@inv_l);

foreach my $type_size (keys %svtype){
    my ($type, $size) = split (/-/, $type_size);
    if (@{$svtype{$type_size}} > 0){
        foreach (@{$svtype{$type_size}}){
            push @{${$tool_pairs{$type}}{$size}}, $_;
        }
    }
}

my $ins_sd = 100;
my $min_overlap_rate = 0.8;

my %vcf;
my %vcf_all;

open (FILE, $vcf) or die "$vcf is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
        print "$line\n";
        next;
    }
    my $type = $1 if ($line =~ /SVTYPE=(.+?);/);
    next if ($sv_type ne 'ALL') and ($type ne $sv_type);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $len = $1 if ($line =~ /SVLEN=(\d+)/);
    ${${$vcf_all{$type}}{$chr}}{$pos} = $len;
    my $size = 'ALL';
    my $size2 = $size;
    if ($type ne 'INS'){
        if ($len <= 1000){
            $size = 'S';
            if (($len <= 100) and ($type eq 'DEL')){
                if (exists ${$tool_pairs{$type}}{'SS'}){
                    $size = 'SS';
                    $size2 = 'SS';
                }
            }
        }
        elsif ($len <= 100000){
            $size = 'M';
        }
        else{
            $size = 'L';
        }
        if ($size eq 'S'){
            if ($len >= 900){
                $size2 = 'M';
            }
            else{
                $size2 = 'S';
            }
        }
        if ($size eq 'M'){
            if ($len <= 1100){
                $size2 = 'S';
            }
            elsif ($len >= 90000){
                $size2 = 'L';
            }
            else{
                $size2 = 'M';
            }
        }
        if ($size eq 'L'){
            if ($len <= 110000){
                $size2 = 'M';
            }
            else{
                $size2 = 'L';
            }
        }
    }
    my $overlap_tools = $1 if ($line =~ /TOOLS=([^\;]+)/);
    my @tool_set = split (/,/, $overlap_tools);
    my %tool_RSS;

    if (@tool_set >= 1){        # select SVs based on the specified algorithm's pairs and RSS
        foreach (@tool_set){
            my $tool = $1 if ($_ =~ /(.+?):/);
            my $rss = $1 if ($_ =~ /:(\d+)$/);
            $tool_RSS{$tool} = $rss;
        }
        my $hit_flag = 0;
        if (($hit_flag == 0) and (exists ${$tool_pairs{$type}}{$size})){
            foreach my $tpair (@{${$tool_pairs{$type}}{$size}}){
                my ($tool_rss1, $tool_rss2) = split (/=/, $tpair) if ($tpair =~ /=/);
                my ($tool1, $rss1) = ('', 0);
                my ($tool2, $rss2) = ('', 0);
                if ($tpair !~ /=/){
                    ($tool1, $rss1) = split (/:/, $tpair);
                    if (exists $tool_RSS{$tool1}){
                        if ($rss1 <= $tool_RSS{$tool1}){
                            $hit_flag = 1;
                            if (($tool1 eq 'CNVnator') and ($type eq 'DEL') and ($size eq 'M') and ($len < 4000)){
                                $hit_flag = 0;
                            }
                            elsif (($tool1 eq 'CNVnator') and ($type eq 'DEL') and ($size eq 'M') and ($rss1 <= 7) and ($len < 10000)){
                                $hit_flag = 0;
                            }
                            else{
                                last;
                            }
                        }
                    }
                }
                else{
                    ($tool1, $rss1) = split (/:/, $tool_rss1);
                    ($tool2, $rss2) = split (/:/, $tool_rss2);
                    if ((exists $tool_RSS{$tool1}) and ($tool2 eq '#')){
                        if (($rss1 <= $tool_RSS{$tool1}) and (scalar keys %tool_RSS > 1)){
                            $hit_flag = 1;
                            last;
                        }
                    }
                    elsif ((exists $tool_RSS{$tool1}) and (exists $tool_RSS{$tool2})){
                        if (($rss1 <= $tool_RSS{$tool1}) and ($rss2 <= $tool_RSS{$tool2})){
                            $hit_flag = 1;
                            last;
                        }
                    }
                }
            }
        }
        if ($hit_flag == 1){
            ${${$vcf{$type}}{$chr}}{$pos} = $line;
        }
        elsif ($size ne $size2){
            if (exists ${$tool_pairs{$type}}{$size2}){
                foreach my $tpair (@{${$tool_pairs{$type}}{$size2}}){
                    my ($tool_rss1, $tool_rss2) = split (/=/, $tpair) if ($tpair =~ /=/);
                    my ($tool1, $rss1) = ('', 0);
                    my ($tool2, $rss2) = ('', 0);
                    if ($tpair !~ /=/){
                        ($tool1, $rss1) = split (/:/, $tpair);
                        if (exists $tool_RSS{$tool1}){
                            if ($rss1 <= $tool_RSS{$tool1}){
                                $hit_flag = 1;
                                last;
                            }
                        }
                    }
                    else{
                        ($tool1, $rss1) = split (/:/, $tool_rss1);
                        ($tool2, $rss2) = split (/:/, $tool_rss2);
                        if ((exists $tool_RSS{$tool1}) and ($tool2 eq '#')){
                            if (($rss1 <= $tool_RSS{$tool1}) and (scalar keys %tool_RSS > 1)){
                                $hit_flag = 1;
                                last;
                            }
                        }
                        elsif ((exists $tool_RSS{$tool1}) and (exists $tool_RSS{$tool2})){
                            if (($rss1 <= $tool_RSS{$tool1}) and ($rss2 <= $tool_RSS{$tool2})){
                                $hit_flag = 1;
                                last;
                            }
                        }
                    }
                }
            }
            if ($hit_flag == 1){
                ${${$vcf{$type}}{$chr}}{$pos} = $line;
            }
        }
    }
}
close (FILE);

foreach my $type (keys %vcf){       # remove redundant SVs
    foreach my $chr (keys %{$vcf{$type}}){
        my $pre_pos = 0;
        my $pre_end = 0;
        my $pre_read = 0;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf{$type}}{$chr}}){
            my $line = ${${$vcf{$type}}{$chr}}{$pos};
            my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            my $read = $1 if ($line =~ /READS=(\d+)/);
            if ($type eq 'INS'){
                if (abs ($pos - $pre_pos) <= $ins_sd){
                    if ($pre_read > $read){
                        delete ${${$vcf{$type}}{$chr}}{$pos};
                        next;
                    }
                    else{
                        delete ${${$vcf{$type}}{$chr}}{$pre_pos};
                    }
                }
            }
            else{
                if ($pos < $pre_end){
                    my $overlap = $pre_end - $pos + 1;
                    $overlap = $len if ($end < $pre_end);
                    my $pre_len = $pre_end - $pre_pos + 1;
                    if (($overlap >= $len * $min_overlap_rate) and ($overlap >= $pre_len * $min_overlap_rate)){
                        if ($pre_read > $read){
                            delete ${${$vcf{$type}}{$chr}}{$pos};
                            next;
                        }
                        else{
                            delete ${${$vcf{$type}}{$chr}}{$pre_pos};
                        }
                    }
                }
            }
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_read = $read;
        }
    }
}

my %vcf2;
my $min_len = 1000;
my $min_distance_rate = 0.01;

foreach my $type (keys %vcf){   #merge adjacent SVs with distance of < 1% of both the SV lengths (with >= 1 Kb in size for both the SVs) and with evidence from the original input SV vcf
    next if ($type eq 'INS');
    foreach my $chr (keys %{$vcf{$type}}){
        my $pre_pos = 0;
        my $pre_len = 0;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf{$type}}{$chr}}){
            my $line = ${${$vcf{$type}}{$chr}}{$pos};
            my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            if ($pre_pos > 0){
                if (($len >= $min_len) and ($pre_len >= $min_len)){
                    my $pre_end = $pre_pos + $pre_len - 1;
                    my $distance = $pos - $pre_end + 1;
                    if (($distance <= $len * $min_distance_rate) and ($distance <= $pre_len * $min_distance_rate) and ($end > $pre_end)){
                        my $new_len = $end - $pre_pos + 1;
                        $new_len = $pre_len if ($pre_end > $end);
                        my $flag = 0;
                        foreach my $pos2 (sort {$a <=> $b} keys %{${$vcf_all{$type}}{$chr}}){
                            last if ($pos2 > $end);
                            next if ($pos2 == $pos) or ($pos2 == $pre_pos);
                            my $len2 = ${${$vcf_all{$type}}{$chr}}{$pos2};
                            my $end2 = $pos2 + $len2 - 1;
                            next if ($end2 < $pre_pos);
                            if (($pos2 <= $pre_pos) and ($end2 >= $end)){
                                if ($new_len >= $len2 * $min_overlap_rate){
                                    $flag = 1;
                                }
                            }
                            elsif (($pos2 >= $pre_pos) and ($pos2 <= $end)){
                                my $overlap = $end - $pos2;
                                $overlap = $len2 if ($end2 < $end);
                                if (($overlap >= $new_len * $min_overlap_rate) and ($overlap >= $len2 * $min_overlap_rate)){
                                    $flag = 1;
                                }
                            }
                            elsif (($end2 >= $pre_pos) and ($end2 <= $end)){
                                my $overlap = $end2 - $pre_pos;
                                $overlap = $len2 if ($pos2 > $pre_pos);
                                if (($overlap >= $new_len * $min_overlap_rate) and ($overlap >= $len2 * $min_overlap_rate)){
                                    $flag = 1;
                                }
                            }
                            last if ($flag == 1);
                        }
                        if ($flag == 1){
                            my $pre_line = ${${$vcf{$type}}{$chr}}{$pre_pos};
                            my $pre_read = $1 if ($pre_line =~ /READS=(\d+)/);
                            my $pre_tools = $1 if ($pre_line =~ /TOOLS=(\S+)/);
                            my $read = $1 if ($line =~ /READS=(\d+)/);
                            my $tools = $1 if ($line =~ /TOOLS=(\S+)/);
                            my $new_read = int (($pre_read + $read) * 0.5 + 0.5);
                            my %tools;
                            my @pre_tools = split (/,/, $pre_tools);
                            my @tools = split (/,/, $tools);
                            foreach (@pre_tools){
                                my ($tool) = split (/:/, $_);
                                $tools{$tool} = 1;
                            }
                            foreach (@tools){
                                my ($tool) = split (/:/, $_);
                                if (!exists $tools{$tool}){
                                    $pre_tools .= ",$_";
                                }
                            }
                            $pre_tools =~ s/,$// if ($pre_tools =~ /,$/);
                            delete ${${$vcf{$type}}{$chr}}{$pos};
                            $pre_line =~ s/SVLEN=\d+/SVLEN=$new_len/ if ($pre_end < $end);
                            $pre_line =~ s/READS=\d+/READS=$new_read/;
                            $pre_line =~ s/TOOLS=.+/TOOLS=$pre_tools/;
                            ${${$vcf{$type}}{$chr}}{$pre_pos} = $pre_line;
                            $pre_len = $new_len;
                            next;
                        }
                    }
                }
            }
            if ($len >= $min_len){
                $pre_pos = $pos;
                $pre_len = $len;
            }
        }
    }
}

foreach my $type (keys %vcf){
    foreach my $chr (keys %{$vcf{$type}}){
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
        foreach my $pos (keys %{${$vcf{$type}}{$chr}}){
            ${${$vcf2{$chr02d}}{$pos}}{$type} = ${${$vcf{$type}}{$chr}}{$pos};
        }
    }
}

foreach my $chr (sort keys %vcf2){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf2{$chr}}){
        foreach my $type (keys %{${$vcf2{$chr}}{$pos}}){
            print "${${$vcf2{$chr}}{$pos}}{$type}\n";
        }
    }
}

