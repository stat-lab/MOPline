#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use File::Basename;

# merge SV calls from multiple SV sets and remove redundant SV calls


my @tools = ();

#@tools = ('inGAP', 'CNVnator', 'MELT.MEI', 'Manta', 'Wham', 'MATCHCLIP');
my @tools_2 = ('inGAP', 'CNVnator', 'MELT', 'Manta', 'SoftSV', 'Wham');

my @min_read = ();

my $min_sv_len = 50;
my $max_sv_len = 30000000;

my $max_inv_size = 200000;

my $sv_type = 'ALL';

my $non_human = 0;

my $read_length = 150;
my $var_sd = 125;
my $ins_sd = 200;
my $ins_sd2 = 100;

my $min_overlap_ratio = 0.5;
my $min_overlap_ratio2 = 0.8;
my $min_overlap_ratio3 = 0.4;

my $min_read = 3;

my $pos_ave = 1;

my $input_dir = '.';

my $sample_id = '';

my $help;

GetOptions(
    'tool|t=s{,}' => \@tools,
    'sv_type|st=s' => \$sv_type,
    'min_read|m=i{,}' => \@min_read,
    'len|l=i' => \$min_sv_len,
    'xlen|xl=i' => \$max_sv_len,
    'var_sd|vsd=i' => \$var_sd,
    'read_len|rl=i' => \$read_length,
    'dir|d=s' => \$input_dir,
    'sample|s=s' => \$sample_id,
    'non_human|nh=i' => \$non_human,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  merge_SV_calls_multi_tools.pl -t <tool name list> -s <sample name> -st <SV type>  > [output vcf]

  Options:
   --tools or -t <STR>      list (space-delimited) of tool names 
   --sample or -s <STR>     sample name/ID
   --dir or -d <STR>        directory of input vcf files [default: ./]
   --sv_type or -st <STR>   SV type, DEL|DUP|INS|INV|ALL [default: ALL]
   --read_len or -rl <INT>  read length [default: 150]
   --len or -l <INT>        minimum length (bp) of SV [default: 50]
   --xlen or -x <INT>       maximum length (bp) of SV [default: 30,000,000]
   --var_sd or -vsd <INT>   standard deviation of break points [default: 150]
   
   --non_human or -nh <INT> species is non-human (human: 0, non-human: 1) [default: 0];
   --help or -h             output help message
   
=cut

@tools = (@tools_2) if (@tools == 0);

$var_sd = $read_length;

my @work_files = <$input_dir/*.vcf>;
my @work_files2;
if (@work_files >= @tools){
    foreach (@work_files){
        push @work_files2, $_ if ($_ !~ /\.TP-FP|\.FP\.vcf$|\.FN\.vcf$|=|~|,|Merge|merge|\.DEL|\.INS|\.DUP|\.INV/);
    }
}
else{
    foreach my $tool (@tools){
        my $tool_dir = $tool;
        $tool_dir =~ s/\./-/ if ($tool_dir =~ /\.MEI|\.NUMT|\.VEI/);
        $tool_dir =~ s/NUMT/MT/ if ($tool_dir =~ /-NUMT/);
        my $vcf = "$input_dir/$tool_dir/$tool.$sample_id.vcf";
        if (!-f $vcf){
            die "$vcf file is not found:\n";
        }
        else{
            push @work_files2, $vcf;
        }
    }
}

my %SRR;
my %SDlen;
my %SDbp;
my %tool_preci;

#die "var_files should be < 4:\n" if (@var_file > 4) or (@var_id > 4);
my $data_dir = "$Bin/../Data";
$data_dir = "$ENV{MOPLINE_DIR}/Data" if (exists $ENV{MOPLINE_DIR});

#my $sv_min_read_file = "$data_dir/SV_min_read_list.txt";    # file showing minimum number of reads supporting SV (RSS) to select a moderate quality of SVs from existing SV detection algorithms

# the data of these files are used for the determination of breakpoints and sizes of SVs overlapping between algorithms

my $SDlen_rank_file = "$data_dir/SV_SDlen_rank_list.txt";   # file showing breakpoint accuracy for existing SV detection algorithms

my $tool_eval_file = "$data_dir/NA78_data1.eval3.txt";      # file showing precision and recall of SV calls detected with NA12878 WGS data, which was generated with evaluate_SV_callers.pl in https://github.com/stat-lab/EvalSVcallers

=pod
open (FILE, $sv_min_read_file) or die "$sv_min_read_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my ($tool, $type, $size, $srr) = split (/\s+/, $line);
    $tool = 'GASV' if ($tool eq 'GASVPro');
    $tool = 'Basil' if ($tool eq 'BASIL-ANISE');
    $tool = '123SV' if ($tool eq '1-2-3-SV');
    $tool = 'GenomeStrip' if ($tool eq 'GenomeSTRiP');
    $tool = 'VH' if ($tool eq 'VariationHunter');
    $tool = 'metaSV' if ($tool eq 'MetaSV');
    $tool = 'inGAP' if ($tool eq 'inGAP-sv');
    $tool = 'MELT' if ($tool eq 'MELT.MEI');
    $tool = 'Mobster' if ($tool eq 'Mobster.MEI');
    ${${$SRR{$type}}{$size}}{$tool} = $srr;
}
close (FILE);
=cut

open (FILE, $SDlen_rank_file) or die ("$SDlen_rank_file is not found: $!\n");
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#/);
    my @line = split (/\s+/, $line);
    my $id = $line[1];
    my $type = $line[0];
    my $sd_len_rank = $line[5];
    my $sd_bp_rank = $line[4];
    my $bprank = 1;
    my $lenrank = 1;
    if ($sd_bp_rank eq 'M'){
        $bprank = 2;
    }
    elsif ($sd_bp_rank eq 'L'){
        $bprank = 3;
    }
    if ($sd_len_rank eq 'M'){
        $lenrank = 2;
    }
    elsif ($sd_len_rank eq 'L'){
        $lenrank = 3;
    }
    elsif ($sd_len_rank eq 'LL'){
        $lenrank = 4;
    }
    ${$SDbp{$type}}{$id} = $bprank;
    ${$SDlen{$type}}{$id} = $lenrank;
}
close (FILE);

my $tool_2 = '';
my $type_2 = '';
open (FILE, $tool_eval_file) or die "$tool_eval_file is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    if ($line !~ /^\s+/){
        $tool_2 = $line[0];
        $type_2 = $line[1];
        $tool_2 = 'GASV' if ($tool_2 eq 'GASVPro');
        $tool_2 = 'Basil' if ($tool_2 eq 'BASIL-ANISE');
        $tool_2 = '123SV' if ($tool_2 eq '1-2-3-SV');
        $tool_2 = 'GenomeStrip' if ($tool_2 eq 'GenomeSTRiP');
        $tool_2 = 'VH' if ($tool_2 eq 'VariationHunter');
        $tool_2 = 'metaSV' if ($tool_2 eq 'MetaSV');
        $tool_2 = 'inGAP' if ($tool_2 eq 'inGAP-sv');
        $tool_2 = 'MELT' if ($tool_2 eq 'MELT.MEI');
        $tool_2 = 'Mobster' if ($tool_2 eq 'Mobster.MEI');
        next;
    }
    if ($line =~ /Precis\s\((\S+)\)/){
        my $size2 = $1;
        $size2 = 'ALL' if ($size2 eq 'A');
        my $count = 0;
        my $read = 2;
        foreach (@line){
            $count ++;
            next if ($count <= 4);
            my $type_size = "$type_2-$size2";
            ${${${$tool_preci{$type_2}}{$size2}}{$tool_2}}{$read} = $_;
            $read ++;
            $read ++ if ($read > 10);
        }
    }
}
close (FILE);

foreach my $type (keys %tool_preci){
    foreach my $size (keys %{$tool_preci{$type}}){
        foreach my $tool (keys %{${$tool_preci{$type}}{$size}}){
            my $pre_preci = 0;
            my $preci_2 = 0;
            for (my $i = 2; $i <= 30; $i++){
                $preci_2 = ${${${$tool_preci{$type}}{$size}}{$tool}}{$i} if ($i == 2);
                if (!exists ${${${$tool_preci{$type}}{$size}}{$tool}}{$i}){
                    ${${${$tool_preci{$type}}{$size}}{$tool}}{$i} = $pre_preci;
                }
                $pre_preci = ${${${$tool_preci{$type}}{$size}}{$tool}}{$i};
            }
            ${${${$tool_preci{$type}}{$size}}{$tool}}{1} = $preci_2;
        }
    }
}

my %var_file;
my %call;
my %call_count;
my %tool_size;
my %tool_sv;

@tools = sort {"\L$a" cmp "\L$b"} @tools;

print STDERR "$sv_type\tNum of tools: ", scalar @tools, "\n";

foreach my $tool (@tools){
    $tool =~ s/:\d+$// if ($tool =~ /:\d+$/);
    my $type = '';
    my $size = '';
    my $tool2 = $tool;
    if ($tool =~ /:/){
        my @tool_str = split (/:/, $tool);
        if (@tool_str == 2){
            $tool2 = $tool_str[0];
            $size = $tool_str[1];
            $type = $sv_type;
        }
        elsif (@tool_str == 3){
            $tool2 = $tool_str[0];
            $type = $tool_str[1];
            $size = $tool_str[2];
        }
    }
=pod
    $tool2 = 'inGAP-sv' if ($tool eq 'inGAP');
    $tool2 = 'PBHoney-NGM' if ($tool =~ /PBHoney/);
    $tool2 = 'MELT.MEI' if ($tool =~ /MELT/);
    $tool2 = 'Mobster.MEI' if ($tool =~ /Mobster/);
    $tool2 = 'DELLY' if ($tool eq 'Delly');
    $tool2 = 'PopIns' if ($tool eq 'Popins');
    $tool2 = '1-2-3-SV' if ($tool eq '123SV');
    $tool2 = 'VariationHunter' if ($tool eq 'VH');
    $tool2 = 'MetaSV' if ($tool eq 'metaSV');
=cut
    my $exist_flag = 0;
    foreach my $file (@work_files2){
        my $file_base = basename ($file);
        if (($file_base =~ /$tool2\./i) and ($file_base !~ /=|~|,|Merge/)){
            $var_file{$tool2} = $file;
            $exist_flag = 1;
            if ($size ne ''){
                ${${$tool_size{$tool2}}{$type}}{$size} = 1;
print STDERR "$tool2\t$type\t$size\n";
            }
print STDERR "$tool:\t$file\n";
            last;
        }
    }
    if ($exist_flag == 0){
        die "vcf file for $tool is not found:\n";
    }
}

foreach my $id (keys %var_file){
    my $var_file = $var_file{$id};
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#/);
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        next if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
        my $pos = $line[1];
        my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        $type = 'INS' if ($type =~ /MEI|NUMT|VEI/);
        if ($sv_type ne 'ALL'){
            next if ($type ne $sv_type);
        }
        next if ($type !~ /^DEL$|^DUP$|^INS$|^INV$/);
        ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos} = $line;
    }
    close (FILE);
    
    foreach my $type (keys %{$tool_sv{$id}}){
        foreach my $chr (keys %{${$tool_sv{$id}}{$type}}){
            my $pre_pos = 0;
            my $pre_len = 0;
            my $pre_read = 0;
            foreach my $pos (sort {$a <=> $b} keys %{${${$tool_sv{$id}}{$type}}{$chr}}){
                my $line = ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                my $len = 0;
                $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
                my $end = $pos + $len - 1;
                my $read = $1 if ($line =~ /READS=(\d+)/);
                if ($pre_pos > 0){
                    if ($type eq 'INS'){
                        if ($pos - $pre_pos <= $ins_sd2){
                            if ($pre_read >= $read){
                                delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                next;
                            }
                            else{
                                delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos};
                            }
                        }
                    }
                    else{
                        my $pre_end = $pre_pos + $pre_len - 1;
                        if ($pos < $pre_end){
                            my $overlap = $pre_end - $pos + 1;
                            $overlap = $len if ($end < $pre_end);
                            my $min_read_fold = 4;
                            $min_read_fold = 3 if ($type eq 'DUP');
                            $min_read_fold = 2 if ($type eq 'INV');
                            if ($type ne 'DEL'){
                                if ($end < $pre_end){
                                    delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                    next;
                                }
                                if ($pre_read >= $read * $min_read_fold){
                                    delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                    next;
                                }
                                elsif ($read >= $pre_read * $min_read_fold){
                                    delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos};
                                }
                                else{
                                    my $pre_line = ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos};
                                    my $new_len = $end - $pre_pos + 1;
                                    $pre_line =~ s/SVLEN=-*\d+/SVLEN=$new_len/;
                                    my $new_read = int (($read + $pre_read) * 0.5 + 0.5);
                                    $pre_line =~ s/READS=\d+/READS=$new_read/;
                                    ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos} = $pre_line;
                                    delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                    $pre_len = $new_len;
                                    $pre_read = $new_read;
                                    next;
                                }
                            }
                            elsif ($type eq 'DEL'){
                                if (($overlap >= $pre_len * $min_overlap_ratio3) or ($overlap >= $len * $min_overlap_ratio3)){
                                    if ($pre_read >= $read * $min_read_fold){
                                        delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                        next;
                                    }
                                    elsif ($read >= $pre_read * $min_read_fold){
                                        delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos};
                                    }
                                    else{
                                        my $pre_line = ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos};
                                        if ($end > $pre_end){
                                            my $new_len = $end - $pre_pos + 1;
                                            $pre_line =~ s/SVLEN=-*\d+/SVLEN=$new_len/;
                                            $pre_len = $new_len;
                                        }
                                        my $new_read = int (($read + $pre_read) * 0.5 + 0.5);
                                        $pre_line =~ s/READS=\d+/READS=$new_read/;
                                        ${${${$tool_sv{$id}}{$type}}{$chr}}{$pre_pos} = $pre_line;
                                        delete ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                                        $pre_read = $new_read;
                                        next;
                                    }
                                }
                            }
                        }
                    }
                }
                $pre_pos = $pos;
                $pre_len = $len;
                $pre_read = $read;
            }
        }
    }
}      

foreach my $id (keys %tool_sv){
    foreach my $type (keys %{$tool_sv{$id}}){
        foreach my $chr (keys %{${$tool_sv{$id}}{$type}}){
            foreach my $pos (sort {$a <=> $b} keys %{${${$tool_sv{$id}}{$type}}{$chr}}){
                my $line = ${${${$tool_sv{$id}}{$type}}{$chr}}{$pos};
                my @line = split (/\t/, $line);
                my $chr = $line[0];
                next if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
                my $pos = $line[1];
                my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
                $type = 'INS' if ($type =~ /MEI|NUMT|VEI/);
                if ($sv_type ne 'ALL'){
                    next if ($type ne $sv_type);
                }
                next if ($type !~ /^DEL$|^DUP$|^INS$|^INV$/);
                my $len = 0;
                $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
                next if ($type eq 'INV') and ($len > $max_inv_size);
                my $end = 0;
                $end = $pos + $len - 1;
                my $size = 'ALL';
                if (($type eq 'DEL') or ($type eq 'DUP')){
                    if ($len <= 1000){
                        $size = 'S';
                        if (($type eq 'DEL') and ($len <= 100)){
                            $size = 'SS';
                        }
                    }
                    elsif ($len <= 100000){
                        $size = 'M';
                    }
                    else{
                        $size = 'L';
                    }
                }
                my %tsize;
                my $tsize_flag = 0;
                my $tsize_all = 0;
                if (exists ${$tool_size{$id}}{$type}){
                    foreach my $tsize (keys %{${$tool_size{$id}}{$type}}){
                        $tsize{$tsize} = 1;
                        $tsize_flag = 1;
                        $tsize_all = 1 if ($tsize eq 'ALL');
                    }
                }
                next if ($tsize_flag == 1) and ($tsize_all == 0) and (!exists $tsize{$size});       
                ${$call_count{$id}}{$type} = 0 if (!exists ${$call_count{$id}}{$type});
                my $reads = 3;
                $reads = $1 if ($line[7] =~ /READS=(\d+)/);
                my $min_read2 = $min_read;
                if ($id =~ /Sniffles|PBHoney/){
                    $min_read2 = 2;
                }
=pod
                if (exists ${${$SRR{$type}}{$size}}{$id}){
                    $min_read2 = ${${$SRR{$type}}{$size}}{$id};
                }
=cut
                next if ($reads < $min_read2);
                my $chr2 = '';
                my $pos2 = 0;
                next if ($len < $min_sv_len) and ($len > 0) and ($type ne 'INS') and ($type ne 'TRA');
                next if ($len > $max_sv_len);
                my $subtype = $line[2];
                $subtype = 'INS' if ($subtype =~ /INS/);
                if ((!exists ${${$call{$type}}{$chr}}{$pos}) or ($type eq 'INS')){
                    push @{${${$call{$type}}{$chr}}{$pos}}, "$id=$pos=$len=$reads" if ($type ne 'INS');
                    push @{${${$call{$type}}{$chr}}{$pos}}, "$id=$pos=$len=$reads=$subtype" if ($type eq 'INS');
                }
                else{
                    my $sumlen = 0;
                    foreach (@{${${$call{$type}}{$chr}}{$pos}}){
                        my ($id2, $pos2, $len2) = split (/=/, $_);
                        $sumlen += $len2;
                    }
                    my $avelen = int ($sumlen / @{${${$call{$type}}{$chr}}{$pos}} + 0.5);
                    if (($len / $avelen > 2) or ($len / $avelen < 0.5)){
                        if (@{${${$call{$type}}{$chr}}{$pos}} == 1){
                            if ($avelen < $len){
                                delete ${${$call{$type}}{$chr}}{$pos};
                                push @{${${$call{$type}}{$chr}}{$pos}}, "$id=$pos=$len=$reads";
                            }
                        }
                    }
                    else{
                        push @{${${$call{$type}}{$chr}}{$pos}}, "$id=$pos=$len=$reads";
                    }
                }
                ${$call_count{$id}}{$type} ++;
            }
        }
    }
}
%tool_sv = ();

foreach my $type (keys %call){
    foreach my $chr (keys %{$call{$type}}){
        my $pre_pos = 0;
        my $pre_len = 0;
        foreach my $pos (sort {$a <=> $b} keys %{${$call{$type}}{$chr}}){
            my $sum_len = 0;
            my @len;
            foreach (@{${${$call{$type}}{$chr}}{$pos}}){
                my ($id, $pos1, $len) = split (/=/, $_);
                push @len, $len;
            }
            map{$sum_len += $_} @len;
            my $avelen = int ($sum_len / @len + 0.5);
            if ($pre_pos > 0){
                if ($type eq 'INS'){
                    if ($pos - $pre_pos <= $var_sd){
                        my $pre_num = scalar @{${${$call{$type}}{$chr}}{$pre_pos}};
                        my $num = scalar @{${${$call{$type}}{$chr}}{$pos}};
                        my $new_pos = int (($pre_pos * $pre_num + $pos * $num) / ($pre_num + $num) + 0.5);
                        my @new_info = (@{${${$call{$type}}{$chr}}{$pre_pos}}, @{${${$call{$type}}{$chr}}{$pos}});
                        delete ${${$call{$type}}{$chr}}{$pre_pos};
                        delete ${${$call{$type}}{$chr}}{$pos};
                        push @{${${$call{$type}}{$chr}}{$new_pos}}, @new_info;
                        $pre_pos = $new_pos;
                        next;
                    }
                }
                else{
                    my $end = $pos + $avelen - 1;
                    my $pre_end = $pre_pos + $pre_len - 1;
                    if ($pos < $pre_end){
                        my $overlap = $pre_end - $pos + 1;
                        $overlap = $avelen if ($end < $pre_end);
                        if (($overlap >= $avelen * $min_overlap_ratio2) and ($overlap >= $pre_len * $min_overlap_ratio2)){
                            my $pre_num = scalar @{${${$call{$type}}{$chr}}{$pre_pos}};
                            my $num = scalar @{${${$call{$type}}{$chr}}{$pos}};
                            my $new_pos = int (($pre_pos * $pre_num + $pos * $num) / ($pre_num + $num) + 0.5);
                            my @new_info = (@{${${$call{$type}}{$chr}}{$pre_pos}}, @{${${$call{$type}}{$chr}}{$pos}});
                            my @new_len;
                            my $sum_len2 = 0;
                            foreach (@new_info){
                                my ($id, $pos1, $len) = split (/=/, $_);
                                push @new_len, $len;
                            }
                            map{$sum_len2 += $_} @new_len;
                            my $new_len = int ($sum_len2 / @new_len + 0.5);
                            delete ${${$call{$type}}{$chr}}{$pre_pos};
                            delete ${${$call{$type}}{$chr}}{$pos};
                            push @{${${$call{$type}}{$chr}}{$new_pos}}, @new_info;
                            $pre_pos = $new_pos;
                            $pre_len = $new_len;
                            next;
                        }
                    }
                }
            }
            $pre_len = $avelen;
            $pre_pos = $pos;
        }
    }
}

foreach my $type (keys %call){
    foreach my $chr (keys %{$call{$type}}){
        my $pre_pos = 0;
        my $pre_len = 0;
        foreach my $pos (sort {$a <=> $b} keys %{${$call{$type}}{$chr}}){
            my $sum_len = 0;
            my @len;
            foreach (@{${${$call{$type}}{$chr}}{$pos}}){
                my ($id, $pos1, $len) = split (/=/, $_);
                push @len, $len;
            }
            map{$sum_len += $_} @len;
            my $avelen = int ($sum_len / @len + 0.5);
            if ($pre_pos > 0){
                if ($type eq 'INS'){
                    if ($pos - $pre_pos <= $var_sd){
                        my $pre_num = scalar @{${${$call{$type}}{$chr}}{$pre_pos}};
                        my $num = scalar @{${${$call{$type}}{$chr}}{$pos}};
                        my $new_pos = int (($pre_pos * $pre_num + $pos * $num) / ($pre_num + $num) + 0.5);
                        my @new_info = (@{${${$call{$type}}{$chr}}{$pre_pos}}, @{${${$call{$type}}{$chr}}{$pos}});
                        delete ${${$call{$type}}{$chr}}{$pre_pos};
                        delete ${${$call{$type}}{$chr}}{$pos};
                        push @{${${$call{$type}}{$chr}}{$new_pos}}, @new_info;
                        $pre_pos = $new_pos;
                        next;
                    }
                }
                else{
                    my $end = $pos + $avelen - 1;
                    my $pre_end = $pre_pos + $pre_len - 1;
                    if ($pos < $pre_end){
                        my $overlap = $pre_end - $pos + 1;
                        $overlap = $avelen if ($end < $pre_end);
                        if (($overlap >= $avelen * $min_overlap_ratio) and ($overlap >= $pre_len * $min_overlap_ratio)){
                            my $pre_num = scalar @{${${$call{$type}}{$chr}}{$pre_pos}};
                            my $num = scalar @{${${$call{$type}}{$chr}}{$pos}};
                            my $new_pos = int (($pre_pos * $pre_num + $pos * $num) / ($pre_num + $num) + 0.5);
                            my @new_info = (@{${${$call{$type}}{$chr}}{$pre_pos}}, @{${${$call{$type}}{$chr}}{$pos}});
                            my @new_len;
                            my $sum_len2 = 0;
                            foreach (@new_info){
                                my ($id, $pos1, $len) = split (/=/, $_);
                                push @new_len, $len;
                            }
                            map{$sum_len2 += $_} @new_len;
                            my $new_len = int ($sum_len2 / @new_len + 0.5);
                            delete ${${$call{$type}}{$chr}}{$pre_pos};
                            delete ${${$call{$type}}{$chr}}{$pos};
                            push @{${${$call{$type}}{$chr}}{$new_pos}}, @new_info;
                            $pre_pos = $new_pos;
                            $pre_len = $new_len;
                            next;
                        }
                    }
                }
            }
            $pre_len = $avelen;
            $pre_pos = $pos;
        }
    }
}

my %call_merge;

foreach my $type (keys %call){
    foreach my $chr (keys %{$call{$type}}){
        foreach my $pos (sort {$a <=> $b} keys %{${$call{$type}}{$chr}}){
            if (@{${${$call{$type}}{$chr}}{$pos}} == 1){
                my ($tool1, $tpos, $tlen, $tread, $subtype) = split (/=/, ${${${$call{$type}}{$chr}}{$pos}}[0]);
                ${${$call_merge{$type}}{$chr}}{$pos} = "$tool1:$tread=$tpos=$tlen=$tread" if ($type ne 'INS');
                ${${$call_merge{$type}}{$chr}}{$pos} = "$tool1:$tread=$tpos=$tlen=$tread=$subtype" if ($type eq 'INS');
            }
            else{
                my %tool_info;
                my %bp_rank;
                my %len_rank;
                my %preci_rank;
                my $ins_subtype = '';
                my @len;
                foreach (@{${${$call{$type}}{$chr}}{$pos}}){
                    my ($tool1, $pos1, $len1, $read1, $subtype) = split (/=/, $_);
                    my $tool2 = $tool1;
                    $tool2 =~ s/\.MEI$// if ($tool2 =~ /\.MEI$/);
                    $tool2 =~ s/\.VEI$// if ($tool2 =~ /\.VEI$/);
                    $tool2 =~ s/\.NUMT$// if ($tool2 =~ /\.NUMT$/);
                    if (!exists $tool_info{$tool2}){
                        $tool_info{$tool2} = $_;
                    }
                    else{
                        my $tread = $1 if ($_ =~ /=(\d+)$/) and ($type ne 'INS');
                        $tread = $1 if ($_ =~ /=(\d+)=\w+?$/) and ($type eq 'INS');
                        my $pre_tread = $1 if ($tool_info{$tool2} =~ /=(\d+)$/) and ($type ne 'INS');
                        $pre_tread = $1 if ($tool_info{$tool2} =~ /=(\d+)=\w+?$/) and ($type eq 'INS');
                        if ($tread > $pre_tread){
                            $tool_info{$tool2} = $_;
                        }
                    }
                    my $size = 'ALL';
                    if ($type ne 'INS'){
                        if ($len1 <= 1000){
                            $size = 'S';
                            if (($len1 <= 100) and ($type eq 'DEL')){
                                $size = 'SS';
                            }
                        }
                        elsif ($len1 <= 100000){
                            $size = 'M';
                        }
                        else{
                            $size = 'L';
                        }
                    }
                    else{
                        $ins_subtype .= "$subtype,";
                        push @len, $len1 if ($len1 > 20);
                    }
                    my $read_2 = $read1;
                    $read_2 = 30 if ($read1 > 30);
print STDERR "Warning No precision data: $tool2\t$type-$size\t$read_2\n" if (!exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
                    my $preci1 = 1;
                    $preci1 = ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2} if (exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
                    $preci_rank{$tool2} = $preci1;
                    my $bp_rank = 3;
                    my $len_rank = 3;
                    $bp_rank = ${$SDbp{$type}}{$tool2} if (exists ${$SDbp{$type}}{$tool2});
                    $len_rank = ${$SDlen{$type}}{$tool2} if (exists ${$SDlen{$type}}{$tool2});
                    $bp_rank{$tool2} = $bp_rank;
                    $len_rank{$tool2} = $len_rank if ($type ne 'INS');
                }
                $ins_subtype =~ s/,$//;
                my $top_tool = '';
                my $sec_tool = '';
                my $top_bp_tool = '';
                my $top_len_tool = '';
                my $top_preci = 0;
                my $sec_preci = 0;
                my $select_pos = 0;
                my $select_len = 0;
                my $select_read = 0;
                foreach my $tool (sort {$preci_rank{$b} <=> $preci_rank{$a}} keys %preci_rank){
                    if ($top_tool eq ''){
                        $top_tool = $tool;
                        $top_preci = $preci_rank{$tool};
                        next;
                    }
                    if ($sec_tool eq ''){
                        $sec_tool = $tool;
                        $sec_preci = $preci_rank{$tool};
                        last;
                    }
                }
                if ($top_preci < $sec_preci * 1.5){
                    if ($bp_rank{$top_tool} > $bp_rank{$sec_tool}){
                        $top_bp_tool = $sec_tool;
                    }
                    else{
                        $top_bp_tool = $top_tool;
                    }
                    if ($type ne 'INS'){
                        if ($len_rank{$top_tool} > $len_rank{$sec_tool}){
                            $top_len_tool = $sec_tool;
                        }
                        else{
                            $top_len_tool = $top_tool;
                        }
                    }
                }
                else{
                    $top_bp_tool = $top_tool;
                    $top_len_tool = $top_tool;
                }
                my ($tool1, $tpos1, $tlen1, $tread1) = split (/=/, $tool_info{$top_bp_tool});
                my ($tool2, $tpos2, $tlen2, $tread2) = split (/=/, $tool_info{$top_len_tool}) if ($type ne 'INS');
                $select_pos = $tpos1;
                $select_len = $tlen2 if ($type ne 'INS');
                $select_read = $tread1;
                if (@len > 0){
                    my $sumlen = 0;
                    map{$sumlen += $_} @len;
                    my $avelen = int ($sumlen / @len + 0.5);
                    $select_len = $avelen;
                }
                my $tool_set = '';
                foreach my $tool (sort {"\L$a" cmp "\L$b"} keys %tool_info){
                    my ($tool1, $tpos, $tlen, $tread) = split (/=/, $tool_info{$tool});
                    $tool_set .= "$tool1:$tpos:$tlen:$tread,";
                }
                $tool_set =~ s/,$//;
                ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read" if ($type ne 'INS');
                ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read=$ins_subtype" if ($type eq 'INS');
            }
        }
    }
}
%call = ();

foreach my $type (keys %call_merge){
    foreach my $chr (keys %{$call_merge{$type}}){
        my $count = 0;
        my %skip;
        while (1){
            my $delete = 0;
            $count ++;
            my $pre_pos = 0;
            my $pre_len = 0;
            my @merge;
            my $match_count = 0;
            foreach my $pos (sort {$a <=> $b} keys %{${$call_merge{$type}}{$chr}}){
                my ($tool, $tpos, $len) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos});
                my $end = $tpos + $len - 1;
                $end = $tpos if ($type eq 'INS');
                if ($pre_pos == 0){
                    push @merge, $pos;
                    $pre_pos = $tpos;
                    $pre_len = $len;
                    next;
                }
                if ($pre_pos > 0){
                    if ($type eq 'INS'){
                        if ($pos - $pre_pos <= $ins_sd){
                            push @merge, $pos;
                            next;
                        }
                        my %toolset;
                        my $select_pos = 0;
                        my $select_len = 0;
                        my $select_read = 0;
                        my $tool_set = '';
                        my $subtype_set = '';
                        foreach my $pos1 (@merge){
                            my ($toolset1, $tpos1, $len1, $read1, $subtype) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos1});
                            my @toolset = split (/,/, $toolset1);
                            my $toolnum = scalar @toolset;
                            $tool_set .= "$toolset1,";
                            $subtype_set .= "$subtype,";
                            push @{$toolset{$toolnum}}, "$toolset1=$pos1=$len1=$read1";
                        }
                        $tool_set =~ s/,$//;
                        $subtype_set =~ s/,$//;
                        foreach my $toolnum (sort {$b <=> $a} keys %toolset){
                            my $setnum = scalar @{$toolset{$toolnum}};
                            if ($setnum == 1){
                                my ($toolset1, $pos1, $len1, $read1) = split (/=/, ${$toolset{$toolnum}}[0]);
                                $select_pos = $pos1;
                                $select_read = $read1;
                                $select_len = $len1;
                            }
                            else{
                                my $count = 0;
                                my %preci_rank;
                                my @len;
                                foreach (@{$toolset{$toolnum}}){
                                    my ($toolset1, $pos1, $len1, $read1) = split (/=/, $_);
                                    my @toolset = split (/,/, $toolset1);
                                    my $size = 'ALL';
                                    if ($type ne 'INS'){
                                        if ($len1 <= 1000){
                                            $size = 'S';
                                            if (($len1 <= 100) and ($type eq 'DEL')){
                                                $size = 'SS';
                                            }
                                        }
                                        elsif ($len1 <= 100000){
                                            $size = 'M';
                                        }
                                        else{
                                            $size = 'L';
                                        }
                                    }
                                    else{
                                        push @len, $len1 if ($len1 > 20);
                                    }
                                    foreach my $tool1 (@toolset){
                                        my $tool2 = $tool1;
                                        my $read2 = 0;
                                        $read2 = $1 if ($tool1 =~ /:(\d+)$/);
                                        $read2 = $read1 if ($read2 == 0);
                                        $tool2 = $1 if ($tool2 =~ /(.+?):/);
                                        $tool2 =~ s/\.MEI$// if ($tool2 =~ /\.MEI$/);
                                        $tool2 =~ s/\.VEI$// if ($tool2 =~ /\.VEI$/);
                                        $tool2 =~ s/\.NUMT$// if ($tool2 =~ /\.NUMT$/);
                                        my $read_2 = $read2;
                                        $read_2 = 30 if ($read2 > 30);
                                        my $preci1 = 1;
                                        $preci1 = ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2} if (exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
                                        if (!exists $preci_rank{$count}){
                                            $preci_rank{$count} = $preci1;
                                        }
                                        else{
                                            if ($preci1 > $preci_rank{$count}){
                                                $preci_rank{$count} = $preci1;
                                            }
                                        }
                                    }
                                    $count ++;
                                }
                                my $select_count = 0;
                                foreach my $count1 (sort {$preci_rank{$b} <=> $preci_rank{$a}} keys %preci_rank){
                                    $select_count = $count1;
                                    last;
                                }
                                my ($stoolset1, $spos1, $slen1, $sread1) = split (/=/, ${$toolset{$toolnum}}[$select_count]);
                                if (@len > 0){
                                    my $sumlen = 0;
                                    map{$sumlen += $_} @len;
                                    my $avelen = int ($sumlen / @len + 0.5);
                                    $select_len = $avelen;
                                }
                                $select_pos = $spos1;
                                $select_read = $sread1;
                            }
                            last;
                        }
                        foreach my $pos1 (@merge){
                            delete ${${$call_merge{$type}}{$chr}}{$pos1};
                        }
                        if ($select_pos > 0){
                            $delete ++ if (@merge > 1);
                            ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read=$subtype_set";
                        }
                        @merge = ();
                        push @merge, $pos;
                    }
                    else{
                        my $pre_end = $pre_pos + $pre_len - 1;
                        if ($pos < $pre_end){
                            my $overlap = $pre_end - $pos + 1;
                            $overlap = $len if ($end < $pre_end);
                            if (($overlap >= $len * $min_overlap_ratio3) and ($overlap >= $pre_len * $min_overlap_ratio3)){
                                push @merge, $pos;
                                $match_count ++;
                                next;
                            }
                            elsif (($overlap >= $len * $min_overlap_ratio2) or ($overlap >= $pre_len * $min_overlap_ratio2)){
                                my ($toolset1) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos});
                                my @toolset1 = split (/,/, $toolset1);
                                my $tool_num = scalar @toolset1;
                                my $pre_tool_num = 0;
                                foreach (@merge){
                                    my ($pre_toolset) = split (/=/, ${${$call_merge{$type}}{$chr}}{$_});
                                    my @pre_toolset = split (/,/, $pre_toolset);
                                    $pre_tool_num += @pre_toolset;
                                }
                                if ($overlap >= $len * $min_overlap_ratio2){
                                    if (($pre_len <= 1000) and ($pre_len / $len <= 3)){
                                        push @merge, $pos;
                                        $match_count ++;
                                        next;
                                    }
                                    elsif (($tool_num == 1) and ($pre_tool_num > 2) and ($len <= 1000)){
                                        delete ${${$call_merge{$type}}{$chr}}{$pos};
                                        $delete ++;
                                        $match_count ++;
                                        next;
                                    }
                                    elsif (($pre_tool_num == 1) and (2 < $tool_num) and ($pre_len <= 1000)){
                                        foreach my $pos1 (@merge){
                                            delete ${${$call_merge{$type}}{$chr}}{$pos1};
                                        }
                                        $delete ++;
                                        @merge = ();
                                        push @merge, $pos;
                                        $pre_pos = $tpos;
                                        $pre_len = $len;
                                        $match_count ++;
                                        next;
                                    }
                                    elsif (exists ${$skip{$type}}{$pos}){
                                    }
                                    else{
                                        ${$skip{$type}}{$pos} = 1;
                                        push @merge, $pos;
                                        next;
                                    }
                                }
                                elsif ($overlap >= $pre_len * $min_overlap_ratio2){
                                    if (($len <= 1000) and ($len / $pre_len <= 3)){
                                        push @merge, $pos;
                                        $match_count ++;
                                        next;
                                    }
                                    elsif (($pre_tool_num == 1) and (2 < $tool_num) and ($pre_len <= 1000)){
                                        foreach my $pos1 (@merge){
                                            delete ${${$call_merge{$type}}{$chr}}{$pos1};
                                        }
                                        $delete ++;
                                        @merge = ();
                                        push @merge, $pos;
                                        $pre_pos = $tpos;
                                        $pre_len = $len;
                                        $match_count ++;
                                        next;
                                    }
                                    elsif (($tool_num == 1) and ($pre_tool_num > 2) and ($len <= 1000)){
                                        delete ${${$call_merge{$type}}{$chr}}{$pos};
                                        $delete ++;
                                        next;
                                    }
                                    elsif (exists ${$skip{$type}}{$pos}){
                                    }
                                    else{
                                        ${$skip{$type}}{$pos} = 1;
                                        push @merge, $pos;
                                        next;
                                    }
                                }
                            }
                        }
                        
                        my %toolset;
                        my $top_tnum = 0;
                        my $total_tools = 0;
                        my $total_top_tnum = 0;
                        my $select_pos = 0;
                        my $select_len = 0;
                        my $select_read = 0;
                        my @pos2;
                        my @len2;
                        my @read2;
                        my $tool_set = '';
                        foreach my $pos1 (@merge){
                            next if (exists ${$skip{$type}}{$pos1});
                            my ($toolset1, $tpos1, $len1, $read1) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos1});
                            my @toolset = split (/,/, $toolset1);
                            my $new_toolset = '';
                            foreach (@toolset){
                                if ($_ !~ /:\d+:/){
                                    my $tool = $1 if ($_ =~ /(.+?):/);
                                    $new_toolset .= "$tool:$tpos1:$len1:$read1,";
                                }
                                else{
                                    $new_toolset .= "$_,";
                                }
                            }
                            $new_toolset =~ s/,$//;
                            my $toolnum = scalar @toolset;
                            $tool_set .= "$new_toolset,";
                            push @{$toolset{$toolnum}}, "$new_toolset=$tpos1=$len1=$read1";
                            $top_tnum = $toolnum if ($toolnum > $top_tnum);
                            $total_tools += $toolnum;
                            push @pos2, $tpos1;
                            push @len2, $len1;
                            push @read2, $read1;
                        }
                        $tool_set =~ s/,$//;
                        foreach (@{$toolset{$top_tnum}}){
                            my @tset = split (/,/, $_);
                            $total_top_tnum += @tset;
                        }
                        
                        if ($total_tools == 1){
                            my ($toolset1, $pos1, $len1, $read1) = split (/=/, ${$toolset{$top_tnum}}[0]);
                            $select_pos = $pos1;
                            $select_len = $len1;
                            $select_read = $read1;
                        }
                        elsif (($top_tnum > 1) and ($total_top_tnum >= $total_tools * 0.5)){
                            my $count = 0;
                            my %preci_rank;
                            my %pos_len;
                            foreach (@{$toolset{$top_tnum}}){
                                my ($toolset1, $pos1, $len1, $read1) = split (/=/, $_);
                                $pos_len{$pos1} = $len1;
                            }
                            my $large_len = 0;
                            my $small_pos = 0;
                            my $small_len = 0;
                            my $delete_pos = 0;
                            
                            foreach my $posL (sort {$pos_len{$a} <=> $pos_len{$b}} keys %pos_len){
                                $small_pos = $posL if ($small_pos == 0);
                                $small_len = $pos_len{$posL} if ($small_len == 0);
                                $large_len = $pos_len{$posL};
                            }
                            if (($small_len > 0) and ($large_len / $small_len > 2)){
                                $delete_pos = $small_pos;
                            }
                            foreach (@{$toolset{$top_tnum}}){
                                my ($toolset1, $pos1, $len1, $read1) = split (/=/, $_);
                                if ($pos1 != $delete_pos){
                                    my @toolset = split (/,/, $toolset1);
                                    my $size = 'ALL';
                                    if ($type ne 'INS'){
                                        if ($len1 <= 1000){
                                            $size = 'S';
                                            if (($len1 <= 100) and ($type eq 'DEL')){
                                                $size = 'SS';
                                            }
                                        }
                                        elsif ($len1 <= 100000){
                                            $size = 'M';
                                        }
                                        else{
                                            $size = 'L';
                                        }
                                    }
                                    foreach my $tool1 (@toolset){
                                        my $tool2 = $tool1;
                                        my $read2 = 0;
                                        $read2 = $1 if ($tool1 =~ /:(\d+)$/);
                                        $read2 = $read1 if ($read2 == 0);
                                        my $read_2 = $read2;
                                        $read_2 = 30 if ($read2 > 30);
                                        $tool2 = $1 if ($tool2 =~ /(.+?):/);
                                        my $preci1 = 1;
                                        $preci1 = ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2} if (exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
print STDERR "Warning No precision data: $type-$size\t$tool1\t$tool2\t$read_2\n" if (!exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
                                        if (!exists $preci_rank{$count}){
                                            $preci_rank{$count} = $preci1;
                                        }
                                        else{
                                            if ($preci1 > $preci_rank{$count}){
                                                $preci_rank{$count} = $preci1;
                                            }
                                        }                                               
                                    }
                                }
                                $count ++;
                            }
                            my $select_count = 0;
                            foreach my $count1 (sort {$preci_rank{$b} <=> $preci_rank{$a}} keys %preci_rank){
                                $select_count = $count1;
                                last;
                            }
                            my ($stoolset1, $spos1, $slen1, $sread1) = split (/=/, ${$toolset{$top_tnum}}[$select_count]);
                            $select_pos = $spos1;
                            $select_len = $slen1;
                            $select_read = $sread1;
                        }
                        else{
                            my %len_pos;
                            my $count1 = 0;
                            foreach (@len2){
                                $len_pos{$_} = "$pos2[$count1]=$read2[$count1]";
                                $count1 ++;
                            }
                            my $med_num = int (scalar (@pos2) / 2);
                            $count1 = 0;
                            foreach my $len1 (sort {$a <=> $b} keys %len_pos){
                                if ($count1 == $med_num){
                                    $select_len = $len1;
                                    my ($pos1, $read1) = split (/=/, $len_pos{$len1});
                                    $select_pos = $pos1;
                                    $select_read = $read1;
                                    last;
                                }
                                $count1 ++;
                            }
                        }
                        my $skip_num = 0;
                        foreach my $pos1 (@merge){
                            if (exists ${$skip{$type}}{$pos1}){
                                $skip_num ++;
                                next;
                            }
                            delete ${${$call_merge{$type}}{$chr}}{$pos1};
                        }
                        if ($select_pos > 0){
                            $delete ++ if (@merge - $skip_num > 1);
                            ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read";
                        }
                        @merge = ();
                        push @merge, $pos;
                        $pre_pos = $select_pos;
                        $pre_len = $select_len;
                        $match_count = 0;
                        %skip = ();
                    }
                }
                $pre_pos = $tpos;
                $pre_len = $len;
            }
            if (@merge > 0){
                my %toolset;
                my $top_tnum = 0;
                my $total_tools = 0;
                my $total_top_tnum = 0;
                my $select_pos = 0;
                my $select_len = 0;
                my $select_read = 0;
                my @pos2;
                my @len2;
                my @read2;
                my $tool_set = '';
                my $subtype_set = '';
                foreach my $pos1 (@merge){
                    next if (exists ${$skip{$type}}{$pos1}) and (@merge > 1);
                    my ($toolset1, $tpos1, $len1, $read1, $subtype) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos1});
                    my @toolset = split (/,/, $toolset1);
                    my $toolnum = scalar @toolset;
                    my $new_toolset = '';
                    foreach (@toolset){
                        if ($_ !~ /:\d+:/){
                            my $tool = $1 if ($_ =~ /(.+?):/);
                            $new_toolset .= "$tool:$tpos1:$len1:$read1,";
                        }
                        else{
                            $new_toolset .= "$_,";
                        }
                    }
                    $new_toolset =~ s/,$//;
                    $tool_set .= "$new_toolset,";
                    $subtype_set .= "$subtype," if ($type eq 'INS');
                    push @{$toolset{$toolnum}}, "$new_toolset=$pos1=$len1=$read1";
                    $top_tnum = $toolnum if ($toolnum > $top_tnum);
                    $total_tools += $toolnum;
                    push @pos2, $tpos1;
                    push @len2, $len1;
                    push @read2, $read1;
                }
                $tool_set =~ s/,$//;
                $subtype_set =~ /,$/ if ($type eq 'INS');
                foreach (@{$toolset{$top_tnum}}){
                    my @tset = split (/,/, $_);
                    $total_top_tnum += @tset;
                }
                
                if ($total_tools == 1){
                    my ($toolset1, $pos1, $len1, $read1) = split (/=/, ${$toolset{$top_tnum}}[0]);
                    $select_pos = $pos1;
                    $select_len = $len1;
                    $select_read = $read1;
                }
                elsif (($top_tnum > 1) and ($total_top_tnum >= $total_tools * 0.5)){
                    my $count = 0;
                    my %preci_rank;
                    my %pos_len;
                    my $delete_pos = 0;
                    if ($type ne 'INS'){
                        foreach (@{$toolset{$top_tnum}}){
                            my ($toolset1, $pos1, $len1, $read1) = split (/=/, $_);
                            $pos_len{$pos1} = $len1;
                        }
                        my $large_len = 0;
                        my $small_pos = 0;
                        my $small_len = 0;
                        foreach my $posL (sort {$pos_len{$a} <=> $pos_len{$b}} keys %pos_len){
                            $small_pos = $posL if ($small_pos == 0);
                            $small_len = $pos_len{$posL} if ($small_len == 0);
                            $large_len = $pos_len{$posL};
                        }
                        if (($small_len > 0) and ($large_len / $small_len > 2)){
                            $delete_pos = $small_pos;
                        }
                    }
                    foreach (@{$toolset{$top_tnum}}){
                        my ($toolset1, $pos1, $len1, $read1) = split (/=/, $_);
                        if ($pos1 != $delete_pos){
                            my @toolset = split (/,/, $toolset1);
                            my $size = 'ALL';
                            if ($type ne 'INS'){
                                if ($len1 <= 1000){
                                    $size = 'S';
                                    if (($len1 <= 100) and ($type eq 'DEL')){
                                        $size = 'SS';
                                    }
                                }
                                elsif ($len1 <= 100000){
                                    $size = 'M';
                                }
                                else{
                                    $size = 'L';
                                }
                            }
                            foreach my $tool1 (@toolset){
                                my $tool2 = $tool1;
                                my $read2 = 0;
                                $read2 = $1 if ($tool1 =~ /:(\d+)$/);
                                $read2 = $read1 if ($read2 == 0);
                                my $read_2 = $read2;
                                $read_2 = 30 if ($read2 > 30);
                                $tool2 = $1 if ($tool2 =~ /(.+?):/);
                                $tool2 =~ s/\.MEI$// if ($tool2 =~ /\.MEI$/);
                                $tool2 =~ s/\.VEI$// if ($tool2 =~ /\.VEI$/);
                                $tool2 =~ s/\.NUMT$// if ($tool2 =~ /\.NUMT$/);
                                my $preci1 = 1;
                                $preci1 = ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2} if (exists ${${${$tool_preci{$type}}{$size}}{$tool2}}{$read_2});
                                if (!exists $preci_rank{$count}){
                                    $preci_rank{$count} = $preci1;
                                }
                                else{
                                    if ($preci1 > $preci_rank{$count}){
                                        $preci_rank{$count} = $preci1;
                                    }
                                } 
                            }
                        }
                        $count ++;
                    }
                    my $select_count = 0;
                    foreach my $count1 (sort {$preci_rank{$b} <=> $preci_rank{$a}} keys %preci_rank){
                        $select_count = $count1;
                        last;
                    }
                    my ($stoolset1, $spos1, $slen1, $sread1) = split (/=/, ${$toolset{$top_tnum}}[$select_count]);
                    $select_pos = $spos1;
                    $select_read = $sread1;
                    $select_len = $slen1 if ($type ne 'INS');
                }
                else{
                    my %len_pos;
                    my $count1 = 0;
                    foreach (@len2){
                        $len_pos{$_} = "$pos2[$count1]=$read2[$count1]";
                        $count1 ++;
                    }
                    my $med_num = int (scalar (@pos2) / 2);
                    $count1 = 0;
                    foreach my $len1 (sort {$a <=> $b} keys %len_pos){
                        if ($count1 == $med_num){
                            $select_len = $len1;
                            my ($pos1, $read1) = split (/=/, $len_pos{$len1});
                            $select_pos = $pos1;
                            $select_read = $read1;
                            last;
                        }
                        $count1 ++;
                    }
                }
                foreach my $pos1 (@merge){
                    delete ${${$call_merge{$type}}{$chr}}{$pos1};
                }
                if ($select_pos > 0){
                    $delete ++ if (@merge > 1);
                    ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read" if ($type ne 'INS');
                    ${${$call_merge{$type}}{$chr}}{$select_pos} = "$tool_set=$select_pos=$select_len=$select_read=$subtype_set" if ($type eq 'INS');
                }
                @merge = ();
            }
            last if ($delete == 0);
        }
    }
}

my %vcf;
my @ins_subtype = ('ALU', 'L1', 'SVA', 'HERVK', 'VEI', 'NUMT');

foreach my $type (keys %call_merge){
    foreach my $chr (keys %{$call_merge{$type}}){
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
        foreach my $pos (keys %{${$call_merge{$type}}{$chr}}){
            my ($toolset, $pos1, $len, $read, $subtypeset) = split (/=/, ${${$call_merge{$type}}{$chr}}{$pos});
            my %toolset;
            my %toolinfo;
            my %subtype;
            my @toolset = split (/,/, $toolset);
            my $new_toolset = '';
            foreach my $tool1 (@toolset){
                my $tool2 = $1 if ($tool1 =~ /(.+?):/);
                my $read = $1 if ($tool1 =~ /:(\d+)$/);
                ${$toolset{$tool2}}{$read} ++;
                ${$toolinfo{$tool2}}{$read} = $tool1;
            }
            foreach my $tool (sort {"\L$a" cmp "\L$b"} keys %toolset){
                foreach my $read (sort {$b <=> $a} keys %{$toolset{$tool}}){
                    $new_toolset .= "${$toolinfo{$tool}}{$read},";
                    last;
                }
            }
            $new_toolset =~ s/,$//;
            my $type2 = '';
            if ($type eq 'INS'){
                my @subtype = split (/,/, $subtypeset);
                my $MEI_flag = 0;
                foreach my $stype (@subtype){
                    next if ($stype eq 'INS');
                    $stype = 'L1' if ($stype eq 'LINE1');
                    $subtype{$stype} = 1;
                }
                if (scalar keys %subtype > 0){
                    foreach my $subtype (sort keys %subtype){
                        $type2 .= "$subtype,";
                    }
                    $type2 =~ s/,$//;
                }
                else{
                    $type2 = $type;
                }
            }
            ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;TOOLS=$new_toolset\n" if ($type ne 'INS');
            ${${$vcf{$chr02d}}{$pos}}{$type} = "$chr\t$pos\t$type2\t.\t.\t.\tPASS\tSVTYPE=$type;SVLEN=$len;READS=$read;TOOLS=$new_toolset\n" if ($type eq 'INS');
        }
    }
}

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (sort keys %{${$vcf{$chr}}{$pos}}){
            print ${${$vcf{$chr}}{$pos}}{$type};
        }
    }
}

