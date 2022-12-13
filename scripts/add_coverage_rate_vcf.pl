#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use File::Basename;

# add SV coverage rate (DPR) and split read rate (SR) to an input vcf file

my $vcf = '';
my $cov_dir = '';
my $flank_len = 1000;
my $flank_rate = 0.1;
my $bin_size = 50;
my $sample_name = '';
my $sample_dir = '';
my $out_prefix = '';
my $non_human = 0;
my $build = 37;
my $ref_index = '';
my $gap_bed = '';

my $min_del_len = 500;
my $min_del_len2 = 10000;
my $min_dup_len = 200;
my $min_dup_len2 = 1000;
my $min_deldup_len = 1000;
my $max_del_dprate = 0.85;
my $max_del_dprate2 = 0.7;
my $max_del_dprate3 = 0.9;
my $min_dup_dprate = 1.35;
my $min_dup_dprate3 = 1.15;
my $max_del_dsr = 0.2;
my $max_dup_dsr = 0.1;

my $data_dir = "$Bin/../Data";
$data_dir = "$ENV{MOPLINE_DIR}/Data" if (exists $ENV{MOPLINE_DIR});

my $min_hom_ins_split_rate = 0.97;
my $min_hom_inv_split_rate = 1.15;

my $help;

GetOptions(
    'vcf|v=s' => \$vcf,
    'sample_name|sn=s' => \$sample_name,
    'sample_dir|sd=s' => \$sample_dir,
    'cov_dname|cd=s' => \$cov_dir,
    'flank|f=i' => \$flank_len,
    'bin=i' => \$bin_size,
    'non_humna|nh=i' => \$non_human,
    'build|b=s' => \$build,
    'ref_index|ri=s' => \$ref_index,
    'gap_bed|gap=s' => \$gap_bed,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  add_coverage_rate_vcf.pl -v <vcf file> -sn <sample name> -sd <sample directory> 
  (-ri <reference index file> -gap <gap bed file> -nh 1 if sample is a non-human species)
  (Add alignment information (DR/DS/SR tags) field of a specified MOP-merged vcf file with cov files that were generated with add_GT_DPR_vcf_thread.pl)

  Options:
   --vcf or -v <STR>        SV vcf file path to be genotyped
   --sample_name or -sn <STR>  sample name [mandatory]
   --sample_dir or -sd <STR>   sample directory (optional, sample_dir is regarded as sample_name if not specified)
   --ref_index or -ri <STR> reference fasta index file [mandatory if non-human species]
   --gap_file or -gap <STR> gap bed file indicating reference gap regions 
                            (1st column: chr-name, 2nd and 3rd columns: start and end positions of a gap) [optional, automatically selected for human]
   --cov_dname or -cd <STR> directory name containing coverage data [default: Cov]
   --flank or -f <INT>      flanking size of SV (DEL and DUP) breakpoints to be analyzed [default: 1000]
   --bin <INT>              bin size of coverage data [default: 50]
   --non_human or -nh <INT> sample is non-human species (human: 0. non-human: 1) [default: 0]
   --build or -b <STR>      human reference build (GRCh37, GRCh38, T2T-CHM13) (37, 38, or T2T) [default: 37]
   --prefix or -p <STR>     output prefix
   --help or -h             output help message
   
=cut

die "SV vcf file to genotype is not specified: \n" if ($vcf eq '');
die "Sample name is not specified: \n" if ($sample_name eq '');
die "Reference fasta index file is not specified:\n" if ($non_human == 1) and ($ref_index eq '');

$sample_dir = $sample_name if ($sample_dir eq '');

if (($non_human == 0) and ($ref_index eq '')){
    $ref_index = "$data_dir/hs37.fa.fai" if ($build eq '37');
    $ref_index = "$data_dir/hs38.fa.fai" if ($build eq '38');
    $ref_index = "$data_dir/chm13v2.0.fa.fai" if ($build eq 'T2T');
}
if (($non_human == 0) and ($gap_bed eq '')){
    $gap_bed = "$data_dir/gap.bed" if ($build eq '37');
    $gap_bed = "$data_dir/gap.b38.bed" if ($build eq '38');
}

my %cov;
my %cov2;
my %split1;
my %split2;
my %cov_sum;
my %cov_sum2;
my %cov_num;
my %del;
my %dup;
my %chrlen;
my %ref_index;
my @chr;
my %gap;
my %hq;

open (FILE, "$ref_index") or die "$ref_index is not fuond: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $len) = split (/\t/, $line);
    next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXY]+$/);
    $chrlen{$chr} = $len - $bin_size;
    push @chr, $chr;
    $ref_index{$chr} = $len;
}
close (FILE);

my $Mbin_size = 1000000;

if ($gap_bed ne ''){
    open (FILE, $gap_bed) or die "$gap_bed is not found: $!\n";
    while (my $line = <FILE>){
        next if ($line =~ /^#/);
        my ($chr, $pos, $end) = split (/\t/, $line);
        my $pos_res = $pos % $bin_size;
        my $pos2 = $pos - $pos_res;
        while (1){
            ${$gap{$chr}}{$pos2} = "$pos=$end";
            $pos2 += $bin_size;
            last if ($pos2 > $end);
        }
    }
    close (FILE);
}

my $pre_chr = '';
my $pre_pos = 0;
my $pre_end = 0;
my $pre_type = '';
my @header;
open (FILE, $vcf) or die "$vcf is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
        next if ($line =~ /^#CHROM/);
        push @header, $line;
    }
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    my $end = $pos + $len - 1;
    if ($type =~ /DEL|DUP/){
        my $pos_res = $pos % $bin_size;
        my $pos2 = $pos - $pos_res;
        while (1){
            ${$del{$chr}}{$pos2} = "$pos=$end" if ($type eq 'DEL');
            ${$dup{$chr}}{$pos2} = "$pos=$end" if ($type eq 'DUP');
            $pos2 += $bin_size;
            last if ($pos2 > $end);
        }
        if ($pre_chr ne $chr){
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_type = $type;
            $pre_chr = $chr; 
        }
        if ($end > $pre_end){
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_type = $type;
        }
    }
}
close (FILE);

my $out2 = "$vcf.cov";
$out2 = "$out_prefix.addDPR.vcf" if ($out_prefix ne '');

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

my %cov_filter;
my %cnv_num;

my @header2;
foreach my $chr (@chr){
    my $len = $ref_index{$chr};
    my $line = "##contig=<ID=$chr,length=$len>";
    push @header2, $line;
}
my $ref = basename ($ref_index);
$ref = $1 if ($ref =~ /(.+)\.fas*t*a*\.fai/);

open (OUT, "> $out2");
foreach my $hline (@header){
    if ($hline =~ /fileDate=/){
        $hline =~ s/fileDate=\d+/fileDate=$time/;
    }
    elsif ($hline =~ /reference=/){
        $hline =~ s/reference=/reference=$ref/;
    }
    print OUT "$hline\n";
}
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
print OUT "##FORMAT=<ID=VP,Number=1,Type=Integer,Description=\"Position of SV for a sample\">\n";
print OUT "##FORMAT=<ID=VL,Number=1,Type=Integer,Description=\"Length of SV for a sample\">\n";
print OUT "##FORMAT=<ID=DR,Number=1,Type=Float,Description=\"Ratio of read depth in CNV region to that in its flanking regions for a sample\">\n";
print OUT "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Ratio of inconsistent DR segments (> 0.8 for DEL, < 1.1 for DUP) in DEL and DUP\">\n";
print OUT "##FORMAT=<ID=SR,Number=1,Type=Float,Description=\"Ratio of breakpoints of soft-clipped reads in the read depth at the 50- or 100-bp window\">\n";
foreach my $hline (@header2){
    print OUT "$hline\n";
}
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";

open (FILE, $vcf) or die "$vcf is not found: $!\n";
$pre_chr = '';
my $cov_ave = 0;
my $cov_ave2 = 0;
my $cov_num = 0;
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    die "No match with $chr between the reference index file and the input vcf file\n" if (!exists $chrlen{$chr});
    if ($pre_chr ne $chr){
        # collect coverage and split read data for a chromosome from a cov file
        my $cov_sum = 0;
        my $cov_sum2 = 0;
        %cov = ();
        %cov2 = ();
        %split1 = ();
        %split2 = ();
        my $cov_chr = "$sample_dir/$cov_dir/$sample_name.$chr.cov.gz";
        $cov_chr = "$sample_dir/$cov_dir/$sample_name.chr$chr.cov.gz" if (!-f $cov_chr);
        open (FILE2, $cov_chr) or die "$cov_chr is not found: $!\n" if ($cov_chr !~ /\.gz$/);
        open (FILE2, "gzip -dc $cov_chr |") or die "$cov_chr is not found: $!\n" if ($cov_chr =~ /\.gz$/);
        while (my $line2 = <FILE2>){
            chomp $line2;
            my ($chr, $pos, $cov, $cov2, $split1, $split2) = split (/\t/, $line2);
            $cov{$pos} = $cov;
            $cov2{$pos} = $cov2 if (defined $cov2);
            $split1{$pos} = $split1 if (defined $split1);
            $split2{$pos} = $split2 if (defined $split2);
            $cov_sum += $cov;
            $cov_sum2 += $cov2 if (defined $cov2);
            $cov_num ++;
        }
        close (FILE2);
        $cov_ave = int ($cov_sum / $cov_num * 10 + 0.5) / 10 if ($cov_num > 0);
        $cov_ave2 = int ($cov_sum2 / $cov_num * 10 + 0.5) / 10 if ($cov_num > 0);
    }
    $pre_chr = $chr;
    my $pos = $line[1];
    if (@line > 8){
        $line[8] = 'GT:GQ:VP:VL:DR:DS:SR';
        my @info = split (/:/, $line[9]);
        if (@info >= 5){
            $line[9] = "$info[0]:$info[1]:$info[2]:$info[3]";
        }
=pod
        if ($type =~ /INS|INV/){
            $line[9] .= ':1:0';
            print OUT join ("\t", @line), "\n";
            next;
        }
=cut
    }
=pod
    elsif ($type =~ /INS|INV/){
        $line[7] .= ";DR=1;DS=0";
        print OUT join ("\t", @line), "\n";
        next;
    }
=cut
    # calculate SR for any type of SV at a called SV site using the split read data in 150 bp (50 bp window * 3) regions around the breakpoints
    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    my $end = $pos + $len - 1;
    $end = $pos if ($type eq 'INS');
    my $start2 = $pos - ($pos % $bin_size);
    my $start1 = $start2 - $bin_size;
    my $start3 = $start2 + $bin_size;
    my $start4 = $start3 + $bin_size;
    my %split_num;
    my %split2_num;
    my $split_rate = 0;
    my $cov_rate = 1;
    my $incons_cov_rate = 0;
    if ($type =~ /INS/){
        $split_num{$start1} = $split1{$start1} if (exists $split1{$start1});
        $split2_num{$start1} = $split2{$start1} if (exists $split2{$start1});
        $split_num{$start2} = $split1{$start2} if (exists $split1{$start2});
        $split2_num{$start2} = $split2{$start2} if (exists $split2{$start2});
        $split_num{$start3} = $split1{$start3} if (exists $split1{$start3});
        $split2_num{$start3} = $split2{$start3} if (exists $split2{$start3});
        $split2_num{$start4} = $split2{$start4} if (exists $split2{$start4});
    }
    elsif ($type =~ /INV/){
        $split_num{$start1} = $split1{$start1} if (exists $split1{$start1});
        $split_num{$start1} += $split2{$start1} if (exists $split2{$start1});
        $split_num{$start2} = $split1{$start2} if (exists $split1{$start2});
        $split_num{$start2} += $split2{$start2} if (exists $split2{$start2});
        $split_num{$start3} = $split1{$start3} if (exists $split1{$start3});
        $split_num{$start3} += $split2{$start3} if (exists $split2{$start3});
    }
    elsif ($type eq 'DEL'){
        $split_num{$start1} = $split2{$start1} if (exists $split2{$start1});
        $split_num{$start2} = $split2{$start2} if (exists $split2{$start2});
        $split_num{$start3} = $split2{$start3} if (exists $split2{$start3});
    }
    elsif ($type eq 'DUP'){
        $split_num{$start1} = $split1{$start1} if (exists $split1{$start1});
        $split_num{$start2} = $split1{$start2} if (exists $split1{$start2});
        $split_num{$start3} = $split1{$start3} if (exists $split1{$start3});
    }
    my $top_num = 0;
    my $top_pos = 0;
    my $sec_num = 0;
    my $sec_pos = 0;
    foreach my $spos (sort {$split_num{$b} <=> $split_num{$a}} keys %split_num){
        if ($top_pos == 0){
            $top_pos = $spos;
            $top_num = $split_num{$spos};
            next;
        }
        elsif (($sec_pos == 0) and (abs ($top_pos - $spos) <= $bin_size)){
            $sec_pos = $spos;
            $sec_num = $split_num{$spos};
        }
        last;
    }
    
    if ($type eq 'INS'){
        my $top_num_2 = 0;
        my $top_pos_2 = 0;
        my $sec_num_2 = 0;
        my $sec_pos_2 = 0;
        my %ins_bp;
        foreach my $spos (sort {$split2_num{$b} <=> $split2_num{$a}} keys %split2_num){
            if ($top_pos_2 == 0){
                $top_pos_2 = $spos;
                $top_num_2 = $split2_num{$spos};
                next;
            }
            elsif (($sec_pos_2 == 0) and (abs ($top_pos_2 - $spos) <= $bin_size)){
                $sec_pos_2 = $spos;
                $sec_num_2 = $split2_num{$spos};
            }
            last;
        }
        
        if (($top_pos > 0) and (exists $cov{$top_pos})){
            $ins_bp{$top_pos} = $top_num;
        }
        if (($sec_num >= $top_num * 0.5) and ($sec_pos > 0) and (exists $cov{$sec_pos})){
            $ins_bp{$sec_pos} += $sec_num;
        }
        if (($top_pos_2 > 0) and (exists $cov{$top_pos_2})){
            $ins_bp{$top_pos_2} += $top_num_2;
        }
        if (($sec_num_2 >= $top_num_2 * 0.5) and ($sec_pos_2 > 0) and (exists $cov{$sec_pos_2})){
            $ins_bp{$sec_pos_2} += $sec_num_2;
        }
        my $total_bin = scalar keys %ins_bp;
        if ($total_bin > 0){
            my $sum = 0;
            my $sum_cov = 0;
            foreach my $spos (keys %ins_bp){
                $sum += $ins_bp{$spos};
                $sum_cov += $cov{$spos};
            }
            my $ave_cov = $sum_cov / $total_bin;
            $split_rate = int ($sum / $ave_cov * 100 + 0.5) / 100 if ($ave_cov > 0);
        }
    }
    else{
        if (($top_pos > 0) and (exists $cov{$top_pos})){
            if (($sec_num >= $top_num * 0.5) and ($sec_pos > 0)){
                my $ave_cov = ($cov{$top_pos} + $cov{$sec_pos}) * 0.5 if (exists $cov{$sec_pos});
                $ave_cov = $cov{$top_pos} if (!exists $cov{$sec_pos});
                $split_rate = int (($top_num + $sec_num) / $ave_cov * 100 + 0.5) / 100 if ($ave_cov > 0);
            }
            else{
                $split_rate = int ($top_num / $cov{$top_pos} * 100 + 0.5) / 100 if ($cov{$top_pos} > 0);
            }
        }
    }
    if ($type ne 'INS'){
        my $end2 = $end - ($end % $bin_size);
        my $end1 = $end2 - $bin_size;
        my $end3 = $end2 + $bin_size;
        my $split_rate2 = 0;
        my %split_num2;
        if ($type eq 'DEL'){
            $split_num2{$end1} = $split1{$end1} if (exists $split1{$end1});
            $split_num2{$end2} = $split1{$end2} if (exists $split1{$end2});
            $split_num2{$end3} = $split1{$end3} if (exists $split1{$end3});
        }
        elsif ($type eq 'DUP'){
            $split_num2{$end1} = $split2{$end1} if (exists $split2{$end1});
            $split_num2{$end2} = $split2{$end2} if (exists $split2{$end2});
            $split_num2{$end3} = $split2{$end3} if (exists $split2{$end3});
        }
        elsif ($type eq 'INV'){
            $split_num2{$end1} = $split1{$end1} if (exists $split1{$end1});
            $split_num2{$end1} += $split2{$end1} if (exists $split2{$end1});
            $split_num2{$end2} = $split1{$end2} if (exists $split1{$end2});
            $split_num2{$end2} += $split2{$end2} if (exists $split2{$end2});
            $split_num2{$end3} = $split1{$end3} if (exists $split1{$end3});
            $split_num2{$end3} += $split2{$end3} if (exists $split2{$end3});
        }
        my $top_num2 = 0;
        my $top_pos2 = 0;
        my $sec_num2 = 0;
        my $sec_pos2 = 0;
        foreach my $spos (sort {$split_num2{$b} <=> $split_num2{$a}} keys %split_num2){
            if ($top_pos2 == 0){
                $top_pos2 = $spos;
                $top_num2 = $split_num2{$spos};
                next;
            }
            elsif (($sec_pos2 == 0) and (abs ($top_pos2 - $spos) <= $bin_size)){
                $sec_pos2 = $spos;
                $sec_num2 = $split_num2{$spos};
            }
            last;
        }
        if (($top_pos2 > 0) and (exists $cov{$top_pos2})){
            if (($sec_num2 >= $top_num2 * 0.5) and ($sec_pos2 > 0)){
                my $ave_cov = ($cov{$top_pos2} + $cov{$sec_pos2}) * 0.5 if (exists $cov{$sec_pos2});
                $ave_cov = $cov{$top_pos2} if (!exists $cov{$sec_pos2});
                $split_rate2 = int (($top_num2 + $sec_num2) / $ave_cov * 100 + 0.5) / 100 if ($ave_cov > 0);
            }
            else{
                $split_rate2 = int ($top_num2 / $cov{$top_pos2} * 100 + 0.5) / 100 if ($cov{$top_pos2} > 0);
            }
        }
        my $ave_split_rate = 0;
        if ($split_rate >= $split_rate2 * 2){
            $ave_split_rate = $split_rate;
        }
        elsif ($split_rate2 >= $split_rate * 2){
            $ave_split_rate = $split_rate2;
        }
        else{
            $ave_split_rate = int (($split_rate + $split_rate2) * 0.5 * 100 + 0.5) / 100;
        }
        $split_rate = $ave_split_rate;
    }
    if ($type =~ /DEL|DUP/){
        # calculate DPR for DEL or DUP at a called SV site using the covergae between the breakpoints and the coverage outsides the breakpoints
        my $CNVnator_flag = 0;
        my $pos2 = $pos - 1;
        $CNVnator_flag = 1 if ($pos2 % 1000 <= 1) and ($len >= 2000);
        my $left_end = $pos - 100;
        $left_end -= 400 if ($CNVnator_flag == 1);
        my $left_res = $left_end % $bin_size;
        $left_end -= $left_res;
        my $left_start = $left_end - $flank_len;
        $left_start = $left_end - $len * $flank_rate if ($len > 10000);
        my $right_start = $end + 100;
        $right_start += 400 if ($CNVnator_flag == 1);
        my $right_res = $right_start % $bin_size;
        $right_start += $bin_size - $right_res;
        my $right_end = $right_start + $flank_len;
        $right_end = $right_start + $len * $flank_rate if ($len > 10000);
        my $start = $pos + $bin_size - ($pos % $bin_size);
        my @left_cov;
        my @right_cov;
        my @sv_cov;
        my @sv_cov2;
        my $ref_cov = \%cov;
        $ref_cov = \%cov2 if ($type eq 'DEL');
        while (1){
            if (exists ${$ref_cov}{$start}){
                push @sv_cov, ${$ref_cov}{$start};
                if (($type eq 'DEL') and ($len >= 1000)){
                    push @sv_cov2, $cov{$start};
                }
            }
            else{
                push @sv_cov, 0;
            }
            $start += $bin_size;
            last if ($start > $end);
        }
        my $sum = 0;
        map{$sum += $_} @sv_cov;
        if ($sum == 0){
            if (@line > 8){
                $line[9] .= ":0:0:0";
            }
            else{
                $line[7] .= ";DR=0;DS=0;SR=0";
            }
            print OUT join ("\t", @line), "\n";
            next;
        }
        my $ave = int ($sum / @sv_cov * 100) / 100;
        my $ave2 = 0;
        if (($type eq 'DEL') and ($len >= 1000)){
            my $sum2 = 0;
            map{$sum2 += $_} @sv_cov2;
            $ave2 = int ($sum2 / @sv_cov2 * 100) / 100 if (@sv_cov2 > 0);
        }
        
        while (1){
            if ((exists ${$ref_cov}{$left_end}) and (!exists ${$gap{$chr}}{$left_end})){
                if (($CNVnator_flag == 0) and ($type eq 'DEL') and (!exists ${$del{$chr}}{$left_end}) ){
                    push @left_cov, ${$ref_cov}{$left_end};
                }
                elsif (($CNVnator_flag == 0) and ($type eq 'DUP') and (!exists ${$dup{$chr}}{$left_end}) ){
                    push @left_cov, ${$ref_cov}{$left_end};
                }
                elsif (($CNVnator_flag == 1) and (!exists ${$del{$chr}}{$left_end}) and (!exists ${$dup{$chr}}{$left_end})){
                    push @left_cov, ${$ref_cov}{$left_end};
                }
                else{
                    $left_start -= $bin_size if ($left_start > $pos - 10000);
                }
            }
            else{
                $left_start -= $bin_size if ($left_start > $pos - 10000);
            }
            $left_end -= $bin_size;
            last if ($left_end < $left_start);
            last if ($left_end < $bin_size);
        }
        while (1){
            if ((exists ${$ref_cov}{$right_start}) and (!exists ${$gap{$chr}}{$right_start})){
                if (($CNVnator_flag == 0) and ($type eq 'DEL') and (!exists ${$del{$chr}}{$right_start})){
                    push @right_cov, ${$ref_cov}{$right_start};
                }
                elsif (($CNVnator_flag == 0) and ($type eq 'DUP') and (!exists ${$dup{$chr}}{$right_start})){
                    push @right_cov, ${$ref_cov}{$right_start};
                }
                elsif (($CNVnator_flag == 1) and (!exists ${$del{$chr}}{$right_start}) and (!exists ${$dup{$chr}}{$right_start})){
                    push @right_cov, ${$ref_cov}{$right_start};
                }
                else{
                    $right_end += $bin_size if ($right_end < $end + 10000);
                }
            }
            else{
                $right_end += $bin_size if ($right_end < $end + 10000);
            }
            $right_start += $bin_size;
            last if ($right_start > $right_end);
            last if ($right_start > $chrlen{$chr});
        }
        
        $cnv_num{$type} ++;
        my $sum_flank = 0;
        my $sum_flank5 = 0;
        my $sum_flank3 = 0;
        my $ave_flank = 0;
        my $ave_flank5 = 0;
        my $ave_flank3 = 0;
        map{$sum_flank += $_} @left_cov;
        map{$sum_flank5 += $_} @left_cov;
        map{$sum_flank += $_} @right_cov;
        map{$sum_flank3 += $_} @right_cov;
        $ave_flank5 = int ($sum_flank5 / @left_cov * 100) / 100 if (@left_cov > 0);
        $ave_flank3 = int ($sum_flank3 / @right_cov * 100) / 100 if (@right_cov > 0);
        my $flanknum = scalar (@left_cov) + scalar (@right_cov);
        if ((@left_cov > 0) and (@right_cov > 0) and ($ave_flank3 > 0) and (($ave_flank5 / $ave_flank3 > 1.5) or ($ave_flank5 / $ave_flank3 < 0.67))){
            my $ave_flank2 = $cov_ave;
            $ave_flank2 = $cov_ave2 if ($type eq 'DEL');
            if (abs ($ave_flank2 - $ave_flank5) <= abs ($ave_flank2 - $ave_flank3)){
                $ave_flank = $ave_flank5;
            }
            else{
                $ave_flank = $ave_flank3;
            }
        }
        elsif ($flanknum >= 10){
            $ave_flank = int ($sum_flank / $flanknum * 100) / 100 if ($flanknum > 0);
        }
        else{
            $ave_flank = $cov_ave;
            $ave_flank = $cov_ave2 if ($type eq 'DEL');
        }
        $cov_rate = int ($ave / $ave_flank * 100 + 0.5) / 100 if ($ave_flank > 0);
        if ($ave2 > 0){
            my $cov_rate2 = 0;
            $cov_rate2 = int ($ave2 / $ave_flank * 100 + 0.5) / 100 if ($ave_flank > 0);
            if ((defined $cov_rate2) and ($cov_rate2 > 1)){
                $cov_rate = $cov_rate2;
            }
        }
        my $incons_cov = 0;
        # inconsistent coverage corresponds to >= 0.91 DPR for DEL and <= 1.1 DPR for DUP in a 50-bp window size
        if ($ave_flank > 0){
            foreach (@sv_cov){
                if ($type eq 'DEL'){
                    $incons_cov ++ if ($_ / $ave_flank >= 0.91);
                }
                elsif ($type eq 'DUP'){
                    $incons_cov ++ if ($_ / $ave_flank <= 1.1);
                }
            }
        }
        # inconsistent coverage rate (DS) is the rate of the window with inconsistent coverage in all the number of window within DEL or DUP
        $incons_cov_rate = int ($incons_cov / @sv_cov * 100 + 0.5) / 100;
    }
    if (@line > 8){
        $line[9] .= ":$cov_rate:$incons_cov_rate:$split_rate";
    }
    else{
        $line[7] .= ";DR=$cov_rate;DS=$incons_cov_rate;SR=$split_rate";
    }
    print OUT join ("\t", @line), "\n";
}
close (FILE);
close (OUT);
