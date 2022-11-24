#!/usr/bin/perl -w
use strict;
use File::Basename;

# script for add_dp_rate_samples_3.pl in add_dp_rate subroutine of genotype_SV_SMC_7.4.pl

my $id = shift @ARGV;
my $pos_file = shift @ARGV;
my $cov_dir = shift @ARGV;
my $min_split_diff = shift @ARGV;

my $flank_rate = 0.1;
my $bin_size = 50;

my %cov;
my %cov2;
my %split;
my %chrlen;
my %pos_info;
my $cov2_flag = 0;
my $split_flag = 0;

# collect the information for SVs that DPR and SR are added to
open (FILE, $pos_file) or die "$pos_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my ($chr, $pos, $len, $type) = split (/\t/, $line);
    my $end = $pos + $len - 1;
    $end = $pos if ($type eq 'INS');
    my $pos_res = $pos % $bin_size;
    my $pos2 = $pos - $pos_res;
    ${$pos_info{$chr}}{$pos} = "$len=$type";
}

foreach my $chr (keys %pos_info){
    # collect coverage and split read data for a chromosome from a cov file
    my $cov_file = "$cov_dir/$id.chr$chr.cov.gz";
    $cov_file = "$cov_dir/$id.$chr.cov.gz" if (!-f $cov_file);
    my %cov_info;
    open (FILE, "gzip -dc $cov_file |") or die "$cov_file is not found:$!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr2, $cpos, $cov1, $cov2, $split1, $split2) = split (/\t/, $line);
        $cov_info{$cpos} = "$cov1=$cov2=$split1=$split2";
    }
    close (FILE);
    # calculate SR for any type of SV at an SV site using the split read data in 150 bp (50 bp window * 3) regions around the breakpoints
    # calculate DPR for DEL or DUP at an SV site using the covergae between the breakpoints and the coverage outsides the breakpoints
    foreach my $pos (sort {$a <=> $b} keys %{$pos_info{$chr}}){
        my ($len, $type) = split (/=/, ${$pos_info{$chr}}{$pos});
        my $end = $pos + $len - 1;
        $end = $pos if ($type eq 'INS');
        my $cov_rate = 1;
        my $incons_cov_rate = 0;
        my $split_rate = 0;
        my $CNVnator_flag = 0;
        my $flank_len = 1000;
        $flank_len = 100 if ($len < 150) and ($type eq 'DUP');
        my $pos2 = $pos - 1;
        $CNVnator_flag = 1 if ($pos2 % 1000 <= 1) and ($len >= 2000);
        my $left_end = $pos - 100;
        $left_end -= 400 if ($CNVnator_flag == 1);
        $left_end = $pos - 10 if ($len < 150);
        my $left_res = $left_end % $bin_size;
        $left_end -= $left_res;
        my $left_start = $left_end - $flank_len;
        $left_start = $left_end - int ($len * $flank_rate / $bin_size) * $bin_size if ($len > 10000);
        my $right_start = $end + 100;
        $right_start += 400 if ($CNVnator_flag == 1);
        $right_start = $end + 10 if ($len < 150);
        my $right_res = $right_start % $bin_size;
        $right_start += $bin_size - $right_res;
        my $right_end = $right_start + $flank_len;
        $right_end = $right_start + int ($len * $flank_rate / $bin_size) * $bin_size if ($len > 10000);
        my $start = $pos + $bin_size - ($pos % $bin_size);
        my $start2 = $pos - ($pos % $bin_size);
        my $start1 = $start2 - $bin_size;
        my $start3 = $start2 + $bin_size;
        my $start4 = $start3 + $bin_size;
        my $end2 = $end - ($end % $bin_size);
        my $end1 = $end2 - $bin_size;
        my $end3 = $end2 + $bin_size;
        my @left_cov;
        my @right_cov;
        my @sv_cov;
        my @sv_cov2;
        my @flankcov;
        my %split_num;
        my %split_num2;
        my %split2_num;
        my $split_diff_fold = 0;
        my %cov;
        
        if ($type =~ /DEL|DUP/){
            for (my $i = $left_start; $i <= $left_end; $i += $bin_size){
                next if (!exists $cov_info{$i});
                my ($cov1, $cov2, $split1, $split2) = split (/=/, $cov_info{$i});
                push @left_cov, $cov1 if ($type eq 'DUP') and ($cov1 > 0);
                push @left_cov, $cov2 if ($type eq 'DEL') and ($cov2 > 0);
            }
            for (my $i = $right_start; $i <= $right_end; $i += $bin_size){
                next if (!exists $cov_info{$i});
                my ($cov1, $cov2, $split1, $split2) = split (/=/, $cov_info{$i});
                push @right_cov, $cov1 if ($type eq 'DUP') and ($cov1 > 0);
                push @right_cov, $cov2 if ($type eq 'DEL') and ($cov2 > 0);
            }
            for (my $i = $start; $i <= $end2; $i += $bin_size){
                next if (!exists $cov_info{$i});
                my ($cov1, $cov2, $split1, $split2) = split (/=/, $cov_info{$i});
                push @sv_cov, $cov1 if ($type eq 'DUP');
                push @sv_cov, $cov2 if ($type eq 'DEL');
                push @sv_cov2, $cov1 if ($type eq 'DEL') and ($len >= 1000);
            }
        }
        for (my $i = $start1; $i <= $start4; $i += $bin_size){
            next if (!exists $cov_info{$i});
            my ($cov1, $cov2, $split1, $split2) = split (/=/, $cov_info{$i});
            $cov{$i} = $cov1;
            next if ($i == $start4) and ($type ne 'INS');
            if ($type eq 'DEL'){
                $split_num{$i} = $split2;
            }
            elsif ($type eq 'DUP'){
                $split_num{$i} = $split1;
            }
            elsif ($type eq 'INV'){
                $split_num{$i} = $split1;
                $split_num{$i} += $split2;
            }
            elsif ($type eq 'INS'){
                $split_num{$i} = $split1;
            }
        }
        for (my $i = $end1; $i <= $end3; $i += $bin_size){
            next if (!exists $cov_info{$i});
            my ($cov1, $cov2, $split1, $split2) = split (/=/, $cov_info{$i});
            if ($type eq 'DEL'){
                $split_num2{$i} = $split1;
            }
            elsif ($type eq 'DUP'){
                $split_num2{$i} = $split2;
            }
            elsif ($type eq 'INV'){
                $split_num2{$i} = $split1;
                $split_num2{$i} = $split2;
            }
            elsif ($type eq 'INS'){
                $split2_num{$i} = $split2;
            }
            $cov{$i} = $cov1;
        }
        if (scalar keys %cov == 0){
            print "$id\t$chr\t$pos\t$cov_rate\t$incons_cov_rate\t$split_rate\t$type\t0\n";
            next;
        }
        
        my $sum = 0;
        my $sum_flank = 0;
        my $ave = 0;
        my $ave2 = 0;
        my $ave_flank = 0;
        my $ave_left_flank = 0;
        my $ave_right_flank = 0;
        map{$sum += $_} @sv_cov;
        $ave = int ($sum / @sv_cov * 100) / 100 if (@sv_cov > 0);
        if (($type eq 'DEL') and ($len >= 1000)){
            my $sum2 = 0;
            map{$sum2 += $_} @sv_cov2;
            $ave2 = int ($sum2 / @sv_cov2 * 100) / 100 if (@sv_cov2 > 0);
        }
        if (@left_cov > 0){
            map{$sum_flank += $_} @left_cov;
            $ave_left_flank = int ($sum_flank / @left_cov * 100) / 100;
        }
        $sum_flank = 0;
        if (@right_cov > 0){
            map{$sum_flank += $_} @right_cov;
            $ave_right_flank = int ($sum_flank / @right_cov * 100) / 100;
        }
        if (($ave_left_flank == 0) and ($ave_right_flank == 0)){
            $ave_flank = 30;
        }
        elsif ($ave_left_flank == 0){
            $ave_flank = $ave_right_flank;
        }
        elsif ($ave_right_flank == 0){
            $ave_flank = $ave_left_flank;
        }
        else{   # for DEL/DUP with imbalanced coverage of upstream and downstream flanking regions, use the lower and higher coverage as flanking coverage for DEL and DUP, respectively
            if (($ave_left_flank / $ave_right_flank >= 1.2) or ($ave_left_flank / $ave_right_flank <= 0.83)){
                if ($type eq 'DEL'){
                    if ($ave_left_flank < $ave_right_flank){
                        $ave_flank = $ave_left_flank;
                    }
                    else{
                        $ave_flank = $ave_right_flank;
                    }
                }
                elsif ($type eq 'DUP'){
                    if ($ave_left_flank < $ave_right_flank){
                        $ave_flank = $ave_right_flank;
                    }
                    else{
                        $ave_flank = $ave_left_flank;
                    }
                }
            }
            else{
                $ave_flank = int (($ave_left_flank + $ave_right_flank) * 0.5);
            }
        }
        $cov_rate = int ($ave / $ave_flank * 100 + 0.5) / 100 if ($ave_flank > 0);
        if ($ave2 > 0){
            my $cov_rate2 = int ($ave2 / $ave_flank * 100 + 0.5) / 100 if ($ave_flank > 0);
            if ($cov_rate2 > 1){
                $cov_rate = $cov_rate2;
            }
        }
        my $incons_cov = 0;
        if (($ave_flank > 0) and (@sv_cov > 0)){
            foreach (@sv_cov){
                if (($type eq 'DEL') and ($_ / $ave_flank >= 0.91)){
                    $incons_cov ++;
                }
                elsif (($type eq 'DUP') and ($_ / $ave_flank <= 1.1)){
                    $incons_cov ++;
                }
            }
            $incons_cov_rate = int ($incons_cov / @sv_cov * 100 + 0.5) / 100 if (@sv_cov > 0);
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
        if (($top_num > 0) and ($sec_num > 0)){
            my $num_diff = abs ($top_num - $sec_num);
            if (($num_diff / $top_num <= 1.2) and ($num_diff / $top_num >= 0.8)){
                if ($top_pos > $sec_pos){
                    my $A = $top_pos;
                    $top_pos = $sec_pos;
                    $sec_pos = $A;
                }
            }
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
                my $split1_num = $top_num + $sec_num;
                my $split2_num = $top_num_2 + $sec_num_2;
                if (($split1_num >= $min_split_diff) or ($split2_num >= $min_split_diff)){
                    if (($split1_num > 0) and ($split2_num > 0)){
                        if ($split1_num >= $split2_num){
                            $split_diff_fold = int ($split1_num / $split2_num * 10 + 0.5) / 10;
                        }
                        else{
                            $split_diff_fold = int ($split2_num / $split1_num * 10 + 0.5) / 10;
                        }
                    }
                    else{
                        $split_diff_fold = 99;
                    }
                }
                else{
                    $split_diff_fold = 1;
                }
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
            my $top_num2 = 0;
            my $top_pos2 = 0;
            my $sec_num2 = 0;
            my $sec_pos2 = 0;
            my $split_rate2 = 0;
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
            if (($top_num2 > 0) and ($sec_num2 > 0)){
                my $num2_diff = abs ($top_num2 - $sec_num2);
                if (($num2_diff / $top_num2 <= 1.2) and ($num2_diff / $top_num2 >= 0.8)){
                    if ($top_pos2 < $sec_pos2){
                        my $A = $top_pos2;
                        $top_pos2 = $sec_pos2;
                        $sec_pos2 = $A;
                    }
                }
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
            # calculate the rate of split reass with 5'-split end and 3'-split end
            my $split1_num = $top_num + $sec_num;
            my $split2_num = $top_num2 + $sec_num2;
            if (($split1_num >= 4) or ($split2_num >= 4)){
                if (($split1_num > 0) and ($split2_num > 0)){
                    if ($split1_num >= $split2_num){
                        $split_diff_fold = int ($split1_num / $split2_num * 10 + 0.5) / 10;
                    }
                    else{
                        $split_diff_fold = int ($split2_num / $split1_num * 10 + 0.5) / 10;
                    }
                }
                else{
                    $split_diff_fold = 99;
                }
            }
            else{
                $split_diff_fold = 1;
            }
            if ($len <= 100){
                if (($top_pos > 0) and ($top_pos2 > 0) and ($top_pos > $top_pos2)){
                    $split_diff_fold = 99;
                }
            }
            else{
                if (($top_pos > 0) and ($top_pos2 > 0) and ($top_pos >= $top_pos2)){
                    $split_diff_fold = 99;
                }
            }
        }
        print "$id\t$chr\t$pos\t$cov_rate\t$incons_cov_rate\t$split_rate\t$type\t$split_diff_fold\n";
    }
    undef %cov_info;
}
