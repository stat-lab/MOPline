#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;

my $vcf = '';

my $non_human = 0;

my $build = 37;

my $segDup_filt = 1;

my $gap_bed = '';

my $segdup_file = '';

my $exclude_bed = '';

my $exclude_cen = '';

my $min_del_len = 10000;
my $min_dup_len1 = 1000;
my $min_dup_len2 = 200;
my $min_dup_len3 = 5000;

my $max_del_dprate1 = 0.85;
my $max_del_dprate2 = 0.75;

my $min_dup_dprate1 = 1.25;
my $min_dup_dprate2= 1.35;


my $max_del_dsr = 0.35;
my $max_del_dsr2 = 0.2;
my $max_dup_dsr = 0.1;

my $min_repeat_overlap = 1;
my $min_repeat_overlap2 = 0.7;

my $total_sample = 0;

my $min_sample_overlap_rate = 0.8;

my $min_overlap_rate = 0.5;
my $ins_sd = 20;

my $min_AF_lowconf = 0.01;

my $help;
 
GetOptions(
    'vcf|v=s' => \$vcf,
    'non_human|nh=i' => \$non_human,
    'build|b=s' => \$build,
    'gap|g=s' => \$gap_bed,
    'segdup|sg=s' => \$segdup_file,
    'exclude|ex=s' => \$exclude_bed,
    'exclude_cen|ec' => \$exclude_cen,
    'overlap_rate|or=f' => \$min_overlap_rate,
    'insdup_sd|ids=i' => \$ins_sd,
    'min_sd_duplen|msd=i' => \$min_dup_len3,
    'min_sd_overlap|mso=f' => \$min_repeat_overlap,
    'dis_sdf|ds' => \$segDup_filt,
    'del_dpr_len|dpl=i' => \$min_del_len,
    'del_dpr1|dp1=f' => \$max_del_dprate1,
    'del_dpr2|dp2=f' => \$max_del_dprate2,
    'dup_dpr_len1|upl1=i' => \$min_dup_len1,
    'dup_dpr_len2|upl2=i' => \$min_dup_len2,
    'dup_dpr1|up1=f' => \$min_dup_dprate1,
    'dup_dpr2|up2=f' => \$min_dup_dprate2,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  filter_MOPline.pl -v <SV vcf file> (-nh 1 -sg <segdup file> -g <gap bed> if sample is a non-human species) > out.filt.vcf
  output is STDOUT

  Options:
   --vcf or -v <STR>        SV vcf file [mandatory]
   --non_human or -nh <INT> sample is non-human species (0: humna, 1: non-human) [default: 0]
   --build or -b <STR>      human reference build (37, 38, or T2T; effective for only human) [default: 37]
   --gap or -g <STR>        gap bed file, indicating gap regions in reference genome [default for human: Data/gap.bed or gap.b38,bed]
   --segdup or -sg <STR>    sgmental duplication file from UCSC [default for human: Data/genomicSuperDups.txt.gz]
   --exclude or -ex <STR>   bed file indicating ambiguous/excludable regions in reference genome. These regions could often generate false positive SV calls
                            If it is prefferable for human data to disable this option, specify any character that is not a file.
                            [default for human: Data/Low_confidence.bed]
   --exclude_cen or -ec     bed file indicating centromere regions to exclude SVs within centromeres
                            If it is prefferable for human data to disable this option, specify any character that is not a file.
                            [default for human build 38 and T2T: Data/hg38.centromere.bed or Data/chm13.v2.0.centromere.bed]
   --overlap_rate or -or <FLOAT> minimum reciprocal overlap rate between SVs and ambiguous/excludable regions [default: 0.5]
   --insdup_sd or -isd <INT>  minimum distance (bp) between INS and DUP to be merged for a single sample vcf [default: 20]

   --dis_sdf or -ds         disable filtering of DUPs that overlap segmental duplicated regions with > 2 copies [default: false]
   --min_sd_duplen or -msd <INT>    minimum length (bp) of DUPs overlapping segmental duplicated regions to be filtered [default: 5000]
   --min_sd_overlap or -mso <FLOAT> minimum overlap rate of an overlap between DUP and segmental duplicated region to DUP length [default: 1.0]
   
   --del_dpr_len or -dpl <INT>   minimum length (bp) of DELs to be filtered based on DPR (combination with -dp1 or -dp2) [default: 10000]
   --del_dpr1 or -dp1 <FLOAT>    maximum DPR of DELs to be filtered (combination with -dpl) [default: 0.85]
   --del_dpr2 or -dp2 <FLOAT>    maximum DPR of DELs to be filtered (combination with -ds and -dpl) [default: 0.75]
   --del_dsr or -ds <FLOAT>      minimum DPS of DELs to be filtered (combination with -dp2 and -dpl) [default: 0.2]
   --dup_dpr_len1 or -upl1 <INT> minimum length (bp) of DUPs to be filtered based on DPR (combination with -up1) [default: 1000]
   --dup_dpr_len2 or -upl2 <INT> minimum length (bp) of DUPs to be filtered based on DPR (combination with -up2 and -us) [default: 200]
   --dup_dpr1 or -up1 <FLOAT>    minimum DPR of DUPs to be filtered (combination with -upl1) [default: 1.25]
   --dup_dpr2 or -up2 <FLOAT>    minimum DPR of DUPs to be filtered (combination with -us and -upl2) [default: 1.35]
   --dup_dsr or -us <FLOAT>      minimum DPR of DUPs to be filtered (combination with -up2 and -upl2) [default: 0.1]

   --help or -h             output help message
   
=cut

my $data_dir = "$Bin/../Data";

if ($non_human == 0){
    if ($gap_bed eq ''){
        $gap_bed = "$data_dir/gap.bed" if ($build eq '37');
        $gap_bed = "$data_dir/gap.b38.bed" if ($build eq '38');
    }
    if ($segdup_file eq ''){
        $segdup_file = "$data_dir/genomicSuperDups.txt.gz";
        $segdup_file = "$data_dir/genomicSuperDups.b38.txt.gz" if ($build eq '38');
        $segdup_file = "$data_dir/genomicSuperDups.T2T-chm13.v2.0.txt.gz" if ($build eq 'T2T');
    }
    if ($exclude_bed eq ''){
        $exclude_bed = "$data_dir/Low_confidence.bed";
        $exclude_bed = "$data_dir/Low_confidence.b38.bed" if ($build eq '38');
        $exclude_bed = "$data_dir/Low_confidence.T2T-chm13.v2.0.bed" if ($build eq 'T2T');
    }
    if ($exclude_cen eq ''){
        $exclude_cen = "$data_dir/hg38.centromere.bed" if ($build eq '38');
        $exclude_cen = "$data_dir/chm13.v2.0.centromere.bed" if ($build eq 'T2T');
    }
}
print STDERR "gap bed: $gap_bed\n" if ($gap_bed ne '');
print STDERR "segDup bed: $segdup_file\n" if ($segdup_file ne '');
print STDERR "centromere bed: $exclude_cen" if ($exclude_cen ne '');
print STDERR "exlude SV bed: $exclude_bed\n" if ($exclude_bed ne '');

my $del_num = 0;
my $dup_num = 0;
my $ins_num = 0;
my $inv_num = 0;
my $del_filt = 0;
my $dup_filt = 0;
my $ins_filt = 0;
my $inv_filt = 0;
my @del_dprate;
my @dup_dprate;


my %gap;
my %repeat;
my %repeat2;
my %segdup;
my %ambiguous;
my %ambiguous_flag;
my %exclude;
my $total_segdup = 0;
my $Mbin_size = 1000000;


my $pre_gchr = '';
my $pre_gend = 0;
my $pre_gpos = 0;

if ($gap_bed ne ''){
    open (FILE, $gap_bed) or die "$gap_bed is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\s+/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $end = $line[2];
        if (($chr eq $pre_gchr) and ($pos - $pre_gend <= 1000)){
            ${$gap{$chr}}{$pre_gpos} = $end;
            $pre_gend = $end;
            next;
        }
        ${$gap{$chr}}{$pos} = $end;
        $pre_gchr = $chr;
        $pre_gend = $end;
        $pre_gpos = $pos;
    }
    close (FILE);
}

if ($segdup_file ne ''){
    open (FILE, "gzip -dc $segdup_file |") or die "$segdup_file is not found: $!\n" if ($segdup_file =~ /\.gz$/);
    open (FILE, $segdup_file) or die "$segdup_file is not found: $!\n" if ($segdup_file !~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my @line = split (/\t/, $line);
        my $chr = $line[1];
        $chr =~ s/^chr// if ($chr =~ /^chr/) and ($non_human == 0) and ($build eq '37');
        next if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
        my $pos = $line[2];
        my $end = $line[3];
        if (!exists ${$repeat{$chr}}{$pos}){
            ${$repeat{$chr}}{$pos} = $end;
        }
        else{
            if ($end > ${$repeat{$chr}}{$pos}){
                ${$repeat{$chr}}{$pos} = $end;
            }
        }
        ${${$segdup{$chr}}{$pos}}{$end} = 1;
    }
    close (FILE);
}

foreach my $chr (keys %repeat){
    my $pre_pos = 0;
    my $pre_end = 0;
    foreach my $pos (sort {$a <=> $b} keys %{$repeat{$chr}}){
        my $end = ${$repeat{$chr}}{$pos};
        if (($pre_pos > 0) and ($pos <= $pre_end + 10)){
            if ($end <= $pre_end){
                delete ${$repeat{$chr}}{$pos};
                next;
            }
            else{
                delete ${$repeat{$chr}}{$pos};
                ${$repeat{$chr}}{$pre_pos} = $end;
                $pre_end = $end;
                next;
            }
        }
        $pre_pos = $pos;
        $pre_end = $end;
    }
}

foreach my $chr (keys %repeat){
    foreach my $pos (sort {$a <=> $b} keys %{$repeat{$chr}}){
        my $end = ${$repeat{$chr}}{$pos};
        my $Mbin_start = int ($pos / $Mbin_size);
        my $Mbin_end = int ($end / $Mbin_size);
        if ($Mbin_start == $Mbin_end){
            ${${$repeat2{$chr}}{$Mbin_start}}{$pos} = $end;
        }
        else{
            my $Mbin = $Mbin_start;
            while ($Mbin <= $Mbin_end){
                ${${$repeat2{$chr}}{$Mbin}}{$pos} = $end;
                $Mbin ++;
            }
        }
        $total_segdup += $end - $pos + 1;
    }
}
%repeat = ();

print STDERR "#Total length of segmental duplication regions:\t$total_segdup\n";

if (($exclude_bed ne '') and (-f $exclude_bed)){
    open (FILE, $exclude_bed) or die "$exclude_bed is not found:$!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $pos, $end, $type) = split (/\t/, $line);
        my $len = $end - $pos + 1;
        if ((defined $type) and ($type =~ /DEL|DUP|INV|INS|ALL/)){
            if ($type =~ /\//){
                my @type = split (/\//, $type);
                foreach (@type){
                    ${${$ambiguous{$_}}{$chr}}{$pos} = $len;
                    ${${$ambiguous{$_}}{$chr}}{$pos} = $len;
                }
            }
            elsif ($type eq 'ALL'){
                ${${$ambiguous{'INS'}}{$chr}}{$pos} = $len;
                ${${$ambiguous{'DUP'}}{$chr}}{$pos} = $len;
                ${${$ambiguous{'DEL'}}{$chr}}{$pos} = $len;
                ${${$ambiguous{'INV'}}{$chr}}{$pos} = $len;
            }
            else{
                ${${$ambiguous{$type}}{$chr}}{$pos} = $len;
                ${${$ambiguous_flag{$type}}{$chr}}{$pos} = 1;
            }
        }
        else{
            ${$exclude{$chr}}{$pos} = $end;
        }
    }
    close (FILE);
}
if (($exclude_cen ne '') and (-f $exclude_cen)){
    open (FILE, $exclude_cen) or die "$exclude_cen is not found:$!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $pos, $end) = split (/\t/, $line);
        ${$exclude{$chr}}{$pos} = $end;
    }
    close (FILE);
}

my %sv;
my %delete;
my %merged_sv;
my %gap_dup;
my $rep_dup = 0;

my @header;

open (FILE, $vcf) or die "$vcf is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my @line = split (/\t/, $line);
    if ($line =~ /^#/){
        push @header, $line;
        if ($line =~ /^#CHROM/){
            my $count = 0;
            $total_sample = 0;
            foreach (@line){
                $count ++;
                next if ($count <= 9);
                $total_sample ++;
            }
        }
        next;
    }
    my $chr = $line[0];
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    ${${$sv{$type}}{$chr}}{$pos} = $line;
    $del_num ++ if ($type eq 'DEL');
    $dup_num ++ if ($type eq 'DUP');
    $ins_num ++ if ($type eq 'INS');
    $inv_num ++ if ($type eq 'INV')   
}
close (FILE);

print STDERR "Total samples: $total_sample\n";

foreach my $type (keys %sv){
    next if ($type eq 'INS');
    foreach my $chr (keys %{$sv{$type}}){
        foreach my $pos (sort {$a <=> $b} keys %{${$sv{$type}}{$chr}}){
            next if (exists ${${$delete{$type}}{$chr}}{$pos});
            my $line = ${${$sv{$type}}{$chr}}{$pos};
            my @line = split (/\t/, $line);
            my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;

            if (($type eq 'DUP') and ($len >= 150)){        # filter DUPs that are located in the flanking regions of gaps
                my $pos_f1 = $pos - $len * 1.2;
                my $pos_f2 = $pos - $len * 0.2;
                my $end_f1 = $end + $len * 1.2;
                my $end_f2 = $end + $len * 0.2;
                my $flank_gap = 0;
                foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
                    my $gend = ${$gap{$chr}}{$gstart};
                    last if ($gstart > $end_f1);
                    if (($gend >= $pos_f2) and ($gend < $pos + $len * 0.1)){
                        $flank_gap += $gend - $gstart + 1 if ($gstart > $pos_f1);
                        $flank_gap += $gend - $pos_f1 + 1 if ($gstart <= $pos_f1);
                    }
                    if (($gstart <= $end_f2) and ($gstart > $end - $len * 0.1)){
                        $flank_gap += $gend - $gstart + 1 if ($gend < $end_f1);
                        $flank_gap += $end_f1 - $gstart + 1 if ($gend >= $end_f1);
                    }
                }
                if (($flank_gap / $len >= 0.5) or ($flank_gap >= 20000)){
                    my $chr_pos = "$chr:$pos";
                    $gap_dup{$chr_pos} = $len;
                    ${${$delete{$type}}{$chr}}{$pos} = 1;
                    next;
                }
            }
            if (($type eq 'DEL') and ($len >= 10000)){      # filter DELs whose breakpoints oveerlap a gap
                my $gap_flag = 0;
                foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
                    my $gend = ${$gap{$chr}}{$gstart};
                    last if ($gstart > $pos + 1000);
                    if (($pos >= $gstart - 200) and ($pos <= $gstart + 1000)){
                        $gap_flag = 1;
                        last;
                    }
                    if (($end >= $gend - 1000) and ($end <= $gend + 200)){
                        $gap_flag = 1;
                        last;
                    }
                }
                if ($gap_flag == 1){
                    ${${$delete{$type}}{$chr}}{$pos} = 1;
                    $del_filt ++;
                }
            }
            if (($type eq 'DUP') and ($len >= $min_dup_len3) and ($segDup_filt == 1)){  # filter DUPs overlapping SegDups with > 2 copies in the genome
                my $Mbin_start = int ($pos / $Mbin_size);
                my $Mbin_end = int ($end / $Mbin_size);
                my $flag = 0;
                my $overlap = 0;
                if ($Mbin_start == $Mbin_end){
                    if (exists ${$repeat2{$chr}}{$Mbin_start}){
                        foreach my $rpos (sort {$a <=> $b} keys %{${$repeat2{$chr}}{$Mbin_start}}){
                            last if ($rpos > $end);
                            my $rend = ${${$repeat2{$chr}}{$Mbin_start}}{$rpos};
                            my $rlen = $rend - $rpos + 1;
                            next if ($rend < $pos);
                            if (($rpos <= $pos) and ($rend >= $end)){
                                $overlap += $len;
                            }
                            elsif (($rpos >= $pos) and ($rpos <= $end)){
                                if ($rend < $end){
                                    $overlap += $rlen ;
                                }
                                else{
                                    $overlap += $end - $rpos + 1;
                                }
                            }
                            elsif (($rend >= $pos) and ($rend <= $end)){
                                if ($rpos > $pos){
                                    $overlap += $rlen;
                                }
                                else{
                                    $overlap += $rend - $pos + 1;
                                }
                            }
                        }
                        if ($overlap >= $len * $min_repeat_overlap){
                            $flag = 1;
                        }
                    }
                }
                else{
                    if (exists ${$repeat2{$chr}}{$Mbin_start}){
                        my %rep = (%{${$repeat2{$chr}}{$Mbin_start}});
                        for (my $i = $Mbin_start + 1; $i <= $Mbin_end; $i++){
                            %rep = (%rep, %{${$repeat2{$chr}}{$i}}) if (exists ${$repeat2{$chr}}{$i});
                        }
                        foreach my $rpos (sort {$a <=> $b} keys %rep){
                            last if ($rpos > $end);
                            my $rend = $rep{$rpos};
                            my $rlen = $rend - $rpos + 1;
                            next if ($rend < $pos);
                            if (($rpos <= $pos) and ($rend >= $end)){
                                $overlap += $len;
                            }
                            elsif (($rpos >= $pos) and ($rpos <= $end)){
                                if ($rend < $end){
                                    $overlap += $rlen ;
                                }
                                else{
                                    $overlap += $end - $rpos + 1;
                                }
                            }
                            elsif (($rend >= $pos) and ($rend <= $end)){
                                if ($rpos > $pos){
                                    $overlap += $rlen;
                                }
                                else{
                                    $overlap += $rend - $pos + 1;
                                }
                            }
                        }
                        if ($overlap >= $len * $min_repeat_overlap){
                            $flag = 1;
                        }
                    }
                }
                if ($flag == 1){
                    my $rcount = 0;
                    foreach my $rpos (sort {$a <=> $b} keys %{$segdup{$chr}}){
                        last if ($rpos > $end);
                        foreach my $rend (sort {$a <=> $b} keys %{${$segdup{$chr}}{$rpos}}){
                            if (($rpos >= $pos) and ($rpos <= $end)){
                                $rcount ++;
                            }
                            elsif (($rend >= $pos) and ($rend <= $end)){
                                $rcount ++;
                            }
                            elsif (($rpos <= $pos) and ($rend >= $end)){
                                $rcount ++;
                            }
                        }
                    }
                    $line[2] = $rcount;
                    $flag = 0 if ($rcount <= 1);
                }
                if ($flag == 1){
                    ${${$delete{$type}}{$chr}}{$pos} = 1;
                    $rep_dup ++;
                }
            }
        }
    }
}

my %vcf;
my %INS;
my %DUP;

foreach my $type (keys %sv){
    foreach my $chr (keys %{$sv{$type}}){
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
        foreach my $pos (keys %{${$sv{$type}}{$chr}}){
            next if (exists ${${$delete{$type}}{$chr}}{$pos});
            ${${$vcf{$chr02d}}{$pos}}{$type} = ${${$sv{$type}}{$chr}}{$pos};
            if ($total_sample == 1){
                ${$INS{$chr02d}}{$pos} = ${${$sv{$type}}{$chr}}{$pos} if ($type eq 'INS');
                ${$DUP{$chr02d}}{$pos} = ${${$sv{$type}}{$chr}}{$pos} if ($type eq 'DUP');
            }
        }
    }
}

if ($total_sample <= 1){
    my $dup_conv = 0;
    foreach my $chr (keys %INS){
        next if (!exists $DUP{$chr});
        foreach my $ipos (sort {$a <=> $b} keys %{$INS{$chr}}){
            foreach my $dpos (sort {$a <=> $b} keys %{$DUP{$chr}}){
                last if ($dpos > $ipos + $ins_sd);
                my $dline = ${$DUP{$chr}}{$dpos};
                my $dlen = $1 if ($dline =~ /SVLEN=(\d+)/);
                my $dend = $dpos + $dlen - 1;
                next if ($dend < $ipos - $ins_sd);
                if ((abs ($ipos - $dpos) <= $ins_sd) or (abs ($ipos = $dend) <= $ins_sd)){  # remove either INS or DUP if they are located in close proximity
                    if ($dlen >= 100){
                        delete ${${$vcf{$chr}}{$ipos}}{'INS'};
                        $ins_filt ++;
                        $dup_conv ++;
                    }
                    else{
                        delete ${${$vcf{$chr}}{$dpos}}{'DUP'};
                        $dup_filt ++;
                        $dup_conv ++;
                    }
                }
            }
        }
    }
    print STDERR "DUP-INS overlaps: $dup_conv\n";
}

foreach my $type (sort keys %merged_sv){
    my $merged_num = scalar @{$merged_sv{$type}};
    print STDERR "Merged $type num: $merged_num\n";
}


my $ambiguous_del = 0;
my $ambiguous_dup = 0;
my $ambiguous_ins = 0;
my @filt_sv;

foreach (@header){
    print "$_\n";
}
foreach my $chr (sort keys %vcf){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
            my $line = ${${$vcf{$chr}}{$pos}}{$type};
            my @line = split(/\s+/, $line);
            my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
            my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            my $af = 1;
            $af = $1 if ($line[7] =~ /AF=([\d\.]+)/);
            my $read = 15;
            $read = $1 if ($line[7] =~ /READS=(\d+)/);
            my $DPR = 0;
            my $DSR = 0;
            if ($line[7] =~ /DPR=([\d\.]+)/){
                $DPR = $1;
            }
            elsif ($line[7] =~ /DR=([\d\.]+)/){
                $DPR = $1;
            }
            elsif (($line[8] =~ /DR/) and (@line > 9)){
                my @format = split (/:/, $line[8]);
                my $fcount = 0;
                foreach (@format){
                    last if ($_ eq 'DR');
                    $fcount ++;
                }
                my @info = split (/:/, $line[9]);
                $DPR = $info[$fcount];
            }
            if ($line[7] =~ /DSR=([\d\.]+)/){
                $DSR = $1;
            }
            elsif ($line[7] =~ /DPS=([\d\.]+)/){
                $DSR = $1;
            }
            elsif (($line[8] =~ /:DS/) and (@line > 9)){
                my @format = split (/:/, $line[8]);
                my $fcount = 0;
                foreach (@format){
                    last if ($_ eq 'DS');
                    $fcount ++;
                }
                my @info = split (/:/, $line[9]);
                $DSR = $info[$fcount];
            }
            if ($type eq 'DEL'){
                my $Mbin_start = int ($pos / $Mbin_size);
                my $Mbin_end = int ($end / $Mbin_size);
                my $flag = 0;
                my $overlap = 0;
                if ($Mbin_start == $Mbin_end){
                    if (exists ${$repeat2{$chr}}{$Mbin_start}){
                        foreach my $rpos (sort {$a <=> $b} keys %{${$repeat2{$chr}}{$Mbin_start}}){
                            last if ($rpos > $end);
                            my $rend = ${${$repeat2{$chr}}{$Mbin_start}}{$rpos};
                            my $rlen = $rend - $rpos + 1;
                            next if ($rend < $pos);
                            if (($rpos <= $pos) and ($rend >= $end)){
                                $overlap += $len;
                            }
                            elsif (($rpos >= $pos) and ($rpos <= $end)){
                                if ($rend < $end){
                                    $overlap += $rlen ;
                                }
                                else{
                                    $overlap += $end - $rpos + 1;
                                }
                            }
                            elsif (($rend >= $pos) and ($rend <= $end)){
                                if ($rpos > $pos){
                                    $overlap += $rlen;
                                }
                                else{
                                    $overlap += $rend - $pos + 1;
                                }
                            }
                        }
                        if ($overlap >= $len * $min_repeat_overlap2){
                            $flag = 1;
                        }
                    }
                }
                else{
                    if (exists ${$repeat2{$chr}}{$Mbin_start}){
                        my %rep = (%{${$repeat2{$chr}}{$Mbin_start}});
                        for (my $i = $Mbin_start + 1; $i <= $Mbin_end; $i++){
                            %rep = (%rep, %{${$repeat2{$chr}}{$i}}) if (exists ${$repeat2{$chr}}{$i});
                        }
                        foreach my $rpos (sort {$a <=> $b} keys %rep){
                            last if ($rpos > $end);
                            my $rend = $rep{$rpos};
                            my $rlen = $rend - $rpos + 1;
                            next if ($rend < $pos);
                            if (($rpos <= $pos) and ($rend >= $end)){
                                $overlap += $len;
                            }
                            elsif (($rpos >= $pos) and ($rpos <= $end)){
                                if ($rend < $end){
                                    $overlap += $rlen ;
                                }
                                else{
                                    $overlap += $end - $rpos + 1;
                                }
                            }
                            elsif (($rend >= $pos) and ($rend <= $end)){
                                if ($rpos > $pos){
                                    $overlap += $rlen;
                                }
                                else{
                                    $overlap += $rend - $pos + 1;
                                }
                            }
                        }
                        if ($overlap >= $len * $min_repeat_overlap2){
                            $flag = 1;
                        }
                    }
                }
                if ($flag == 0){
                    if ($len >= $min_del_len){
                        if (($DPR > $max_del_dprate2) and ($DSR > $max_del_dsr2)){      # filter DELs based on DPR, DPS, and size in non repeat regions
                            push @filt_sv, $line;
                            $del_filt ++;
                            next;
                        }
                    }
                }
                if ($len >= $min_del_len){
                    if (($DPR > $max_del_dprate1) or ($DSR > $max_del_dsr)){            # filter DELs based on DPR, DPS, and size
                        push @filt_sv, $line;
                        $del_filt ++;
                        next;
                    }
                }
            }
=pod
            if (($type eq 'DEL') and ($len >= $min_del_len) and ($len < $min_del_len2)){
                if (($DPR > $max_del_dprate) and ($DSR > $max_del_dsr) and ($read < 15)){
                    $del_filt ++;
                    push @filt_sv, $line;
                    next;
                }
            }
            if (($type eq 'DEL') and ($len >= $min_del_len2)){
                if (($DPR > $max_del_dprate2) and ($DSR > $max_del_dsr)){
                    push @filt_sv, $line;
                    $del_filt ++;
                    next;
                }
            }
            if (($type eq 'DEL') and ($len >= $min_len) and ($DPR > $max_del_dprate)){
                push @filt_sv, $line;
                $del_filt ++;
                next;
            }
=cut
            if (($type eq 'DUP') and ($len >= $min_dup_len2)){
                if (($DPR < $min_dup_dprate2) and ($DSR > $max_dup_dsr)){               # filter DUPs based on DPR, DPS, and size
                    push @filt_sv, $line;
                    $dup_filt ++;
                    next;
                }
            }
            if (($type eq 'DUP') and ($len >= $min_dup_len1) and ($DPR < $min_dup_dprate1)){    # filter DUPs based on DPR and size in non repeat regions
                push @filt_sv, $line;
                $dup_filt ++;
                next;
            }
            if ((exists ${$ambiguous{$type}}{$chr2}) and ($af >= $min_AF_lowconf)){     # filter DELs/DUPs/INSs in low confidence regions specific to SV type, based on AF
                my $match = 0;
                foreach my $apos (sort {$a <=> $b} keys %{${$ambiguous{$type}}{$chr2}}){
                    last if ($apos > $end);
                    my $alen = ${${$ambiguous{$type}}{$chr2}}{$apos};
                    my $aend = $apos + $alen - 1;
                    next if ($aend < $pos);
                    if ($type eq 'INS'){
                        $match = 1;
                        last;
                    }
                    else{
                        my $overlap = 0;
                        if (($apos <= $pos) and ($aend >= $end)){
                            $overlap = $len;
                        }
                        elsif (($apos >= $pos) and ($apos <= $end)){
                            $overlap = $end - $apos + 1;
                            $overlap = $alen if ($aend < $end);
                        }
                        elsif (($aend >= $pos) and ($aend <= $end)){
                            $overlap = $aend - $pos + 1;
                            $overlap = $alen if ($apos > $pos);
                        }
                        if ($overlap >= $len * $min_repeat_overlap2){
                            if (exists ${${$ambiguous_flag{$type}}{$chr2}}{$apos}){
                                if ($len >= 8000){
                                    $match = 1;
                                    last;
                                }
                            }
                            else{
                                $match = 1;
                                last;
                            }
                        }
                    }
                }
                if ($match == 1){
                    if ($type eq 'DEL'){
                        $del_filt ++;
                        $ambiguous_del ++;
                    }
                    elsif ($type eq 'DUP'){
                        $dup_filt ++;
                        $ambiguous_dup ++;
                    }
                    elsif ($type eq 'INS'){
                        $ins_filt ++;
                        $ambiguous_ins ++;
                    }
                    next;
                }
            }
            if (exists $exclude{$chr}){
                my $match = 0;
                foreach my $cpos (sort {$a <=> $b} keys %{$exclude{$chr}}){
                    last if ($cpos > $end);
                    my $cend = ${$exclude{$chr}}{$cpos};
                    next if ($cend < $pos);
                    if ($type eq 'INS'){
                        $match = 1;
                        last;
                    }
                    else{
                        my $clen = $cend - $cpos + 1;
                        my $overlap = 0;
                        if (($cpos <= $pos) and ($cend >= $end)){
                            $match = 1;
                            last;
                        }
                        elsif (($cpos >= $pos) and ($cpos <= $end)){
                            $overlap = $end - $cpos + 1;
                            $overlap = $clen if ($cend < $end);
                        }
                        elsif (($cend >= $pos) and ($cend <= $end)){
                            $overlap = $cend - $pos + 1;
                            $overlap = $clen if ($cpos > $pos);
                        }
                        if ($overlap >= $len *  $min_overlap_rate){
                            $match = 1;
                            last;
                        }
                    }
                }
                if ($match == 1){
                    if ($type eq 'DEL'){
                        $del_filt ++;
                    }
                    elsif ($type eq 'DUP'){
                        $dup_filt ++;
                    }
                    elsif ($type eq 'INS'){
                        $ins_filt ++;
                    }
                    elsif ($type eq 'INV'){
                        $inv_filt ++;
                    }
                    next;
                }
            }
            print "$line\n";
        }
    }
}

my $gap_dup_num = scalar keys %gap_dup;
print STDERR "Filterd DUPs with flanking Gaps: $gap_dup_num\n";
print STDERR "Filterd DUPs overlapping seg-dup: $rep_dup\n";
print STDERR "Filterd DUPs overlapping ambiguous DUPs: $ambiguous_dup\n";
print STDERR "Filterd DELs overlapping ambiguous DELs: $ambiguous_del\n";


$dup_filt += $gap_dup_num;
$dup_filt += $rep_dup;

my $del_filt_rate = 0;
my $dup_filt_rate = 0;
my $ins_filt_rate = 0;
my $inv_filt_rate = 0;

$del_filt_rate = int ($del_filt / $del_num * 1000) / 10 if ($del_num > 0);
$dup_filt_rate = int ($dup_filt / $dup_num * 1000) / 10 if ($dup_num > 0);
$ins_filt_rate = int ($ins_filt / $ins_num * 1000) / 10 if ($ins_num > 0);
$inv_filt_rate = int ($inv_filt / $inv_num * 1000) / 10 if ($inv_num > 0);

print STDERR "\nTotal DEL num: $del_num\tFiltered: $del_filt ($del_filt_rate%)\n";
print STDERR "Total DUP num: $dup_num\tFiltered: $dup_filt ($dup_filt_rate%)\n";
print STDERR "Total INS num: $ins_num\tFiltered: $ins_filt ($ins_filt_rate%)\n";
print STDERR "Total INV num: $inv_num\tFiltered: $inv_filt ($inv_filt_rate%)\n";

