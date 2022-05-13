#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);
use threads;

my $sv_vcf = '';

my $sample_dir = '.';

my $target_chr = 'ALL';

my $target_type = 'ALL';

my $tool_set = '7tools';

my $tool_set_all = 'TRUE';

my $non_human = 0;

my $build = 37;

my $ins_sd = 100;
my $BPsd = $ins_sd;

#my $min_overlap_rate = 0.7;
my $min_overlap_rate = 0.5;
my $min_overlap_rate2 = 0.5;

my $min_repeat_overlap = 0.3;

my $min_sv_len = 50;

my $thread = 1;

my $out_prefix = 'MOPline';

my $out_dir = 'MOPline_SMC';

my $disable_SMC_2nd = 0;

my $DP_filter = 0;

my $min_prob = 0.9;
my $min_prob3 = 0.9;
my $min_prob2 = 0.99;

my $min_SF = 0.01;
my $min_SF2 = 0.8;
my $min_SC = 5;
my $min_SC_DEL = 3;
my $min_SMC_DEL = 50;

my $min_split_diff = 8;
my $min_split_diff_dup = 4;
my $min_split_diff_del = 4;
my $min_ref_sample = 10;
my $min_ref_sf = 0.01;
my $min_split_filt_rate = 0.9; # minimum rate of samples with ref allele, that had been assigned due to biased direction (> $min_split_diff-fold difference) of split-alignment reads

my $data_dir = "$Bin/../Data";

my $gender_list = '';

my $simple_repeat = '';

my $segdup_file = '';

my $help;


GetOptions(
    'vcf|v=s' => \$sv_vcf,
    'target_chr|c=s' => \$target_chr,
    'target_type|tt=s' => \$target_type,
    'out_prefix|p=s' => \$out_prefix,
    'sample_dir|sd=s' => \$sample_dir,
    'out_dir|od=s' => \$out_dir,
    'toolset|ts=s' => \$tool_set,
    'toolset_all|ta=s' => \$tool_set_all,
    'nthread|n=i' => \$thread,
    'gender_list|g=s' => \$gender_list,
    'simple_repeat|sr=s' => \$simple_repeat,
    'segdup|sg=s' => \$segdup_file,
    'non_human|nh=i' => \$non_human,
    'build|b=i' => \$build,
    'ins_sd|is=i' => \$ins_sd,
    'overlap_rate|or=f' => \$min_overlap_rate,
    'repeat_overlap|ro=f' => \$min_repeat_overlap,
    'min_prob_mc|mpc=f' => \$min_prob,
    'min_prob|mp=f' => \$min_prob3,
    'min_prob_mb|mpb=f' => \$min_prob2,
    'min_sf|ms=f' => \$min_SF,
    'min_sf2|ms2=f' => \$min_SF2,
    'min_sc|msc=i' => \$min_SC,
    'min_sc2|msc2=i' => \$min_SC_DEL,
    'min_dl|mdl=i' => \$min_SMC_DEL,
    'min_split_diff_ins|msdi=f' => \$min_split_diff,
    'min_split_diff_del|msdd=f' => \$min_split_diff_del,
    'min_split_diff_dup|msdu=f' => \$min_split_diff_dup,
    'min_filt_rate_split|mfs=f' => \$min_split_filt_rate,
    'dp_filt|df' => \$DP_filter,
    'dis_smc2|ds' => \$disable_SMC_2nd,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;


=head1 SYNOPSIS

  genotype_SV_SMC_7.4.pl -v <vcf file> -ts <preset name of tool set> -od <output directory> -p <output prefix> -n <num of threads> (-nh 1 -sr <simple repeat file> -sd <segmental dup file> if sample is a non-human species)

  Options:
   --vcf or -v <STR>         input vcf file of SVs [mandatory]
   --sample_dir or -sd <STR> directory containing sample directories or group_directories/sample_directories [default: ./]
   --out_dir or -od <STR>    output directory [default: MOPline_SMC]
   --out_prefix pr -p <STR>  prefix name of an output vcf file [default: MOPline]
   --target_chr or -c <STR>  target chromosome [default: ALL]
   --target_type or -tt <STR> target SV type [default: ALL]
   --toolset or -ts <STR>    a preset toolset of SV detection tools used (6tools_1, 6tools_2, 7tools, 9tools, or 11tools) or a list file describing tool names in each line [default: 7tools]
   --toolset_all or -ta <STR> consider all the SVs calls (except for < 100 bp DEL/DUP) from the original single tools to execute SMC (if FLASE, consider minimum supporting reads indicated in a tool parameter file, specified with -ts or -ct option) [default: TRUE]
   --non_human or -nh <INT>  samples are non-human species (0: human, 1: non-human) [default: 0]
   --build or -b <INT>       reference build (GRCh37, GRCh38) number (37 or 38) [default: 37]
   --gender_list or -g <STR> sample-gender list file to adjust genotypes for sex chromosomes (only for X, Y, chrX, or chrY) (1st column: sample ID, 2nd column: M|F) [optional]
   --simple_repeat or -sr <STR> simple repeat file from UCSC (simpleRepeat.txt.gz)
   --segdup or -sg <STR> segmental duplication file from UCSC (genomicSuperDups.txt.gz)
   
   --ins_sd or -is <INT>     maximum distance (bp) between proximal INS breakpoints to consider identical INSs [default: 100]
   --overlap_rate or -or <FLOAT>   minimum reciprocal overlap rate to consider identical DELs, DUPs, or INVs [default: 0.5]
   --repeat_overlap or -ro <FLOAT> minimum overlap rate to cinsider SVs overlapping repeats [default: 0.3]

   --min_prob_mc or -mpc <FLOAT> minimum probability of converting undef genotype to non-ref genotype for SMC genotype prediction [default: 0.9]
   --min_prob_mb or -mpb <FLOAT> minimum probability of converting ref genotype to non-ref genotype for SMC genotype prediction [default: 0.99]
   --min_prob or -mp <FLOAT> minimum probability of converting non-ref genotype to ref genotype for SMC genotype prediction [default: 0.9]
   --min_sf or -ms <FLOAT>   minimum site frequency of INSs subject to the 2nd SMC [default: 0.01]
   --min_sf2 or -ms2 <FLOAT> minimum site frequency of INVs subject to the 2nd SMC [default: 0.8]
   --min_sc or -msc <INT>    minimum site count of INSs subject to the 2nd SMC [default: 5]
   --min_sc2 or -msc2 <INT>  minimum site count of DELs/DUPs subject to the 2nd SMC [default: 3]
   --min_dl or -mdl <INT>    minimum length of DEL/DUP subject to the 2nd SMC [default: 50]

   --min_split_diff_ins or -msdi <FLOAT>  minimum fold of the difference between the numbers of forward and reverse soft-clipped reads supporting INS (INS non-ref genotype that exceeds this value is converted to ref genotype) [default: 8.0]
   --min_split_diff_del or -msdd <FLOAT>  minimum fold of the difference between the numbers of forward and reverse soft-clipped reads supporting <= 500 bp DEL (DEL non-ref genotype that exceeds this value is converted to ref genotype) [default: 4.0]
   --min_split_diff_dup or -msdu <FLOAT>  minimum fold of the difference between the numbers of forward and reverse soft-clipped reads supporting <= 200 bp DUP (DUP non-ref genotype that exceeds this value is converted to ref genotype) [default: 4.0]
   --min_filt_rate_split or -mfs <FLOAT>  minimum rate of samples that exceed the value specified with --min_split_diff_* at an SV site (SV sites exceeding this value is removed) [default: 0.9]
   
   --nthread or -n <INT>     number of threads to use [default: 1]
   --dp_filter or -df        allow depth-based filter for DEL/DUP [default: FALSE]
   --dis_smc2 or -ds         cancel 2nd SMC step [default: FALSE]
   --help or -h              output help message
   
=cut

die "input vcf file is not specified: \n" if ($sv_vcf eq '');
die "inappropriate toolset option:\n" if ($tool_set !~ /^[479]tools|^11tools|^6tools_[12]/);

my $tool_set_list = "$data_dir/SVtool_param.txt";

my @tools = ();

@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham') if ($tool_set =~ /4tools/);
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'MELT') if ($tool_set =~ /6tools_1/);
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'GRIDSS') if ($tool_set =~ /6tools_2/);
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'GRIDSS', 'MELT') if ($tool_set =~ /7tools/);
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'DELLY', 'MELT', 'Lumpy', 'SoftSV') if ($tool_set =~ /9tools/);
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'DELLY', 'MELT', 'Lumpy', 'SoftSV', 'forestSV', 'Mobster') if ($tool_set =~ /11tools/);
if ((@tools == 0) and (-f $tool_set)){
    open (FILE, $tool_set);
    while (my $line = <FILE>){
        chomp $line;
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        push @tools, $line;
    }
    close (FILE);
}

my $training_data_dir = "$data_dir/R.nnet.models";
my $rscript = "$Bin//multiLR.R";

if ($non_human == 0){
	$simple_repeat = "$data_dir/simpleRepeat.txt.gz";
	$simple_repeat = "$data_dir/simpleRepeat.b38.txt.gz" if ($build == 38);

	$segdup_file = "$data_dir/genomicSuperDups.txt.gz";
	$segdup_file = "$data_dir/genomicSuperDups.b38.txt.gz" if ($build == 38);
}

my $Gid = '';
$Gid = $1 if ($sv_vcf =~ /^([WZ]\d+)/);

my $sv_dir = '.';
$sv_dir = $1 if ($sv_vcf =~ /(.+)\//);

my $sv_vcf_base = basename ($sv_vcf);

my @sample_list;
open (FILE, $sv_vcf) or die "$sv_vcf is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
		if ($line =~ /^#CHROM/){
		    my @header2 = split (/\s+/, $line);
		    my $count = 0;
		    foreach (@header2){
				$count ++;
				next if ($count <= 9);
				push @sample_list, $_;
		    }
		}
		next;
    }
    else{
		last;
    }
}
close (FILE);

my $total_samples = scalar @sample_list;
my $div_num = int ($total_samples / $thread) if ($total_samples % $thread == 0);
$div_num = int ($total_samples / $thread) + 1 if ($total_samples % $thread > 0);

print STDERR "Total sample: $total_samples\n";

my %div_samples;
my $sn = 0;
my $div = 1;
foreach (@sample_list){
    $sn ++;
    if ($sn <= $div_num){
		push @{$div_samples{$div}}, $_;
		if ($sn == $div_num){
		    $div ++;
		    $sn = 0;
		}
    }
}

my %toolset;
my %toolset2;
my %gender;

open (FILE, $tool_set_list) or die "$tool_set_list is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($tool, $type, $size, $read, $read2) = split (/\s+/, $line);
    ${${$toolset{$tool}}{$type}}{$size} = $read;
}
close (FILE);

foreach my $tool (@tools){
    if (!exists $toolset{$tool}){
        die "$tool parameter info is lacking in $data_dir/SVtool_param.txt\n";
    }
}

if ($gender_list ne ''){
    open (FILE, $gender_list) or die "$gender_list is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#/);
        my ($ID, $sex) = split (/\s+/, $line);
        $gender{$ID} = $sex;
    }
    close (FILE);
}

system ("mkdir $out_dir") if (!-d $out_dir);
my $temp_dir = "$out_dir/$out_prefix.tmp";
system ("rm -r $temp_dir") if (-d $temp_dir);
system ("mkdir $temp_dir");

my %repeat;
my %repeat2;
my %simple_repeat;
my %simple_repeat2;
my $Mbin_size = 1000000;
my $ins_sd2 = 50;
my $max_STR = 500;

my $total_repeat = 0;

if ($simple_repeat ne ''){
	open (FILE, "gzip -dc $simple_repeat |") or die "$simple_repeat is not found: $!\n" if ($simple_repeat =~ /\.gz$/);
    open (FILE, $simple_repeat) or die "$simple_repeat is not found: $!\n" if ($simple_repeat !~ /\.gz$/);
	while (my $line = <FILE>){
	    chomp $line;
	    next if ($line =~ /^$|^#/);
	    my @line = split (/\s+/, $line);
	    my $chr = $line[1];
	    $chr =~ s/^chr// if ($chr =~ /^chr/) and ($non_human == 0) and ($build == 37);
	    next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
	    next if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
	    my $start = $line[2];
	    my $end = $line[3];
        if (!exists ${$repeat{$chr}}{$start}){
    	    ${$repeat{$chr}}{$start} = $end;
    	    ${$simple_repeat{$chr}}{$start} = $end if ($end - $start + 1 <= $max_STR);
        }
        else{
            my $end2 = ${$repeat{$chr}}{$start};
            if ($end2 < $end){
                ${$repeat{$chr}}{$start} = $end;
                ${$simple_repeat{$chr}}{$start} = $end if ($end - $start + 1 <= $max_STR);
            }
        }
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
	    $chr =~ s/^chr// if ($chr =~ /^chr/) and ($non_human == 0) and ($build == 37);
	    next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
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
	$total_repeat += $end - $pos + 1;
    }
}
%repeat = ();

foreach my $chr (keys %simple_repeat){
    foreach my $pos (sort {$a <=> $b} keys %{$simple_repeat{$chr}}){
        my $end = ${$simple_repeat{$chr}}{$pos};
        my $Mbin_start = int ($pos / $Mbin_size);
        my $Mbin_end = int ($end / $Mbin_size);
        if ($Mbin_start == $Mbin_end){
            ${${$simple_repeat2{$chr}}{$Mbin_start}}{$pos} = $end;
        }
        else{
            my $Mbin = $Mbin_start;
            while ($Mbin <= $Mbin_end){
                ${${$simple_repeat2{$chr}}{$Mbin}}{$pos} = $end;
                $Mbin ++;
            }
        }
    }
}
%simple_repeat = ();

print STDERR "#Total length of repeats:\t$total_repeat\n";

my %MC;
my %MC2;
my %select;
my %select_dup;
my %vcf;
my %vcf2;
my %sv;
my %sv_dup;
my %INS;
my %DUP;
my %sample_order;
my %order_sample;
my %group;
my %sv_repeat;
my @header;
my $Missed_calls = 0;
my $correct_calls = 0;
my $ref_correct_calls = 0;
my $total_calls = 0;

my @MC_chr_outfiles;
my $pre_chr = '';
open (FILE, $sv_vcf) or die "$sv_vcf is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
		if ($line =~ /^#CHROM/){
			my @header2 = split (/\s+/, $line);
			my $count = 0;
			foreach (@header2){
				$count ++;
				next if ($count <= 9);
				$sample_order{$count} = $_;
				my ($group, $id) = split (/:/, $_) if ($_ =~ /:/);
				$id = $_ if ($_ !~ /:/);
				$group = '' if ($_ !~ /:/);
#print STDERR "NA12878: $count\n" if ($id eq 'NA12878');
				$group{$id} = $group if ($group ne '');
				$order_sample{$id} = $count;
			}
		}
		push @header, $line;
		next;
    }
    my @line = split (/\s+/, $line);
    my $chr = $line[0];
    next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
    my $pos = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    next if ($target_type ne 'ALL') and ($type ne $target_type);
    $type = 'INS' if ($type eq 'ALU') or ($type eq 'LINE1') or ($type eq 'SVA') or ($type eq 'HERVK') or ($type eq 'VEI') or ($type eq 'NUMT');
    my $subtype = $1 if ($line[4] =~ /<(.+)>/);
    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    my $end = $pos + $len - 1;
    $end = $pos if ($type eq 'INS');
    my $type2 = $type;
    $type2 = 'INS' if ($type eq 'DUP');
    if ($type eq 'DUP'){
    	${${$sv{$type2}}{$chr}}{$pos} = $pos;
    	${${$sv{$type2}}{$chr}}{$end} = $end;
    	${$sv_dup{$chr}}{$pos} = $end;
    	${$sv_dup{$chr}}{$end} = $pos;
    }
    else{
    	${${$sv{$type2}}{$chr}}{$pos} = $end;
    }
    my $Mbin_start = int ($pos / $Mbin_size);
    my $Mbin_end = int ($end / $Mbin_size);
    my $flag = 0;
    if ($Mbin_start == $Mbin_end){
		if (exists ${$repeat2{$chr}}{$Mbin_start}){
		    foreach my $rpos (sort {$a <=> $b} keys %{${$repeat2{$chr}}{$Mbin_start}}){
				last if ($rpos > $end + $ins_sd);
				my $rend = ${${$repeat2{$chr}}{$Mbin_start}}{$rpos};
				my $rlen = $rend - $rpos + 1;
				next if ($rend < $pos - $ins_sd);
				if ($type eq 'INS'){
				    $flag = 1;
				    last;
				}
				else{
					my $overlap = 0;
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
				    if ($overlap >= $len * $min_repeat_overlap){
						$flag = 1;
						last;
				    }
				}
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
				last if ($rpos > $end + $ins_sd);
				my $rend = $rep{$rpos};
				my $rlen = $rend - $rpos + 1;
				next if ($rend < $pos - $ins_sd);
				if ($type eq 'INS'){
				    $flag = 1;
				    last;
				}
				else{
					my $overlap = 0;
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
				    if ($overlap >= $len * $min_repeat_overlap){
						$flag = 1;
						last;
				    }
				}
		    }
		}
    }
    if ($flag == 1){
		${${$sv_repeat{$chr}}{$pos}}{$type} = 1;
    }
    my $sample_num = $1 if ($line[7] =~ /SC=(\d+)/);
    my $chr02d = $chr;
    $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
    ${${$vcf{$chr02d}}{$pos}}{$type} = $line;
    my @missed_call = ();
    my $count = 0;
    if ($total_samples >= 200){
		if (($pre_chr ne $chr) and ($pre_chr ne '')){
			my $MC_chr_out = "$temp_dir/MC.chr$pre_chr.txt";
			open (OUT, "> $MC_chr_out");
			foreach my $GID (sort keys %MC){
				foreach my $type2 (keys %{$MC{$GID}}){
					foreach my $chr2 (keys %{${$MC{$GID}}{$type2}}){
						foreach my $pos2 (sort {$a <=> $b} keys %{${${$MC{$GID}}{$type2}}{$chr2}}){
							my $len2 = ${${${$MC{$GID}}{$type2}}{$chr2}}{$pos2};
							print OUT "$GID\t$type2\t$pos2\t$len2\n";
						}
					}
				}
			}
			close (OUT);
			push @MC_chr_outfiles, $MC_chr_out;
            undef %MC;
		}
    }
    foreach (@line){
		$count ++;
		next if ($count <= 9);
		my $gid = $sample_order{$count};
		my ($group, $id) = split (/:/, $gid) if ($gid =~ /:/);
		$id = $gid if ($gid !~ /:/);
		if ($_ =~ /^0\/0/){
			push @missed_call, $gid;
		}
	}
	if (@missed_call > 0){
		foreach my $gid (@missed_call){
			${${${$MC{$gid}}{$type}}{$chr}}{$pos} = $len;
            $MC2{$gid} = 1;
			$Missed_calls ++;
		}
    }
    $pre_chr = $chr;
}
close (FILE);
if ($total_samples >= 200){
    if (scalar keys %MC > 0){
		my $MC_chr_out = "$temp_dir/MC.chr$pre_chr.txt";
		open (OUT, "> $MC_chr_out");
		foreach my $GID (sort keys %MC){
			foreach my $type2 (keys %{$MC{$GID}}){
				foreach my $chr2 (keys %{${$MC{$GID}}{$type2}}){
					foreach my $pos2 (sort {$a <=> $b} keys %{${${$MC{$GID}}{$type2}}{$chr2}}){
						my $len2 = ${${${$MC{$GID}}{$type2}}{$chr2}}{$pos2};
						print OUT "$GID\t$type2\t$pos2\t$len2\n";
					}
				}
			}
		}
		close (OUT);
		push @MC_chr_outfiles, $MC_chr_out;
        undef %MC;
    }
    foreach my $chr_out (@MC_chr_outfiles){
    	my $chr_out_base = basename ($chr_out);
		my $chr = $1 if ($chr_out_base =~ /chr(.+)\.txt$/);
		my $pre_gid = '';
		my $MC_gid_out = '';
		open (FILE, $chr_out) or die "$chr_out is not found: $!\n";
		while (my $line = <FILE>){
			chomp $line;
			my ($gid, $type, $pos, $len) = split (/\t/, $line);
			if ($pre_gid ne $gid){
				$MC_gid_out = "$temp_dir/$gid.MC.txt";
				close (OUT) if ($pre_gid ne '');
				open (OUT, ">> $MC_gid_out");
			}
			print OUT "$type\t$chr\t$pos\t$len\n";
			$pre_gid = $gid;
		}
		close (OUT);
		close (FILE);
    }
}

foreach my $n (keys %sample_order){
    my $gid = $sample_order{$n};
    my $MC_temp = "$temp_dir/$gid.MC.txt";
    if (!-f $MC_temp){
        open (OUT, "> $MC_temp");
        close (OUT);
    }
}

print STDERR "1st step completed:\n";

my @jobs;
my @out_files;
foreach my $dnum (sort {$a <=> $b} keys %div_samples){
    my ($thread_t) = threads->new(\&collect_call, $dnum, \@{$div_samples{$dnum}});
    push @jobs, $thread_t;
}
foreach (@jobs){
    my ($file1) = $_->join;
    push @out_files, $file1 if (-f $file1);
}
undef @jobs;

undef %MC;

print STDERR "2nd step completed:\n";

foreach my $out_file (@out_files){
    open (FILE, $out_file) or die "$out_file is not found: $!\n";
    while (my $line = <FILE>){
		chomp $line;
		my ($chr, $epos, $type, $id, $info) = split (/\t/, $line);
		push @{${${${$select{$chr}}{$epos}}{$type}}{$id}}, $info;
    }
    close (FILE);
}

foreach my $chr (keys %select){
    foreach my $pos (keys %{$select{$chr}}){
		foreach my $type (keys %{${$select{$chr}}{$pos}}){
			foreach my $id (keys %{${${$select{$chr}}{$pos}}{$type}}){
				my $CNVnator_flag = 0;
				my $large_CNV_flag = 0;
				foreach (@{${${${$select{$chr}}{$pos}}{$type}}{$id}}){
					my ($spos, $len, $tool) = split (/=/, $_);
					$CNVnator_flag = 1 if ($tool eq 'CNVnator');
					$large_CNV_flag = 1 if ($len >= 50000);
				}
				if (($large_CNV_flag == 1) and ($CNVnator_flag == 0)){
					delete ${${${$select{$chr}}{$pos}}{$type}}{$id};
					if (scalar keys %{${${$select{$chr}}{$pos}}{$type}} == 0){
						delete ${${$select{$chr}}{$pos}}{$type};
					}
					if (scalar keys %{${$select{$chr}}{$pos}} == 0){
						delete ${$select{$chr}}{$pos};
					}
					next;
				}
				$correct_calls ++;
			}
		}
    }
}

print STDERR "MC: $Missed_calls\tCorrected MC: $correct_calls\n";

print STDERR "3rd step completed:\n";

my @MCcor_chr_out;
my %vcf_dr;
my $min_ld = 0.8;
my $min_ld2 = 0.1;

foreach my $chr (sort keys %vcf){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (sort keys %{${$vcf{$chr}}{$pos}}){
		    my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
		    my $sample_num = $1 if ($line[7] =~ /SC=(\d+)/);
		    my $sample_num2 = $sample_num;
		    my $SF = int ($sample_num / $total_samples * 1000 + 0.5) / 1000;
		    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		    my $subtype = $line[4];
		    my $count = 0;
            my $maxlen = 0;
            my $maxend = 0;
            my $minend = 999999999999;
            my $maxpos = 0;
            my $minpos = 0;
            my $medpos = 0;
			my @sumlen = ();
			my @sumpos = ();
			foreach (@line){
			    $count ++;
			    next if ($count <= 9);
			    if ($_ !~ /^0\/0/){
					my ($idgt, $idgq, $idpos, $idlen) = split (/:/, $_);
					push @sumlen, $idlen;
					push @sumpos, $idpos;
                    $maxlen = $idlen if ($maxlen < $idlen);
                    my $idend = $idpos + $idlen - 1;
                    $maxend = $idend if ($idend > $maxend);
                    $minend = $idend if ($idend < $minend);
			    }
			}
			my $sumlen = 0;
			map{$sumlen += $_} @sumlen;
			my $avelen = int ($sumlen / @sumlen);
			my $newlen = $avelen;
			my $halfnum = int (@sumpos * 0.5);
			$medpos = $sumpos[$halfnum];
			if ($type eq 'INS'){
                if ($maxlen >= $min_sv_len){
                    $newlen = $maxlen;
                }
                else{
                    $newlen = 0;
                }
            }
			else{
			    if (@sumlen > 2){
					$newlen = &cons_len (\@sumlen, $avelen);
			    }
			}

		    if (exists ${${$select{$chr2}}{$pos}}{$type}){
				my $sample_count = 0;
				my %id_info;
				my $id_num = scalar keys %{${${$select{$chr2}}{$pos}}{$type}};
				foreach my $id (sort keys %{${${$select{$chr2}}{$pos}}{$type}}){
				    my $cons_pos = 0;
				    my $cons_len = 0;
				    my $cons_read = 0;
				    my $cons_gt = './.';
				    my $tool_str = '';
				    if (@{${${${$select{$chr2}}{$pos}}{$type}}{$id}} == 1){
						($cons_pos, $cons_len, $tool_str, $cons_gt, $cons_read) = split (/=/, ${${${${$select{$chr2}}{$pos}}{$type}}{$id}}[0]);
						$tool_str .= "-$cons_pos-$cons_len";
				    }
				    else{
						my @cons_pos = ();
						my @cons_len = ();
						my @cnvnator_pos = ();
						my @cnvnator_len = ();
						my @cons_read;
						my $sum_pos = 0;
						my $sum_len = 0;
						my $sum_read = 0;
						my %gt;
						foreach my $info (@{${${${$select{$chr2}}{$pos}}{$type}}{$id}}){
						    my ($tpos, $tlen, $tool, $gt, $read) = split (/=/, $info);
						    if ($tool =~ /CNVnator/i){
								push @cnvnator_pos, $tpos;
								push @cnvnator_len, $tlen;
						    }
						    else{
								push @cons_pos, $tpos;
								push @cons_len, $tlen;
						    }
						    push @cons_read, $read;
						    $gt{$gt} ++;
						    $tool_str .= "$tool-$tpos-$tlen,";
						}
						if (@cons_pos == 0){
						    @cons_pos = (@cnvnator_pos);
						    @cons_len = (@cnvnator_len);
						}
						map{$sum_pos += $_} @cons_pos;
						map{$sum_len += $_} @cons_len;
						map{$sum_read += $_} @cons_read;
						$cons_pos = int ($sum_pos / @cons_pos + 0.5);
						$cons_len = int ($sum_len / @cons_len + 0.5);
						$cons_read = int ($sum_read / @cons_read + 0.5);
						if (($type ne 'INS') and (@cons_len > 2)){
						    my ($new_len) = &cons_len (\@cons_len, $cons_len);
						    $cons_len = $new_len if ($new_len > 0);
						}
						foreach my $gt (sort{$gt{$b} <=> $gt{$a}} keys %gt){
						    next if (scalar keys %gt > 1) and ($gt eq './.');
						    $cons_gt = $gt;
						    last;
						}
				    }
				    $tool_str =~ s/,$//;
				    $sample_count ++;
                    my $dr = 1;
                    my $ds = 0;
				    my $sr = 0;
                    ${${$vcf_dr{$id}}{$chr2}}{$cons_pos} = "$cons_len\t$type";
				    $id_info{$id} = "$cons_gt:100:$cons_pos:$cons_len:$dr:$ds:$sr:1";
				}
                $sample_num += $sample_count;
                my $SF2 = int ($sample_num / $total_samples * 1000 + 0.5) / 1000;
				my $count = 0;
				foreach (@line){
				    $count ++;
				    next if ($count <= 9);
				    my $gid = $sample_order{$count};
				    my ($group, $id2) = split (/:/, $gid) if ($gid =~ /:/);
				    $id2 = $gid if ($gid !~ /:/);
				    if (exists $id_info{$id2}){
						$line[$count - 1] = $id_info{$id2};
				    }
				    else{
#						if (($_ =~ /^0\/0/) and ((($type eq 'DUP') and ($SF >= $min_SF)) or (($type ne 'DUP') and ($SF >= $min_SF2)))){
#						if (($_ =~ /^0\/0/) and ($type eq 'DUP') and ($SF >= $min_SF)){
						if (($disable_SMC_2nd == 0)  and ($_ =~ /^0\/0/) and ($sample_num2 >= $min_SC_DEL) and ($newlen >= $min_SMC_DEL) and ($type =~ /DEL|DUP/)){
							${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
							$line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
					    }
					    elsif (($disable_SMC_2nd == 0) and ($_ =~ /^0\/0/) and ($SF >= $min_SF2) and ($sample_num2 >= $min_SC) and ($type eq 'INV') and ($newlen >= 1000)){
					    	${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
							$line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
					    }
						elsif (($disable_SMC_2nd == 0) and ($_ =~ /^0\/0/) and ($SF >= $min_SF) and ($sample_num2 >= $min_SC) and ($type eq 'INS')){
						    ${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
						    $line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
						}
						else{
						    $line[$count - 1] .= ':0';
						}
				    }
				}
				$line[7] =~ s/SC=\d+/SC=$sample_num/;
				$count = 0;
                my $max_len = 0;
				my @sum_len = ();
				foreach (@line){
				    $count ++;
				    next if ($count <= 9);
				    if ($_ !~ /^0\/0/){
						my @idinfo = split (/:/, $_);
						push @sum_len, $idinfo[3];
	                    $max_len = $idinfo[3] if ($max_len < $idinfo[3]);
				    }
				}
				my $sum_len = 0;
				map{$sum_len += $_} @sum_len;
				my $ave_len = int ($sum_len / @sum_len);
				my $new_len = $ave_len;
				if ($type eq 'INS'){
                    if ($maxlen >= $min_sv_len){
                        $new_len = $maxlen;
                    }
                    else{
                        $new_len = 0;
                    }
                }
				else{
				    if (@sum_len > 2){
						$new_len = &cons_len (\@sum_len, $ave_len);
				    }
				}
				$new_len = 0 - $new_len if ($type eq 'DEL');       
				$line[7] =~ s/SVLEN=-*\d+/SVLEN=$new_len/;
		    }
		    else{
				my $count = 0;
				foreach (@line){
				    $count ++;
				    next if ($count <= 9);
				    my $gid = $sample_order{$count};
				    my ($group, $id2) = split (/:/, $gid) if ($gid =~ /:/);
				    $id2 = $gid if ($gid !~ /:/);
#				    if (($_ =~ /^0\/0/) and ((($type eq 'DUP') and ($SF >= $min_SF)) or (($type ne 'DUP') and ($SF >= $min_SF2)))){
#				    if (($_ =~ /^0\/0/) and ($type eq 'DUP') and ($SF >= $min_SF)){
					if (($disable_SMC_2nd == 0)  and ($_ =~ /^0\/0/) and ($sample_num2 >= $min_SC_DEL) and ($newlen >= $min_SMC_DEL) and ($type =~ /DEL|DUP/)){
						${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
						$line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
				    }
				    elsif (($disable_SMC_2nd == 0) and ($_ =~ /^0\/0/) and ($SF >= $min_SF2) and ($sample_num2 >= $min_SC) and ($type eq 'INV') and ($newlen >= 1000)){
				    	${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
						$line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
				    }
					elsif (($disable_SMC_2nd == 0) and ($_ =~ /^0\/0/) and ($SF >= $min_SF) and ($sample_num2 >= $min_SC) and ($type eq 'INS')){
					    ${${$vcf_dr{$id2}}{$chr2}}{$medpos} = "$newlen\t$type";
					    $line[$count - 1] = "0/0:0:$medpos:$newlen:1:0:0:2";
					}
				    else{
						$line[$count - 1] .= ':0';
				    }
				}
		    }
			$line[8] = 'GT:GQ:VP:VL:DR:DS:SR:MC';     
		    $line[7] =~ s/;$// if ($line[7] =~ /;$/);
		    my $new_line = join ("\t", @line);
		    ${${$vcf{$chr}}{$pos}}{$type} = $new_line;
		}
    }
}
close (OUT);
undef %select;

my %dp_rate_input;
my %dp_rate;

foreach my $id (sort keys %vcf_dr){
    my $gr = '';
    $gr = $group{$id} if (exists $group{$id});
    my $out_pos = "$temp_dir/$id.dr-input.txt";
    open (OUT1, "> $out_pos");
    foreach my $chr (sort keys %{$vcf_dr{$id}}){
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf_dr{$id}}{$chr}}){
        	my ($len, $type) = split (/\t/, ${${$vcf_dr{$id}}{$chr}}{$pos});
	        print OUT1 "$chr\t$pos\t$len\t$type\n";
        }
    }
    close (OUT1);
    my $cov_dir = "$sample_dir/$gr/$id/Cov" if ($gr ne '');
    $cov_dir = "$sample_dir/$id/Cov" if ($gr eq '') or (!-d $cov_dir);
    $cov_dir = "$sample_dir/Cov" if (!-d $cov_dir);
    if (!-d $cov_dir){
        print STDERR "Warning $sample_dir/$gr/$id/Cov directory is not found:\n";
        next;
    }
    $dp_rate_input{$id} = "$out_pos=$cov_dir";
}

#print STDERR "MC keys: ", scalar keys %MC, "\n";

my @out_files2;
my @jobs2;
foreach my $dnum (sort {$a <=> $b} keys %div_samples){
    my ($thread_t) = threads->new(\&add_dp_rate, $dnum, \@{$div_samples{$dnum}});
    push @jobs2, $thread_t;
}
foreach (@jobs2){
    my ($file1) = $_->join;
    push @out_files2, $file1 if (-f $file1);
}
undef @jobs2;

foreach my $file (@out_files2){
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($id, $chr, $pos, $cov_rate, $incons_covr, $split_rate, $type, $split_diff_fold) = split (/\t/, $line);
        ${${${$dp_rate{$id}}{$chr}}{$pos}}{$type} = "$cov_rate=$incons_covr=$split_rate=$split_diff_fold";
    }
    close (FILE);
}

undef %MC;
undef %MC2;

$pre_chr = '';

my %rm_sv;
my %DR_boundary = ('DEL' => 0.3, 'DELr' => 0.35, 'DUP' => 1.65, 'DUPr' => 1.8);
my %SR_boundary = ('DEL' => 0.9, 'DELr' => 0.9, 'DUP' => 0.5, 'DUPr' => 0.55, 'INS' => 0.85, 'INSr' => 0.85, 'INV' => 1.1, 'INVr' => 1.1);
my %DR_het = ('DEL' => 0.54, 'DELr' => 0.6, 'DUP' => 1.33, 'DUPr' => 1.5);
my %SR_het = ('DEL' => 0.43, 'DELr' => 0.41, 'DUP' => 0.3, 'DUPr' => 0.35, 'INS' => 0.5, 'INSr' => 0.5, 'INV' => 0.45, 'INVr' => 0.4);
=pod
my %DR_boundary = ('DEL' => 0.35, 'DUP' => 1.9);
my %SR_boundary = ('DEL' => 1, 'DUP' => 1, 'INS' => 0.97, 'INV' => 1.15);
my %DR_het = ('DEL' => 0.56, 'DUP' => 1.48);
my %SR_het = ('DEL' => 0.5, 'DUP' => 0.58, 'INS' => 0.47, 'INV' => 0.57);
=cut
my %mLR = ('DEL1' => "model_del.r", 'DEL2' => "model_del1k.r", 'DEL3' => "model_del150bp.r", 'DUP1' => "model_dup.r", 'DUP2' => "model_dup1k.r", 'DUP3' => "model_dup100bp.r", 'INS' => "model_ins.r", 'INV' => "model_inv.r");
my %mLR_r = ('DEL1' => "model_delR.r", 'DEL2' => "model_del1kR.r", 'DEL3' => "model_del150bpR.r", 'DUP1' => "model_dupR.r", 'DUP2' => "model_dup1kR.r", 'DUP3' => "model_dup100bpR.r", 'INS' => "model_insR.r", 'INV' => "model_invR.r");

my %GQ_filt = ('DEL' => 14, 'DUP' => 12, 'INS' => 16, 'INV' => 14);
my %GQr_filt = ('DEL' => 6, 'DUP' => 0, 'INS' => 0, 'INV' => 0);
my %GQDEV_filt = ('DEL' => 1.5, 'DUP' => 2, 'INS' => 2, 'INV' => 3);
my %GQDEVr_filt = ('DEL' => 2.5, 'DUP' => 2, 'INS' => 2, 'INV' => 3);
my %MSC_filt = ('DEL' => 0.63, 'DUP' => 0.54, 'INS' => 0.6, 'INV' => 1);
my %MSCr_filt = ('DEL' => 0.67, 'DUP' => 1, 'INS' => 0.88, 'INV' => 1);

my $min_len = 1000;
my $min_del_len = 500;
my $min_del_len2 = 10000;
my $min_dup_len = 200;
my $max_del_dprate = 0.85;
my $max_del_dprate2 = 0.7;
my $max_del_dprate3 = 0.9;
my $min_dup_dprate = 1.35;
my $min_dup_dprate3 = 1.25;
my $max_del_dsr = 0.2;
my $max_dup_dsr = 0.1;
my $del_filt = 0;
my $dup_filt = 0;
my $del_num = 0;
my $dup_num = 0;

foreach my $chr (sort keys %vcf){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    my %new_line;
    my %gt_pred;
    my %DR_ave;
    my %SR_ave;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (sort keys %{${$vcf{$chr}}{$pos}}){
		    my $rep_flag = 0;
		    $rep_flag = 1 if (exists ${${$sv_repeat{$chr2}}{$pos}}{$type});
		    my $type2 = $type;
		    $type2 .= 'r' if ($rep_flag == 1);
            my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
            my $SN = $1 if ($line[7] =~ /SC=(\d+)/);
            my $SF = int ($SN / $total_samples * 1000 + 0.5) / 1000;
            my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		    my $type3 = $type;
		    if ($type eq 'DEL'){
				if ($len <= 150){
				    $type3 = $type . '3';
				}
				elsif ($len <= 1000){
				    $type3 = $type . '1';
				}
				else{
				    $type3 = $type . '2';
				}
		    }
		    elsif ($type eq 'DUP'){
				if ($len <= 100){
				    $type3 = $type . '3';
				}
				elsif ($len <= 2000){
				    $type3 = $type . '1';
				}
				else{
				    $type3 = $type . '2';
				}
		    }
		    $type3 .= 'r' if ($rep_flag == 1);
		    
            $del_num ++ if ($type eq 'DEL');
            $dup_num ++ if ($type eq 'DUP');
		    my $het_sr = $SR_het{$type2};
		    my $het_dr = 0;
		    $het_dr = $DR_het{$type2} if ($type =~ /DEL|DUP/);
		    my $sr_boundary = $SR_boundary{$type2};
		    my $dr_boundary = 0;
		    $dr_boundary = $DR_boundary{$type2} if ($type =~ /DEL|DUP/);
		    if ($type3 =~ /DUP1/){
				$dr_boundary += 0.15;
				$sr_boundary += 0.15;
				$het_dr += 0.15;
				$het_sr += 0.15;
		    }
		    my $dr_ave = 1;
		    my $sr_ave = 0;
		    my @dr;
		    my @sr;
		    my @dr2;
		    my @sr2;
		    my $count = 0;
		    foreach (@line){
				$count ++;
				next if ($count <= 9);
				next if ($_ =~ /^0\/0/);
				my $gid = $sample_order{$count};
				my ($group, $id2) = split (/:/, $gid) if ($gid =~ /:/);
				$id2 = $gid if ($gid !~ /:/);
				my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
				($dr, $ds, $sr) = split (/=/, ${${${$dp_rate{$id2}}{$chr2}}{$vp}}{$type}) if (exists ${${${$dp_rate{$id2}}{$chr2}}{$vp}}{$type});
				next if ($dr >= 4);
				next if ($dr >= 1) and ($type eq 'DEL');
				next if ($type eq 'DUP') and ($vl >= 200) and ($dr <= 1);
				next if ($type eq 'DUP') and ($vl <= 1000) and ($sr < 0.1) and ($mc >= 1);
				next if ($type eq 'INS') and ($sr < 0.02);
				next if ($type eq 'INV') and ($sr < 0.1);
				push @dr, $dr if ($type =~ /DEL|DUP/);
				push @sr, $sr;
				if ($mc == 0){
				    push @dr2, $dr if ($type =~ /DEL|DUP/);
				    push @sr2, $sr;
				}
		    }
		    my $sum = 0;
		    map{$sum += $_} @dr;
		    $dr_ave = int ($sum / @dr * 100 + 0.5) / 100 if (@dr > 0);
		    $sum = 0;
		    map{$sum += $_} @sr;
		    $sr_ave = int ($sum / @sr * 100 + 0.5) / 100 if (@sr > 0);
		    if (@dr2 > 0){
				$sum = 0;
				map{$sum += $_} @dr2;
				my $dr_ave2 = int ($sum / @dr2 * 100 + 0.5) / 100;
				${$DR_ave{$type}}{$pos} = $dr_ave2;
		    }
		    if (@sr2 > 0){
				$sum = 0;
				map{$sum += $_} @sr2;
				my $sr_ave2 = int ($sum / @sr2 * 100 + 0.5) / 100;
				${$SR_ave{$type}}{$pos} = $sr_ave2;
		    }
		    $count = 0;
		    my $sc = 0;
		    foreach (@line){
				$count ++;
				next if ($count <= 9);
				my $gid = $sample_order{$count};
				my ($group, $id2) = split (/:/, $gid) if ($gid =~ /:/);
				$id2 = $gid if ($gid !~ /:/);
				my ($gt, $gq, $tpos, $tlen, $dr, $ds, $sr, $MCtag) = split (/:/, $_);
				my $sr_diff = 0;
				if (($_ =~ /^0\/0:/) and (!exists ${${${$dp_rate{$id2}}{$chr2}}{$tpos}}{$type})){
				    $line[$count - 1] = "0/0:0:0:0:0:0:0:0";
				}
				else{
				    $gq = 99 if ($tpos > 0);
				    ($dr, $ds, $sr, $sr_diff) = split (/=/, ${${${$dp_rate{$id2}}{$chr2}}{$tpos}}{$type}) if (exists ${${${$dp_rate{$id2}}{$chr2}}{$tpos}}{$type});
                    if (($type eq 'INS') and ($sr_diff >= $min_split_diff)){
                        $line[$count - 1] = "0/0:0:0:0:$dr:$ds:0:3";
						next;
                    }
                    elsif (($type eq 'DUP') and ($len <= 200) and ($sr_diff >= $min_split_diff_dup)){
                    	$line[$count - 1] = "0/0:0:0:0:$dr:$ds:0:3";
						next;
                    }
                    elsif (($type eq 'DEL') and ($len <= 500) and ($sr_diff >= $min_split_diff_del)){
                    	$line[$count - 1] = "0/0:0:0:0:$dr:$ds:0:3";
						next;
                    }
#print STDERR "$count\t$id2\t$gt\t$tpos\t$dr\t$sr\t$MCtag\n";
				    if ($tpos == 0){
						$tpos = $pos;
						$tlen = $len;
				    }
				    if ($type eq 'DEL'){
						if (($MCtag >= 1) and ($dr >= 1) and ($len >= 200)){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($MCtag >= 2) and (($dr > 0.65) or ($ds >= 0.5))){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($MCtag >= 2) and ($len <= 1000) and ($sr < 0.03)){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif ($gt eq './.'){
						    if ($DP_filter == 1){
								if ($dr >= 0.8){
								    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
								    next;
								}
								elsif (($len >= 1000) and ($ds > 0.25)){
								    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
								    next;
								}
						    }
						    if (($dr >= $dr_boundary + 0.05) and ($dr <= $het_dr + 0.1)){
								$gt = '0/1';
						    }
						    elsif ($dr <= $dr_boundary - 0.05){
								$gt = '1/1';
						    }
						}
						${${$gt_pred{$type3}}{$pos}}{$count} = "$dr=$sr=$dr_ave=$sr_ave";
				    }
				    elsif ($type eq 'DUP'){
						if (($tlen >= 200) and ($dr <= 1)){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($tlen <= 1000) and ($sr < 0.1) and ($MCtag >= 1)){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						if (($MCtag == 2) and ($dr < 1.3)){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif ($gt ne '0/0'){
						    if ($DP_filter == 1){
								if (($dr <= 1) or (($len >= $min_dup_len) and ($dr < $min_dup_dprate) and ($ds > $max_dup_dsr))){
								    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
								    next;
								}
								elsif (($len >= $min_len) and ($dr < $min_dup_dprate3)){
								    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
								    next;
								}
						    }
						    if (($dr >= $dr_boundary + 0.05) and ($len >= 2000)){
								$gt = '1/1';
						    }
						    elsif (($dr <= $dr_boundary - 0.05) and ($dr >= $het_dr - 0.2) and ($len >= 2000)){
								$gt = '0/1';
						    }
						    else{
								$gt = './.';
						    }
						}
						${${$gt_pred{$type3}}{$pos}}{$count} = "$dr=$sr=$dr_ave=$sr_ave";
				    }
				    elsif ($type =~ /INS|INV/){
						if (($sr < 0.1) and ($type eq 'INV')){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($sr < 0.1) and ($MCtag <= 1) and ($type eq 'INS')){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($sr < 0.2) and ($MCtag >= 2) and ($type eq 'INS')){
						    $line[$count - 1] = "0/0:0:0:0:$dr:$ds:$sr:0";
						    next;
						}
						elsif (($gt eq './.') and ($type eq 'INS')){
						    if (($sr <= $sr_boundary - 0.1) and ($sr >= $het_sr - 0.1)){
								$gt = '0/1';
						    }
						}
						elsif (($gt ne '0/0') and ($type eq 'INV')){
						    if ($sr >= $sr_boundary + 0.1){
								$gt = '1/1';
						    }
						    elsif (($sr <= $sr_boundary) and ($sr >= $het_sr - 0.2)){
								$gt = '0/1';
						    }
						    else{
								$gt = './.';
						    }
						}
						${${$gt_pred{$type3}}{$pos}}{$count} = "$dr=$sr=$dr_ave=$sr_ave";
				    }
				    $line[$count - 1] = "$gt:$gq:$tpos:$tlen:$dr:$ds:$sr:$MCtag";
				    $sc ++ if ($gt ne '0/0');
				}
		    }
		    next if ($sc == 0);
		    my $new_line = join ("\t", @line);
		    ${$new_line{$pos}}{$type} = $new_line;
		}
    }
    my %pred;
    my $del_pred_input = "$temp_dir/chr$chr.del.pred.input.txt";
    my $del1k_pred_input = "$temp_dir/chr$chr.del1k.pred.input.txt";
    my $del150_pred_input = "$temp_dir/chr$chr.del150.pred.input.txt";
    my $dup_pred_input = "$temp_dir/chr$chr.dup.pred.input.txt";
    my $dup1k_pred_input = "$temp_dir/chr$chr.dup1k.pred.input.txt";
    my $dup100_pred_input = "$temp_dir/chr$chr.dup100.pred.input.txt";
    my $ins_pred_input = "$temp_dir/chr$chr.ins.pred.input.txt";
    my $inv_pred_input = "$temp_dir/chr$chr.inv.pred.input.txt";
    
    my $del_pred_inputR = "$temp_dir/chr$chr.del.pred.inputR.txt";
    my $del1k_pred_inputR = "$temp_dir/chr$chr.del1k.pred.inputR.txt";
    my $del150_pred_inputR = "$temp_dir/chr$chr.del150.pred.inputR.txt";
    my $dup_pred_inputR = "$temp_dir/chr$chr.dup.pred.inputR.txt";
    my $dup1k_pred_inputR = "$temp_dir/chr$chr.dup1k.pred.inputR.txt";
    my $dup100_pred_inputR = "$temp_dir/chr$chr.dup100.pred.inputR.txt";
    my $ins_pred_inputR = "$temp_dir/chr$chr.ins.pred.inputR.txt";
    my $inv_pred_inputR = "$temp_dir/chr$chr.inv.pred.inputR.txt";
    
    my $del_pred_out = "$temp_dir/chr$chr.del.pred.out.txt";
    my $del1k_pred_out = "$temp_dir/chr$chr.del1k.pred.out.txt";
    my $del150_pred_out = "$temp_dir/chr$chr.del150.pred.out.txt";
    my $dup_pred_out = "$temp_dir/chr$chr.dup.pred.out.txt";
    my $dup1k_pred_out = "$temp_dir/chr$chr.dup1k.pred.out.txt";
    my $dup100_pred_out = "$temp_dir/chr$chr.dup100.pred.out.txt";
    my $ins_pred_out = "$temp_dir/chr$chr.ins.pred.out.txt";
    my $inv_pred_out = "$temp_dir/chr$chr.inv.pred.out.txt";
    
    my $del_pred_outR = "$temp_dir/chr$chr.del.pred.outR.txt";
    my $del1k_pred_outR = "$temp_dir/chr$chr.del1k.pred.outR.txt";
    my $del150_pred_outR = "$temp_dir/chr$chr.del150.pred.outR.txt";
    my $dup_pred_outR = "$temp_dir/chr$chr.dup.pred.outR.txt";
    my $dup1k_pred_outR = "$temp_dir/chr$chr.dup1k.pred.outR.txt";
    my $dup100_pred_outR = "$temp_dir/chr$chr.dup100.pred.outR.txt";
    my $ins_pred_outR = "$temp_dir/chr$chr.ins.pred.outR.txt";
    my $inv_pred_outR = "$temp_dir/chr$chr.inv.pred.outR.txt";
    
    open (OUT1, "> $del_pred_input");
    open (OUT2, "> $del1k_pred_input");
    open (OUT3, "> $dup_pred_input");
    open (OUT4, "> $dup1k_pred_input");
    open (OUT5, "> $ins_pred_input");
    open (OUT6, "> $inv_pred_input");
    open (OUT7, "> $del150_pred_input");
    open (OUT8, "> $dup100_pred_input");
    open (OUT1R, "> $del_pred_inputR");
    open (OUT2R, "> $del1k_pred_inputR");
    open (OUT3R, "> $dup_pred_inputR");
    open (OUT4R, "> $dup1k_pred_inputR");
    open (OUT5R, "> $ins_pred_inputR");
    open (OUT6R, "> $inv_pred_inputR");
    open (OUT7R, "> $del150_pred_inputR");
    open (OUT8R, "> $dup100_pred_inputR");
    print OUT1 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT2 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT3 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT4 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT5 "GT\tSR\tMSR\tPOS\n";
    print OUT6 "GT\tSR\tMSR\tPOS\n";
    print OUT7 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT8 "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT1R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT2R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT3R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT4R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT5R "GT\tSR\tMSR\tPOS\n";
    print OUT6R "GT\tSR\tMSR\tPOS\n";
    print OUT7R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    print OUT8R "GT\tDR\tSR\tMDR\tMSR\tPOS\n";
    foreach my $type2 (keys %gt_pred){
		my %test_order;
		my $count = 0;
		foreach my $pos2 (sort {$a <=> $b} keys %{$gt_pred{$type2}}){
		    foreach my $count2 (sort {$a <=> $b} keys %{${$gt_pred{$type2}}{$pos2}}){
				$count ++;
				$test_order{$count} = "$pos2=$count2";
				my ($dr, $sr, $mdr, $msr) = split (/=/, ${${$gt_pred{$type2}}{$pos2}}{$count2});
				if ($type2 eq 'DEL1'){
				    print OUT1 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DEL1r'){
				    print OUT1R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DEL2'){
				    print OUT2 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DEL2r'){
				    print OUT2R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DEL3'){
				    print OUT7 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DEL3r'){
				    print OUT7R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP1'){
				    print OUT3 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP1r'){
				    print OUT3R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP2'){
				    print OUT4 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP2r'){
				    print OUT4R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP3'){
				    print OUT8 "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'DUP3r'){
				    print OUT8R "0\t$dr\t$sr\t$mdr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'INS'){
				    print OUT5 "0\t$sr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'INSr'){
				    print OUT5R "0\t$sr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'INV'){
				    print OUT6 "0\t$sr\t$msr\t$pos2:$type2:$count2\n";
				}
				elsif ($type2 eq 'INVr'){
				    print OUT6R "0\t$sr\t$msr\t$pos2:$type2:$count2\n";
				}
		    }
		}
		my $type3 = $type2;
		$type3 =~ s/r$// if ($type3 =~ /r$/);
		my $model_file = "$training_data_dir/$mLR{$type3}";
		$model_file = "$training_data_dir/$mLR_r{$type3}" if ($type2 =~ /r$/);
		$type3 =~ s/\d$// if ($type3 =~ /\d$/);
		my @out_pred;
		my $test_file = '';
		my $out_file = '';
		if ($type2 eq 'DEL1'){
		    $test_file = $del_pred_input;
		    $out_file = $del_pred_out;
		}
		elsif ($type2 eq 'DEL1r'){
		    $test_file = $del_pred_inputR;
		    $out_file = $del_pred_outR;
		}
		elsif ($type2 eq 'DEL2'){
		    $test_file = $del1k_pred_input;
		    $out_file = $del1k_pred_out;
		}
		elsif ($type2 eq 'DEL2r'){
		    $test_file = $del1k_pred_inputR;
		    $out_file = $del1k_pred_outR;
		}
		elsif ($type2 eq 'DEL3'){
		    $test_file = $del150_pred_input;
		    $out_file = $del150_pred_out;
		}
		elsif ($type2 eq 'DEL3r'){
		    $test_file = $del150_pred_inputR;
		    $out_file = $del150_pred_outR;
		}
		elsif ($type2 eq 'DUP1'){
		    $test_file = $dup_pred_input;
		    $out_file = $dup_pred_out;
		}
		elsif ($type2 eq 'DUP1r'){
		    $test_file = $dup_pred_inputR;
		    $out_file = $dup_pred_outR;
		}
		elsif ($type2 eq 'DUP2'){
		    $test_file = $dup1k_pred_input;
		    $out_file = $dup1k_pred_out;
		}
		elsif ($type2 eq 'DUP2r'){
		    $test_file = $dup1k_pred_inputR;
		    $out_file = $dup1k_pred_outR;
		}
		elsif ($type2 eq 'DUP3'){
		    $test_file = $dup100_pred_input;
		    $out_file = $dup100_pred_out;
		}
		elsif ($type2 eq 'DUP3r'){
		    $test_file = $dup100_pred_inputR;
		    $out_file = $dup100_pred_outR;
		}
		elsif ($type2 eq 'INS'){
		    $test_file = $ins_pred_input;
		    $out_file = $ins_pred_out;
		}
		elsif ($type2 eq 'INSr'){
		    $test_file = $ins_pred_inputR;
		    $out_file = $ins_pred_outR;
		}
		elsif ($type2 eq 'INV'){
		    $test_file = $inv_pred_input;
		    $out_file = $inv_pred_out;
		}
		elsif ($type2 eq 'INVr'){
		    $test_file = $inv_pred_inputR;
		    $out_file = $inv_pred_outR;
		}
	#print STDERR "chr$chr: $type2\t$model_file $test_file $out_file\n";
		system ("Rscript --vanilla $rscript $model_file $test_file $out_file");
		if (-f $out_file){
		    if (scalar keys %test_order > 1){
				open (FILE2, $out_file);
				while (my $line2 = <FILE2>){
				    chomp $line2;
				    my ($num, $ref, $het, $hom) = split (/\s+/, $line2);
				    my %prob;
				    $prob{0} = $ref;
				    $prob{1} = $het;
				    $prob{2} = $hom;
				    my $top_gt = 0;
				    my $sec_gt = 0;
				    my $third_gt = 0;
				    my $top_prob = 0;
				    my $sec_prob = 0;
				    my $third_prob = 0;
				    foreach my $gt (sort {$prob{$b} <=> $prob{$a}} keys %prob){
						if ($top_prob == 0){
						    $top_prob = $prob{$gt};
						    $top_gt = $gt;
						    next;
						}
						if ($sec_prob == 0){
						    $sec_prob = $prob{$gt};
						    $sec_gt = $gt;
						}
						if ($third_prob == 0){
						    $third_prob = $prob{$gt};
						    $third_gt = $gt;
						}
				    }
				    my $top_pvalue = 1 - $top_prob;
				    my $sec_pvalue = 1 - $sec_prob;
				    my $third_pvalue = 1 - $third_prob;
				    my $GQ = 99;
				    my $GQ2 = 99;
				    if ($top_pvalue > 0){
						my $top_gl = int (-10 * (log ($top_pvalue) / log(10)) + 0.5);
						my $sec_gl = 0;
						$sec_gl = int (-10 * (log ($sec_pvalue) / log(10)) + 0.5) if ($sec_pvalue > 0);
						$GQ = $top_gl - $sec_gl;
						$GQ = 99 if ($GQ > 99);
				    }
				    if ($sec_pvalue > 0){
						my $sec_gl = int (-10 * (log ($sec_pvalue) / log(10)) + 0.5);
						my $third_gl = 0;
						$third_gl = int (-10 * (log ($third_pvalue) / log(10)) + 0.5) if ($third_pvalue > 0);
						$GQ2 = $sec_gl - $third_gl;
						$GQ2 = 99 if ($GQ2 > 99);
				    }
$GQ = 999;
$GQ = int (-10 * (log ($prob{0}) / log(10)) + 0.7) if ($prob{0} > 0);
$GQ = 999 if ($GQ > 999);
$GQ2 = $GQ;
				    my ($test_pos, $test_order) = split (/=/, $test_order{$num});
				    ${${$pred{$test_pos}}{$type3}}{$test_order} = "$top_gt=$top_prob=$sec_gt=$GQ=$GQ2";
				}
				close (FILE2);
		    }
		    elsif (scalar keys %test_order == 1){
				my %prob;
				open (FILE2, $out_file);
				while (my $line2 = <FILE2>){
				    chomp $line2;
				    my ($n, $prob) = split (/\s+/, $line2);
				    $prob{$n} = $prob;
				}
				close (FILE2);
				my $top_gt = 0;
				my $sec_gt = 0;
				my $third_gt = 0;
				my $top_prob = 0;
				my $sec_prob = 0;
				my $third_prob = 0;
				foreach my $gt (sort {$prob{$b} <=> $prob{$a}} keys %prob){
				    if ($top_prob == 0){
						$top_prob = $prob{$gt};
						$top_gt = $gt;
						next;
				    }
				    if ($sec_prob == 0){
						$sec_prob = $prob{$gt};
						$sec_gt = $gt;
				    }
				    if ($third_prob == 0){
						$third_prob = $prob{$gt};
						$third_gt = $gt;
				    }
				}
				my $top_pvalue = 1 - $top_prob;
				my $sec_pvalue = 1 - $sec_prob;
				my $third_pvalue = 1 - $third_prob;
				my $GQ = 99;
				my $GQ2 = 99;
				if ($top_pvalue > 0){
				    my $top_gl = int (-10 * (log ($top_pvalue) / log(10)) + 0.5);
				    my $sec_gl = 0;
				    $sec_gl = int (-10 * (log ($sec_pvalue) / log(10)) + 0.5) if ($sec_pvalue > 0);
				    $GQ = $top_gl - $sec_gl;
				    $GQ = 99 if ($GQ > 99);
				}
				if ($sec_pvalue > 0){
				    my $sec_gl = int (-10 * (log ($sec_pvalue) / log(10)) + 0.5);
				    my $third_gl = 0;
				    $third_gl = int (-10 * (log ($third_pvalue) / log(10)) + 0.5) if ($third_pvalue > 0);
				    $GQ2 = $sec_gl - $third_gl;
				    $GQ2 = 99 if ($GQ2 > 99);
				}
$GQ = 999;
$GQ = int (-10 * (log ($prob{0}) / log(10)) + 0.7) if ($prob{0} > 0);
$GQ = 999 if ($GQ > 999);
$GQ2 = $GQ;
				my ($test_pos, $test_order) = split (/=/, $test_order{1});
				${${$pred{$test_pos}}{$type3}}{$test_order} = "$top_gt=$top_prob=$sec_gt=$GQ=$GQ2";
		    }
		}
    }
    close (OUT1);
    close (OUT2);
    close (OUT3);
    close (OUT4);
    close (OUT5);
    close (OUT6);
    close (OUT7);
    close (OUT8);
    close (OUT1R);
    close (OUT2R);
    close (OUT3R);
    close (OUT4R);
    close (OUT5R);
    close (OUT6R);
    close (OUT7R);
    close (OUT8R);
    
    foreach my $pos (sort {$a <=> $b} keys %new_line){
		foreach my $type (keys %{$new_line{$pos}}){
		    my @line = split (/\t/, ${$new_line{$pos}}{$type});
		    my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		    my $sn = 0;
		    $sn = $1 if ($line[7] =~ /SC=(\d+)/);
		    my $rep_flag = 0;
		    $rep_flag = 1 if (exists ${${$sv_repeat{$chr2}}{$pos}}{$type});
		    my $GQ_min = 0;
		    my $MC_max = 1;
		    my $GQDEV_max = 3;
		    my $GQDEV_min = -3;
		    $MC_max = $MSC_filt{$type} if ($rep_flag == 0) and (exists $MSC_filt{$type});
		    $MC_max = $MSCr_filt{$type} if ($rep_flag == 1) and (exists $MSCr_filt{$type});
		    $GQ_min = $GQ_filt{$type} if ($rep_flag == 0) and (exists $GQ_filt{$type});
		    $GQ_min = $GQr_filt{$type} if ($rep_flag == 1) and (exists $GQr_filt{$type});
		    $GQDEV_max = $GQDEV_filt{$type} if ($rep_flag == 0) and (exists $GQDEV_filt{$type});
		    $GQDEV_max = $GQDEVr_filt{$type} if ($rep_flag == 1) and (exists $GQDEVr_filt{$type});
		    $GQDEV_min = -5 if ($type =~ /INS|INV/) or ($rep_flag == 1);
		    my $ave_dr = 1;
		    my $ave_sr = 0;
		    $ave_dr = ${$DR_ave{$type}}{$pos} if (exists ${$DR_ave{$type}}{$pos});
		    $ave_sr = ${$SR_ave{$type}}{$pos} if (exists ${$SR_ave{$type}}{$pos});
	#print STDERR "$pos\t$type\tDR_av: $ave_dr\tSR_ave: $ave_sr\n";
		    my $mc = 0;
		    my $mc_gq_sum = 0;
		    my @gqM0;
		    my $count = 0;
		    foreach (@line){
				$count ++;
				next if ($count <= 9);
				my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $MC) = split (/:/, $_);
				if (exists ${${$pred{$pos}}{$type}}{$count}){
				    my ($gt2, $prob, $gt3, $qual, $qual2) = split (/=/, ${${$pred{$pos}}{$type}}{$count});
				    $gq = $qual;
				    if ($gt eq '0/0'){
		#print STDERR "$count\t$_\n";
		#print STDERR "${${$pred{$pos}}{$type}}{$count}\n";
						if ($MC == 0){
						    $line[$count - 1] = "0/0:$gq:0:0:$dr:$ds:$sr:0";
						}
						elsif (($gt2 == 1) and ($prob >= $min_prob2) and ($gq > 0)){
						    $gt = '0/1';
						    $line[$count - 1] = "$gt:$gq:$vp:$vl:$dr:$ds:$sr:$MC";
						    $correct_calls ++;
						    $ref_correct_calls ++;
						}
						elsif (($gt2 == 2) and ($gq > 0)){
						    $gt = '1/1';
						    $line[$count - 1] = "$gt:$gq:$vp:$vl:$dr:$ds:$sr:$MC";
						    $correct_calls ++;
						    $ref_correct_calls ++;
						}
						else{
						    $line[$count - 1] = "0/0:$gq:0:0:$dr:$ds:$sr:0";
						}
				    }
				    elsif ($gt eq './.'){
						if (($gt2 == 1) and ($MC == 0)){
						    $gt = '0/1';
						}
						elsif (($gt2 == 1) and ($MC >= 1)){
						    if ($gt3 == 0){
								if ($prob >= $min_prob){
								    $gt = '0/1';
								}
								elsif (($type =~ /INS|INV/) and ($sr >= $ave_sr * 0.9)){
								    $gt = '0/1';
								}
								elsif (($type eq 'DEL') and ($sr >= $ave_sr * 0.9) and ($dr <= $ave_dr * 1.1)){
								    $gt = '0/1';
								}
								elsif (($type eq 'DUP') and ($sr >= $ave_sr * 0.9) and ($dr >= $ave_dr * 0.9)){
								    $gt = '0/1';
								}
								else{
								    $gt = '0/0';
								    $correct_calls --;
								}
						    }
						    else{
								$gt = '0/1';
						    }
						}
						elsif ($gt2 == 2){
						    $gt = '1/1';
						}
                        elsif (($gt2 == 0) and ($MC == 0)){
                            $gt = '0/1';
                            $gq = $qual2;			
                        }
						elsif (($gt2 == 0) and ($prob >= $min_prob3) and ($MC >= 1)){
						     $gt = '0/0';
						     $correct_calls --;
						}
						elsif (($gt2 == 0) and ($type =~ /INS|INV/) and ($sr >= $ave_sr)){
						    $gt = '0/1';
						}
						elsif (($gt2 == 0) and ($type eq 'DEL') and ($sr >= $ave_sr * 0.9) and ($dr <= $ave_dr * 1.1)){
						    $gt = '0/1';
						}
						elsif (($gt2 == 0) and ($type eq 'DUP') and ($sr >= $ave_sr * 0.9) and ($dr >= $ave_dr * 0.9)){
						    $gt = '0/1';
						}
						elsif ($gt2 == 0) {
						    $gt = '0/0';
						    $correct_calls --;
						}
						
						if (($gt ne '0/0') and ($gq == 0) and ($type eq 'DUP')){
						    $gt = './.';
						}
						$line[$count - 1] = "$gt:$gq:$vp:$vl:$dr:$ds:$sr:$MC" if ($gt ne '0/0');
						$line[$count - 1] = "0/0:$gq:0:0:$dr:$ds:$sr:0" if ($gt eq '0/0');
				    }
				    else{
						if ($MC == 1){
						    if (($gt2 == 0) and ($prob >= $min_prob3)){
								$gt = '0/0';
								$correct_calls -- ;
						    }
						    elsif (($gt2 == 0) and ($type =~ /INS|INV/) and ($sr >= $ave_sr)){
								$gt = '0/1';
						    }
						    elsif (($gt2 == 0) and ($type eq 'DEL') and ($sr >= $ave_sr * 0.9) and ($dr <= $ave_dr * 1.1)){
								$gt = '0/1';
						    }
						    elsif (($gt2 == 0) and ($type eq 'DUP') and ($sr >= $ave_sr * 0.9) and ($dr >= $ave_dr * 0.9)){
								$gt = '0/1';
						    }
						    elsif ($gt2 == 0){
								$gt = '0/0';
								$correct_calls -- ;
						    }
						    elsif (($gt eq '1/1') and ($gt2 == 1)){
								$gt = '0/1';
						    }
						}
						if (($gt ne '0/0') and ($gq == 0) and ($type eq 'DUP')){
						    $gt = './.';
						}
						$line[$count - 1] = "$gt:$gq:$vp:$vl:$dr:$ds:$sr:$MC" if ($gt ne '0/0');
						$line[$count - 1] = "0/0:$gq:0:0:$dr:$ds:$sr:0" if ($gt eq '0/0');
				    }
		#print STDERR "GT: $gt\n";
				    if (($MC >= 1) and ($gt ne '0/0')){
						$mc ++;
						$mc_gq_sum += $gq;
				    }
				    elsif (($MC == 0) and ($gt ne '0/0')){
						push @gqM0, $gq;
				    }
				}
		    }
		    my $MCR_filt_flag = 0;
		    my $msc_rate = 0;
		    $msc_rate = int ($mc / $total_samples * 1000 + 0.5) / 1000;
		    my $mcgq_ave = 100;
		    $mcgq_ave = int ($mc_gq_sum / $mc + 0.5) if ($mc > 0);
		    if ((@gqM0 > 0) and ($mc_gq_sum > 0)){
				my $sumgq = 0;
				map{$sumgq += $_} @gqM0;
				my $gqave = int ($sumgq / @gqM0 * 10 + 0.5) / 10;
				my $sum_diff = 0;
				map{$sum_diff += ($gqave - $_) ** 2} @gqM0;
				my $sd = int (($sum_diff / @gqM0) ** 0.5 * 10 + 0.5) / 10;
				if ($mcgq_ave > $gqave + $sd * $GQDEV_max){
				    $MCR_filt_flag = 1;
				}
				elsif ($mcgq_ave < $gqave - $sd * $GQDEV_min){
				    $MCR_filt_flag = 1;
				}
		    }
		    $MCR_filt_flag = 1 if ($msc_rate > $MC_max) and ($total_samples >= 20);
		    $MCR_filt_flag = 1 if ($mcgq_ave < $GQ_min);
		    my $ac = 0;
		    my $sc = 0;
		    my @DR;
		    my @DS;
		    my @SR;
		    $count = 0;
		    foreach (@line){
				$count ++;
				next if ($count <= 9);
				my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $MC) = split (/:/, $_);
				next if ($gt eq '0/0');
				push @DR, $dr;
				push @DS, $ds;
				push @SR, $sr;
				if ($gt ne './.'){
				    $ac ++;
				    $ac ++ if ($gt eq '1/1');
				}
				$total_calls ++;
				$sc ++;
		    }
		    if ($sc == 0){
				print STDERR "Sample count: 0: $chr2\t$pos\t$type\n";
				next;
		    }
		    my $sum_dr = 0;
		    my $ave_dr2 = 1;
		    my $sum_ds = 0;
		    my $ave_ds = 0;
		    my $sum_sr = 0;
		    my $ave_sr2 = 0;
		    if (@SR > 0){
				map{$sum_sr += $_} @SR;
				$ave_sr2 = int ($sum_sr / @SR * 100 + 0.5) / 100;
		    }
		    if ($type =~ /DEL|DUP/){
				if (@DR > 0){
				    map{$sum_dr += $_} @DR;
				    $ave_dr2 = int ($sum_dr / @DR * 100 + 0.5) / 100;
				}
				if (@DS > 0){
				    map{$sum_ds += $_} @DS;
				    $ave_ds = int ($sum_ds / @DS * 100 + 0.5) / 100;
				}
				my $filt_flag = 0;
				if ($DP_filter == 1){
				    if (($type eq 'DEL') and ($len >= $min_del_len) and ($len < $min_del_len2)){
						if (($ave_dr > $max_del_dprate) and ($ave_ds > $max_del_dsr)){
						    $del_filt ++;
						    $filt_flag = 1;
						}
				    }
				    elsif (($type eq 'DEL') and ($len >= $min_del_len2)){
						if (($ave_dr > $max_del_dprate2) and ($ave_ds > $max_del_dsr)){
						    $del_filt ++;
						    $filt_flag = 1;
						}
				    }
				    if (($type eq 'DEL') and ($len >= $min_len) and ($ave_dr > $max_del_dprate3)){
						$del_filt ++;
						$filt_flag = 1;
				    }
				    if (($type eq 'DUP') and ($len >= $min_dup_len)){
						if (($ave_dr < $min_dup_dprate) and ($ave_ds > $max_dup_dsr)){
						    $dup_filt ++;
						    $filt_flag = 1;
						}
				    }
				    if (($type eq 'DUP') and ($len >= $min_len) and ($ave_dr < $min_dup_dprate3)){
						$dup_filt ++;
						$filt_flag = 1;
				    }
				    if ($filt_flag == 1){
						${${$rm_sv{$type}}{$chr}}{$pos} = 1;
						next;
				    }
				}
		    }
		    my $af = int ($ac / $total_samples * 0.5 * 10000 + 0.5) / 10000;
		    if (($af >= 0.99) and ($total_samples >= 20)){
				print STDERR "High-AF-SV\t$chr2\t$pos\t$type\t$len\t$sc\t$af\n";
		#		next;
		    }
		    if (($af >= 0.2) and ($type eq 'DUP') and ($len <= 300)){
				$count = 0;
				my $ac2 = 0;
				my $flag = 0;
				foreach (@line){
				    $count ++;
				    next if ($count <= 9);
				    my @info = split (/:/, $_);
				    next if ($info[0] eq '0/0');
				    my $gt = $info[0];
				    my $dr = $info[4];
				    my $sr = $info[6];
				    if (($dr >= 1.2) and ($sr >= 0.35) and (($gt eq '0/1') or ($gt eq './.'))){
						$gt = '1/1';
						$flag = 1;
						$info[0] = $gt;
						$line[$count - 1] = join (':', @info);
				    }
				    $ac2 ++ if ($gt eq '0/1');
				    $ac2 += 2 if ($gt eq '1/1');
				}
				if ($flag == 1){
				    $af = int ($ac2 / $total_samples * 0.5 * 10000 + 0.5) / 10000;
				    $ac = $ac2;
				}
		    }
		    
		    $line[7] =~ s/DPR=[\d\.]+/DPR=$ave_dr2/ if ($line[7] =~ /DPR=/);
		    $line[7] =~ s/DPRATE=[\d\.]+/DPR=$ave_dr2/ if ($line[7] =~ /DPRATE=/);
		    $line[7] =~ s/SR=[\d\.]+/SR=$ave_sr2/ if ($line[7] =~ /SR=/);
		    $line[7] .= ";DPR=$ave_dr2" if ($line[7] !~ /DPR=/);
		    $line[7] .= ";DPS=$ave_ds" if ($line[7] !~ /DPS=/);
		    $line[7] .= ";SR=$ave_sr2" if ($line[7] !~ /SR=/);
	        if ($line[7] =~ /AF=/){
	            $line[7] =~ s/SC=\d+/SC=$sc/;
	            $line[7] =~ s/AF=[\d\.\-Ee]+/AF=$af/;
	            $line[7] =~ s/AC=\d+/AC=$ac/;
	        }
	        else{
	            $line[7] =~ s/SC=\d+/SC=$sc;AC=$ac;AF=$af/;
	        }
		    my $new_line = join ("\t", @line);
		    ${${$vcf2{$type}}{$chr}}{$pos} = $new_line;
		    ${$INS{$chr}}{$pos} = $new_line if ($type eq 'INS');
		    ${$DUP{$chr}}{$pos} = $new_line if ($type eq 'DUP');
        }
    }   
    $pre_chr = $chr2;
}
%vcf = ();

print STDERR "4th step completed: \n";

my $min_merge_len = 1000;
my $min_distance_rate = 0.01;
my $ins_dup_merge = 0;

if ($total_samples >= 5){      # merge overlapping SVs in LD (>= 0.7 r2) 
    foreach my $type (keys %vcf2){
		next if ($type eq 'INS');
        foreach my $chr (sort keys %{$vcf2{$type}}){
            my $pre_pos = 0;
            my $pre_end = 0;
            my $pre_sn = 0;
            my $pre_line = '';
            foreach my $pos (sort {$a <=> $b} keys %{${$vcf2{$type}}{$chr}}){
                my $line = ${${$vcf2{$type}}{$chr}}{$pos};
                my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
                my $end = $pos + $len - 1;
                my $sn = $1 if ($line =~ /SC=(\d+)/);
                my $pre_len = $pre_end - $pre_pos + 1;
                my $distance = $pos - $pre_end + 1;
				my $overlap_flag = 0;
				$overlap_flag = 1 if ($pre_len >= $min_merge_len) and ($len >= $min_merge_len) and ($distance <= $len * $min_distance_rate) and ($distance <= $pre_len * $min_distance_rate) and ($sn /$pre_sn >= 0.8) and ($sn / $pre_sn <= 1.2);
				if ($pos < $pre_end){
				    my $overlap = $pre_end - $pos + 1;
				    if (($overlap >= $len * $min_overlap_rate2) or ($overlap >= $pre_len * $min_overlap_rate2)){
						$overlap_flag = 1 if ($len / $pre_len < 3) and ($len / $pre_len > 0.33);
				    }
				}
                if ($overlap_flag == 1){
				    my $temp_vcf = "$temp_dir/$out_prefix.chr$chr.$pos.vcf";
				    open (OUT1, "> $temp_vcf");
				    map{print OUT1 "$_\n"} @header;
				    print OUT1 "$pre_line\n";
				    print OUT1 "$line\n";
				    close (OUT);
				    my $ld_outprefix = "$temp_dir/$out_prefix.chr$chr.$pos";
				    system ("vcftools --vcf $temp_vcf --geno-r2 --ld-window-bp 1000000 --min-r2 0 --out $ld_outprefix 2>>$out_prefix.vcftool.log");
				    my $ld_out = "$ld_outprefix.geno.ld";
				    my $hit_ld = 0;
				    if (-f $ld_out){
						my $hit_sv = 0;
						open (FILE2, $ld_out);
						while (my $line2 = <FILE2>){
						    chomp $line2;
						    next if ($line2 =~ /^CHR|^$/);
						    my @line2 = split (/\t/, $line2);
						    my $r2 = $line2[4];
						    $hit_ld = $r2;
						}
						close (FILE2);   
				    }
				    if ($hit_ld >= $min_ld){
						my $new_len = $end - $pre_pos + 1;
						my %info;
						my @pre_line = split (/\t/, $pre_line);
						my @line = split (/\t/, $line);
						my $count1 = 0;
						my $sc = 0;
						my $ac = 0;
						my $sum_dpr = 0;
						my $sum_dps = 0;
						my $sum_sr = 0;
						foreach my $preinf (@pre_line){
						    $count1 ++;
						    next if ($count1 <= 9);
						    my $inf = $line[$count1 - 1];
						    if (($inf !~ /^0\/0/) and ($preinf !~ /^0\/0/)){
								my @preinf = split (/:/, $preinf);
								my @inf = split (/:/, $inf);
								$preinf[3] = $new_len;
								$preinf[0] = $inf[0] if ($preinf[0] eq './.');
								$preinf[6] = 1;
								my $dpr2 = int (($preinf[4] + $inf[4]) * 0.5 * 100 + 0.5) / 100;
								my $dps2 = int (($preinf[5] + $inf[5]) * 0.5 * 100 + 0.5) / 100;
								my $sr2 = int (($preinf[6] + $inf[6]) * 0.5 * 100 + 0.5) / 100;
								$preinf[4] = $dpr2;
								$preinf[5] = $dps2;
								$preinf[6] = $sr2;
								$sum_dpr += $dpr2;
								$sum_dps += $dps2;
								$sum_sr += $sr2;
								$ac ++ if ($preinf[0] eq '0/1');
								$ac += 2 if ($preinf[0] eq '1/1');
								my $new_inf = join (':', @preinf);
								$pre_line[$count1 - 1] = $new_inf;
						    }
						    elsif (($inf !~ /^0\/0/) and ($preinf =~ /^0\/0/)){
								my @inf = split (/:/, $inf);
								$inf[2] = $pre_pos;
								$inf[3] = $new_len;
								my $new_inf = join (':', @inf);
								$sum_dpr += $inf[4];
								$sum_dps += $inf[5];
								$sum_sr += $inf[6];
								$ac ++ if ($inf[0] eq '0/1');
								$ac += 2 if ($inf[0] eq '1/1');
								$pre_line[$count1 - 1] = $new_inf;
						    }
						    elsif (($inf =~ /^0\/0/) and ($preinf !~ /^0\/0/)){
								my @preinf = split (/:/, $preinf);
								$preinf[3] = $new_len;
								$ac ++ if ($preinf[0] eq '0/1');
								$ac += 2 if ($preinf[0] eq '1/1');
								my $new_inf = join (':', @preinf);
								$pre_line[$count1 - 1] = $new_inf;
								$sum_dpr += $preinf[4];
								$sum_dps += $preinf[5];
								$sum_sr += $preinf[6];
						    }
						    $sc ++ if ($inf !~ /^0\/0/) or ($preinf !~ /^0\/0/);
						}
						my $ave_dpr = int ($sum_dpr / $sc * 100 + 0.5) / 100;
						my $ave_dps = int ($sum_dps / $sc * 100 + 0.5) / 100;
						my $ave_sr = int ($sum_sr / $sc * 100 + 0.5) / 100;
						my $af = int ($ac / $total_samples * 0.5 * 10000 + 0.5) / 10000;
						$new_len = 0 - $new_len if ($type eq 'DEL');
						$pre_line[7] =~ s/SVLEN=-*\d+/SVLEN=$new_len/;
						$pre_line[7] =~ s/SC=\d+/SC=$sc/ if ($pre_line[7] =~ /SC=/);
						$pre_line[7] =~ s/AC=\d+/AC=$ac/ if ($pre_line[7] =~ /AC=/);
						$pre_line[7] =~ s/AF=[\d\.]+/AF=$af/ if ($pre_line[7] =~ /AF=/);
						$pre_line[7] =~ s/END=\d+/END=$end/ if ($pre_line[7] =~ /END=/);
						$pre_line[7] =~ s/DPR=[\d\.]+/DPR=$ave_dpr/ if ($pre_line[7] =~ /DPR=/);
						$pre_line[7] =~ s/DPS=[\d\.]+/DPS=$ave_dps/ if ($pre_line[7] =~ /DPS=/);
						$pre_line[7] =~ s/SR=[\d\.]+/SR=$ave_sr/ if ($pre_line[7] =~ /SR=/);
						$pre_line = join ("\t", @pre_line);
						delete ${${$vcf2{$type}}{$chr}}{$pos};
						${${$vcf2{$type}}{$chr}}{$pre_pos} = $pre_line;
						$pre_end = $end;
print STDERR "Merged $type: $chr $pre_pos-$pos ($distance)\t$pre_len + $len ($new_len)\tld: $hit_ld\tSN: $sn\n";
						next;
				    }
                }
                $pre_pos = $pos;
                $pre_end = $end;
                $pre_sn = $sn;
                $pre_line = $line;
            }
        }
    }
}
    
foreach my $chr (sort keys %INS){		# mergeing overlapped DUP-INS
    next if (!exists $DUP{$chr});
    my $pre_pos = 0;
    foreach my $ipos (sort {$a <=> $b} keys %{$INS{$chr}}){
		my @overlap;
		foreach my $dpos (sort {$a <=> $b} keys %{$DUP{$chr}}){
		    last if ($dpos > $ipos + $ins_sd);
		    my $ins_line = ${$INS{$chr}}{$ipos};
		    my $dup_line = ${$DUP{$chr}}{$dpos};
		    my $duplen = $1 if ($dup_line =~ /SVLEN=(\d+)/);
		    my $dend = $dpos + $duplen - 1;
		    if ((abs ($ipos - $dpos) <= $ins_sd2) or (abs ($ipos - $dend) <= $ins_sd2)){
				my $hit_flag = 0;
				my $hit_ld = 0;
				if ((abs ($ipos - $dpos) <= 20) and ($duplen < 100)){
				    $hit_flag = 1;
				}
				elsif ((abs ($ipos - $dend) <= 20) and ($duplen < 100)){
				    $hit_flag = 1;
				}
				elsif ($total_samples == 1){
				    $hit_flag = 1;
				}
				else{
				    my $temp_vcf = "$temp_dir/$out_prefix.chr$chr.$ipos.vcf";
				    open (OUT1, "> $temp_vcf");
				    map{print OUT1 "$_\n"} @header;
				    print OUT1 "$ins_line\n";
				    print OUT1 "$dup_line\n";
				    close (OUT);
				    my $ld_outprefix = "$temp_dir/$out_prefix.chr$chr.$ipos";
				    system ("vcftools --vcf $temp_vcf --geno-r2 --ld-window-bp 1000000 --min-r2 0 --out $ld_outprefix 2>>$out_prefix.vcftool.log");
				    my $ld_out = "$ld_outprefix.geno.ld";
				    if (-f $ld_out){
						my $hit_sv = 0;
						open (FILE2, $ld_out);
						while (my $line2 = <FILE2>){
						    chomp $line2;
						    next if ($line2 =~ /^CHR|^$/);
						    my @line2 = split (/\t/, $line2);
						    my $r2 = $line2[4];
						    $hit_ld = $r2;
						}
						close (FILE2);
						$hit_flag = 1 if ($hit_ld >= $min_ld2);
				    }
				}
				if ($hit_flag == 1){
				    my @ins_line = split (/\t/, $ins_line);
				    my @dup_line = split (/\t/, $dup_line);
				    my $count = 0;
				    my %ins_gt;
				    my %dup_gt;
				    foreach (@ins_line){
						$count ++;
						next if ($count <= 9);
						next if ($_ =~ /^0\/0/);
						my $gt = $1 if ($_ =~ /^(.+?):/);
						$ins_gt{$count} = $_;
				    }
				    $count = 0;
				    foreach (@dup_line){
						$count ++;
						next if ($count <= 9);
						next if ($_ =~ /^0\/0/);
						my $gt = $1 if ($_ =~ /^(.+?):/);
						$dup_gt{$count} = $_;
				    }
				    my $sn = 0;
				    my $ac = 0;
				    my $count1 = 0;
				    if ($duplen >= 100){
						foreach (@dup_line){
						    $count1 ++;
						    next if ($count1 <= 9);
						    my @info = split (/:/, $_);
						    my $gt = $info[0];
						    if (exists $ins_gt{$count1}){
								if ($gt eq '0/0'){
								    my @ins_info = split (/:/, $ins_gt{$count1});
								    $ins_info[3] = $duplen;
								    my $new_info = join (':', @ins_info);
								    $dup_line[$count1 - 1] = $new_info;
								    $sn ++;
								    $ac ++ if ($ins_info[0] eq '0/1');
								    $ac += 2 if ($ins_info[0] eq '1/1');
								}
								else{
								    my $ins_gt = $1 if ($ins_gt{$count1} =~ /^(.+?):/);
								    $info[0] = $ins_gt;
								    my $new_info = join (':', @info);
								    $dup_line[$count1 - 1] = $new_info;
								    $ac ++ if ($ins_gt eq '0/1');
								    $ac += 2 if ($ins_gt eq '1/1');
								}
						    }
						    else{
								$ac ++ if ($gt eq '0/1');
								$ac += 2 if ($gt eq '1/1');
						    }
						    $sn ++ if ($gt ne '0/0');
						}
						delete ${${$vcf2{'INS'}}{$chr}}{$ipos};
						my $af = int ($ac / $total_samples * 0.5 * 10000 + 0.5) / 10000;
						$dup_line[7] =~ s/SC=\d+/SC=$sn/ if ($dup_line[7] =~ /SC=/);
						$dup_line[7] =~ s/AC=\d+/AC=$ac/ if ($dup_line[7] =~ /AC=/);
						$dup_line[7] =~ s/AF=[\d\.]+/AF=$af/ if ($dup_line[7] =~ /AF=/);
						my $new_dupline = join ("\t", @dup_line);
						${${$vcf2{'DUP'}}{$chr}}{$dpos} = $new_dupline;
						$ins_dup_merge ++;
#print STDERR "Merged INS to DUP: $chr $ipos-$dpos\t$duplen\tld: $hit_ld\tSC: $sn\n";
				    }
				    else{
						foreach (@ins_line){
						    $count1 ++;
						    next if ($count1 <= 9);
						    my @info = split (/:/, $_);
						    my $gt = $info[0];
						    if ((exists $dup_gt{$count1}) and ($gt eq '0/0')){
								$ins_line[$count1 - 1] = $dup_gt{$count1};
								$sn ++;
								$ac ++ if ($dup_gt{$count1} =~ /^0\/1/);
								$ac += 2 if ($dup_gt{$count1} =~ /^1\/1/);
						    }
						    else{
								$ac ++ if ($gt eq '0/1');
								$ac += 2 if ($gt eq '1/1');
						    }
						    $sn ++ if ($gt ne '0/0');
						}
						delete ${${$vcf2{'DUP'}}{$chr}}{$dpos};
						my $af = int ($ac / $total_samples * 0.5 * 10000 + 0.5) / 10000;
						$ins_line[7] =~ s/SC=\d+/SC=$sn/ if ($ins_line[7] =~ /SC=/);
						$ins_line[7] =~ s/AC=\d+/AC=$ac/ if ($ins_line[7] =~ /AC=/);
						$ins_line[7] =~ s/AF=[\d\.]+/AF=$af/ if ($ins_line[7] =~ /AF=/);
						my $new_insline = join ("\t", @ins_line);
						${${$vcf2{'INS'}}{$chr}}{$ipos} = $new_insline;
						$ins_dup_merge ++;
#print STDERR "Merged DUP to INS: $chr $dpos-$ipos\t$duplen\tld: $hit_ld\tSC: $sn\n";
				    }
				}
		    }
		} 
    }
}
%INS = ();
%DUP = ();

print STDERR "Overlapped INS/DUP in LD with >= $min_ld2 r2 or <= 20 bp BP difference for < 100 bp DUPs: $ins_dup_merge\n";

my %delete;

foreach my $type (sort keys %vcf2){		# delete either SVs of overlapping sample SVs
    foreach my $chr (keys %{$vcf2{$type}}){
        my %pre_info;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf2{$type}}{$chr}}){
            my $line = ${${$vcf2{$type}}{$chr}}{$pos};
            my @line = split (/\t/, $line);
            my $end = $1 if ($line[7] =~ /END=(\d+)/);
            my $count = 0;
            foreach (@line){
                $count ++;
                next if ($count <= 9);
                next if ($_ =~ /^0\/0/);
                my ($gt, $gq, $vpos, $vlen) = split (/:/, $_);
                my $vend = $vpos + $vlen - 1;
                $vend = $vpos if ($type eq 'INS');
                if (exists $pre_info{$count}){
                    my ($pre_pos, $pre_end, $pre_vpos, $pre_vend) = split (/=/, $pre_info{$count});
                    if ($type eq 'INS'){
                        if (abs ($vpos - $pre_vpos) <= $ins_sd){
                            my $diff = abs ($pos - $vpos);
                            my $pre_diff = abs ($pre_pos - $pre_vpos);
                            if ($diff <= $pre_diff){
                                ${${${$delete{$type}}{$chr}}{$pre_pos}}{$count} = 1;
                            }
                            else{
                                ${${${$delete{$type}}{$chr}}{$pos}}{$count} = 1;
                                next;
                            }
                        }
                    }
                    else{
                        my $overlap = 0;
                        my $flag = 0;
                        my $pre_vlen = $pre_vend - $pre_vpos + 1;
                        if (($vpos <= $pre_vpos) and ($vend >= $pre_vend)){
                            $overlap = $pre_vlen;
                        }
                        elsif (($vpos >= $pre_vpos) and ($vpos <= $pre_vend)){
                            $overlap = $pre_vend - $vpos + 1;
                            $overlap = $vlen if ($vend < $pre_vend);
                        }
                        elsif (($vend >= $pre_vpos) and ($vend <= $pre_vend)){
                            $overlap = $vend - $pre_vpos + 1;
                            $overlap = $vlen if ($vpos > $pre_vpos);
                        }
                        if (($overlap >= $vlen * $min_overlap_rate2) and ($overlap >= $pre_vlen * $min_overlap_rate2)){
                            $flag = 1;
                        }
                        if ($flag == 1){
                            my $diff = abs ($pos - $vpos) + abs ($end - $vend);
                            my $pre_diff = abs ($pos - $pre_vpos) + abs ($end - $pre_vend);
                            if ($diff <= $pre_diff){
                                ${${${$delete{$type}}{$chr}}{$pre_pos}}{$count} = 1;
                            }
                            else{
                                ${${${$delete{$type}}{$chr}}{$pos}}{$count} = 1;
                                next;
                            }
                        }
                    }
                }
                $pre_info{$count} = "$pos=$end=$vpos=$vend";
            }
        }
    }
}

my %sv_len;
foreach my $type (keys %vcf2){
    foreach my $chr (keys %{$vcf2{$type}}){
    	my $chr2 = $chr;
    	$chr2 =~ s/^0*//;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf2{$type}}{$chr}}){
            my $line = ${${$vcf2{$type}}{$chr}}{$pos};
            my @line = split (/\t/, $line);
            my $ori_len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
            my $count = 0;
			my $sn = 0;
            my @len;
            my @sr_het2;
            my @sr_het_mc0;
            my @sr_hom;
            my @dr_het;
            my @drsr_het_del;
            my @dr_hom;
            my $flag = 0;
            $flag = 1 if (exists ${${$delete{$type}}{$chr}}{$pos});
            foreach (@line){
                $count ++;
                next if ($count <= 9);
                next if ($_ =~ /^0\/0/);
                if (exists ${${${$delete{$type}}{$chr}}{$pos}}{$count}){
                    $line[$count - 1] = '0/0:0:0:0:0:0:0:0';
                }
			    else{
					my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
                    
                    if ($type ne 'INV'){
                    	if ($type eq 'INS'){
                    		if ($sr <= 0.1){
                    			$line[$count - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
                    			next;
                    		}
                    	}
                    	push @sr_hom, $sr if ($gt eq '1/1');
                    	push @sr_het2, $sr if ($gt eq '0/1') and ($mc <= 1);
                    	push @sr_het_mc0, $sr if (($gt eq '0/1') or ($gt eq './.')) and ($mc == 0);
                    }
                    if ($type =~ /DEL|DUP/){
                    	if (($type eq 'DEL') and ($ori_len >= 200) and ($dr > 0.9) and ($mc >= 1)){
                    		$line[$count - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
                    		next;
                    	}
                    	elsif (($type eq 'DUP') and ($ori_len >= 100) and ($dr <= 1.1)){
                    		$line[$count - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
                    		next;
                    	}
                    	push @dr_het, $dr if ($gt eq '0/1') and ($mc <= 1);
                    	push @dr_hom, $dr if ($gt eq '1/1') and ($mc <= 1);
                    	if (($type eq 'DEL') and ($mc == 0)){
                    		my $drsr = $sr + (1 - $dr);
	                    	push @drsr_het_del, $drsr;
	                    }
                    }
                    push @len, $vl;
                    $sn ++;
			    }
            }
			next if ($sn == 0);
			my $sumlen = 0;
			map{$sumlen += $_} @len;
			my $len = int ($sumlen / @len + 0.5);
			${${$sv_len{$chr}}{$pos}}{$type} = $len;
			my $sf = int ($sn / $total_samples * 1000 + 0.5) / 1000;
			my $min_sr_mc0 = 100;
			my $med_sr_mc0 = 0;
			my $ave_sr_mc0 = 0;
			my $sr_mc0_sd = 0;
			if ($type ne 'INV'){
				push @sr_het_mc0, @sr_het2 if (@sr_het_mc0 == 0);
				push @sr_het_mc0, @sr_hom if (@sr_het_mc0 == 0);
				my $med_num = int (scalar (@sr_het_mc0) * 0.5) - 1;
				$med_num = 0 if ($med_num < 0);
				@sr_het_mc0 = sort {$a <=> $b} @sr_het_mc0;
				$med_sr_mc0 = $sr_het_mc0[$med_num];
				my @sr_mc0_2;
				my $sum_sr_mc0 = 0;
				foreach (@sr_het_mc0){
					next if ($med_sr_mc0 > 0) and ($_ >= $med_sr_mc0 * 2);
					$min_sr_mc0 = $_ if ($_ < $min_sr_mc0);
					push @sr_mc0_2, $_;
					$sum_sr_mc0 += $_;
				}
				$ave_sr_mc0 = int ($sum_sr_mc0 / @sr_mc0_2 * 100 + 0.5) / 100 if (@sr_mc0_2 > 0);
				my $sr_mc0_diff = 0;
				map{$sr_mc0_diff += ($_ - $ave_sr_mc0) ** 2} @sr_mc0_2;
				$sr_mc0_sd = int (($sr_mc0_diff / @sr_mc0_2) ** 0.5 * 100 + 0.5) / 100 if (@sr_mc0_2 > 0);
			}
			my $min_drsr_mc0 = 0;
			foreach (@drsr_het_del){
				$min_drsr_mc0 = $_ if ($_ < $min_drsr_mc0);
			}
			if (((@sr_het2 > 0) or (@sr_hom > 0)) and ($sf >= $min_SF) and ($sn >= $min_SC) and ($type eq 'INS')){
				my $sum_sr_het = 0;
				my $sum_sr_het_diff = 0;
				my $min_sr_het = 100;
				my $min_sr_hom = 200;
				my $med_sr_hom = 0;
				my $sr_hom_sd = 0;
				my @sr_het_2;
				push @sr_het2, @sr_hom if (@sr_het2 == 0);
				my $med_het_num = int (scalar (@sr_het2) * 0.5) - 1;
				$med_het_num = 0 if ($med_het_num < 0);
				@sr_het2 = sort {$a <=> $b} @sr_het2;
				my $med_sr_het = $sr_het2[$med_het_num];
				foreach (@sr_het2){
					next if ($med_sr_het > 0) and ($_ >= $med_sr_het * 2);
					push @sr_het_2, $_;
					$sum_sr_het += $_;
					$min_sr_het = $_ if ($_ < $min_sr_het);
				}
				if (@sr_hom > 0){
					my $med_hom_num = int (scalar (@sr_hom) * 0.5) - 1;
					$med_hom_num = 0 if ($med_hom_num < 0);
					@sr_hom = sort {$a <=> $b} @sr_hom;
					$med_sr_hom = $sr_hom[$med_hom_num];
					my $sum_sr_hom = 0;
					my $sum_sr_hom_diff = 0;
					my @sr_hom_2;
					foreach (@sr_hom){
						next if ($med_sr_hom > 0) and ($_ >= $med_sr_hom * 2);
						push @sr_hom_2, $_;
						$sum_sr_hom += $_;
						$min_sr_hom = $_ if ($_ < $min_sr_hom);
					}
					if (@sr_hom_2 > 2){
						my $ave_sr_hom = int ($sum_sr_hom / @sr_hom_2 * 100 + 0.5) / 100;
						map{$sum_sr_hom_diff += ($ave_sr_hom - $_) ** 2} @sr_hom_2;
						$sr_hom_sd = int (($sum_sr_hom_diff / @sr_hom_2) ** 0.5 * 100 + 0.5) / 100;
					}
				}
				my $ave_sr_het = 0;
				$ave_sr_het = int ($sum_sr_het / @sr_het_2 * 100 + 0.5) / 100 if (@sr_het_2 > 0);
				map{$sum_sr_het_diff += ($ave_sr_het - $_) ** 2} @sr_het_2;
				my $sr_sd = int (($sum_sr_het_diff / @sr_het_2) ** 0.5 * 100 + 0.5) / 100;
				my $count1 = 0;
#print STDERR "$sf\tHet: $ave_sr_het\tHom: $ave_sr_hom\tSD: $sr_sd\tMin-Het: $min_sr_het\n";
				foreach (@line){
                    $count1 ++;
                    next if ($count1 <= 9);
                    my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);

					if (($gt eq '0/0') and ($sr >= 0.3) and ($mc == 0) and (($sr >= $min_sr_mc0) or (($ave_sr_het > 0) and ($sr >= $ave_sr_het - $sr_sd)))){
                    	$gt = '0/1';
                    	$mc = 3;
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:3";
						$sn ++;
                    }
                    if ($gt ne '0/0'){
#	                    if (($mc >= 2) and ($sr < $min_sr_mc0) and ($min_sr_het < 1) and (($sr < $min_sr_het) or ($sr < $ave_sr_het - $sr_sd * 1.5) or ($sr < $ave_sr_hom))){
						if (($mc >= 1) and ($sr < $min_sr_mc0) and ($min_sr_het < 1) and ($sr < $med_sr_het - $sr_sd * 2.5)){
	                    	$gt = '0/0';
	                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
							$sn --;
	                    }
	                    elsif (($mc >= 2) and ($sr <= 0.2) and ($min_sr_het < 1) and (($sr < $min_sr_het) or ($sr < $med_sr_het - $sr_sd * 2.5))){
	                    	$gt = '0/0';
	                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
							$sn --;
	                    }
	                    elsif (($mc >= 2) and ($gt eq '0/1') and ($med_sr_hom > 0) and ($sr > $min_sr_hom) and ($sr > $med_sr_hom - $sr_hom_sd * 2)){
	                    	$gt = '1/1';
	                    	$line[$count1 - 1] = "1/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
	                    }
	                    elsif (($mc >= 2) and ($gt eq '1/1') and ($med_sr_hom > 0) and (($sr < $ave_sr_het) or ($sr < $med_sr_hom - $sr_hom_sd * 2))){
	                    	$gt = '0/1';
	                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
	                    }
	                }
                }
	        }
=pod
	        if ((@sr_hom >= 5) and ($sf >= 0.8) and ($type eq 'INS')){
	        	my $sum_sr_het = 0;
				my $sum_sr_hom = 0;
				my $sum_sr_hom_diff = 0;
				my $max_sr_het = 0;
				my $min_sr_hom = 1;
				foreach (@sr_het){
					$max_sr_het = $_ if ($_ > $max_sr_het);
					$sum_sr_het += $_;
				}
				foreach (@sr_het){
					$min_sr_hom = $_ if ($_ < $min_sr_hom);
					$sum_sr_hom += $_;
				}
				my $ave_sr_het = 0;
				$ave_sr_het = int ($sum_sr_het / @sr_het * 100 + 0.5) / 100 if (@sr_het > 0);
				my $ave_sr_hom = int ($sum_sr_hom / @sr_hom * 100 + 0.5) / 100;
				map{$sum_sr_hom_diff += ($ave_sr_hom - $_) ** 2} @sr_hom;
				my $sr_sd = int (($sum_sr_hom_diff / @sr_hom) ** 0.5 * 100 + 0.5) / 100;
				my $count1 = 0;
				foreach (@line){
                    $count1 ++;
                    next if ($count1 <= 9);
                    my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
                    next if ($gt eq '0/0') and ($sr == 0);
                    next if ($gt eq '0/1') or ($gt eq '1/1');
                    if (($gt eq '0/0') and ($sr >= 0.3) and ($sr >= $min_sr_hom) and ($sr >= $ave_sr_hom - $sr_sd * 2.5)){
                    	$mc = 2;
                    	$mc = 3 if (exists ${${$LR_overlap{$chr2}}{$pos}}{$type});
						$line[$count1 - 1] = "1/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }
                    elsif (($gt eq '0/0') and ($sr >= 0.2) and ($sr < $min_sr_hom) and ($sr >= $ave_sr_het)){
                    	$mc = 2;
                    	$mc = 3 if (exists ${${$LR_overlap{$chr2}}{$pos}}{$type});
						$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }
                }
			}
=cut
			if (((@dr_het > 0) or (@dr_hom > 0)) and ($sn >= $min_SC_DEL) and ($type =~ /DEL|DUP/)){
				my $sum_dr_het = 0;
				my $sum_dr_hom = 0;
				my $sum_dr_het_diff = 0;
				my $sum_dr_hom_diff = 0;
				my $min_dr_het = 2;
				my $max_dr_het = 0;
				my $min_dr_hom = 4;
				my $max_dr_hom = 0;
				my $ave_dr_het = -1;
				my $ave_dr_hom = -1;

				my $dr_sd = 0;
				my $dr_hom_sd = 0;
				push @dr_het, @dr_hom if (@dr_het == 0);
				foreach (@dr_het){
					$min_dr_het = $_ if ($_ < $min_dr_het);
					$max_dr_het = $_ if ($_ > $max_dr_het);
					$sum_dr_het += $_;
				}
				foreach (@dr_hom){
					$min_dr_hom = $_ if ($_ < $min_dr_hom);
					$max_dr_hom = $_ if ($_ > $max_dr_hom);
					$sum_dr_hom += $_;
				}

				$ave_dr_het = int ($sum_dr_het / @dr_het * 100 + 0.5) / 100 if (@dr_het > 0);
				$ave_dr_hom = int ($sum_dr_hom / @dr_hom * 100 + 0.5) / 100 if (@dr_hom > 0);
				map{$sum_dr_het_diff += ($ave_dr_het - $_) ** 2} @dr_het;
				map{$sum_dr_hom_diff += ($ave_dr_hom - $_) ** 2} @dr_hom;
				$dr_sd = int (($sum_dr_het_diff / @dr_het) ** 0.5 * 100 + 0.5) / 100 if (@dr_het > 0);
				$dr_hom_sd = int (($sum_dr_hom_diff / @dr_hom) ** 0.5 * 100 + 0.5) / 100 if (@dr_hom > 0);
				my $count1 = 0;
				foreach (@line){
                    $count1 ++;
                    next if ($count1 <= 9);
                    my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_); 

                    if (($gt eq '0/0') and ($type eq 'DEL') and ($len < 1000) and (@dr_het > 0) and ($dr < 0.7) and ($dr < $max_dr_het) and ($dr < $ave_dr_het + $dr_sd) and (($sr > $min_sr_mc0) or ($sr > $ave_sr_mc0 - $sr_mc0_sd * 2))){
                    	$gt = '0/1';
                    	$mc = 3;
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }
                    elsif (($gt eq '0/0') and ($type eq 'DEL') and ($len >= 1000) and (@dr_het > 0) and ($dr < 0.7) and ($dr < $max_dr_het) and ($dr < $ave_dr_het + $dr_sd) and ($sr > $min_sr_mc0) and ($sr > $ave_sr_mc0 - $sr_mc0_sd * 2)){
                    	$gt = '0/1';
                    	$mc = 3;
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }
                    elsif (($gt eq '0/0') and ($type eq 'DUP') and ($mc == 0) and (@dr_het > 0) and ($dr > 1.3) and ($dr > $min_dr_het) and ($dr > $ave_dr_het - $dr_sd) and (($sr > $min_sr_mc0) or ($sr > $ave_sr_mc0 - $sr_mc0_sd * 2))){
                    	$gt = '0/1';
                    	$mc = 3;
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }
                    elsif (($gt eq '0/0') and ($type eq 'DUP') and ($mc == 0) and ($len <= 100) and (@sr_het_mc0 > 0) and ($sr >= 0.2) and ($sr >= $med_sr_mc0 - $sr_mc0_sd * 2)){
                    	$gt = '0/1';
                    	$mc = 3;
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
						$sn ++;
                    }

                    next if ($gt eq '0/0');
                    next if ($mc == 0);
                    my $drsr = $sr + (1 - $dr);
    #                if (($type eq 'DEL') and ($len > 1000) and (@dr_het > 0) and ($drsr < $min_drsr_mc0)){
                    if (($type eq 'DEL') and ($len > 100) and (@dr_het > 0) and ($drsr < $min_drsr_mc0)){
                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
						$sn --;
                    }
                    elsif (($type eq 'DEL') and ($len > 1000) and (@dr_het > 0) and (($dr > $max_dr_het) or ($dr > $ave_dr_het + $dr_sd * 2))){
                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
						$sn --;
                    }
                    elsif (($type eq 'DEL') and ($len <= 1000) and (@dr_het > 0) and ($sr < $min_sr_mc0) and ($sr < $ave_sr_mc0 - $sr_mc0_sd * 2)){
                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
						$sn --;
                    }
                    elsif (($type eq 'DUP') and ($len > 1000) and (@dr_het > 0) and (($dr < $min_dr_het) or ($dr < $ave_dr_het - $dr_sd * 2))){
                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
						$sn --;
                    }
                    elsif (($type eq 'DUP') and ($len <= 1000) and (@dr_het > 0) and ($sr < $min_sr_mc0) and ($sr < $ave_sr_mc0 - $sr_mc0_sd * 2)){
                    	$line[$count1 - 1] = "0/0:$gq:$vp:$vl:$dr:$ds:$sr:0";
						$sn --;
                    }
                    elsif (($type eq 'DEL') and ($gt eq '0/1') and (@dr_hom > 0) and ($dr <= $max_dr_hom) and ($dr < $ave_dr_hom + $dr_hom_sd * 2)){
                    	$line[$count1 - 1] = "1/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                    }
                    elsif (($type eq 'DEL') and ($gt eq '1/1') and (@dr_het > 0) and ($dr >= $min_dr_het) and ($dr > $ave_dr_het - $dr_sd * 2)){
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                    }
                    elsif (($type eq 'DUP') and ($gt eq '0/1') and (@dr_hom > 0) and ($dr >= $min_dr_hom) and ($dr > $ave_dr_hom - $dr_hom_sd * 2)){
                    	$line[$count1 - 1] = "1/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                    }
                    elsif (($type eq 'DUP') and ($gt eq '1/1') and (@dr_het > 0) and ($dr <= $max_dr_het) and ($dr < $ave_dr_het + $dr_hom_sd * 2)){
                    	$line[$count1 - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                    }
                }
			}
			my $new_line = join ("\t", @line);
	        ${${$vcf{$chr}}{$pos}}{$type} = $new_line;
        }
    }
}

my %merge_str_sv;
my %str_overlap;

foreach my $chr (sort keys %vcf){
	my $chr2 = $chr;
	$chr2 =~ s/^0*//;
	foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
			next if ($type eq 'INV');
			my $len = $1 if (${${$vcf{$chr}}{$pos}}{$type} =~ /SVLEN=-*(\d+)/);
			next if ($len > $max_STR) and ($type ne 'INS');
			my $end = $pos + $len - 1;
			$end = $pos if ($type eq 'INS');
			my $Mbin_start = int ($pos / $Mbin_size);
		    my $Mbin_end = int ($end / $Mbin_size);
		    my $flag = 0;
		    my $hit_rpos = 0;
		    if ($Mbin_start == $Mbin_end){
				if (exists ${$simple_repeat2{$chr2}}{$Mbin_start}){
				    foreach my $rpos (sort {$a <=> $b} keys %{${$simple_repeat2{$chr2}}{$Mbin_start}}){
						last if ($rpos > $end + 5);
						my $rend = ${${$simple_repeat2{$chr2}}{$Mbin_start}}{$rpos};
						my $rlen = $rend - $rpos + 1;
						next if ($rend < $pos - 5);
						if ($type eq 'INS'){
						    $flag = 1;
						    $hit_rpos = $rpos;
						    last;
						}
						else{
							my $overlap = 0;
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
						    if ($overlap >= $len * $min_repeat_overlap){
								$flag = 1;
								$hit_rpos = $rpos;
								last;
						    }
						}
				    }
				}
		    }
		    else{
				if (exists ${$simple_repeat2{$chr2}}{$Mbin_start}){
				    my %rep = (%{${$simple_repeat2{$chr2}}{$Mbin_start}});
				    for (my $i = $Mbin_start + 1; $i <= $Mbin_end; $i++){
						%rep = (%rep, %{${$simple_repeat2{$chr2}}{$i}}) if (exists ${$simple_repeat2{$chr2}}{$i});
				    }
				    foreach my $rpos (sort {$a <=> $b} keys %rep){
						last if ($rpos > $end + $ins_sd);
						my $rend = $rep{$rpos};
						my $rlen = $rend - $rpos + 1;
						next if ($rend < $pos - $ins_sd);
						if ($type eq 'INS'){
						    $flag = 1;
						    $hit_rpos = $rpos;
						    last;
						}
						else{
							my $overlap = 0;
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
						    if ($overlap >= $len * $min_repeat_overlap){
								$flag = 1;
								$hit_rpos = $rpos;
								last;
						    }
						}
				    }
				}
		    }
		    if ($flag == 1){
				 ${${${$str_overlap{$chr}}{$hit_rpos}}{$type}}{$pos} = $len;
		    }
		}
	}
}

foreach my $chr (keys %str_overlap){
	foreach my $str_pos (keys %{$str_overlap{$chr}}){
        foreach my $type (keys %{${$str_overlap{$chr}}{$str_pos}}){
        	next if (scalar keys %{${${$str_overlap{$chr}}{$str_pos}}{$type}} <= 1);
        	my %pos_item;
        	my $top_pos = 0;
            if ($type eq 'INS'){
            	foreach my $pos (sort {$a <=> $b} keys %{${${$str_overlap{$chr}}{$str_pos}}{$type}}){
            		next if (!exists ${${$vcf{$chr}}{$pos}}{$type});
            		my $ac = $1 if (${${$vcf{$chr}}{$pos}}{$type} =~ /AC=(\d+)/);
					$pos_item{$pos} = $ac;
				}
				if (scalar keys %pos_item > 1){
					foreach my $pos (sort {$pos_item{$b} <=> $pos_item{$a}} keys %pos_item){
						$top_pos = $pos;
						last;
					}
				}
			}
			else{
            	my $pre_pos = 0;
            	my $pre_len = 0;
            	my $pre_ac = 0;
            	foreach my $pos (sort {${${${$str_overlap{$chr}}{$str_pos}}{$type}}{$b} <=> ${${${$str_overlap{$chr}}{$str_pos}}{$type}}{$a}} keys %{${${$str_overlap{$chr}}{$str_pos}}{$type}}){
            		my $len = ${${${$str_overlap{$chr}}{$str_pos}}{$type}}{$pos};
            		next if (!exists ${${$vcf{$chr}}{$pos}}{$type});
            		my $ac = $1 if (${${$vcf{$chr}}{$pos}}{$type} =~ /AC=(\d+)/);
            		if ($pre_pos > 0){
            			if ($pre_len / $len <= 2){
            				$pos_item{$pre_pos} = $pre_ac;
            				$pos_item{$pos} = $ac;
            			}
            		}
            		$pre_pos = $pos;
            		$pre_len = $len;
            		$pre_ac = $ac;
            	}
            	if (scalar keys %pos_item > 1){
					foreach my $pos (sort {$pos_item{$b} <=> $pos_item{$a}} keys %pos_item){
						$top_pos = $pos;
						last;
					}
            	}
            }
            if ($top_pos > 0){
				my %info;
				foreach my $pos (sort {$pos_item{$b} <=> $pos_item{$a}} keys %pos_item){
					next if ($pos == $top_pos);
					next if (!exists ${${$vcf{$chr}}{$pos}}{$type});
					my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
					my $count = 0;
					foreach (@line){
						$count ++;
						next if ($count <= 9);
						next if ($_ =~ /^0\/0/);
						$info{$count} = $_ if (!exists $info{$count});
					}
					delete ${${$vcf{$chr}}{$pos}}{$type};
					delete ${$vcf{$chr}}{$pos} if (scalar keys %{${$vcf{$chr}}{$pos}} == 0);
					$merge_str_sv{$type} ++;
				}
				if (scalar keys %info > 0){
					my @line = split (/\t/, ${${$vcf{$chr}}{$top_pos}}{$type});
					my $count = 0;
					foreach (@line){
						$count ++;
						next if ($count <= 9);
						if (($_ =~ /^0\/0/) and (exists $info{$count})){
							$line[$count - 1] = $info{$count};
						}
					}
					my $new_line = join ("\t", @line);
					${${$vcf{$chr}}{$top_pos}}{$type} = $new_line;
				}
            }
        }
    }
}

print STDERR "#Merged SVs within STRs (<= $max_STR bp)\n";
foreach my $type (sort keys %merge_str_sv){
	print STDERR "$type\t$merge_str_sv{$type}\n";
}

=pod
foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
		    next if ($type ne 'INS');
		    my $line = ${${$vcf{$chr}}{$pos}}{$type};
            my @line = split (/\t/, $line);
            my $count = 0;
            my @sr_het;
            my @sr_hom;
            my $ac = 0;
            my $sn = 0;
            my $flag = 0;
		    foreach (@line){
                $count ++;
                next if ($count <= 9);
                next if ($_ =~ /^0\/0/);
				my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
                $sn ++;
                $ac ++ if ($gt eq '0/1');
                $ac += 2 if ($gt eq '1/1');
                push @sr_het, $sr if ($gt eq '0/1');
                push @sr_hom, $sr if ($gt eq '1/1');
            }
            my $af = int ($ac / $total_samples * 0.5 * 1000 + 0.5) / 1000;
            next if ($af < $min_SF);
            next if (@sr_het == 0) or (@sr_hom == 0);
            my $sum_sr_het = 0;
            my $sum_sr_hom = 0;
            my $sum_het_diff = 0;
            my $sum_hom_diff = 0;
            my $max_het = 0;
            my $min_hom = 2;
            foreach (@sr_het){
            	$sum_sr_het += $_;
            	$max_het = $_ if ($_ > $max_het);
            }
            foreach (@sr_hom){
            	$sum_sr_hom += $_;
            	$min_hom = $_ if ($_ < $min_hom);
            }
            my $ave_het = int ($sum_sr_het / @sr_het * 100 + 0.5) / 100;
            my $ave_hom = int ($sum_sr_hom / @sr_hom * 100 + 0.5) / 100;
            map{$sum_het_diff += ($_ - $ave_het) ** 2} @sr_het;
            map{$sum_hom_diff += ($_ - $ave_hom) ** 2} @sr_hom;
            my $sd_het = int (($sum_het_diff / @sr_het) ** 0.5 * 100 + 0.5) / 100;
            my $sd_hom = int (($sum_hom_diff / @sr_hom) ** 0.5 * 100 + 0.5) / 100;
            $count = 0;
            foreach (@line){
                $count ++;
                next if ($count <= 9);
                next if ($_ =~ /^0\/0/);
				my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
				if (($gt eq '0/1') and ($sr >= $min_hom)){
					$line[$count - 1] = "1/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
					$ac ++;
					$flag = 1;
				}
				elsif (($gt eq '1/1') and ($sr <= $max_het) and ($sr <= $ave_hom - $sd_hom * 2.5)){
					$line[$count - 1] = "0/1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
					$ac --;
					$flag = 1;
				}
			}
			if ($flag == 1){
				$af = int ($ac / $total_samples * 0.5 * 1000 + 0.5) / 1000;
				$line[7] =~ s/AF=[\d\.]+/AF=$af/ if ($line[7] =~ /AF=/);
				my $new_line = join ("\t", @line);
	            ${${$vcf{$chr}}{$pos}}{$type} = $new_line;
			}
		}
    }
}
=cut
my $out_vcf = "$out_dir/$out_prefix.vcf";

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

my $version = 1.7;
$version = $1 if ($Bin =~ /MOPline[_\-\.](v[\d\.]+)/i);

open (OUT, "> $out_vcf");
my $alt_flag = 0;
my $contig_flag = 0;
foreach (@header){
	if ($_ =~ /##fileDate=/){
		print OUT "##fileDate=$time\n";
		next;
	}
	elsif ($_ =~ /source=/){
		print OUT "##source=MOPline-$version\n";
		next;
	}
    elsif (($_ =~ /^##ALT/) and ($alt_flag == 0)){
        print OUT "##INFO=<ID=DPS,Number=1,Type=Float,Description=\"Mean ratio of 50-bp-window-fractions with inconsistent DPR\">\n";
        print OUT "$_\n";
        $alt_flag = 1;
        next;
    }
    elsif (($_ =~ /^##contig/) and ($contig_flag == 0)){
        print OUT "##FORMAT=<ID=MC,Number=1,Type=Integer,Description=\"Site corrected for missing calls\">\n";
        print OUT "$_\n";
        $contig_flag = 1;
        next;
    }
    print OUT "$_\n";
}
foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
		foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
			my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
			my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
			my $count = 0;
			my $ac = 0;
			my $sc = 0;
			my $sum_dpr = 0;
			my $sum_dps = 0;
			my $sum_sr = 0;
			my $ref = 0;
			my $ref_mc3 = 0;
			my @len;
			foreach (@line){
                $count ++;
                next if ($count <= 9);
                my ($gt, $gq, $vp, $vl, $dr, $ds, $sr, $mc) = split (/:/, $_);
                if ($_ =~ /^0\/0/){
                	$ref ++;
                	$ref_mc3 ++ if ($mc == 3);
                    if (($chr eq 'Y') or ($chr eq 'chrY')){
                        my $info = "0:$gq:$vp:$vl:$dr:$ds:$sr:0";
                        $line[$count - 1] = $info;
                    }
                	next;
                }
                $ac ++;
                $ac ++ if ($gt eq '1/1');
                $sc ++;
                push @len, $vl;
                $sum_dpr += $dr;
                $sum_dps += $ds;
                $sum_sr += $sr;
                if ($vp == 0){
                	my $info = "$gt:$gq:$pos:$len:$dr:$ds:$sr:$mc";
                	$line[$count - 1] = $info;
                }
                if ((($chr eq 'Y') or ($chr eq 'chrY')) and ($gender_list ne '')){
                	my $id = $sample_order{$count};
                	$id = $1 if ($id =~ /:(.+)/);
                	if (exists $gender{$id}){
                		my $sex = $gender{$id};
                		if ($sex eq 'F'){
                			my $info = "0:$gq:$vp:$vl:$dr:$ds:$sr:0";
                			$line[$count - 1] = $info;
                		}
                		else{
                			my $info = "1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                			$line[$count - 1] = $info;
                		}
                	}
                }
                elsif ((($chr eq 'X') or ($chr eq 'chrX')) and ($gender_list ne '')){
                	my $id = $sample_order{$count};
                	$id = $1 if ($id =~ /:(.+)/);
                	if (exists $gender{$id}){
                		my $sex = $gender{$id};
                		if ($sex eq 'M'){
                			my $info = "1:$gq:$vp:$vl:$dr:$ds:$sr:$mc";
                			$line[$count - 1] = $info;
                		}
                	}
                }
            }
            next if ($sc == 0);
            if ($type eq 'INS'){
	            my $ref_sf = int ($ref / $total_samples * 1000 + 0.5) / 1000;
	            if (($ref >= $min_ref_sample) and ($ref_sf >= $min_ref_sf)){
	            	my $ref_mc3_rate = int ($ref_mc3 / $ref * 100 + 0.5) / 100;
	        		if ($ref_mc3_rate >= $min_split_filt_rate){
	        			next;   		# >= 90% of ref alleles are derived from biased direction of split reads for INSs
	        		}
	            }
	        }

            my $ave_dpr = int ($sum_dpr / $sc * 100 + 0.5) / 100;
            my $ave_dps = int ($sum_dps / $sc * 100 + 0.5) / 100;
            my $ave_sr = int ($sum_sr / $sc * 100 + 0.5) / 100;
            my $af = int ($ac / $total_samples * 0.5 * 10000 + 0.5) / 10000;
            my $half = int (@len * 0.5);
            @len = sort {$a <=> $b} @len;
            my $med_len = $len[$half];
            $med_len = 0 - $med_len if ($type eq 'DEL');
            $line[7] =~ s/SVLEN=-*\d+/SVLEN=$med_len/ if ($med_len > 0);
            $line[7] =~ s/SC=\d+/SC=$sc/ if ($line[7] =~ /SC=/);
            $line[7] =~ s/AC=\d+/AC=$ac/ if ($line[7] =~ /AC=/);
            $line[7] =~ s/AF=[\d\.]+/AF=$af/ if ($line[7] =~ /AF=/);
            $line[7] =~ s/DPR=[\d\.]+/DPR=$ave_dpr/ if ($line[7] =~ /DPR=/);
            $line[7] =~ s/DPS=[\d\.]+/DPS=$ave_dps/ if ($line[7] =~ /DPS=/);
            $line[7] =~ s/SR=[\d\.]+/SR=$ave_sr/ if ($line[7] =~ /SR=/);
            my $new_line = join ("\t", @line);
		    print OUT $new_line, "\n";
		}
    }
}
close (OUT);

if ($DP_filter == 1){
    my $filt_delrate = 0;
    my $filt_duprate = 0;
    $filt_delrate = int ($del_filt / $del_num * 1000 + 0.5) / 100 if ($del_num > 0);
    $filt_duprate = int ($dup_filt / $dup_num * 1000 + 0.5) / 100 if ($dup_num > 0);
    
    print STDERR "Filtered DELs: $del_filt ($filt_delrate%)\n";
    print STDERR "Filtered DUPs: $dup_filt ($filt_duprate%)\n";
}

system ("rm -r $temp_dir");

my $correct_rate = int ($correct_calls / $Missed_calls * 1000) / 10;
print STDERR "Total Sample Calls:    $total_calls\n";
print STDERR "Total Missing Calls:   $Missed_calls\n";
print STDERR "Total Corrected Calls from reference alleles: $ref_correct_calls\n";
print STDERR "Total Corrected Calls: $correct_calls ($correct_rate%)\n";


sub cons_len{
    my ($ref_len, $len) = @_;
    my $num = scalar @{$ref_len};
    my %len;
    my $len_sd = 3;
    if ($len > 300){
		$len_sd = int ($len * 0.01);
    }
    map{$len{$_} ++} @{$ref_len};
    my $pre_len = 0;
    foreach my $slen (sort {$a <=> $b} keys %len){
		if (($pre_len > 0) and ($slen - $pre_len <= $len_sd)){
		    my $num1 = $len{$pre_len};
		    my $num2 = $len{$slen};
		    if ($num2 > $num1){
				$len{$slen} += $num1;
		    }
		    elsif ($num2 < $num1){
				$len{$pre_len} += $num2;
		    }
		    else{
				$len{$slen} += $num1;
				$len{$pre_len} += $num2;
		    }
		}
		$pre_len = $slen;
    }
    my $high_freq_len = 0;
    my $high_freq_len2 = 0;
    my $high_freq = 0;
    my $high_freq2 = 0;
    my $new_len = 0;
    foreach my $slen (sort {$len{$b} <=> $len{$a}} keys %len){
		next if ($len{$slen} == 1);
		if ($high_freq_len == 0){
		    $high_freq_len = $slen;
		    $high_freq = $len{$slen};
		    next;
		}
		if ($high_freq_len2 == 0){
		    $high_freq_len2 = $slen;
		    $high_freq2 = $len{$slen};
		    last;
		}
    }
    if (($high_freq > 1) and ($high_freq >= $num * 0.3)){
		if ($high_freq == $high_freq2){
		    $new_len = int (($high_freq_len + $high_freq_len2) * 0.5 + 0.5);
		}
		else{
		    $new_len = $high_freq_len;
		}
    }
    else{
		my $half = int ($num * 0.5);
		my $half2 = $half - 1;
		if ($num % 2 == 0){
		    $new_len = int ((${$ref_len}[$half2] + ${$ref_len}[$half]) * 0.5 + 0.5);
		}
		else{
		    $new_len = ${$ref_len}[$half];
		}
    }
    return ($new_len);
}

sub collect_call{
    my ($dnum, $ref_samples) = @_;
    my $select_out = "$temp_dir/select.$dnum.txt";
    my %MC3;
    open (OUT1, "> $select_out");
    foreach my $GID (@{$ref_samples}){
		next if (($total_samples < 200)) and (!exists $MC{$GID});
		my ($Gr, $id) = split (/:/, $GID) if ($GID =~ /:/);
		$id = $GID if ($GID !~ /:/);
		$Gr = '' if ($GID !~ /:/);
		my $Gdir = '.';
		if (-d "$sample_dir/$Gr"){
		    $Gdir = $Gr;
		}
		if (($Gr ne '') and (!-d "$sample_dir/$Gdir/$id")){
		    die "$sample_dir/$Gdir/$id directory is not found:\n";
		}
		if ($total_samples >= 200){
		    my $MC_gid_out = "$temp_dir/$GID.MC.txt";
		    open (FILE, $MC_gid_out) or die "$MC_gid_out is not found: $!\n";
		    while (my $iline = <FILE>){
				chomp $iline;
				my ($itype, $ichr, $ipos, $ilen) = split (/\t/, $iline);
				${${${$MC3{$GID}}{$itype}}{$ichr}}{$ipos} = $ilen;
		    }
		    close (FILE);
		}
        foreach my $tool (@tools){
		    my $tool2 = $tool;
		    $tool2 =~ s/-/./ if ($tool2 =~ /Mobster|MELT/);
		    $tool2 = 'MELT.NUMT' if ($tool2 eq 'MELT.MT');
		    my %epos;
		    my $tool_vcf = "$sample_dir/$Gdir/$id/$tool/$tool2.$id.vcf";
		    $tool_vcf = "$sample_dir/$id/$tool/$tool2.$id.vcf" if (!-f $tool_vcf);
		    $tool_vcf = "$sample_dir/$id/MELT-MEI/MELT.MEI.$id.vcf" if (!-f $tool_vcf) and ($tool eq 'MELT');
		    $tool_vcf = "$sample_dir/$Gdir/$id/MELT-MEI/MELT.MEI.$id.vcf" if (!-f $tool_vcf) and ($tool eq 'MELT');
		    $tool_vcf = "$sample_dir/$id/$tool2.$id.vcf" if (!-f $tool_vcf);
		    $tool_vcf = "$sample_dir/$tool/$tool2.$id.vcf" if (!-f $tool_vcf);
		    $tool_vcf = "$sample_dir/$tool2.$id.vcf" if (!-f $tool_vcf);
print STDERR "$GID\t$tool_vcf\n";		    
		    open (FILE, $tool_vcf) or die "$tool_vcf is not found: $!\n";
		    while (my $line = <FILE>){
				chomp $line;
				next if ($line =~ /^#|^$/);
				my @line = split (/\t/, $line);
				my $chr = $line[0];
				next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
				my $pos = $line[1];
				my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
				my $subtype = $type;
				$type = 'INS' if ($type =~ /ALU|LINE1|SVA|HERVK|VEI|NUMT/);
				next if ($tool_set_list ne '') and (!exists ${$toolset{$tool}}{$type});
				if ($total_samples < 200){
				    next if (!exists ${$MC{$GID}}{$type});
				    next if (!exists ${${$MC{$GID}}{$type}}{$chr});
				}
				else{
				    next if (!exists ${$MC3{$GID}}{$type});
				    next if (!exists ${${$MC3{$GID}}{$type}}{$chr});
				}
				next if ($tool =~ /Delly/i) and ($line[6] ne 'PASS');
				next if ($tool =~ /metaSV/i) and ($line[6] ne 'PASS');
				my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
				my $size = 'A';
				if ($type ne 'INS'){
				    if ($len < 100){
						$size = 'SS';
				    }
				    elsif ($len <= 1000){
						$size = 'S';
				    }
				    elsif ($len <= 100000){
						$size = 'M';
				    }
				    else{
						$size = 'L';
				    }
				}
				my $read = 0;
				$read = $1 if ($line[7] =~ /READS=(\d+)/);
				if ($tool_set_all =~ /FALSE/i){
				    next if (!exists ${${$toolset{$tool}}{$type}}{$size});
				    my $min_read = ${${$toolset{$tool}}{$type}}{$size};
				    next if ($read < $min_read);
				}
				else{
				    if (($len < 100) and ($type ne 'INS')){
				    	next if (!exists ${${$toolset{$tool}}{$type}}{$size});
				    	my $min_read = ${${$toolset{$tool}}{$type}}{$size};
				    	next if ($read < $min_read);
				    }
				    next if ($len < 1000) and ($read < 10) and ($type eq 'INV');
				    next if ($read < 5) and ($type eq 'INV');
				}
				my $end = $pos + $len - 1 if ($type ne 'INS');
				$end = $pos if ($type eq 'INS');
				my $gt = './.';
				$gt = $1 if ($line[7] =~ /GT=([^;]+)/);
				$gt = '0/1' if ($gt eq '1/0');
				$gt = './.' if ($gt eq '0/0');
				if ($total_samples < 200){
				    foreach my $epos (sort {$a <=> $b} keys %{${${$MC{$GID}}{$type}}{$chr}}){
						my $elen = ${${${$MC{$GID}}{$type}}{$chr}}{$epos};
						my $eend = $epos + $elen - 1 if ($type ne 'INS');
						$eend = $epos if ($type eq 'INS');
						next if ($eend + $BPsd < $pos);
						last if ($epos - $BPsd > $end);
						next if (exists $epos{$epos});
						my $flag = 0;
						if ((abs ($pos - $epos) <= $BPsd) and (abs ($end - $eend) <= $BPsd)){
						    $flag = 1;
						}
						if (($flag == 0) and ($type ne 'INS')){
						    if (($epos <= $pos) and ($eend >= $end)){
								if ($len >= $elen * $min_overlap_rate){
								    $flag = 1;
								}
						    }
						    elsif (($epos >= $pos) and ($eend <= $end)){
								if ($elen >= $len * $min_overlap_rate){
								    $flag = 1;
								}
						    }
						    elsif (($epos >= $pos) and ($epos <= $end)){
								if (($end - $epos + 1 >= $len * $min_overlap_rate) and ($end - $epos + 1 >= $elen * $min_overlap_rate)){
								    $flag = 1;
								}
						    }
						    elsif (($eend >= $pos) and ($eend <= $end)){
								if (($eend - $pos + 1 >= $len * $min_overlap_rate) and ($eend - $pos + 1 >= $elen * $min_overlap_rate)){
								    $flag = 1;
								}
						    }
						}
						if ($flag == 1){
						    print OUT1 "$chr\t$epos\t$type\t$id\t$pos=$len=$tool=$gt=$read\n";
						    $epos{$epos} = 1;
						    last;
						}
				    }
				}
				else{
				    foreach my $epos (sort {$a <=> $b} keys %{${${$MC3{$GID}}{$type}}{$chr}}){
						my $elen = ${${${$MC3{$GID}}{$type}}{$chr}}{$epos};
						my $eend = $epos + $elen - 1 if ($type ne 'INS');
						$eend = $epos if ($type eq 'INS');
						next if ($eend + $BPsd < $pos);
						last if ($epos - $BPsd > $end);
						next if (exists $epos{$epos});
						my $flag = 0;
						if ((abs ($pos - $epos) <= $BPsd) and (abs ($end - $eend) <= $BPsd)){
						    $flag = 1;
						}
						if (($flag == 0) and ($type ne 'INS')){
						    if (($epos <= $pos) and ($eend >= $end)){
								if ($len >= $elen * $min_overlap_rate){
								    $flag = 1;
								}
						    }
						    elsif (($epos >= $pos) and ($eend <= $end)){
								if ($elen >= $len * $min_overlap_rate){
								    $flag = 1;
								}
						    }
						    elsif (($epos >= $pos) and ($epos <= $end)){
								if (($end - $epos + 1 >= $len * $min_overlap_rate) and ($end - $epos + 1 >= $elen * $min_overlap_rate)){
								    $flag = 1;
								}
						    }
						    elsif (($eend >= $pos) and ($eend <= $end)){
								if (($eend - $pos + 1 >= $len * $min_overlap_rate) and ($eend - $pos + 1 >= $elen * $min_overlap_rate)){
								    $flag = 1;
								}
						    }
						}
						if ($flag == 1){
						    print OUT1 "$chr\t$epos\t$type\t$id\t$pos=$len=$tool=$gt=$read\n";
						    $epos{$epos} = 1;
						    last;
						}
				    }
				}
		    }
		    close (FILE);
		}
print STDERR "$dnum: $id re-searching completed:\n";
    }
    close (OUT1);
    undef %MC3;
    
    threads->yield();
    sleep 1;
    return ($select_out);
}

sub add_dp_rate{
    my ($dnum, $ref_samples) = @_;
    my $sample_out = "$temp_dir/Sample_list.$dnum.txt";
    open (OUT2, "> $sample_out");
    foreach my $GID (@{$ref_samples}){
        next if (!exists $MC2{$GID});
        my ($Gr, $id) = split (/:/, $GID) if ($GID =~ /:/);
		$id = $GID if ($GID !~ /:/);
		$Gr = '' if ($GID !~ /:/);
        next if (!exists $dp_rate_input{$id});
        my ($pos_file, $cov_file) = split (/=/, $dp_rate_input{$id});
        print OUT2 "$Gr\t$id\t$pos_file\t$cov_file\n";
    }
    close (OUT2);
    
    my $dr_out = "$temp_dir/dprate.$dnum.txt";
    my $command = "add_dp_rate_samples_3.pl $sample_out $dr_out $temp_dir $min_split_diff 2>$temp_dir/calc_dprate.$dnum.log";
    system ($command);
    
    threads->yield();
    sleep 1;
    return ($dr_out);
}
