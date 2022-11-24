#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use File::Basename;

# Joint calling SVs from multiple sample vcf files

my $sample_list = '';

my $out_prefix = 'MOPline';

my $target_chr = 'ALL';

my $min_sv_len = 50;

my $ins_sd = 200;

my $min_overlap_ratio = 0.5;

my $max_inv_size = 200000;

my $pos_ave = 1;

my $non_human = 0;

my $build = 37;

my $merge_dir = "Merge_7tools";

my $output_dir = $merge_dir;

my $data_dir = "$Bin/../Data";

my $gap_bed = '';

my $gender_list = '';

my $help;

GetOptions(
    'sample_list|s=s' => \$sample_list,
    'prefix|p=s' => \$out_prefix,
    'merge_dir|md=s' => \$merge_dir,
    'outdir|od=s' => \$output_dir,
    'target_chr|c=s' => \$target_chr,
    'gender|g=s' => \$gender_list,
    'gap_bed|gap=s' => \$gap_bed,
    'non_human|nh=i' => \$non_human,
    'build|b=s' => \$build,
    'min_size|ms=i' => \$min_sv_len,
    'max_inv|mv=i' => \$max_inv_size,
    'ins_sd|is=i' => \$ins_sd,
    'overlap_rate|or=f' => \$min_overlap_ratio,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  merge_SV_calls_ALLsamples.pl -s <sample list> -md <merged vcf directory> -od <output directory> 
  (-gap <gap bed -nh 1 if samples are non-human sepecies)

  Options:
   --sample_list or -s <STR> sample list file, which indicates sample_directory,sample_name in each line (e.g., SampleA,sampleA). [mandatory]
                            If sample name and sample directory name are identical, only sample name can be indicted.
                            MOPline expects that the working directory contains ${sample_directory}/${merge_dir}/${sample_name}.Merge.ALL.vcf, 
                            where ${merge_dir} is a directory name specified with the -md option (e.g., Merge_7tools).
                            If sample directories are located under a group directory and/or sample name has a group name at the head (e.g., ${Group}.${sample_name}), 
                            the group name # tag (i.e., #GROUP:${Group}) should be added at the lines before a group of sample lines.
							
   --merge_dir or -md <STR> directory name containing merged vcf files from multiple SV callers for every sample [default: Merge_7tools]
   --outdir or -od <STR>    output directory [default: directory name specified with -md]
   --prefix or -p <STR>     prefix name of an output joint called vcf file [default: MOPline]
   --target or -t <STR>     target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --build or -b <STR>      human reference build [GRCh37, GRCh38, T2T-CHM13] (37, 38, or T2T: only effective for human) [default: 37]
   --gap or -gap <STR>      gap bed file, indicating gap regions in reference genome. 
                            This may be specified for non-human species [automatically selected for human: Data/gap.bed or gap.b38,bed]
   --gender or -g <STR>     sample name-gender table file to exclude chrY for female (sample_name and M/F, separated with tab in each line) [optional]
   --min_size or -ms <INT>  minimum size (bp) of SV, except for INS [default: 50]
   --max_inv or -mv <INT>   maximum size (bp) of INV [default: 200000]
   --ins_sd or -is <INT>    maximum distance (bp) between proximal INS breakpoints to be merged [default: 200]
   --overlap_rate or -or <FLOAT>  minimum reciprocal overlap rate for merging overlapped DELs, DUPs, or INVs [default: 0.5]
   --help or -h             output help message
   
=cut

die "sample list file is not specified: \n" if ($sample_list eq '');

my $var_sd = 150;
my $mei_sd = 150;

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

my @depth_callers = ('CNVnator', 'readDepth');
my $depth_callers = join ('=', @depth_callers);

if ($non_human == 0){
	$gap_bed = "$data_dir/gap.bed" if ($build eq '37');
	$gap_bed = "$data_dir/gap.b38.bed" if ($build eq '38');
}

my %gap;

if ($gap_bed ne ''){
	my $pre_gchr = '';
	my $pre_gend = 0;
	my $pre_gpos = 0;
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

my %sample_id;
my %Gid;
my %group;
my %group_id;
my %gender;
my @var_file;
my @sample_id;

if ($gender_list ne ''){
    open (FILE, $gender_list) or die "$gender_list is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#/);
        my ($ID, $sex) = split (/\t/, $line);
        $gender{$ID} = $sex;
    }
    close (FILE);
}

my $group = '';
open (FILE, $sample_list) or die "$sample_list is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#GROUP:(\S+)/){
		$group = $1;
	}
	next if ($line =~ /^#|^$/);
	my @line = split (/\t/, $line);
	my $gid = $line[0];
	$gid = $1 if ($gid =~ /(.+)\.bam/);
	my $sample_dir = '';
	my $ID = '';
	if ($gid !~ /,/){
		$sample_dir = $gid;
		$ID = $gid;
	}
	else{
		(my $gr_sdir, $ID) = split (/,/, $gid);
		$sample_dir = $gr_sdir;
	}
	if (!-d $sample_dir){
		$sample_dir = $1 if ($sample_dir =~ /(.+?)\./);
		$ID = $1 if ($ID =~ /(.+?)\./);
	}
	my $vcf_name = "$ID.Merge.ALL.vcf";
	my $vcf = "$sample_dir/$merge_dir/$vcf_name";
	if ($group ne ''){
		$vcf = "$group/$sample_dir/$merge_dir/$vcf_name";
		$vcf = "$group/$sample_dir/$merge_dir/$group.$vcf_name" if (!-f $vcf);
		$vcf = "$sample_dir/$merge_dir/$vcf_name";
		$vcf = "$sample_dir/$merge_dir/$group.$vcf_name";
	}
	if (!-f $vcf){
		print STDERR "Group-directory: $group\tSample-directory: $sample_dir\tMerge-directory: $merge_dir\tvcf-name: $vcf_name\n" if ($group ne '');
		print STDERR "Sample-directory: $sample_dir\tMergedirectory: $merge_dir\tvcf-name: $vcf_name\n" if ($group eq '');
		print STDERR "Expected vcf path: $vcf\n";
		die "The specified vcf file is not found: Check the folders/files under your working directory and the specified options (-s and -id)\n";
	}
	push @sample_id, $ID;
	push @var_file, $vcf;
	$sample_id{$ID} = $vcf;
	if ($group ne ''){
		$group{$ID} = $group;
		$group_id{$group} = 1;
		$Gid{"$group:$ID"} = 1;
	}
	else{
		$Gid{$ID} = 1;
	}
}
close (FILE);

print STDERR "Merging vcf files of sample names: @sample_id\n";

my $total_sample_num = scalar @sample_id;

my %call;
my %call_subtype;
my %call_cons;
my %call_cons_len;
my %used_pos;
my %call_ave_len;
my %call_diverged;
my %call_diverged_ave;
my $count = 0;
my @header;
my %meitype;

foreach my $id (@sample_id){
    my $var_file = $sample_id{$id};
#print STDERR "$var_file\n";
	$count ++;
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
		chomp $line;
		if ($line =~ /^#/){
			if ($count == 1){
				push @header, $line if ($line !~ /^#CHROM/);
			}
			if ($line =~ /ID=INS:ME:(.+?),/){
				$meitype{$1} = 1;
			}
			elsif ($line =~ /ID=INS:NUMT/){
				$meitype{'NUMT'} = 1;
			}
			elsif ($line =~ /ID=INS:VEI/){
				$meitype{'VEI'} = 1;
			}
		    next;
		}
		my @line = split (/\t/, $line);
		my $chr = $line[0];
		if ($target_chr eq 'ALL'){
	        next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXY]+$/);
	    }
	    elsif (!exists $target_chr{$chr}){
	        next;
	    }
		next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
		next if ($gender{$id}) and ($gender{$id} eq 'F') and ($chr =~ /c*h*r*Y/i);
		my $pos = $line[1];
		my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
		my $subtype = '';
		$subtype = $1 if ($line[4] =~ /<(.+)>/);
		$subtype = $1 if ($subtype =~ /:([^:]+)$/);
		$subtype = $type if ($subtype eq '');
		$type = 'INS' if ($type =~ /ALU|LINE1|L1|SVA|HERVK|NUMT|VEI|MEI/);
		my $len = 0;
		$len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		next if ($type eq 'INV') and ($len > $max_inv_size);
		my $end = 0;
		$end = $pos + $len - 1 if ($type ne 'INS');
		my $chr2 = '';
		my $pos2 = 0;
		if ($type eq 'TRA'){
		    $chr2 = $1 if ($line[7] =~ /CHR2=(.+?);/);
		    $pos2 = $1 if ($line[7] =~ /POS2=(\d+)/);
		    $end = 0;
		    $len = 0;
		}
		my $tools = '';
		$tools = $1 if ($line[7] =~ /TOOLS=(.+)/);
		$tools = '' if ($tools =~ /,/);
		my $gt_gq = './.:0';
		$gt_gq = $1 if (@line > 9) and ($line[8] =~ /GT:GQ/) and ($line[9] =~ /([^:]+:[^:]+)/);
		my $htgt = '0/1';
		$gt_gq =~ s/^1\/0/$htgt/ if ($gt_gq =~ /^1\/0/);
		my $dr = 1;
		my $ds = 0;
		my $sr = 0;
		if (($line[8] =~ /DR:DS:SR$/) and ($line[9] =~ /:([\d\.]+):([\d\.]+):([\d\.]+)$/)){
		    $dr = $1;
		    $ds = $2;
		    $sr = $3;
		    if (($type eq 'DEL') and ($len >= 1000)){	# filter >= 1 Kb DELs with >= 0.9 DPR or >= 0.5 DPS
		    	next if ($dr >= 0.9);
		    	next if ($ds >= 0.5);
		    }
		    if (($type eq 'DEL') and ($len >= 100000)){	# filter >= 100 Kb DELs with >= 0.8 DPR or >= 0.35 DPS
		    	next if ($dr >= 0.8);
		    	next if ($ds >= 0.35);
		    }
		    if (($type eq 'DUP') and ($len >= 1000)){	# filter >= 1 Kb DUPs with <= 1.1 DPR or < 10 Kb DUPs with >= 0.5 DPS
		    	next if ($dr <= 1.1);
		    	next if ($len < 10000) and ($ds >= 0.5);
		    }
		    if (($type eq 'DUP') and ($len >= 10000)){	# filter >= 10 Kb DUPs with >= 0.3 DPS
		    	next if ($ds >= 0.3);
		    }
		}
		${${${$call{$type}}{$chr}}{$pos}}{$id} = $len;
		${${${$call_subtype{$type}}{$chr}}{$pos}}{$id} = "$subtype==$gt_gq==$dr==$ds==$sr";
    }
    close (FILE);
}

foreach my $type (keys %call){			# clustering DEL/DUP/INV with 0.5- to 2.0-fold difference in size at the same called positions
    foreach my $chr (keys %{$call{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call{$type}}{$chr}}){
		    my @len = ();
		    foreach my $id (keys %{${${$call{$type}}{$chr}}{$pos}}){
				my $len = ${${${$call{$type}}{$chr}}{$pos}}{$id};
				push @len, $len if ($type ne 'INS');
				push @{${${$call_cons{$type}}{$chr}}{$pos}}, "$id=$pos=$len";
			}
			if (@len == 1){
				${${$call_ave_len{$type}}{$chr}}{$pos} = $len[0];
			}
			elsif (@len > 1){
				my %clust_len;
				my %clust_ave;
				my %clust_cont;
				my $clust_num = 1;
				my $clust_len = 0;
				my $count_len = 0;
				foreach (sort {$a <=> $b} @len){
				    if ($clust_len == 0){
						push @{$clust_len{$clust_num}}, $_;
						$clust_len = $_;
				    }
				    elsif ($clust_len >= $_ * 0.5){
						push @{$clust_len{$clust_num}}, $_;
				    }
				    else{
						$clust_num ++;
						$clust_len = $_;
						push @{$clust_len{$clust_num}}, $_;
				    }
				}
				foreach my $clustn (sort {$a <=> $b} keys %clust_len){
				    my $sum_len = 0;
				    my $clust_cont = scalar @{$clust_len{$clustn}};
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$sum_len += $len1;
				    }
				    my $ave_len = int ($sum_len / $clust_cont);
				    $clust_ave{$clustn} = $ave_len;
				    $clust_cont{$clustn} = $clust_cont;
				}
				my $pre_ave = 0;
				my $pre_clustn = 0;
				foreach my $clustn (sort {$a <=> $b} keys %clust_ave){
				    my $avelen = $clust_ave{$clustn};
				    if ($pre_ave == 0){
						$pre_ave = $avelen;
						$pre_clustn = $clustn;
						next;
				    }
				    if ($avelen <= $pre_ave * 2){
						foreach my $len1 (@{$clust_len{$clustn}}){
						    push @{$clust_len{$pre_clustn}}, $len1;
						}
						my $new_avelen = int (($pre_ave * $clust_cont{$pre_clustn} + $avelen * $clust_cont{$clustn}) / ($clust_cont{$pre_clustn} + $clust_cont{$clustn}));
						$clust_cont{$pre_clustn} += $clust_cont{$clustn};
						delete $clust_len{$clustn};
						delete $clust_cont{$clustn};
						$clust_ave{$pre_clustn} = $new_avelen;
						$pre_ave = $new_avelen;
				    }
				    else{
						$pre_ave = $avelen;
						$pre_clustn = $clustn;
				    }
				}
				my $top_clust = 0;
				my %top_clust;
				my %low_clust;
				foreach my $clustn (sort {$clust_cont{$b} <=> $clust_cont{$a}} keys %clust_cont){
				    $top_clust = $clustn;
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$top_clust{$len1} = $clustn;
				    }
				    last;
				}
				my $new_clust = 0;
				foreach my $clustn (sort {$a <=> $b} keys %clust_len){
				    next if ($clustn == $top_clust);
				    $new_clust ++;
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$low_clust{$len1} = $new_clust;
				    }
				    ${${${$call_diverged_ave{$type}}{$chr}}{$pos}}{$new_clust} = $clust_ave{$clustn};
				}
				if (scalar keys %clust_len > 1){
				    my @info = (@{${${$call_cons{$type}}{$chr}}{$pos}});
				    delete ${${$call_cons{$type}}{$chr}}{$pos};
				    foreach (@info){
						my ($id1, $pos1, $len1) = split (/=/, $_);
						if (exists $top_clust{$len1}){
						    push @{${${$call_cons{$type}}{$chr}}{$pos}}, $_;
						    ${${$call_ave_len{$type}}{$chr}}{$pos} = $clust_ave{$top_clust};
						}
						else{
						    my $clustn = $low_clust{$len1};
						    push @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}}, $_;
						}
				    }
				}
				else{
				    my $sum_len = 0;
				    map{$sum_len += $_} @len;
				    my $ave_len2 = int ($sum_len / @len);
				    ${${$call_ave_len{$type}}{$chr}}{$pos} = $ave_len2;
				}
		    }
		}
    }
}
%call = ();

foreach my $type (keys %call_diverged){			# reassign the above extracted diverged sites to the most appropriate clusetr
    foreach my $chr (keys %{$call_diverged{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_diverged{$type}}{$chr}}){
		    foreach my $clustn (sort {$a <=> $b} keys %{${${$call_diverged{$type}}{$chr}}{$pos}}){
				my $avelen = ${${${$call_diverged_ave{$type}}{$chr}}{$pos}}{$clustn};
				my $end = $pos + $avelen - 1;
				my %len_diff;
				my $added_pos = 0;
				foreach my $pos2 (sort {$a <=> $b} keys %{${$call_ave_len{$type}}{$chr}}){
				    my $len2 = ${${$call_ave_len{$type}}{$chr}}{$pos2};
				    my $end2 = $pos2 + $len2 - 1;
				    next if ($end2 < $pos);
				    last if ($pos2 > $end);
				    my $diff = abs ($avelen - $len2);
				    if (($pos <= $pos2) and ($end >= $end2)){
						if ($len2 > $avelen * 0.5){
						    $len_diff{$pos2} = $diff;
						}
				    }
				    elsif (($pos >= $pos2) and ($end <= $end2)){
						if ($avelen > $len2 * 0.5){
						    $len_diff{$pos2} = $diff;
						}
				    }
				    elsif (($pos >= $pos2) and ($pos <= $end2)){
						if (($end2 - $pos > $avelen * 0.5) and ($end2 - $pos > $len2 * 0.5)){
						    $len_diff{$pos2} = $diff;
						}
				    }
				}
				if (scalar keys %len_diff > 0){
				    foreach my $pos3 (sort {$len_diff{$a} <=> $len_diff{$b}} keys %len_diff){
						push @{${${$call_cons{$type}}{$chr}}{$pos3}}, @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}};
						$added_pos = $pos3;
						last;
				    }
				}
				else{
				    my $new_pos = $pos;
				    my $flag = 0;
				    while ($flag == 0){
						if (!exists ${${$call_cons{$type}}{$chr}}{$new_pos}){
						    push @{${${$call_cons{$type}}{$chr}}{$new_pos}}, @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}};
						    $added_pos = $new_pos;
						    $flag = 1;
						}
						else{
						    $new_pos ++;
						}
				    }
				}
				if ($added_pos != 0){
				    my $sum_pos = 0;
				    my $sum_len = 0;
				    my $item_num = 0;
				    my @info = (@{${${$call_cons{$type}}{$chr}}{$added_pos}});
				    foreach (@{${${$call_cons{$type}}{$chr}}{$added_pos}}){
						my ($id1, $pos1, $len1) = split (/=/, $_);
						$sum_pos += $pos1;
						$sum_len += $len1;
						$item_num ++;
				    }
				    my $ave_pos = int ($sum_pos / $item_num);
				    my $ave_len = int ($sum_len / $item_num);
				    delete ${${$call_cons{$type}}{$chr}}{$added_pos};
				    delete ${${$call_ave_len{$type}}{$chr}}{$added_pos};
				    if ((!exists ${${$call_cons{$type}}{$chr}}{$ave_pos}) or ($added_pos == $ave_pos)){
						push @{${${$call_cons{$type}}{$chr}}{$ave_pos}}, @info;
						${${$call_ave_len{$type}}{$chr}}{$ave_pos} = $ave_len;
				    }
				    else{
						my $new_pos = $ave_pos;
						my $flag = 0;
						while ($flag == 0){
						    if (!exists ${${$call_cons{$type}}{$chr}}{$new_pos}){
								push @{${${$call_cons{$type}}{$chr}}{$new_pos}}, @info;
								${${$call_ave_len{$type}}{$chr}}{$new_pos} = $ave_len;
								$flag = 1;
						    }
						    else{
								$new_pos ++;
						    }
						}
				    }
				}
		    }
		}
    }
}
%call_ave_len = ();
%call_diverged = ();
%call_diverged_ave = ();
print STDERR "1st step completed:\n";

foreach my $type (keys %call_cons){	# merge positions present within 30-bp of consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $id_num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $sum_len = 0;
		    my $ave_len = 0;
		    if ($type ne 'INS'){
				foreach my $info (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				    my ($id, $pos1, $len1) = split (/=/, $info);
				    $sum_len += $len1;
				}
				$ave_len = int ($sum_len / $id_num);
		    }
			for (my $i = $pos + 1; $i <= $pos + 30; $i++){
				next if (exists ${${$used_pos{$type}}{$chr}}{$i});
				if (exists ${${$call_cons{$type}}{$chr}}{$i}){
				    if ($type ne 'INS'){
						my $sum_leni = 0;
						foreach my $infoi (@{${${$call_cons{$type}}{$chr}}{$i}}){
						    my ($idi, $posi, $leni) = split (/=/, $infoi);
						    $sum_leni += $leni;
						}
						my $ave_leni = int ($sum_leni / scalar @{${${$call_cons{$type}}{$chr}}{$i}});
						next if ($ave_len > $ave_leni * 2) or ($ave_leni > $ave_len * 2);
				    }
				    push @{${${$call_cons{$type}}{$chr}}{$pos}}, @{${${$call_cons{$type}}{$chr}}{$i}};
				    delete ${${$call_cons{$type}}{$chr}}{$i};
				    ${${$used_pos{$type}}{$chr}}{$i} = 1;
				}
		    }
		}
    }
}

foreach my $type (keys %call_cons){	# merge positions present within 100-bp of consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $id_num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $sum_len = 0;
		    my $ave_len = 0;
		    if ($type ne 'INS'){
				foreach my $info (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				    my ($id, $pos1, $len1) = split (/=/, $info);
				    $sum_len += $len1;
				}
				$ave_len = int ($sum_len / $id_num);
		    }
		    my %ids;
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id) = split (/=/, $item);
				$ids{$id} = 1;
		    }
		    for (my $i = $pos + 1; $i <= $pos + 100; $i++){
				next if (exists ${${$used_pos{$type}}{$chr}}{$i});
				if (exists ${${$call_cons{$type}}{$chr}}{$i}){
				    my $flag = 0;
				    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$i}}){		# next if samples at pos and pre_pos share a common sample id
						my ($id) = split (/=/, $item);
						if (exists $ids{$id}){
						    $flag = 1;
						    last;
						}
					}
					if ($type ne 'INS'){
						my $sum_leni = 0;
						foreach my $infoi (@{${${$call_cons{$type}}{$chr}}{$i}}){
						    my ($idi, $posi, $leni) = split (/=/, $infoi);
						    $sum_leni += $leni;
						}
						my $ave_leni = int ($sum_leni / scalar @{${${$call_cons{$type}}{$chr}}{$i}});
						next if ($ave_len > $ave_leni * 2) or ($ave_leni > $ave_len * 2);
					}
					if ($flag == 0){
						push @{${${$call_cons{$type}}{$chr}}{$pos}}, @{${${$call_cons{$type}}{$chr}}{$i}};
						delete ${${$call_cons{$type}}{$chr}}{$i};
						${${$used_pos{$type}}{$chr}}{$i} = 1;
					}
					else{
						if (($type eq 'INS') and ($chr eq '1')){

						}
				    }
				}
		    }
		}
    }
}

foreach my $type (keys %call_cons){		# assign a median position from multiple sample-positions as a consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    my $sample_num = @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $median_pos = $pos;
		    my %id_pos;
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id, $pos2, $len2) = split (/=/, $item);
				my $id_pos = "$id=$pos2";
				$id_pos{$id_pos} = $pos2;
		    }
		    my $sample_num_2 = scalar keys %id_pos;
		    my $half_num = int ($sample_num_2 * 0.5 + 0.5);
		    my $count = 0;
		    if ($sample_num_2 > 2){
				foreach my $id_pos (sort {$id_pos{$a} <=> $id_pos{$b}} keys %id_pos){
				    $count ++;
				    if ($count == $half_num){
						my ($id, $pos2) = split (/=/, $id_pos);
						$median_pos = $pos2;
				    }
				}
		    }
		    if ($pos != $median_pos){
				push @{${${$call_cons{$type}}{$chr}}{$median_pos}}, @{${${$call_cons{$type}}{$chr}}{$pos}};
				delete ${${$call_cons{$type}}{$chr}}{$pos};
				${${$used_pos{$type}}{$chr}}{$pos} = 1;
				delete ${${$used_pos{$type}}{$chr}}{$median_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pos});
		    }
		}
    }
}

foreach my $type (keys %call_cons){		# assign median SV length for each consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my @len;
		    my $sum_len = 0;
		    my $new_len = 0;
	        my $maxlen = 0;
		    my $num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    foreach (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my $len = $1 if ($_ =~ /=(\d+)$/);
				if (($len == 0) or ($len == 1)){
				    $num --;
				    next;
				}
				$sum_len += $len;
				push @len, $len;
	            $maxlen = $len if ($maxlen < $len);
		    }
		    if ($type eq 'INS'){
                if ($maxlen >= $min_sv_len){
                    $new_len = $maxlen;
                }
                else{
                    $new_len = 0;
                }
            }
		    else{
				if ($num == 1){
				    $new_len = $sum_len;
				}
				elsif ($num == 2){
				    $new_len = int ($sum_len / $num + 0.5);
				}
				elsif ($num >= 3){
				    my $num_half = int ($num * 0.5 + 0.5);
				    my $count = 0;
				    my $median = 0;
				    foreach my $len (sort {$a <=> $b} @len){
						$count ++;
						if ($count == $num_half){
						    $median = $len;
						}
				    }
				    ($new_len) = &cons_len (\@len, $median);
				}
		    }
		    ${${$call_cons_len{$type}}{$chr}}{$pos} = $new_len;
		}
    }
}

foreach my $type (keys %call_cons){		# merge neighboring consensus pos with < 200 bp distance for INS or with 50%-reciprocal overlap for the other types of SVs
    foreach my $chr (keys %{$call_cons{$type}}){
		my %pre_info;
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $len = ${${$call_cons_len{$type}}{$chr}}{$pos};
		    my $end = $pos + $len - 1;
		    my %select_pos;
		    my @match;
		    foreach my $pos2 (sort {$b <=> $a} keys %pre_info){
				next if ($pos == $pos2);
				next if (exists ${${$used_pos{$type}}{$chr}}{$pos2});
				my $len2 = $pre_info{$pos2};
				my $end2 = $pos2 + $len2 - 1;
				last if ($pos - $pos2 > 20000000);
				next if ($end2 + $ins_sd < $pos);
				my $ovlrate = 0;
				if ($type eq 'INS'){
				    $end = $pos;
				    if (abs ($pos - $pos2) <= $ins_sd){
						$ovlrate = 1 / abs ($pos - $pos2) if (abs ($pos - $pos2) > 0);
						$ovlrate = 1 if (abs ($pos - $pos2) == 0);
				    }
				}
				else{
				    if ((abs ($pos - $pos2) <= $var_sd) and (abs ($end - $end2) <= $var_sd)){
						if ($len >= $len2){
						    $ovlrate = $len2 / $len;
						}
						else{
						    $ovlrate = $len / $len2;
						}
				    }
				    elsif (($pos >= $pos2) and ($pos <= $end2)){
						if ($end <= $end2){
						    if ($len >= $len2 * $min_overlap_ratio){
								$ovlrate = $len / $len2;
						    }
						}
						else{
						    my $overlen = $end2 - $pos + 1;
						    if (($overlen >= $len * $min_overlap_ratio) and ($overlen >= $len2 * $min_overlap_ratio)){
								my $ovlrate1 = $overlen / $len;
								my $ovlrate2 = $overlen / $len2;
								$ovlrate = ($ovlrate1 + $ovlrate2) / 2;
						    }
						}
				    }
				}
				if ($ovlrate > 0){
				    push @match, "$pos2=$len2=$ovlrate";
				}
		    }
		    if (@match > 0){
				if (@match == 1){
				    my ($mpos, $mlen) = split (/=/, $match[0]);
				    $select_pos{$mpos} = $mlen if (exists ${${$call_cons{$type}}{$chr}}{$mpos});
				}
				else{
				    my $max_rate = 0;
				    my $best_pos = 0;
				    my $best_len = 0;
				    foreach (@match){
						my @info = split (/=/, $_);
						my $rate = $info[2];
						if ($rate > $max_rate){
						    $max_rate = $rate;
						    $best_pos = $info[0];
						    $best_len = $info[1];
						}
				    }
				    $select_pos{$best_pos} = $best_len if (exists ${${$call_cons{$type}}{$chr}}{$best_pos});
				}
		    }
		    if (scalar keys %select_pos > 0){		# check if samples at pos and pre_pos share a common sample id
				my @cur_info = ();
				my @pre_info = ();
				my @pre_pos = ();
				@cur_info = (@{${${$call_cons{$type}}{$chr}}{$pos}});
				foreach my $prepos (sort {$a <=> $b} keys %select_pos){
				    push @pre_info, @{${${$call_cons{$type}}{$chr}}{$prepos}};
				    push @pre_pos, $prepos;
				}
				my $overlap_flag = 0;
				my %ids;
				my %overlap_ids;
				my $overlap_num = 0;
				my $small_sample = 0;
				foreach my $item (@cur_info){
				    my ($id) = split (/=/, $item);
				    $ids{$id} = 1;
				}
				foreach my $item (@pre_info){
				    my ($pre_id) = split (/=/, $item);
				    if ((exists $ids{$pre_id}) and ($ids{$pre_id} == 1)){
						$ids{$pre_id} ++;
						$overlap_num ++;
						$overlap_ids{$pre_id} = 1;
				    }
				}
				if ($overlap_num <= 1){
				    $overlap_flag = 0;
				}
				elsif ($overlap_num > 1){
				    my @pos = ();
				    my @pre_pos2 = ();
				    my @overlap_pos = ();
				    my @overlap_pre_pos = ();
				    my @len = ();
				    my @pre_len = ();
				    my @overlap_len = ();
				    my @overlap_pre_len = ();
				    my @item;
				    my @pre_item;
				    my @overlap_item = ();
				    my @overlap_pre_item = ();
				    my $sum_pos = 0;
				    my $sum_pre_pos = 0;
				    my $sum_overlap_pos = 0;
				    my $sum_overlap_pre_pos = 0;
				    my $sum_len = 0;
				    my $sum_pre_len = 0;
				    my $sum_overlap_len = 0;
				    my $sum_overlap_pre_len = 0;
				    my %used_item;
				    my %used_pre_item;
				    foreach my $item (@cur_info){
						next if (exists $used_item{$item});
						my ($id, $pos1, $len1) = split (/=/, $item);
						if (exists $overlap_ids{$id}){
						    push @overlap_item, $item;
						    push @overlap_pos, $pos1;
						    push @overlap_len, $len1 if ($len1 >= $min_sv_len);
						}
						else{
						    push @item, $item;
						    push @pos, $pos1;
						    push @len, $len1 if ($len1 >= $min_sv_len);
						}
						$used_item{$item} = 1;
				    }
				    foreach my $item (@pre_info){
						next if (exists $used_pre_item{$item});
						my ($pre_id, $pre_pos1, $pre_len1) = split (/=/, $item);
						if (exists $ids{$pre_id}){
						    push @overlap_pre_item, $item;
						    push @overlap_pre_pos, $pre_pos1;
						    push @overlap_pre_len, $pre_len1 if ($pre_len1 >= $min_sv_len);
						}
						else{
						    push @pre_item, $item;
						    push @pre_pos2, $pre_pos1;
						    push @pre_len, $pre_len1 if ($pre_len1 >= $min_sv_len);
						}
						$used_pre_item{$item} = 1;
				    }
				    map{$sum_pos += $_} @pos;
				    map{$sum_pre_pos += $_} @pre_pos2;
				    map{$sum_overlap_pos += $_} @overlap_pos;
				    map{$sum_overlap_pre_pos += $_} @overlap_pre_pos;
				    map{$sum_len += $_} @len;
				    map{$sum_pre_len += $_} @pre_len;
				    map{$sum_overlap_len += $_} @overlap_len;
				    map{$sum_overlap_pre_len += $_} @overlap_pre_len;
				    my $ave_pos = 0;
				    my $ave_pre_pos = 0;
				    my $ave_overlap_pos = 0;
				    my $ave_overlap_pre_pos = 0;
				    my $ave_len = 0;
				    my $ave_pre_len = 0;
				    my $ave_overlap_len = 0;
				    my $ave_overlap_pre_len = 0;
				    $ave_pos = int ($sum_pos / @pos + 0.5) if (@pos > 0);
				    $ave_pre_pos = int ($sum_pre_pos / @pre_pos2 + 0.5) if (@pre_pos2 > 0);
				    $ave_overlap_pos = int ($sum_overlap_pos / @overlap_pos + 0.5) if (@overlap_pos > 0);
				    $ave_overlap_pre_pos = int ($sum_overlap_pre_pos / @overlap_pre_pos + 0.5) if (@overlap_pre_pos > 0);
				    $ave_len = int ($sum_len / @pos + 0.5) if (@len > 0);
				    $ave_pre_len = int ($sum_pre_len / @pre_len + 0.5) if (@pre_len > 0);
				    $ave_overlap_len = int ($sum_overlap_len / @overlap_len + 0.5) if (@overlap_len > 0);
				    $ave_overlap_pre_len = int ($sum_overlap_pre_len / @overlap_pre_len + 0.5) if (@overlap_pre_len > 0);
				    $ave_pos = $ave_overlap_pos if ($ave_pos == 0);
				    $ave_pre_pos = $ave_overlap_pre_pos if ($ave_pre_pos == 0);
				    $ave_len = $ave_overlap_len if ($ave_len == 0);
				    $ave_pre_len = $ave_overlap_pre_len if ($ave_pre_len == 0);
				    my $AB = abs ($ave_pos - $ave_pre_pos);
				    my $ovAovB = abs ($ave_overlap_pos - $ave_overlap_pre_pos);
				    my $AB_len = abs ($ave_len - $ave_pre_len);
				    my $ovAovB_len = abs ($ave_overlap_len - $ave_overlap_pre_len);
				    if ($type eq 'INS'){	# compare absolute distances between SV positions/lengths for overlapped and non-overlapped sample IDs
						if ($AB < $ovAovB){
						    $overlap_flag = 1;
						}
						else{
						    $overlap_flag = 0;
						}
				    }
				    else{
						if (($AB < $ovAovB) and ($AB_len < $ovAovB_len)){
						    $overlap_flag = 1;
						}
						else{
						    $overlap_flag = 0;
						}
				    }
				    if ($overlap_flag == 1){		# re-assign non-overlapped items of pre-pos and pos to overlapped pre-pos and overlapped pos
						foreach my $item (@item){
						    my ($id, $pos1, $len1) = split (/=/, $item);
						    if (abs ($pos1 - $ave_overlap_pos) <= abs ($pos1 - $ave_overlap_pre_pos)){
							push @overlap_item, $item;
						    }
						    else{
							push @overlap_pre_item, $item;
						    }
						}
						foreach my $item (@pre_item){
						    my ($id, $pos1, $len1) = split (/=/, $item);
						    if (abs ($pos1 - $ave_overlap_pos) <= abs ($pos1 - $ave_overlap_pre_pos)){
							push @overlap_item, $item;
						    }
						    else{
							push @overlap_pre_item, $item;
						    }
						}
						my $median_pos = 0;
						my $median_pre_pos = 0;
						my $median_len = 0;
						my $median_pre_len = 0;
						my @sum_pos = ();
						my @sum_pre_pos = ();
						my @sum_len = ();
						my @sum_pre_len = ();
						foreach my $item (@overlap_item){
						    my ($id, $pos1, $len1) = split (/=/, $item);
						    push @sum_pos, $pos1;
						    push @sum_len, $len1 if ($len1 >= $min_sv_len);
						}
						foreach my $item (@overlap_pre_item){
						    my ($id, $pos1, $len1) = split (/=/, $item);
						    push @sum_pre_pos, $pos1;
						    push @sum_pre_len, $len1 if ($len1 >= $min_sv_len);
						}
						my $half_num = int (@sum_pos * 0.5 + 0.5);
						my $half_pre_num = int (@sum_pre_pos * 0.5 + 0.5);
						my $count = 0;
						foreach my $pos1 (sort {$a <=> $b} @sum_pos){
						    $count ++;
						    if ($count == $half_num){
								$median_pos = $pos1;
								last;
						    }
						}
						$count = 0;
						foreach my $pos1 (sort {$a <=> $b} @sum_pre_pos){
						    $count ++;
						    if ($count == $half_pre_num){
								$median_pre_pos = $pos1;
								last;
						    }
						}
						$count = 0;
						foreach my $len1 (sort {$a <=> $b} @sum_len){
						    $count ++;
						    if ($count == $half_num){
								$median_len = $len1;
								last;
						    }
						}
						$count = 0;
						foreach my $len1 (sort {$a <=> $b} @sum_pre_len){
						    $count ++;
						    if ($count == $half_pre_num){
								$median_pre_len = $len1;
								last;
						    }
						}
						$median_pre_pos -- if ($median_pos == $median_pre_pos);
						delete ${${$call_cons{$type}}{$chr}}{$pos};
						map {delete ${${$call_cons{$type}}{$chr}}{$_}} @pre_pos;
						map {delete $pre_info{$_}} @pre_pos;
                        if (exists ${${$call_cons{$type}}{$chr}}{$median_pre_pos}){
                            while (1){
                                $median_pre_pos += 10;
                                last if (!exists ${${$call_cons{$type}}{$chr}}{$median_pre_pos});
                            }
                        }
                        if (exists ${${$call_cons{$type}}{$chr}}{$median_pos}){
                            while (1){
                                $median_pos += 10;
                                last if (!exists ${${$call_cons{$type}}{$chr}}{$median_pos});
                            }
                        }
						push @{${${$call_cons{$type}}{$chr}}{$median_pre_pos}}, @overlap_pre_item;
						${${$call_cons_len{$type}}{$chr}}{$median_pre_pos} = $median_pre_len;
						$pre_info{$median_pre_pos} = $median_pre_len;
						push @{${${$call_cons{$type}}{$chr}}{$median_pos}}, @overlap_item;
						${${$call_cons_len{$type}}{$chr}}{$median_pos} = $median_len;
						$pre_info{$median_pos} = $median_len;
						delete ${${$used_pos{$type}}{$chr}}{$median_pre_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pre_pos});
						delete ${${$used_pos{$type}}{$chr}}{$median_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pos});
						next;
				    }
				}
				if ($overlap_flag == 0){
				    my %pos;
				    my $sum_pos = 0;;
				    my $top_pos1 = 0;
				    my $top_pos2 = 0;
				    my $pos1_freq = 0;
				    my $pos2_freq = 0;
				    my $select_pos = 0;
				    my $sum_len2 = 0;
				    my $sum_len_num = 0;
				    my $ave_len2 = 0;
				    my $total_num = @cur_info + @pre_info;
				    delete ${${$call_cons{$type}}{$chr}}{$pos};
				    delete ${${$call_cons_len{$type}}{$chr}}{$pos};
				    map {delete ${${$call_cons{$type}}{$chr}}{$_}} @pre_pos;
				    map {delete ${${$call_cons_len{$type}}{$chr}}{$_}} @pre_pos;
				    map {delete $pre_info{$_}} @pre_pos;
					
				    foreach (@cur_info){
						my ($id, $bp, $len) = split (/=/, $_);
						$pos{$bp} ++;
						$sum_pos += $bp;
						$sum_len2 += $len if ($len >= $min_sv_len);
						$sum_len_num ++ if ($len >= $min_sv_len);
				    }
				    foreach (@pre_info){
						my ($id, $bp, $len) = split (/=/, $_);
						$pos{$bp} ++;
						$sum_pos += $bp;
						$sum_len2 += $len if ($len >= $min_sv_len);
						$sum_len_num ++ if ($len >= $min_sv_len);
				    }
				    $ave_len2 = int ($sum_len2 / $sum_len_num) if ($sum_len_num > 0);
				    
				    foreach my $pos2 (sort {$pos{$b} <=> $pos{$a}} keys %pos){
						$top_pos2 = $pos2 if ($top_pos1 > 0) and ($top_pos2 == 0);
						$pos2_freq = $pos{$pos2} if ($pos1_freq > 0) and ($pos2_freq == 0);
						$top_pos1 = $pos2 if ($top_pos1 == 0);
						$pos1_freq = $pos{$pos2} if ($pos1_freq == 0);
						last if ($top_pos2 > 0);
				    }
				    if ($total_num == 2){
						$select_pos = int (($top_pos1 + $top_pos2) / 2 + 0.5);
				    }
				    elsif ($total_num == 3){
						if ($pos1_freq >= 2){
						    $select_pos = $top_pos1;
						}
						else{
						    $select_pos = int ($sum_pos / 3 + 0.5);
						}
				    }
				    else{
						if (($pos1_freq >= 2) and ($pos1_freq >= $pos2_freq * 2)){
						    $select_pos = $top_pos1;
						}
						elsif ($pos2_freq >= 2){
						    $select_pos = int (($top_pos1 + $top_pos2) / 2 + 0.5);
						}
						else{
						    $select_pos = int ($sum_pos / $total_num + 0.5);
						}
				    }
                    if (exists ${${$call_cons{$type}}{$chr}}{$select_pos}){
                        while (1){
                            $select_pos += 10;
                            last if (!exists ${${$call_cons{$type}}{$chr}}{$select_pos});
                        }
                    }
				    push @{${${$call_cons{$type}}{$chr}}{$select_pos}}, @pre_info;
				    push @{${${$call_cons{$type}}{$chr}}{$select_pos}}, @cur_info;
				    ${${$call_cons_len{$type}}{$chr}}{$select_pos} = $ave_len2;
				    map {${${$used_pos{$type}}{$chr}}{$_} = 1} @pre_pos;
				    ${${$used_pos{$type}}{$chr}}{$pos} = 1;
				    delete ${${$used_pos{$type}}{$chr}}{$select_pos} if (exists ${${$used_pos{$type}}{$chr}}{$select_pos});
				    $pre_info{$select_pos} = $ave_len2;
				}
		    }
		    else{
				$pre_info{$pos} = $len;
		    }
		}
    }
}

print STDERR "2nd step completed:\n";

system ("mkdir $output_dir") if (!-d $output_dir);
my $temp_dir = "$output_dir/tmp";
system ("rm -r $temp_dir") if (-d $temp_dir);
system ("mkdir $temp_dir") if ($total_sample_num >= 100);

my %vcf_type;
my %vcf_cons;

foreach my $type (keys %call_cons){
    foreach my $chr (keys %{$call_cons{$type}}){
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    next if (scalar @{${${$call_cons{$type}}{$chr}}{$pos}} == 0);
		    my %id_pos;
		    my %id_pos_2;
		    my %ins_subtype;
		    my $id_pos = '';
		    my $ins_subtype = '';
		    my $subtype = '';
		    my $len = ${${$call_cons_len{$type}}{$chr}}{$pos};
		    my $end = $pos + $len - 1 if ($type ne 'INS');
		    $end = $pos if ($type eq 'INS');
		    my $gap_overlap = 0;
		    foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
				my $gend = ${$gap{$chr}}{$gstart};
				if (($pos >= $gstart) and ($pos <= $gend)){
				    if ($type eq 'INS'){
						$gap_overlap = 1;
						last;
				    }
				    else{
						my $overlap = $gend - $pos + 1 if ($end >= $gend);
						$overlap = $end - $pos + 1 if ($end < $gend);
						if ($overlap >= $len * 0.5){
						    $gap_overlap = 1;
						    last;
						}
				    }
				}
				elsif (($type ne 'INS') and ($gstart >= $pos) and ($gstart <= $end)){
				    my $overlap = $end - $gstart + 1 if ($end <= $gend);
				    $overlap = $gend - $gstart + 1 if ($end > $gend);
				    if ($overlap >= $len * 0.5){
						$gap_overlap = 1;
						last;
				    }
				}
				last if ($gstart > $end);
		    }
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id, $pos2, $len2) = split (/=/, $item);
				my $group = '';
				$group = $group{$id} if (exists $group{$id});
				my $Gid = $id;
                $Gid = "$group:$id" if ($group ne '') and (scalar keys %group > 1);
                ${$id_pos{$Gid}}{$pos2} = $len2;
                $id_pos_2{$Gid} = $pos2;
				if ($type eq 'INS'){
				    my ($subtype) = split (/==/, ${${${$call_subtype{$type}}{$chr}}{$pos2}}{$id});
				    ${$ins_subtype{$subtype}}{$Gid} = 1;
				}
		    }
            foreach my $Gid (sort keys %id_pos){
                foreach my $pos2 (keys %{$id_pos{$Gid}}){
                    my $len2 = ${$id_pos{$Gid}}{$pos2};
                    my $Gid_pos = "$Gid==$pos2==$len2";
                    $id_pos .= $Gid_pos . ',';
                }
            }
            $id_pos =~ s/,$//;
            my $sample_num = scalar keys %id_pos;
            my $half_num = int ($sample_num * 0.5 + 0.5);
            my $count = 0;
            my $new_pos = $pos;
            if ($sample_num >= 3){
                foreach my $Gid (sort {$id_pos_2{$a} <=> $id_pos_2{$b}} keys %id_pos_2){		# determine median position of multiple sample-positions
                    $count ++;
                    if ($count == $half_num){
                        $new_pos = $id_pos_2{$Gid};
                        last;
                    }
                }
            }
		    
		    next if ($gap_overlap == 1);
		    
		    if ($type eq 'INS'){
				if (scalar keys %ins_subtype == 0){
				    $ins_subtype = 'INS';
				}
				elsif (scalar keys %ins_subtype == 1){
				    foreach my $subtype (sort keys %ins_subtype){
						$ins_subtype = $subtype;
				    }
				}
				elsif (scalar keys %ins_subtype > 1){
				    my %subtype_rank;
				    my $subtype_1st = '';
				    my $subtype_2nd = '';
				    foreach my $subtype (sort keys %ins_subtype){
						$subtype_rank{$subtype} = scalar keys %{$ins_subtype{$subtype}};
				    }
				    foreach my $subtype (sort {$subtype_rank{$b} <=> $subtype_rank{$a}} keys %subtype_rank){
						$subtype_2nd = $subtype if ($subtype_2nd eq '') and ($subtype_1st ne '');
						$subtype_1st = $subtype if ($subtype_1st eq '');
						last if ($subtype_2nd ne '');
				    }
				    if ($subtype_1st ne 'INS'){
						$ins_subtype = $subtype_1st;
				    }
				    else{
						if ($subtype_rank{$subtype_2nd} >= $subtype_rank{$subtype_1st} * 0.2){
						    $ins_subtype = $subtype_2nd;
						}
						else{
						    $ins_subtype = $subtype_1st;
						}
				    }
				}
		    }
		    $subtype = $type;
		    $subtype = $ins_subtype if ($type eq 'INS');
		    my $alt = "<$type>" if ($type eq 'DEL') or ($type eq 'DUP') or ($type eq 'INV');
		    $alt = '<INS>' if ($type eq 'INS') and ($subtype eq 'INS');
		    $alt = "<INS:ME:$subtype>" if ($subtype ne $type);
		    $alt = '<INS:NUMT>' if ($subtype eq 'NUMT');
		    $alt = '<INS:VEI>' if ($subtype eq 'VEI');
		    $len = 0 - $len if ($type eq 'DEL');
		    if (exists ${${$vcf_cons{$chr02d}}{$new_pos}}{$type}){
				while (1){
				    $new_pos ++;
				    last if (!exists ${${$vcf_cons{$chr02d}}{$new_pos}}{$type});
				}
		    }
            ${${$vcf_type{$type}}{$chr02d}}{$new_pos} = "$chr\t$new_pos\t.\t.\t$alt\t.\tPASS\tSVTYPE=$type;SVLEN=$len;SC=$sample_num;SAMPLES=$id_pos";
		}
    }
}
%call_cons = ();
%call_cons_len = ();

print STDERR "3rd step completed:\n";

foreach my $type (keys %vcf_type){              # adjust variants for each sample between neighboring positions
	next if ($type eq 'INS');
    foreach my $chr (keys %{$vcf_type{$type}}){
        my $pre_pos = 0;
        my $pre_end = 0;
        my $pre_len = 0;
        my $pre_line = '';
        my $chr2 = $chr;
        $chr2 =~ s/^0*//;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
            my $line = ${${$vcf_type{$type}}{$chr}}{$pos};
            my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            if (($pre_pos > 0) and ($pos < $pre_end)){
				my $overlap = $pre_end - $pos + 1;
				$overlap = $len if ($end < $pre_end);
				if (($overlap >= $pre_len * $min_overlap_ratio) and ($overlap >= $len * $min_overlap_ratio)){
				    my $idpos = $1 if ($line =~ /SAMPLES=(.+)/);
				    my $pre_idpos = $1 if ($pre_line =~ /SAMPLES=(.+)/);
				    my @idpos = split (/,/, $idpos);
				    my @pre_idpos = split (/,/, $pre_idpos);
				    my $sn = scalar @idpos;
				    my $pre_sn = scalar @pre_idpos;
				    my %add_idpos;
				    my %add_preidpos;
				    my @add_preidpos;
				    my @add_idpos;
				    foreach (@pre_idpos){
						my ($gid, $pos2, $len2) = split (/==/, $_);
						my $diff_rate1 = abs ($pos2 - $pre_pos) / $pre_pos;
						$diff_rate1 = $diff_rate1 * (abs ($len2 - $pre_len) / $pre_len);
						my $diff_rate2 = abs ($pos2 - $pos) / $pos;
						$diff_rate2 = $diff_rate2 * (abs ($len2 - $len) / $len);
						if ($diff_rate1 > $diff_rate2){
						    $add_idpos{$_} = 1;
						    push @add_idpos, "$_|$diff_rate1==$diff_rate2";
						}
				    }
				    foreach (@idpos){
						my ($gid, $pos2, $len2) = split (/==/, $_);
						my $diff_rate1 = abs ($pos2 - $pre_pos) / $pre_pos;
						$diff_rate1 = $diff_rate1 * (abs ($len2 - $pre_len) / $pre_len);
						my $diff_rate2 = abs ($pos2 - $pos) / $pos;
						$diff_rate2 = $diff_rate2 * (abs ($len2 - $len) / $len);
						if ($diff_rate1 < $diff_rate2){
						    $add_preidpos{$_} = 1;
						    push @add_preidpos, "$_|$diff_rate1==$diff_rate2";
						}
				    }
				    foreach (@pre_idpos){
						$add_preidpos{$_} = 1 if (!exists $add_idpos{$_});
				    }
				    foreach (@idpos){
						$add_idpos{$_} = 1 if (!exists $add_preidpos{$_});
				    }
				    if (@pre_idpos != scalar keys %add_preidpos){
						my $new_sn = 0;
						my $new_presn = 0;
						my $new_pos = 0;
						my $new_prepos = 0;
						my $new_len = 0;
						my $new_prelen = 0;
						my $sum_pos = 0;
						my $sum_prepos = 0;
						my $sum_len = 0;
						my $sum_prelen = 0;
						my $new_idpos = '';
						my $new_preidpos = '';
						$new_presn = scalar keys %add_preidpos;
						$new_sn = scalar keys %add_idpos;
						foreach my $idpos2 (sort keys %add_preidpos){
						    my ($gid, $pos2, $len2) = split (/==/, $idpos2);
						    $sum_prepos += $pos2;
						    $sum_prelen += $len2;
						    $new_preidpos .= "$idpos2,"
						}
						$new_prepos = int ($sum_prepos / $new_presn + 0.5) if ($new_presn > 0);
						$new_prelen = int ($sum_prelen / $new_presn + 0.5) if ($new_presn > 0);
						$new_preidpos =~ s/,$//;
						foreach my $idpos2 (sort keys %add_idpos){
						    my ($gid, $pos2, $len2) = split (/==/, $idpos2);
						    $sum_pos += $pos2;
						    $sum_len += $len2;
						    $new_idpos .= "$idpos2,"
						}
						$new_pos = int ($sum_pos / $new_sn + 0.5) if ($new_sn > 0);
						$new_len = int ($sum_len / $new_sn + 0.5) if ($new_sn > 0);
						$new_idpos =~ s/,$//;
						delete ${${$vcf_type{$type}}{$chr}}{$pre_pos};
						delete ${${$vcf_type{$type}}{$chr}}{$pos};
						${${$vcf_type{$type}}{$chr}}{$new_prepos} = "$chr2\t$new_prepos\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$new_prelen;SC=$new_presn;SAMPLES=$new_preidpos" if ($new_presn > 0);
						${${$vcf_type{$type}}{$chr}}{$new_pos} = "$chr2\t$new_pos\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$new_len;SC=$new_sn;SAMPLES=$new_idpos" if ($new_sn > 0);
						if ($new_sn > 0){
						    $pre_pos = $new_pos;
						    $pre_end = $new_pos + $new_len - 1;
						    $pre_len = $new_len;
						    $pre_line = ${${$vcf_type{$type}}{$chr}}{$new_pos};
						}
						else{
						    $pre_pos = $new_prepos;
						    $pre_end = $new_prepos + $new_prelen - 1;
						    $pre_len = $new_prelen;
						    $pre_line = ${${$vcf_type{$type}}{$chr}}{$new_prepos};
						}
						next;
				    }
				}
            }
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_len = $len;
            $pre_line = $line;
        }
    }
}

foreach my $type (keys %vcf_type){                  # merge overlapping variants again
	next if ($type eq 'INS');
    foreach my $chr (keys %{$vcf_type{$type}}){
        my $pre_pos = 0;
        my $pre_end = 0;
        my $pre_len = 0;
        my $pre_line = '';
        my $chr2 = $chr;
        $chr2 =~ s/^0*//;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
            my $line = ${${$vcf_type{$type}}{$chr}}{$pos};
            my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            if (($pre_pos > 0) and ($pos < $pre_end)){
                my $idpos = $1 if ($line =~ /SAMPLES=(.+)/);
                my $pre_idpos = $1 if ($pre_line =~ /SAMPLES=(.+)/);
                my @idpos = split (/,/, $idpos);
                my @pre_idpos = split (/,/, $pre_idpos);
                my $sn = scalar @idpos;
                my $pre_sn = scalar @pre_idpos;
                my $overlap = $pre_end - $pos + 1;
                $overlap = $len if ($end < $pre_end);
                if (($overlap >= $len * $min_overlap_ratio) and ($overlap >= $pre_len * $min_overlap_ratio)){
                    my %ids;
                    my $new_sn = 0;
                    my $new_idpos = '';
                    my $new_pos = 0;
                    my $new_len = 0;
                    my $sum_pos = 0;
                    my $sum_len = 0;
                    my @len;
                    if ($sn >= $pre_sn){
                        foreach (@pre_idpos){
                            my ($id) = split (/==/, $_);
                            $ids{$id} = $_;
                        }
                        foreach (@idpos){
                            my ($id) = split (/==/, $_);
                            $ids{$id} = $_;
                        }
                    }
                    else{
                        foreach (@idpos){
                            my ($id) = split (/==/, $_);
                            $ids{$id} = $_;
                        }
                        foreach (@pre_idpos){
                            my ($id) = split (/==/, $_);
                            $ids{$id} = $_;
                        }
                    }
                    $new_sn = scalar keys %ids;
                    foreach my $gid (sort keys %ids){
                        my ($gid2, $pos2, $len2) = split (/==/, $ids{$gid});
                        $new_idpos .= "$ids{$gid},";
                        $sum_pos += $pos2;
                        push @len, $len2 if ($len2 > 1);
                    }
                    if (@len > 0){
                    	map{$sum_len += $_} @len;
                    	$new_len = int ($sum_len / @len + 0.5);
                    }
                    $new_pos = int ($sum_pos / $new_sn + 0.5);
                    $new_idpos =~ s/,$//;
                    delete ${${$vcf_type{$type}}{$chr}}{$pre_pos};
                    delete ${${$vcf_type{$type}}{$chr}}{$pos};
                    if (exists ${${$vcf_type{$type}}{$chr}}{$new_pos}){
                        while (1){
                            $new_pos += 10;
                            last if (!exists ${${$vcf_type{$type}}{$chr}}{$new_pos});
                        }
                    }
                    ${${$vcf_type{$type}}{$chr}}{$new_pos} = "$chr2\t$new_pos\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$new_len;SC=$new_sn;SAMPLES=$new_idpos";
                    $pre_pos = $new_pos;
                    $pre_end = $new_pos + $new_len - 1;
                    $pre_len = $new_len;
                    $pre_line = ${${$vcf_type{$type}}{$chr}}{$new_pos};
                    next;
                }
            }
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_len = $len;
            $pre_line = $line;
        }
    }
}

foreach my $chr (keys %{$vcf_type{'INS'}}){				# merge proximal INS and DUP to DUP again
	my $chr2 = $chr;
    $chr2 =~ s/^0*//;
	foreach my $pos (keys %{${$vcf_type{'INS'}}{$chr}}){
		my $match_flag = 0;
		foreach my $dpos (keys %{${$vcf_type{'DUP'}}{$chr}}){
			last if ($dpos > $pos + 10);
			my $dlen = $1 if (${${$vcf_type{'DUP'}}{$chr}}{$dpos} =~ /SVLEN=(\d+)/);
			my $dend = $dpos + $dlen - 1;
			next if ($dend < $pos - 10);
			my $flag = 0;
			if (abs ($dpos - $pos) <= 10){
				$flag = 1;
			}
			elsif (abs ($dend - $pos) <= 10){
				$flag = 2;
			}
			if ($flag > 0){
				my $ins_idpos = $1 if (${${$vcf_type{'INS'}}{$chr}}{$pos} =~ /SAMPLES=(.+)/);
                my $dup_idpos = $1 if (${${$vcf_type{'DUP'}}{$chr}}{$dpos} =~ /SAMPLES=(.+)/);
                my @ins_idpos = split (/,/, $ins_idpos);
                my @dup_idpos = split (/,/, $dup_idpos);
                my %ids;
                foreach (@dup_idpos){
                	my ($gid) = split (/==/, $_);
                	$ids{$gid} = $_;
                }
                foreach (@ins_idpos){
                	my ($gid, $pos2, $len2) = split (/==/, $_);
                	if (!exists $ids{$gid}){
                		my $new_dpos = $pos2;
                		$new_dpos = $dpos if ($flag == 2);
                		my $info = "$gid==$new_dpos==$dlen";
                		$ids{$gid} = $info;
                		my $id = $gid;
                		$id = $1 if ($gid =~ /[^:]+:(.+)/);
                		${${${$call_subtype{'DUP'}}{$chr2}}{$new_dpos}}{$id} = ${${${$call_subtype{'INS'}}{$chr2}}{$pos2}}{$id};
                	}
                }
                my $new_idpos = '';
                my $new_sn = scalar keys %ids;
                foreach my $gid (sort keys %ids){
                    my ($gid2, $pos2, $len2) = split (/==/, $ids{$gid});
                    $new_idpos .= "$ids{$gid},";
                }
                $new_idpos =~ s/,$//;
                ${${$vcf_type{'DUP'}}{$chr}}{$dpos} = "$chr2\t$dpos\t.\t.\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=$dlen;SC=$new_sn;SAMPLES=$new_idpos";
                $match_flag = 1;
                last;
			}
			if ($match_flag == 1){
				delete ${${$vcf_type{'INS'}}{$chr}}{$pos};
				next;
			}
		}
	}
}

foreach my $type (keys %vcf_type){
    foreach my $chr (keys %{$vcf_type{$type}}){
        foreach my $pos (keys %{${$vcf_type{$type}}{$chr}}){
            ${${$vcf_cons{$chr}}{$pos}}{$type} = ${${$vcf_type{$type}}{$chr}}{$pos};
        }
    }
}

print STDERR "4th step completed:\n";

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

my $version = 1.7;
$version = $1 if ($Bin =~ /MOPline[_\-\.](v[\d\.]+)/i);

my $GID_str = '';

foreach my $gid (sort keys %Gid){
    $GID_str .= "$gid\t";
}
$GID_str =~ s/\t$//;

my $AN = $total_sample_num * 2;

my $merge_out = "$output_dir/$out_prefix.All-samples.vcf";

open (OUT, "> $merge_out");
my %mei_h;
foreach my $hline (@header){
	if ($hline =~ /fileDate=/){
		$hline =~ s/fileDate=\d+/fileDate=$time/;
	}
	elsif ($hline =~ /source=MOPline/){
		$hline =~ s/source=MOPline.+/source=MOPline-$version/;
	}
	next if ($hline =~ /ID=READS|ID=TOOLS/);
	if ($hline =~ /ID=INS:ME:(.+?),/){
		$mei_h{$1} = 1;
	}
	elsif ($hline =~ /ID=INS:NUMT/){
		$mei_h{'NUMT'} = 1;
	}
	elsif ($hline =~ /ID=INS:VEI/){
		$mei_h{'VEI'} = 1;
	}
	if ($hline =~ /ID=INV/){
		if (scalar keys %meitype > scalar keys %mei_h){
			foreach my $mei (sort keys %meitype){
				next if (exists $mei_h{$mei});
				print OUT "##ALT=<ID=INS:ME:$mei,Description=\"Insertion of $mei element\">\n";
			}
		}
	}
	print OUT "$hline\n";
	if ($hline =~ /ID=END/){
		print OUT "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n";
		print OUT "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">\n";
		print OUT "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n";
		print OUT "##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Number of samples carrying the SV allele\">\n";
		print OUT "##INFO=<ID=DPR,Number=1,Type=Float,Description=\"Mean ratio of read depth in CNV region to that in its flanking regions\">\n";
		print OUT "##INFO=<ID=SR,Number=1,Type=Float,Description=\"Mean ratio of breakpoints from soft-clipped reads in the read depth at the window\">\n";
	}
}
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$GID_str\n";

foreach my $chr (sort keys %vcf_cons){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf_cons{$chr}}){
        foreach my $type (sort keys %{${$vcf_cons{$chr}}{$pos}}){
            my @line = split (/\t/, ${${$vcf_cons{$chr}}{$pos}}{$type});
            my $id_set = $1 if ($line[7] =~ /SAMPLES=([^;]+)/);
            my @id_set = split (/,/, $id_set);
            my %id_info;
            my $sample_info = '';
            my @DR = ();
            my @SR = ();
            my @len;
            my $ac = 0;
            foreach (@id_set){
                my ($gid, $ipos, $ilen) = split (/==/, $_);
                my $id = $gid;
                $id = $1 if ($gid =~ /[^:]+:(.+)/);
print STDERR "$chr2\t$type\t$pos\t$ipos\t$id\n" if (!exists ${${${$call_subtype{$type}}{$chr2}}{$ipos}}{$id});
                my ($subtype, $gt_gq, $dr, $ds, $sr) = split (/==/, ${${${$call_subtype{$type}}{$chr2}}{$ipos}}{$id});
                $id_info{$gid} = "$gt_gq:$ipos:$ilen:$dr:$ds:$sr";
                push @DR, $dr;
                push @SR, $sr;
                push @len, $ilen;
            }
            foreach my $gid (sort keys %Gid){
                if (exists $id_info{$gid}){
                    $sample_info .= "$id_info{$gid}" . "\t";
                    my ($gt) = split (/:/, $id_info{$gid});
				    if (($gt eq '0/1') or ($gt eq '1/0') or ($gt eq './.')){
				    	$ac ++;
				    }
				    elsif ($gt eq '1/1'){
				    	$ac += 2;
				    }
                }
                else{
                    $sample_info .= "0/0:0:0:0:0:0:0" . "\t";
                }
            }
            $sample_info =~ s/\t$//;
            my $line = ${${$vcf_cons{$chr}}{$pos}}{$type};
            $line =~ s/;SAMPLES=[^;]+/\tGT:GQ:VP:VL:DR:DS:SR\t$sample_info/;
            
            my $sample_num = $1 if ($line[7] =~ /SC=(\d+)/);
            $line =~ s/SC=\d+;/SC=$sample_num/ if ($line[7] =~ /SC=$sample_num;/);
            
            if (($type ne 'INS') and (@len > 2)){
                my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
                my ($new_len) = &cons_len (\@len, $len);
                $new_len = 1 - $new_len if ($type eq 'DEL');
                $line =~ s/SVLEN=-\d+/SVLEN=$new_len/;
            }
            if ($line !~ /END=/){
                my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
                my $end = $pos + $len - 1;
                $end = $pos if ($type eq 'INS');
                $line =~ s/;SC=/;END=$end;SC=/;
            }
            my $sum_dr = 0;
            my $ave_dr = 0;
            my $sum_sr = 0;
            my $ave_sr = 0;
            map{$sum_dr += $_} @DR;
            map{$sum_sr += $_} @SR;
            $ave_dr = int ($sum_dr / @DR * 100 + 0.5) / 100;
            $ave_sr = int ($sum_sr / @SR * 100 + 0.5) / 100;
            my @line2 = split (/\t/, $line);
            my $af = int ($ac / $total_sample_num * 0.5 * 10000 + 0.5) / 10000;
		    $line2[7] .= ";AF=$af;AC=$ac;AN=$AN;DPR=$ave_dr;SR=$ave_sr";
            $line = join ("\t", @line2);
            print OUT $line, "\n";
        }
    }
}
close (OUT);

%vcf_cons = ();


sub cons_len{	# determine consensus length of SV
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
