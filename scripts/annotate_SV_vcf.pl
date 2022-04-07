#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $data_dir = "$Bin/../Data";

my $sv_vcf = '';

my $max_flank5 = 50000;
my $max_flank3 = 50000;
my $flank = $max_flank5;
$flank = $max_flank3 if ($max_flank3 > $max_flank5);

my $core_flank5 = 5000;
my $core_flank3 = 5000;

my $out_prefix = 'out';

my $target_chr = 'ALL';

my $non_human = 0;

my $build = 37;

my $help;

my $temp_dir = 'temp';

my $ref_gff = '';

GetOptions(
    'sv|v=s' => \$sv_vcf,
    'ref|r=s' => \$ref_gff,
    'flank5|f5=i' => \$max_flank5,
    'flank3|f3=i' => \$max_flank3,
    'core5|c5=i' => \$core_flank5,
    'core3|c3=i' => \$core_flank3,
    'target|t=s' => \$target_chr,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'build|b=i' => \$build,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  annotate_SV_vcf.pl -v <vcf file> -p <prefix of output vcf>
  output vcf files: ${out_prefix}.annot.vcf and ${out_prefix}.AS.annot.vcf

  Options:
   --sv or -v <STR>         vcf file of SVs
   --ref or -r <STR>        ref gff3 annotation file from Ensembl [default for human: Data/Homo_sapiens.GRCh37.87.gff3.gz or Data/Homo_sapiens.GRCh38.104.gff3.gz]
   --flank5 or -f5 <INT>    maximum size (bp) of 5' distal flanking gene region, that overlaps SV [default: 50000]
   --flank3 or -f3 <INT>    maximum size (bp) of 3' distal flanking gene region, that overlaps SV [default: 50000]
   --core5 or -c5 <INT>     maximum size (bp) of 5' proximal flanking gene region, that overlaps SV [default: 5000]
   --core3 or -c3 <INT>     maximum size (bp) of 3' proximal flanking gene region, that overlaps SV [default: 5000]
   --non_human or -nh <INT>  samples are non-human species (0: human, 1: non-human) [default: 0]
   --build <INT>            reference build (GRCh37, GRCh38) number (37 or 38) [default: 37]
   --target or -t <STR>     target chromosome [default: ALL]
   --prefix or -p <STR>     outpout prefix
   --help or -h             output help message
   
=cut

if ($non_human == 0){
	$ref_gff = "$data_dir/Homo_sapiens.GRCh37.87.gff3.gz";
	$ref_gff = "$data_dir/Homo_sapiens.GRCh38.104.gff3.gz" if ($build == 38);
}

die "input vcf file is not specified: \n" if ($sv_vcf eq '');

system ("mkdir $temp_dir") if (!-d $temp_dir) and ($target_chr ne 'ALL');

my $flank5_dist_str = '';
my $flank3_dist_str = '';
my $flank5_prox_str = '';
my $flank3_prox_str = '';
if ($max_flank5 >= 1000000){
	my $scale = int ($max_flank5 / 1000000);
	$flank5_dist_str = 'flank5_' . $scale . 'Mb';
}
elsif ($max_flank5 >= 1000){
	my $scale = int ($max_flank5 / 1000);
	$flank5_dist_str = 'flank5_' . $scale . 'Kb';
}
else{
	$flank5_dist_str = 'flank5_' . $max_flank5 . 'bp';
}
if ($max_flank3 >= 1000000){
	my $scale = int ($max_flank3 / 1000000);
	$flank3_dist_str = 'flank3_' . $scale . 'Mb';
}
elsif ($max_flank3 >= 1000){
	my $scale = int ($max_flank3 / 1000);
	$flank3_dist_str = 'flank3_' . $scale . 'Kb';
}
else{
	$flank3_dist_str = 'flank3_' . $max_flank3 . 'bp';
}
if ($core_flank5 >= 1000){
	my $scale = int ($core_flank5 / 1000);
	$flank5_prox_str = 'flank5_' . $scale . 'Kb';
}
else{
	$flank5_prox_str = 'flank5_' . $core_flank5 . 'bp';
}
if ($core_flank3 >= 1000){
	my $scale = int ($core_flank3 / 1000);
	$flank3_prox_str = 'flank3_' . $scale . 'Kb';
}
else{
	$flank3_prox_str = 'flank3_' . $core_flank3 . 'bp';
}

my %gene;
my %ORF;
my %CDS;
my %exon;
my %utr5;
my %utr3;
my %annot;
my $total_samples = 0;;
my $flag = 0;
my $gene_id = '';
my %ref_chr;
my $chr_prefix = '';

open (FILE, $sv_vcf) or die "$sv_vcf is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^##contig=<ID=(.+),/){
		my $chr = $1;
		$ref_chr{$chr} = 1;
		$chr_prefix = $1 if ($chr =~ /^(.+?)[\dIVX]+/);
	}
	last if ($line !~/^#/);
}
close (FILE);
#print STDERR "chr-prefix: $chr_prefix\n";


open (FILE, "gzip -dc $ref_gff |") or die "$ref_gff file is not found: $!\n" if ($ref_gff =~ /\.gz$/);
open (FILE, $ref_gff) or die "$ref_gff file is not found: $!\n" if ($ref_gff !~ /\.gz$/);
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#{1,2}[^#]/);
	if ($line =~ /^###/){
		if (exists $CDS{$gene_id}){
			delete $exon{$gene_id};
		}
		$flag = 1;
		$gene_id = '';
		next;
	}
	my @line = split (/\s+/, $line);
	my $chr = $line[0];
	if (!exists $ref_chr{$chr}){
		if ($chr_prefix ne ''){
			$chr = $chr_prefix . $chr;
		}
		else{
			$chr = 'chr' . $chr;
			if (!exists $ref_chr{$chr}){
				$chr =~ s/^chr//;
				$chr = 'Chr' . $chr;
			}
		}
	}
	next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
	my $code = $line[2];
	next if ($line[1] eq 'GRCh37') or ($line[1] eq 'GRCh38');
	my $start = $line[3];
	my $end = $line[4];
	my $strand = $line[6];
	my $annot = $line[8];
	if (($flag == 1) and ($annot =~ /ID=gene/)){
		if ($code eq 'pseudogene'){
			$flag = 2;
			next;
		}
		$gene_id = $1 if ($annot =~ /ID=gene:(.+?);/);
		my $gene_name = 'NA';
		$gene_name = $1 if ($annot =~ /Name=(.+?);/);
		$gene_name =~ s/\s+/-/g if ($gene_name =~ /-/);
		if ($gene_name eq 'NA'){
			if ($annot =~ /biotype=(.+?);/){
				$gene_name = $1;
			}
		}
		${$gene{$chr}}{$start} = "$end==$gene_id";
		$ORF{$gene_id} = "$chr=$start=$end=$strand";
		$annot{$gene_id} = "$strand==$gene_name";
		$flag = 0;
	}
	elsif ($flag == 0){
		if ($code eq 'CDS'){
			my $start_ext = $start - 1;
			my $end_ext = $end + 1;
			if (!exists ${$CDS{$gene_id}}{$start_ext}){
				${$CDS{$gene_id}}{$start_ext} = $end_ext;
			}
			else{
				if ($end_ext > ${$CDS{$gene_id}}{$start_ext}){
					${$CDS{$gene_id}}{$start_ext} = $end_ext;
				}
			}
		}
		elsif ($code eq 'exon'){
			my $start_ext = $start - 1;
			my $end_ext = $end + 1;
			if (!exists ${$exon{$gene_id}}{$start_ext}){
				${$exon{$gene_id}}{$start_ext} = $end_ext;
			}
			else{
				if ($end_ext > ${$exon{$gene_id}}{$start_ext}){
					${$exon{$gene_id}}{$start_ext} = $end_ext;
				}
			}
		}
		elsif ($code eq 'five_prime_UTR'){
			if (!exists ${$utr5{$gene_id}}{$start}){
				${$utr5{$gene_id}}{$start} = $end;
			}
			else{
				if ($end > ${$utr5{$gene_id}}{$start}){
					${$utr5{$gene_id}}{$start} = $end;
				}
			}
		}
		elsif ($code eq 'three_prime_UTR'){
			if (!exists ${$utr3{$gene_id}}{$start}){
				${$utr3{$gene_id}}{$start} = $end;
			}
			else{
				if ($end > ${$utr3{$gene_id}}{$start}){
					${$utr3{$gene_id}}{$start} = $end;
				}
			}
		}
	}
}
close (FILE);

print STDERR "Total genes: ", scalar keys %annot, "\n";

my @header;
my $total_hit_genes = 0;
my $total_hit_genes2 = 0;
my $flag_1st = 0;

my $vcf_out1 = "$out_prefix.annot.vcf" if ($target_chr eq 'ALL');
$vcf_out1 = "$temp_dir/$out_prefix.annot.$target_chr.vcf" if ($target_chr ne 'ALL');
my $vcf_out2 = "$out_prefix.AS.annot.vcf" if ($target_chr eq 'ALL');
$vcf_out2 = "$temp_dir/$out_prefix.AS.annot.$target_chr.vcf" if ($target_chr ne 'ALL');
open (OUT1, "> $vcf_out1");
open (OUT2, "> $vcf_out2");
open (FILE, $sv_vcf) or die "$sv_vcf is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#/){
		push @header, $line;
		next;
	}
	my %ANN;
	my %ANN2;
	my @line = split (/\t/, $line);
	my @line2 = split (/\t/, $line);
	my $chr = $line[0];
	next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
	my $pos = $line[1];
	my $pos2 = $pos;
	my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
	$type = 'INS' if ($type eq 'ALU') or ($type eq 'LINE1') or ($type eq 'L1') or ($type eq 'SVA') or ($type eq 'HERVK') or ($type eq 'VEI') or ($type eq 'NUMT');
	$type = 'DUP' if ($type eq 'CNV') and ($line =~ /DUP/);
	$type = 'DEL' if ($type eq 'CNV') and ($line =~ /DEL/);
	next if ($type !~ /DEL|DUP|INS|INV|TRA/);
	my $len = 0;
	my $len2 = 0;
	my $count = 0;
	if (($type ne 'INS') and ($line[8] =~ /VP:VL/)){
		$len2 = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		foreach (@line){
			$count ++;
			next if ($count <= 9);
			next if ($_ =~ /^0\/0/);
			my @info = split (/:/, $_);
			if (($info[2] > 0) and ($info[2] < $pos)){
				$pos = $info[2];
			}
			if ($info[3] > $len){
				$len = $info[3];
			}
		}
	}
	my $end = 0;
	my $end2 = 0;
	$end = $pos + $len - 1;
	$end2 = $pos2 + $len2 - 1;
	if ($type eq 'INS'){
		$len = 1;
		$len2 = 1;
		$end = $pos2;
		$end2 = $pos2;
	}
	$line[7] =~ s/;$// if ($line[7] =~ /;$/);
	$line2[8] .= ':AN';
	foreach my $gstart (sort {$a <=> $b} keys %{$gene{$chr}}){
		my $hit_gene = '';
		my $hit_gene2 = '';
		my ($gend, $gene) = split (/==/, ${$gene{$chr}}{$gstart});
		my ($strand, $gname) = split (/==/, $annot{$gene});
		if (($pos <= $gstart) and ($end >= $gend)){
			$hit_gene = $gene;
		}
		elsif ($strand eq '+'){
			if (($end >= $gstart - $max_flank5) and ($end <= $gend + $max_flank3)){
				$hit_gene = $gene;
			}
			elsif (($pos >= $gstart - $max_flank5) and ($pos <= $gend + $max_flank3)){
				$hit_gene = $gene;
			}
		}
		elsif ($strand eq '-'){
			if (($end >= $gstart - $max_flank3) and ($end <= $gend + $max_flank5)){
				$hit_gene = $gene;
			}
			elsif (($pos >= $gstart - $max_flank3) and ($pos <= $gend + $max_flank5)){
				$hit_gene = $gene;
			}
		}
		if (($pos2 <= $gstart) and ($end2 >= $gend)){
			$hit_gene2 = $gene;
		}
		elsif ($strand eq '+'){
			if (($end2 >= $gstart - $max_flank5) and ($end2 <= $gend + $max_flank3)){
				$hit_gene2 = $gene;
			}
			elsif (($pos2 >= $gstart - $max_flank5) and ($pos2 <= $gend + $max_flank3)){
				$hit_gene2 = $gene;
			}
		}
		elsif ($strand eq '-'){
			if (($end2 >= $gstart - $max_flank3) and ($end2 <= $gend + $max_flank5)){
				$hit_gene2 = $gene;
			}
			elsif (($pos2 >= $gstart - $max_flank3) and ($pos2 <= $gend + $max_flank5)){
				$hit_gene2 = $gene;
			}
		}
		if ($hit_gene ne ''){
			$total_hit_genes ++;
			$ANN{$hit_gene} = $gname;
		}
		if ($hit_gene2 ne ''){
			$total_hit_genes2 ++;
			$ANN2{$hit_gene2} = $gname;
		}
		last if ($gstart - $flank > $end);
	}
	if ($flag_1st == 0){
		my $info_flag = 0;
		my $format_flag = 0;
		foreach my $header (@header){
			if ($header =~ /^##INFO=/){
				$info_flag = 1;
			}
			if ($header =~ /^##FORMAT/){
				$format_flag = 1;
			}
			if (($info_flag == 1) and ($header !~ /^##INFO=/)){
				print OUT1 "##INFO=<ID=SVANN,Number=1,Type=String,Description=\"Affected_Gene,GeneID|GeneName|Region\">\n";
				print OUT1 "$header\n";
				print OUT2 "##INFO=<ID=SVANN,Number=1,Type=String,Description=\"Affected_Gene,GeneID|GeneName|Region\">\n";
				print OUT2 "$header\n";
				$info_flag = 0;
			}
			else{
				if (($format_flag == 1) and ($header !~ /^##FORMAT/)){
					print OUT2 "##FORMAT=<ID=AN,Number=1,Type=String,Description=\"Annotation for each sample, E:exon, AE:all exons, I:intron, 5U:5UTR, 3U:3UTR, 5F:flank5, 3F:flank3\">\n";
					$format_flag = 0;
				}
				print OUT1 "$header\n";
				print OUT2 "$header\n";
			}
		}
		$flag_1st = 1;
	}
	if (scalar keys %ANN > 0){
		my %SVANN;
		$count = 0;
		if ($line[8] =~ /VP:VL/){
			foreach (@line){
				$count ++;
				next if ($count <= 9);
				my @info = split (/:/, $_);
				next if ($info[0] eq '0/0');
				my $spos = $info[2];
				my $slen = $info[3];
				my $send = $spos + $slen - 1 if ($type ne 'INS');
				$send = $spos if ($type eq 'INS');
				my $ann_str = '';
				foreach my $gid (sort keys %ANN){
					my ($gchr, $gstart, $gend, $strand) = split (/=/, $ORF{$gid});
					my $gname = $ANN{$gid};
					$gname = $gid if ($gname =~ /protein_coding|miRNA|tRNA|rRNA|lncRNA/);
					my $hit_region = '';
					next if ($gchr ne $chr) and ("chr$gchr" ne $chr);
					if ($send < $gstart){
						if ($strand eq '+'){
							if ($send >= $gstart - $max_flank5){
								$hit_region = '5F';
							}
						}
						else{
							if ($send >= $gstart - $max_flank3){
								$hit_region = '3F';
							}
						}
					}
					elsif ($spos > $gend){
						if ($strand eq '+'){
							if ($send >= $gstart - $max_flank3){
								$hit_region = '3F';
							}
						}
						else{
							if ($send >= $gstart - $max_flank5){
								$hit_region = '5F';
							}
						}
					}
					elsif (($spos <= $gstart) and ($send >= $gend)){
						$hit_region = 'AE';
					}
					else{
						if (exists $CDS{$gid}){
							foreach my $cstart (sort {$a <=> $b} keys %{$CDS{$gid}}){
								my $cend = ${$CDS{$gid}}{$cstart};
								if ((($spos <= $cstart) and ($send >= $cend)) or (($spos >= $cstart) and ($spos <= $cend)) or (($send >= $cstart) and ($send <= $cend))){
									$hit_region = 'E';
									last;
								}
							}
						}
						elsif (exists $exon{$gid}){
							foreach my $estart (sort {$a <=> $b} keys %{$exon{$gid}}){
								my $eend = ${$exon{$gid}}{$estart};
								if ((($spos <= $estart) and ($send >= $eend)) or (($spos >= $estart) and ($spos <= $eend)) or (($send >= $estart) and ($send <= $eend))){
									$hit_region = 'E';
									last;
								}
							}
						}
						if ($hit_region eq ''){
							if (exists $utr5{$gid}){
								foreach my $ustart (sort {$a <=> $b} keys %{$utr5{$gid}}){
									my $uend = ${$utr5{$gid}}{$ustart};
									if ((($spos <= $ustart) and ($send >= $uend)) or (($spos >= $ustart) and ($spos <= $uend)) or (($send >= $ustart) and ($send <= $uend))){
										$hit_region = '5U';
										last;
									}
								}
							}
							if (($hit_region eq '') and (exists $utr3{$gid})){
								foreach my $ustart (sort {$a <=> $b} keys %{$utr3{$gid}}){
									my $uend = ${$utr3{$gid}}{$ustart};
									if ((($spos <= $ustart) and ($send >= $uend)) or (($spos >= $ustart) and ($spos <= $uend)) or (($send >= $ustart) and ($send <= $uend))){
										$hit_region = '3U';
										last;
									}
								}
							}
						}
						if ($hit_region eq ''){
							if ((($spos >= $gstart) and ($spos <= $gend)) or (($send >= $gstart) and ($send <= $gend))){
								$hit_region = 'I';
							}
						}
					}
					$ann_str .= "$gname-$hit_region," if ($hit_region ne '');
	#print STDERR "$count\t$spos\t$send\t$gstart\t$gend\t$gname\t$hit_region\n";
				}
				$ann_str =~ s/,$// if ($ann_str =~ /,$/);
				$line2[$count - 1] .= ":$ann_str" if ($ann_str ne '');
			}
		}
		foreach my $gid (sort keys %ANN2){
print STDERR "Missing ORF $gid\n" if (!exists $ORF{$gid});
			my ($gchr, $gstart, $gend, $strand) = split (/=/, $ORF{$gid});
			my $gname = $ANN2{$gid};
			my $hit_region = '';
			next if ($gchr ne $chr) and ("chr$gchr ne $chr");
			if ($end2 < $gstart){
				if ($strand eq '+'){
					if ($end2 >= $gstart - $max_flank5){
						$hit_region = '5F';
						${$SVANN{$gid}}{6} = $gname if ($end2 < $gstart - $core_flank5);
						${$SVANN{$gid}}{5} = $gname if ($end2 >= $gstart - $core_flank5);
					}
				}
				else{
					if ($end2 >= $gstart - $max_flank3){
						$hit_region = '3F';
						${$SVANN{$gid}}{6.1} = $gname if ($end2 < $gstart - $core_flank3);
						${$SVANN{$gid}}{5.1} = $gname if ($end2 >= $gstart - $core_flank3);
					}
				}
			}
			elsif ($pos2 > $gend){
				if ($strand eq '+'){
					if ($pos2 <= $gend + $max_flank3){
						$hit_region = '3F';
						${$SVANN{$gid}}{6.1} = $gname if ($pos2 > $gend + $core_flank3);
						${$SVANN{$gid}}{5.1} = $gname if ($pos2 <= $gend + $core_flank3);
					}
				}
				else{
					if ($pos2 <= $gend + $max_flank5){
						$hit_region = '5F';
						${$SVANN{$gid}}{6} = $gname if ($pos2 > $gend + $core_flank5);
						${$SVANN{$gid}}{5} = $gname if ($pos2 <= $gend + $core_flank5);
					}
				}
			}
			elsif (($pos2 <= $gstart) and ($end2 >= $gend)){
				$hit_region = 'AE';
				${$SVANN{$gid}}{1} = $gname;
			}
			else{
				if (exists $CDS{$gid}){
					foreach my $cstart (sort {$a <=> $b} keys %{$CDS{$gid}}){
						my $cend = ${$CDS{$gid}}{$cstart};
						if ((($pos2 <= $cstart) and ($end2 >= $cend)) or (($pos2 >= $cstart) and ($pos2 <= $cend)) or (($end2 >= $cstart) and ($end2 <= $cend))){
							$hit_region = 'E';
							${$SVANN{$gid}}{2} = $gname;
							last;
						}
					}
				}
				elsif (exists $exon{$gid}){
					foreach my $estart (sort {$a <=> $b} keys %{$exon{$gid}}){
						my $eend = ${$exon{$gid}}{$estart};
						if ((($pos2 <= $estart) and ($end2 >= $eend)) or (($pos2 >= $estart) and ($pos2 <= $eend)) or (($end2 >= $estart) and ($end2 <= $eend))){
							$hit_region = 'E';
							${$SVANN{$gid}}{2.1} = $gname;
							last;
						}
					}
				}
				if ($hit_region eq ''){
					if (exists $utr5{$gid}){
						foreach my $ustart (sort {$a <=> $b} keys %{$utr5{$gid}}){
							my $uend = ${$utr5{$gid}}{$ustart};
							if ((($pos2 <= $ustart) and ($end2 >= $uend)) or (($pos2 >= $ustart) and ($pos2 <= $uend)) or (($end2 >= $ustart) and ($end2 <= $uend))){
								$hit_region = '5U';
								${$SVANN{$gid}}{3} = $gname;
								last;
							}
						}
					}
					if (($hit_region eq '') and (exists $utr3{$gid})){
						foreach my $ustart (sort {$a <=> $b} keys %{$utr3{$gid}}){
							my $uend = ${$utr3{$gid}}{$ustart};
							if ((($pos2 <= $ustart) and ($end2 >= $uend)) or (($pos2 >= $ustart) and ($pos2 <= $uend)) or (($end2 >= $ustart) and ($end2 <= $uend))){
								$hit_region = '3U';
								${$SVANN{$gid}}{3.1} = $gname;
								last;
							}
						}
					}
				}
				if ($hit_region eq ''){
					if ((($pos2 >= $gstart) and ($pos2 <= $gend)) or (($end2 >= $gstart) and ($end2 <= $gend))){
						$hit_region = 'I';
						${$SVANN{$gid}}{4} = $gname;
					}
				}
			}
		}
		my $svann_str = '';
		foreach my $gid (sort keys %SVANN){
			foreach my $grade (sort {$a <=> $b} keys %{$SVANN{$gid}}){
				my $region = '';
				if ($grade == 1){
					$region = 'All_exons';
				}
				elsif ($grade == 2){
					$region = 'CDS';
				}
				elsif ($grade == 2.1){
					$region = 'exon';
				}
				elsif ($grade == 3){
					$region = '5UTR';
				}
				elsif ($grade == 3.1){
					$region = '3UTR';
				}
				elsif ($grade == 4){
					$region = 'intron';
				}
				elsif ($grade == 5){
					$region = $flank5_prox_str;
				}
				elsif ($grade == 5.1){
					$region = $flank3_prox_str;
				}
				elsif ($grade == 6){
					$region = $flank5_dist_str;
				}
				elsif ($grade == 6.1){
					$region = $flank3_dist_str;
				}
				my $name = ${$SVANN{$gid}}{$grade};
				$svann_str .= "$gid|$name|$region,";
				last;
			}
		}
		$svann_str =~ s/,$// if ($svann_str =~ /,$/);
		$line[7] .= ';SVANN=' . $svann_str if ($svann_str ne '');
		$line2[7] .= ';SVANN=' . $svann_str if ($svann_str ne '');
	}
	print OUT1 join ("\t", @line), "\n";
	print OUT2 join ("\t", @line2), "\n";
}
close (FILE);
close (OUT1);
close (OUT2);

print STDERR "Total hit genes: $total_hit_genes\n";
