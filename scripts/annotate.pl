#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use threads;

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

my $cores = 1;

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
    'build|b=s' => \$build,
    'num_threads|n=i' => \$cores,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  annotate.pl -v <vcf file> -p <prefix of output vcf> -n <num of threads> 
  (-nh 1 -r <reference gff3 annotation file from Ensemble> if sample is a non-human sepecies)
  output vcf files: ${out_prefix}.annot.vcf and ${out_prefix}.AS.annot.vcf

  Options:
   --sv or -v <STR>         vcf file of SVs
   --ref or -r <STR>        ref gff3 annotation file from Ensembl 
                            [default for human: Data/Homo_sapiens.GRCh37.87.gff3.gz or Data/Homo_sapiens.GRCh38.104.gff3.gz]
   --flank5 or -f5 <INT>    maximum size (bp) of 5' distal flanking gene region, that overlaps SV [default: 50000]
   --flank3 or -f3 <INT>    maximum size (bp) of 3' distal flanking gene region, that overlaps SV [default: 50000]
   --core5 or -c5 <INT>     maximum size (bp) of 5' proximal flanking gene region, that overlaps SV [default: 5000]
   --core3 or -c3 <INT>     maximum size (bp) of 3' proximal flanking gene region, that overlaps SV [default: 5000]
   --non_human or -nh <INT>  samples are non-human species (0: human, 1: non-human) [default: 0]
   --build <STR>            reference build [GRCh37, GRCh38, T2T-CHM13] (37, 38, or T2T) [default: 37]
   --target or -t <STR>     target chromosome [default: ALL]
   --prefix or -p <STR>     outpout prefix
   --num_threads or -n <INT> number of threads to be used [default: 1]
   --help or -h             output help message
   
=cut

die "-v option shold be specified:\n" if ($sv_vcf eq '');

system ("mkdir $temp_dir") if (!-d $temp_dir);

my $arg = "-v $sv_vcf -p $out_prefix ";
if ($ref_gff ne ''){
	$arg .= "-r $ref_gff ";
}
if ($non_human != 0){
	$arg .= "-nh $non_human ";
}
if ($build ne '37'){
	$arg .= "-b $build ";
}
if ($target_chr ne 'ALL'){
	$arg .= "-t $target_chr ";
}
if ($core_flank5 != 5000){
	$arg .= "-c5 $core_flank5 ";
}
if ($core_flank3 != 5000){
	$arg .= "-c3 $core_flank3 ";
}
if ($max_flank5 != 50000){
	$arg .= "-f5 $max_flank5 ";
}
if ($max_flank3 != 50000){
	$arg .= "-f3 $max_flank3 ";
}

if (($cores == 1) or ($target_chr ne 'ALL')){
	system ("$Bin/annotate_SV_vcf.pl $arg");
	if ($target_chr ne 'ALL'){
		system ("mv $temp_dir/$out_prefix.annot.$target_chr.vcf $out_prefix.annot.$target_chr.vcf") if (-f "$temp_dir/$out_prefix.annot.$target_chr.vcf");
		system ("mv $temp_dir/$out_prefix.AS.annot.$target_chr.vcf $out_prefix.AS.annot.$target_chr.vcf") if (-f "$temp_dir/$out_prefix.AS.annot.$target_chr.vcf");
	}
}
else{
	my %ref_chr;
	my @ref_chr;
	open (FILE, $sv_vcf) or die "$sv_vcf is not found:$!\n" if ($sv_vcf !~ /\.gz$/);
	open (FILE, "gzip -dc $sv_vcf |") or die "$sv_vcf is not found:$!\n" if ($sv_vcf =~ /\.gz$/);
	while (my $line = <FILE>){
		chomp $line;
		next if ($line =~ /^#/);
		my ($chr) = split (/\t/, $line);
		if (!exists $ref_chr{$chr}){
			push @ref_chr, $chr;
		}
		$ref_chr{$chr} = 1;
	}
	close (FILE);

	my @jobs;
	my %out_files1;
	my %out_files2;
	my $count = 0;
	foreach my $chr (@ref_chr){
		$count ++;
		my ($thread_t) = threads->new(\&annotate, $chr);
		push @jobs, $thread_t;
		if ($count == $cores){
			foreach (@jobs){
				my ($chr2, $out1_chr, $out2_chr) = $_->join;
				$out_files1{$chr2} = $out1_chr;
				$out_files2{$chr2} = $out2_chr if (-f $out2_chr);
			}
			$count = 0;
			undef @jobs;
		}
	}
	if (@jobs > 0){
		foreach (@jobs){
			my ($chr2, $out1_chr, $out2_chr) = $_->join;
			$out_files1{$chr2} = $out1_chr;
			$out_files2{$chr2} = $out2_chr if (-f $out2_chr);
		}
		$count = 0;
		undef @jobs;
	}
	my $out_vcf1 = "$out_prefix.annot.vcf";
	my $count2 = 0;
	open (OUT1, "> $out_vcf1");
	foreach my $chr (@ref_chr){
		$count2 ++;
		next if (!exists $out_files1{$chr});
		my $out1_chr = $out_files1{$chr};
		open (FILE, $out1_chr) or die "$out1_chr is not found:$!\n";
		while (my $line = <FILE>){
			chomp $line;
			if ($line =~ /^#/){
				if ($count2 == 1){
					print OUT1 "$line\n";
				}
				next;
			}
			print OUT1 $line, "\n";
		}
		close (FILE);
	} 
	close (OUT1);
	if (scalar keys %out_files2 > 0){
		my $out_vcf2 = "$out_prefix.AS.annot.vcf";
		open (OUT2, "> $out_vcf2");
		$count2 = 0;
		foreach my $chr (@ref_chr){
			$count2 ++;
			my $out2_chr = $out_files2{$chr};
			open (FILE, $out2_chr) or die "$out2_chr is not found:$!\n";
			while (my $line = <FILE>){
				chomp $line;
				if ($line =~ /^#/){
					if ($count2 == 1){
						print OUT2 "$line\n";
					}
					next;
				}
				print OUT2 $line, "\n";
			}
			close (FILE);
		}
		close (OUT2);
	}
}

print STDERR "SV vcf annotation finished:\n";


sub annotate{
	my ($chr) = @_;
	my $arg_chr = "$arg -t $chr";
	system ("$Bin/annotate_SV_vcf.pl $arg_chr");

	my $out1 = "$temp_dir/$out_prefix.annot.$chr.vcf";
	my $out2 = "$temp_dir/$out_prefix.AS.annot.$chr.vcf";

	threads->yield();
	sleep 1;
	return ($chr, $out1, $out2);
}

