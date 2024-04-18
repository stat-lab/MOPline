#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use FindBin qw($Bin);

my $sif_file = '';
my $ref = '';
my $sample_list = '';
my $config = '';
my $melt_config = '';
my $temp_dir = '';
my $bind_dir = '';
my $no_home = 0;
my $add_MCMQ_tag = 0;
my $help;
 
GetOptions(
		'sif_file|sf=s' => \$sif_file,
    'ref|r=s' => \$ref,
    'sample_list|sl=s' => \$sample_list,
    'config|c=s' => \$config,
    'config2|c2=s' => \$melt_config,
    'temp_dir|td=s' => \$temp_dir,
    'bind_dir|bd=s' => \$bind_dir,
    'no_home|noh' => \$no_home,
    'add_MC|am' => \$add_MCMQ_tag,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_SVcallers_batch.pl -sl <sample list file> -r <reference fasta> -c <config file> -sf <mopline.sif file> -td <temp directory> -c2 <config file for MELT and/or INSurVeyor> (-am)

  Options:
   --sample_list or -sl <STR> sample list file, which must indicate sample name and the path of the corresponding bam/cram file seprated with tab in each line [mandatory]
   --sif_file or -sf <STR>  absolute path of mopline.sif singularity image file [mandatory if the -c option is specified]
   --ref or -r <STR>        reference fasta
   --config or -c <STR>     config file for singularity run
   --config2 or -c2 <STR>   config file for MELT and/or INSurVeyor
   --temp_dir or -td <STR>  absolute path of tmp directory on the host [mandatory]
   --bind_dir or -bd <STR>  comma-separated list of absolute path (except for the working directory) on the host to be added to singularity container (optional)
   --no_home or -noh <BOOLEAN>  do not add $HOME on the host to singularity container [default: false]
   --add_MC or -am <BOOLEAN> add MC/MQ tag in bam file [default: false]
   --help or -h             output help message
   
=cut

die "sif file is not specified: or does not exist\n" if ($config ne '') and (($sif_file eq '') or (!-f $sif_file));
die "sample list file is not specified:\n" if ($sample_list eq '') or (!-f $sample_list);
die "config and/or config2 file is not specified:\n" if ($config eq '') and ($melt_config eq '');
die "tmp directory not specified:\n" if ($temp_dir eq '');

my %sample_bam;
my @samples;

$config = File::Spec->rel2abs($config);
$melt_config = File::Spec->rel2abs($melt_config);

my $target_chr = 'ALL';
my $non_human = 0;
my $run_svcaller_dir = '/opt/local/tools/MOPline/scripts/run_SVcallers';
my $run_svcaller_dir2 = '';

my $bind_dir2 = $temp_dir;
my $bind_dir3 = `pwd`;
chomp $bind_dir3;
$bind_dir2 .= ",$bind_dir3";
$bind_dir2 .= ",$bind_dir" if ($bind_dir ne '');

open (FILE, $sample_list) or die "$sample_list is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my ($sample_name, $bam) = split (/\s+/, $line);
	my $abs_bam = File::Spec->rel2abs($bam);
	$bam = $abs_bam;
	my $sample_dir_name = $sample_name;
	system ("mkdir $sample_dir_name") if (!-d $sample_dir_name);

	$bind_dir2 .= ",$bind_dir3/$sample_dir_name";

	chdir $sample_dir_name;

	if ((!-f "$sample_name.bam.bai") or (-z "$sample_name.bam")){
		if ($bam =~ /\.bam$/){
			if ($add_MCMQ_tag == 0){
				system ("ln -s $bam $sample_name.bam");
				if (-f "$bam.bai"){
					system ("ln -s $bam.bai $sample_name.bam.bai");
				}
				else{
					system ("samtools index $sample_name.bam");
				}
				if ($abs_bam =~ /(.+)\//){
					my $bam_dir = $1;
					$bind_dir2 .= ",$bam_dir";
				}
			}
			else{
				if ($bam ne "$sample_name.bam"){
					system ("samtools view -hb $bam | samtools sort -n - | samtools fixmate -@ 2 -m - $sample_name.fm.bam");
					die "$sample_name.fm.bam not created:\n" if (!-f "$sample_name.fm.bam") or (-z "$sample_name.fm.bam");
					system ("samtools sort -@ 4 $sample_name.fm.bam | samtools markdup -r -u - $sample_name.fm.rmdup.bam");
					system ("rm -f $sample_name.fm.bam");
					die "$sample_name.fm.rmdup.bam not created:\n" if (!-f "$sample_name.fm.rmdup.bam") or (-z "$sample_name.fm.rmdup.bam");
					my $header_sam = "$sample_name.header.sam";
					my @header;
					system ("samtools view -H $sample_name.fm.rmdup.bam > $header_sam");
					open (FILE, $header_sam) or die "$header_sam is not found:$!\n";
					while (my $line = <FILE>){
					    chomp $line;
					    next if ($line =~ /^\@PG/);
					    push @header, $line;
					}
					close (FILE);

					open (OUT, "> $header_sam");
					foreach (@header){
					    print OUT "$_\n";
					}
					close (OUT);
					system ("samtools reheader $header_sam $sample_name.fm.rmdup.bam > $sample_name.bam");
					system ("samtools index $sample_name.bam");
					system ("rm -f $sample_name.fm.rmdup.bam");
					$bam = "$bind_dir3/$sample_dir_name/$sample_name.bam";
				}
				else{
					die "$sample_name.bam is already present in $sample_dir_name\n";
				}
			}
		}
		elsif ($bam =~ /\.cram$/){
			die "Ref fasta not present:\n" if (!-f $ref);
			if ($add_MCMQ_tag == 0){
				system ("samtools view -@ 2 -T $ref -hb $bam -o $sample_name.bam");
				system ("samtools index $sample_name.bam");
				$bam = "$bind_dir3/$sample_dir_name/$sample_name.bam";
			}
			else{
				system ("samtools view -@ 2 -T $ref -hb $bam | samtools sort -n -o bam - | samtools fixmate -@ 2 -m - $sample_name.fm bam");
				die "$sample_name.fm.bam not created:\n" if (!-f "$sample_name.fm.bam") or (-z "$sample_name.fm.bam");
				system ("samtools sort -@ 4 $sample_name.fm.bam | samtools markdup -r -u - $sample_name.fm.rmdup.bam");
				system ("rm -f $sample_name.fm.bam");
				die "$sample_name.fm.rmdup.bam not created:\n" if (!-f "$sample_name.fm.rmdup.bam") or (-z "$sample_name.fm.rmdup.bam");
				my $header_sam = "$sample_name.header.sam";
				my @header;
				system ("samtools view -H $sample_name.fm.rmdup.bam > $header_sam");
				open (FILE, $header_sam) or die "$header_sam is not found:$!\n";
				while (my $line = <FILE>){
				    chomp $line;
				    next if ($line =~ /^\@PG/);
				    push @header, $line;
				}
				close (FILE);

				open (OUT, "> $header_sam");
				foreach (@header){
				    print OUT "$_\n";
				}
				close (OUT);
				system ("samtools reheader $header_sam $sample_name.fm.rmdup.bam > $sample_name.bam");
				system ("samtools index $sample_name.bam");
				system ("rm -f $sample_name.fm.rmdup.bam");
				$bam = "$bind_dir3/$sample_dir_name/$sample_name.bam";
			}
		}
	}
	$sample_bam{$sample_name} = $bam;
	push @samples, $sample_name;
	chdir '..';
}
close (FILE);

my @tools;

my %tool_opt;
my %tool_script;
my $tool = '';

if ($config ne ''){
	open (FILE, $config) or die "$config is not found: $!\n";
	while (my $line = <FILE>){
		chomp $line;
		if ($line =~ /^#/){
			$tool = '';
		}
		next if ($line =~ /^#|^$/);
		if ($line =~ /^<(.+)>/){
			$tool = $1;
			push @tools, $tool;
			next;
		}
		if ($line =~ /^(\S+)\s*=\s*(\S+)\s*;/){
			my $opt1 = $1;
			my $arg = $2;
			$opt1 =~ s/[\'\"]//g if ($opt1 =~ /[\'\"']/);
			$arg =~ s/[\'\"]//g if ($arg =~ /[\'\"']/);
			if ($tool eq 'general'){
				if (($opt1 eq 'ref') and ($ref eq '')){
					$ref = $arg;
				}
				elsif ($opt1 eq 'target'){
					$target_chr = $arg;
				}
				elsif ($opt1 eq 'non_human'){
					$non_human = $arg;
				}
				elsif ($opt1 eq 'run_svcaller_dir'){
					$run_svcaller_dir = $arg;
				}
				if (($tool eq 'CNVnator') and ($opt1 eq 'ref_dir')){
					$ref = $arg if ($arg ne '');
				}
			}
			elsif ($tool ne ''){
				if ($line =~ /\[\s*(\S+)\s+(\S+)\s*\]/){
					my $run_script = $1;
					my $opt2 = $2;
					$run_script =~ s/[\'\"]//g if ($run_script =~ /[\'\"']/);
					$opt2 =~ s/[\'\"]//g if ($opt2 =~ /[\'\"']/);
					${$tool_opt{$tool}}{$opt1} = "$opt2 $arg";
					$tool_script{$tool} = "$run_svcaller_dir/$run_script" if ($run_svcaller_dir ne '');
				}
			}
		}
	}
	close (FILE);
}

if ($melt_config ne ''){
	open (FILE, $melt_config) or die "$melt_config is not found: $!\n";
	while (my $line = <FILE>){
		chomp $line;
		if ($line =~ /^#/){
			$tool = '';
		}
		next if ($line =~ /^#|^$/);
		if ($line =~ /^<(.+)>/){
			$tool = $1;
			push @tools, $tool if ($tool eq 'MELT') or ($tool eq 'INSurVeyor');
			next;
		}
		if ($line =~ /^(\S+)\s*=\s*(\S+)\s*;/){
			my $opt1 = $1;
			my $arg = $2;
			$opt1 =~ s/[\'\"]// if ($opt1 =~ /[\'\"']/);
			$arg =~ s/[\'\"]// if ($arg =~ /[\'\"']/);
			if ($tool eq 'general'){
				if ($opt1 eq 'run_svcaller_dir'){
					$run_svcaller_dir2 = $arg;
				}
			}
			elsif ($tool ne ''){
				if ($line =~ /\[\s*(\S+)\s+(\S+)\s*\]/){
					my $run_script = $1;
					my $opt2 = $2;
					$run_script =~ s/[\'\"]// if ($run_script =~ /[\'\"']/);
					$opt2 =~ s/[\'\"]// if ($opt2 =~ /[\'\"']/);
					${$tool_opt{$tool}}{$opt1} = "$opt2 $arg";
					$run_svcaller_dir2 = $Bin if ($run_svcaller_dir2 eq '');
					$tool_script{$tool} = "$run_svcaller_dir2/$run_script" if ($run_svcaller_dir2 ne '');
				}
			}
		}
	}
	close (FILE);
}

foreach my $sample_name (@samples){
	print STDERR "Sample: $sample_name\n";
	chdir $sample_name;
	foreach my $tool_name (@tools){
		next if ($tool_name eq 'general');
		system ("mkdir $tool_name") if (!-d $tool_name);
		next if (-f "$tool_name/$tool_name.$sample_name.vcf") and (!-z "$tool_name/$tool_name.$sample_name.vcf");
		chdir $tool_name;
		my $run_script = $tool_script{$tool_name};
		my $bam = $sample_bam{$sample_name};
		my $opt_str = '';
		foreach my $opt (keys %{$tool_opt{$tool_name}}){
			$opt_str .= "${$tool_opt{$tool_name}}{$opt} ";
		}
		$opt_str .= "-nh $non_human ";
		$opt_str =~ s/\s$//;
		$opt_str = "-b $bam -p $sample_name " . $opt_str if ($tool_name eq 'CNVnator');
		$opt_str = "-b $bam -p $sample_name -r $ref " . $opt_str if ($tool_name ne 'CNVnator');
		if ($target_chr ne 'ALL'){
			$opt_str .= " -c $target_chr";
		}
		print STDERR "$tool_name\n";
		my $command = "singularity exec --bind $bind_dir2 $sif_file $run_script $opt_str";
		$command = "singularity exec --bind $bind_dir2 --no-home $sif_file $run_script $opt_str" if ($no_home == 1);
		if ($tool_name =~ /MELT|INSurVeyor/){
			$command = "$run_script $opt_str";
		}
		my $comman_log = "$sample_name.command.log";
		my $error_log = "$sample_name.error.log";
		open (OUT, "> $comman_log");
		print OUT "$tool_name run: $command\n";
		close (OUT);
		my $run = `$command 2>$error_log`;
		if ($tool_name eq 'inGAP'){
			if ((-f "inGAP.$sample_name.vcf") and (!-z "inGAP.$sample_name.vcf")){
				system ("rm -f *.out");
				system ("rm -f *.fa");
			}
		}
		elsif ($tool_name eq 'MELT'){
			if ((-f "MELT.$sample_name.vcf") and (!-z "MELT.$sample_name.vcf")){
				system ("rm -f *.bam");
				system ("rm -f *.bam.fq");
				system ("rm -f *.bam.disc");
			}
		}
		chdir '..';
	}
	chdir '..';
}
