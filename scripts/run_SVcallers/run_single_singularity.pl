#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

my $sif_file = '';
my $bam = '';
my $sample_dir_name = '';
my $config = '';
my $temp_dir = '';
my $bind_dir = '';
my $no_home = 0;
my $help;
 
GetOptions(
		'sif|s=s' => \$sif_file,
    'bam|b=s' => \$bam,
    'sample_dir|sd=s' => \$sample_dir_name,
    'config|c=s' => \$config,
    'temp_dir|td=s' => \$temp_dir,
    'bind_dir|bd=s' => \$bind_dir,
    'no_home|noh' => \$no_home,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_single.pl -b <input bam> -c <config file> (-sn <sample name>)

  Options:
   --sif or -s <STR>        absolute path of mopline sif file generated with MOPline-Definition.txt [mandatory]
   --bam or -b <STR>        input bam file [mandatory]
   --sample_dir or -sd <STR> sample directory name where tool directories are created (if not specified, sample directory is not created and tools directories are created in the working directory)
   --config or -c <STR>     config file [mandatory]
   --temp_dir or -td <STR>  absolute path of tmp directory on the host [mandatory]
   --bind_dir or -bd <STR>  comma-separated list of absolute path (except for the working directory) on the host to be added to singularity container (optional)
   --no_home or -noh <BOOLEAN>  do not add $HOME on the host to singularity container [default: false]
   --help or -h             output help message
   
=cut

die "sif file is not specified: or does not exist\n" if ($sif_file eq '') or (!-f $sif_file);
die "bam is not specified:\n" if ($bam eq '');
die "config file is not specified:\n" if ($config eq '');
die "tmp directory not specified:\n" if ($temp_dir eq '');

my $bam_base = basename ($bam);
my $sample_name = $sample_dir_name if ($sample_dir_name ne '');
$sample_name = $1 if ($sample_dir_name eq '') and ($bam_base =~ /(.+?)\./);

my $abs_bam = File::Spec->rel2abs($bam);
$abs_bam = readlink ($abs_bam) if (-l $abs_bam);
$bam = $abs_bam;

my $ref = '';
my $target_chr = 'ALL';
my $non_human = 0;
my $run_svcaller_dir = '/opt/local/tools/MOPline/scripts/run_SVcallers';

my $bind_dir2 = $temp_dir;
my $bind_dir3 = `pwd`;
chomp $bind_dir3;
$bind_dir2 .= ",$bind_dir3";
$bind_dir2 .= ",$bind_dir" if ($bind_dir ne '');

my @tools;

my %tool_opt;
my %tool_script;
my $tool = '';

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
			if ($opt1 eq 'ref'){
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

my $abs_ref = File::Spec->rel2abs($ref);
$abs_ref = readlink ($abs_ref) if (-l $abs_ref);
$ref = $abs_ref;

my $bam_dir = '';
my $ref_dir = '';
$bam_dir = $1 if ($bam =~ /(.+)\//);
$ref_dir = $1 if ($ref =~ /(.+)\//);

$bind_dir2 .= ",$bam_dir" if ($bam_dir ne '') and ($bam_dir ne $bind_dir3);
$bind_dir2 .= ",$ref_dir" if ($ref_dir ne '') and ($ref_dir ne $bind_dir3) and ($ref_dir ne $bam_dir);

if ($sample_dir_name ne ''){
	system ("mkdir $sample_dir_name") if (!-d $sample_dir_name);
	chdir $sample_dir_name;
}
print STDERR "Sample: $sample_name\n";
foreach my $tool_name (@tools){
	next if ($tool_name eq 'general');
	system ("mkdir $tool_name") if (!-d $tool_name);
	chdir $tool_name;
	my $run_script = $tool_script{$tool_name};
	my $opt_str = '';
	foreach my $opt (keys %{$tool_opt{$tool_name}}){
		$opt_str .= "${$tool_opt{$tool_name}}{$opt} ";
	}
	$opt_str .= "-nh $non_human ";
	$opt_str =~ s/\s$//;
	my $bind_dir = $bind_dir2;
	if ($tool_name eq 'CNVnator'){
		$opt_str = "-b $bam -p $sample_name " . $opt_str ;
		my $ref_dir2 = '';
		$ref_dir2 = $1 if ($opt_str =~ /-r\s+(\S+)/);
		$bind_dir .= ",$ref_dir2" if ($ref_dir2 ne '');
	}
	else{
		$opt_str = "-b $bam -p $sample_name -r $ref " . $opt_str;
	}
	if ($target_chr ne 'ALL'){
		$opt_str .= " -c $target_chr";
	}
	print STDERR "$tool_name\n";
	my $command = "singularity exec --bind $bind_dir $sif_file $run_script $opt_str";
	$command = "singularity exec --bind $bind_dir --no-home $sif_file $run_script $opt_str" if ($no_home == 1);
	if ($tool_name =~ /MELT/){
		$command = "$run_script $opt_str";
	}
	my $comman_log = "$sample_name.command.log";
	my $error_log = "$sample_name.error.log";
	open (OUT, "> $comman_log");
	print OUT "$tool_name run: $command\n";
	close (OUT);
	my $run = `$command 2>$error_log`;
	chdir '..';
}
chdir '..';
