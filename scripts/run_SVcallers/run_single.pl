#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

my $bam = '';
my $sample_dir_name = '';
my $config = '';
my $help;
 
GetOptions(
    'bam|b=s' => \$bam,
    'sample_dir|sd=s' => \$sample_dir_name,
    'config|c=s' => \$config,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_single.pl -b <input bam> -c <config file> (-sn <sample name>)

  Options:
   --bam or -b <STR>        input bam file [mandatory]
   --sample_dir or -sd <STR> sample directory name where tool directories are created (if not specified, sample directory is not created and tools directories are created in the working directory)
   --config or -c <STR>     config file [mandatory]
   --help or -h             output help message
   
=cut

die "bam is not specified:\n" if ($bam eq '');
die "config file is not specified:\n" if ($config eq '');

my $script_dir = '';
$script_dir = $1 if ($0 =~ /(.+)\//);

my $bam_base = basename ($bam);
my $sample_name = $sample_dir_name if ($sample_dir_name ne '');
$sample_name = $1 if ($sample_dir_name eq '') and ($bam_base =~ /(.+?)\./);

my $abs_bam = File::Spec->rel2abs($bam);
$bam = $abs_bam;

my $ref = '';
my $target_chr = 'ALL';
my $non_human = 0;
my $run_svcaller_dir = $script_dir;

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
		$opt1 =~ s/[\'\"]// if ($opt1 =~ /[\'\"']/);
		$arg =~ s/[\'\"]// if ($arg =~ /[\'\"']/);
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
				$run_script =~ s/[\'\"]// if ($run_script =~ /[\'\"']/);
				$opt2 =~ s/[\'\"]// if ($opt2 =~ /[\'\"']/);
				${$tool_opt{$tool}}{$opt1} = "$opt2 $arg";
				$tool_script{$tool} = "$run_svcaller_dir/$run_script" if ($run_svcaller_dir ne '');
			}
		}
	}
}
close (FILE);

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
	$opt_str .= "-nh $non_human " if ($tool_name =~ /CNVnator|inGAP|MELT|Wham|DELLY|Lumpy|SoftSV|Manta|MATCHCLIP|INSurVeylor/);
	$opt_str =~ s/\s$//;
	$opt_str = "-b $bam -p $sample_name " . $opt_str if ($tool_name eq 'CNVnator');
	$opt_str = "-b $bam -p $sample_name -r $ref " . $opt_str if ($tool_name ne 'CNVnator');
	if ($target_chr ne 'ALL'){
		$opt_str .= " -c $target_chr";
	}
	print STDERR "$tool_name\n";
	my $command = "$run_script $opt_str";
	my $comman_log = "$sample_name.command.log";
	my $error_log = "$sample_name.error.log";
	open (OUT, "> $comman_log");
	print OUT "$tool_name run: $command\n";
	close (OUT);
	system ("$command 2>$error_log");
	chdir '..';
}
chdir '..';
