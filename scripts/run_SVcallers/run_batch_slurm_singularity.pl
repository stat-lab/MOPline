#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $sif_file = '';
my $bam_list = '';
my $config = '';
my $temp_dir = '';
my $bind_dir = '';
my $no_home = 0;
my $account = '';
my $quos = '';
my $partition = '';
my $memory = 30000;
my $time = 0;
my $node_list = '';
my $help;
 
GetOptions(
	  'sif|s=s' => \$sif_file,
    'bam_list|b=s' => \$bam_list,
    'config|c=s' => \$config,
    'temp_dir|td=s' => \$temp_dir,
    'bind_dir|bd=s' => \$bind_dir,
    'no_home|noh' => \$no_home,
    'account|a=s' => \$account,
    'quos|q=s' => \$quos,
    'partition|p=s' => \$partition,
    'mem|m=i' => \$memory,
    'time|t=i' => \$time,
    'node|w=s' => \$node_list,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_batch_slurm.pl -b <bam_list> -c <config_file> -a <account> -m <memory(MB)> (-q <quos> -p <partition>)

  Options:
   --sif or -s <STR>        absolute path of mopline sif file generated with MOPline-Definition.txt [mandatory]
   --bam_list or -b <STR>   bam list file [mandatory]
   --config or -c <STR>     config file [mandatory]
   --temp_dir or -td <STR>  tmp directory on the host [mandatory]
   --bind_dir or -bd <STR>  comma-separated list of path (except for the working directory: sample directory) on the host to be added to singularity container (optional)
   --no_home or -noh <BOOLEAN>  do not add $HOME on the host to singularity container [default: false]
   --account or -a <STR>    account name for slurm [mandatory]
   --quos or -q             quality of service level if specified
   --partition or -p <INT>  partition name if specified
   --memory or -m <INT>     minimum amount of real memory in MB [default: 30000]
   --time or -t <INT>  		  time limit in min (optional)
   --node or -w <STR>       node list (optional)
   --help or -h             output help message
   
=cut

die "sif file is not specified: or does not exist\n" if ($sif_file eq '') or (!-f $sif_file);
die "bam list is not specified:\n" if ($bam_list eq '');
die "config file is not specified:\n" if ($config eq '');
die "tmp directory not specified:\n" if ($temp_dir eq '');
die "account is not specified:\n" if ($account eq '');

my $sbatch_opt = "sbatch -A $account";
if ($quos ne ''){
	$sbatch_opt .= " --qos $quos";
}
if ($partition ne ''){
	$sbatch_opt .= " -p $partition";
}
if ($memory ne ''){
	$sbatch_opt .= " --mem $memory";
}
if ($time > 0){
	$sbatch_opt .= " -t $time";
}
if ($node_list ne ''){
	$sbatch_opt .= " -w $node_list";
}

my $ref = '';
my $target_chr = 'ALL';
my $non_human = 0;
my $run_svcaller_dir = '/opt/local/tools/MOPline/scripts/run_SVcallers';

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

my $cur_dir = `pwd`;
chomp $cur_dir;

my $bind_dir2 = $temp_dir;
$bind_dir2 .= ",$cur_dir";
$bind_dir2 .= ",$bind_dir" if ($bind_dir ne '');

my $bam_path = '';

my $sbatch_flag = 0;
my $return = `which sbatch`;
chomp $return;
$sbatch_flag = 1 if ($return !~ /no|\s/);


open (FILE, $bam_list) or die "$bam_list is not found:$!\n";
while (my $line = <FILE>){
	chomp $line;
	if ($line =~ /^#\s*(\S+)/){
		my $path = $1;
		if ($path =~ /\//){
			$bam_path = $path;
		}
	}
	next if ($line =~ /^#|^$/);
	my $bam = $line;
	my $bam_base = basename ($line);
	my $ID = $bam_base;
	$ID = $1 if ($bam_base =~ /(.+?)\./);
	if ($bam_base =~ /\.cram$/){
		die "input alignment file should be in bam format (not cram)\n";
	}
	if ($bam_base !~ /\.bam$/){
		$bam_base = "$ID.bam";
		$bam .= '.bam';
	}
	system ("mkdir $ID") if (!-d $ID);
	chdir $ID;
	if (!-f $bam_base){
		$bam = "$bam_path/$bam" if ($bam_path ne '');
		system ("ln -s $bam");
		system ("ln -s $bam.bai");
		$bam = $bam_base;
	}
	print STDERR "Sample: $ID\n";
	foreach my $tool_name (@tools){
		next if ($tool_name eq 'general');
		system ("mkdir $tool_name") if (!-d $tool_name);
		chdir $tool_name;
		my $run_script = $tool_script{$tool_name};
		my $opt_str = '';
		my $thread = 1;
		foreach my $opt (keys %{$tool_opt{$tool_name}}){
			$opt_str .= "${$tool_opt{$tool_name}}{$opt} ";
			if ($opt eq 'threads'){
				$thread = $1 if (${$tool_opt{$tool_name}}{$opt} =~ /\s(\d+)/);
			}
		}
#print STDERR "$tool_name\t$run_script $opt_str\n";
		$bam = "../$bam" if ($bam !~ /\//);
		$opt_str .= "-nh $non_human " if ($tool_name =~ /CNVnator|inGAP|MELT|Wham|DELLY|Lumpy|SoftSV|Manta|GRIDSS/);
		$opt_str =~ s/\s$//;
		if (($tool_name eq 'CNVnator') and ($opt_str =~ /-r\s/)){
			$opt_str = "-b $bam -p $ID " . $opt_str;
		}
		else{
			$opt_str = "-b $bam -p $ID -r $ref " . $opt_str;
		}
		my $sbatch_opt2 = $sbatch_opt;
		$sbatch_opt2 .= " -c $thread";
		my $error_log = "$ID.error.log";
		my $out_log = "$ID.out.log";
		my $command = "singularity exec --bind $bind_dir2 $sif_file $run_script $opt_str";
		$command = "singularity exec --bind $bind_dir2 --no-home $sif_file $run_script $opt_str" if ($no_home == 1);
		if ($tool_name =~ /MELT/){
			$command = "$run_script $opt_str";
		}
		open (OUT, "> run.sh");
		print OUT "#! /usr/bin/bash\n\n";
		print OUT "$command\n";
		close (OUT);
		open (OUT, "> $ID.command.log");
		print OUT "slurm command: $sbatch_opt2 -o $out_log -e $error_log run.sh\n";
		close (OUT);
		if ($sbatch_flag == 1){
			my $jobid = `$sbatch_opt2 -o $out_log -e $error_log run.sh`;
			$jobid = $1 if ($jobid =~ /(\d+)/);
			
      print STDERR "$tool_name:$jobid\n";
    }
    else{
    	print STDERR "No sbatch command: $sbatch_opt2\n";
    }
    chdir '..';
	}
	print STDERR "\n";
	chdir $cur_dir;
}
close (FILE);

