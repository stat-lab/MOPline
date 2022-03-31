#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;

my $tool_set = '';

my $tool_pair_config = '';

my $data_cov = 'H';

my $out_prefix = '';

my $help;

GetOptions(
    'tool_set|t=s' => \$tool_set,
    'tool_config|tc=s' => \$tool_pair_config,
    'data_cov|dc=s' => \$data_cov,
    'prefix|p=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  make_merge_SV_vcf_script.pl -t <comma-separated tool-names or a list file of tool names> -p <output prefix>
  output a script for overlap call selection: ${out_prefix}.merge_sv_vcf.pl
  (This script creates a custom perl script, such as merge_SV_vcf.7tools.pl, for selecting overlap calls from tool pairs using a user-specified tool set)

  Options:
   --tool_set or -t <STR>       Comma-sepated tool-names or a list file with tool name per line. Tool names should be ones contained in pre-defined tool set, specified in ${tool_config} file (e.g., CNVnator,DELLY,inGAP,Manta,Lumpy) [mandatory]
                                If tools of interest are not specified in ${tool_config} file, the tool names and their minimum RSSs in each SV category should be added to the ${tool_config} file.
                                All the tool names should exactly match the tool directory name and the tool names contained in the SV vcf file in the tool directory.
   --tool_config or -tc <STR>   configuration file for tool pairs and minimum RSS to select overlap calls [default: MOPline_xx/Data/SVtool_pairs_config.txt]
   --data_cov -dc <STR>         data coverage (H or L, specify 'L' for ~20x coverage) [default: H]
   --prefix or -p <STR>         Output prefix [mandatory]
   --help or -h                 Output help message
   
=cut

if ($tool_pair_config eq ''){
    if ($data_cov eq 'H'){
        $tool_pair_config = "$Bin/../Data/SVtool_pairs_config.txt";
    }
    else{
        $tool_pair_config = "$Bin/../Data/SVtool_pairs_config_20x.txt";
    }
}

my $script_template = "$Bin/merge_SV_vcf_template.pl";

my %tools;
my %tool_conf;
my %tool_pairs;
my %type_tools;
my %used;

if (-f $tool_set){
    open (FILE, $tool_set) or die "$tool_set is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        $tools{$line} = 1;
    }
    close (FILE);
}
else{
    my @tools = split (/,/, $tool_set);
    map{$tools{$_} = 1} @tools;
}

my $type = '';

open (FILE, $tool_pair_config) or die "$tool_pair_config is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    if ($line =~ /^<(.+)>/){
        $type = $1;
        next;
    }
    my ($tool, $pair) = split (/\s+/, $line);
    my @pairs = split (/,/, $pair);
    foreach my $tool2 (@pairs){
        ${${$tool_conf{$type}}{$tool}}{$tool2} = 1;
    }
}
close (FILE);

foreach my $type (keys %tool_conf){
    my $type_base = $type;
    $type_base = $1 if ($type =~ /(.+)-/);
    foreach my $tool1 (sort {lc $a cmp lc $b} keys %{$tool_conf{$type}}){
        my ($tool1_base, $rss1) = split (/:/, $tool1);
        next if (!exists $tools{$tool1_base});
        ${$type_tools{$type_base}}{$tool1_base} = 1;
        my $tool1a = $tool1;
        foreach my $tool2 (sort {lc $a cmp lc $b} keys %{${$tool_conf{$type}}{$tool1}}){
            my ($tool2_base, $rss2) = split (/:/, $tool2);
            next if ($tool2 ne 'ANY') and ($tool2 ne 'SINGLE') and (!exists $tools{$tool2_base});
            my $tool2a = $tool2;
            if (($tool1_base eq 'inGAP') and ($type eq 'INS')){
                if ($tool2 eq 'SINGLE'){
                    $tool1a = 'inGAP:$inGAP_ins_mr2';
                }
                else{
                    $tool1a = 'inGAP:$inGAP_ins_mr1';
                }
            }
            elsif (($tool2_base eq 'inGAP') and ($type eq 'INS')){
                $tool2a = 'inGAP:$inGAP_ins_mr1';
            }
            elsif (($tool1_base eq 'inGAP') and ($type eq 'DEL-L') and ($tool2 eq 'SINGLE')){
                $tool1a = 'inGAP:$inGAP_del_mr1';
            }
            my $pair = "$tool1a=$tool2a";
            $pair = "$tool1a=#" if ($tool2 eq 'ANY');
            $pair = $tool1a if ($tool2 eq 'SINGLE');
            my $pair2 = "$tool2a=$tool1a";
            next if (exists ${$used{$type}}{$pair2}) and ($tool2 ne 'ANY') and ($tool2 ne 'SINGLE');
            push @{$tool_pairs{$type}}, $pair;
            ${$used{$type}}{$pair} = 1;
        }
    }
}

my $all_tools = '';
my $ins_tools = '';
my $del_tools = '';
my $dup_tools = '';
my $inv_tools = '';
my $inGAP_flag = 0;
$inGAP_flag = 1 if (exists $tools{'inGAP'});

foreach my $tool (sort keys %tools){
    $all_tools .= "$tool ";
}
$all_tools =~ s/\s$//;

if (exists $type_tools{'INS'}){
    foreach my $tool (sort keys %{$type_tools{'INS'}}){
        $ins_tools .= "$tool ";
    }
    $ins_tools =~ s/\s$//;
}
if (exists $type_tools{'DEL'}){
    foreach my $tool (sort keys %{$type_tools{'DEL'}}){
        $del_tools .= "$tool ";
    }
    $del_tools =~ s/\s$//;
}
if (exists $type_tools{'DUP'}){
    foreach my $tool (sort keys %{$type_tools{'DUP'}}){
        $dup_tools .= "$tool ";
    }
    $dup_tools =~ s/\s$//;
}
if (exists $type_tools{'INV'}){
    foreach my $tool (sort keys %{$type_tools{'INV'}}){
        $inv_tools .= "$tool ";
    }
    $inv_tools =~ s/\s$//;
}

my $ins_pairs = '';
my $del_ss_pairs = '';
my $del_s_pairs = '';
my $del_m_pairs = '';
my $del_l_pairs = '';
my $dup_s_pairs = '';
my $dup_m_pairs = '';
my $dup_l_pairs = '';
my $inv_s_pairs = '';
my $inv_m_pairs = '';
my $inv_l_pairs = '';

if (exists $tool_pairs{'INS'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'INS'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'INS'}};
    $ins_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $ins_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DEL-SS'}){
    $del_ss_pairs = join (' ', @{$tool_pairs{'DEL-SS'}});
}
if (exists $tool_pairs{'DEL-S'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DEL-S'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DEL-S'}};
    $del_s_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $del_s_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DEL-M'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DEL-M'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DEL-M'}};
    $del_m_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $del_m_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DEL-L'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DEL-L'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DEL-L'}};
    $del_l_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $del_l_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DUP-S'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DUP-S'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DUP-S'}};
    $dup_s_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $dup_s_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DUP-M'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DUP-M'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DUP-M'}};
    $dup_m_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $dup_m_pairs .= " $single";
    }
}
if (exists $tool_pairs{'DUP-L'}){
    my @single = grep {$_ !~ /=/} @{$tool_pairs{'DUP-L'}};
    my @pairs = grep {$_ =~ /=/} @{$tool_pairs{'DUP-L'}};
    $dup_l_pairs = join (' ', @pairs);
    if (@single > 0){
        my $single = join (' ', @single);
        $dup_l_pairs .= " $single";
    }
}
if (exists $tool_pairs{'INV-S'}){
    $inv_s_pairs = join (' ', @{$tool_pairs{'INV-S'}});
}
if (exists $tool_pairs{'INV-M'}){
    $inv_m_pairs = join (' ', @{$tool_pairs{'INV-M'}});
}
if (exists $tool_pairs{'INV-L'}){
    $inv_l_pairs = join (' ', @{$tool_pairs{'INV-L'}});
}

my $flag = 0;
my $out_script = "$out_prefix.merge_sv_vcf.pl";

open (OUT, "> $out_script");
open (FILE, $script_template) or die "$script_template is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^##/){
        $flag = 1;
    }
    elsif ($line =~ /^#INS-start/){
        $flag = 2;
    }
    if (($flag == 1) and ($inGAP_flag == 0)){
        next;
    }
    if ($line =~ /Tools\sused:/){
        $line .= $all_tools;
    }
    elsif ($line =~ /my\s\$ins_tools\s=\s\'STR/){
        $line =~ s/STR/$ins_tools/;
    }
    elsif ($line =~ /my\s\$del_tools\s=\s\'STR/){
        $line =~ s/STR/$del_tools/;
    }
    elsif ($line =~ /my\s\$dup_tools\s=\s\'STR/){
        $line =~ s/STR/$dup_tools/;
    }
    elsif ($line =~ /my\s\$inv_tools\s=\s\'STR/){
        $line =~ s/STR/$inv_tools/;
    }
    elsif ($line =~ /my\s\$ins_set\s=\s\'STR/){
        $line =~ s/STR/$ins_pairs/ if ($inGAP_flag == 0);
        $line =~ s/\'STR\'/\"$ins_pairs\"/ if ($inGAP_flag == 1);
    }
    elsif ($line =~ /my\s\$del_set_SS\s=\s\'STR/){
        $line =~ s/STR/$del_ss_pairs/;
    }
    elsif ($line =~ /my\s\$del_set_S\s=\s\'STR/){
        $line =~ s/STR/$del_s_pairs/;
    }
    elsif ($line =~ /my\s\$del_set_M\s=\s\'STR/){
        $line =~ s/STR/$del_m_pairs/;
    }
    elsif ($line =~ /my\s\$del_set_L\s=\s\'STR/){
        $line =~ s/STR/$del_l_pairs/ if ($inGAP_flag == 0);
        $line =~ s/\'STR\'/\"$del_l_pairs\"/ if ($inGAP_flag == 1);
    }
    elsif ($line =~ /my\s\$dup_set_S\s=\s\'STR/){
        $line =~ s/STR/$dup_s_pairs/;
    }
    elsif ($line =~ /my\s\$dup_set_M\s=\s\'STR/){
        $line =~ s/STR/$dup_m_pairs/;
    }
    elsif ($line =~ /my\s\$dup_set_L\s=\s\'STR/){
        $line =~ s/STR/$dup_l_pairs/;
    }
    elsif ($line =~ /my\s\$inv_set_S\s=\s\'STR/){
        $line =~ s/STR/$inv_s_pairs/;
    }
    elsif ($line =~ /my\s\$inv_set_M\s=\s\'STR/){
        $line =~ s/STR/$inv_m_pairs/;
    }
    elsif ($line =~ /my\s\$inv_set_L\s=\s\'STR/){
        $line =~ s/STR/$inv_l_pairs/;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools\.pl\s-st\sINS/) and ($ins_tools eq '')){
        $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_eachSize\.pl\s-st\sDEL/) and ($del_tools eq '')){
        $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_eachSize\.pl\s-st\sDUP/) and ($dup_tools eq '')){
        $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_eachSize\.pl\s-st\sINV/) and ($inv_tools eq '')){
        $line = '#' . $line;
    }
    if (($line =~ /merge_SV_calls_multi_tools_filter\.pl\s-t\sDEL/) and ($del_s_pairs eq '') and ($del_m_pairs eq '') and ($del_l_pairs eq '')){
        $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_filter\.pl\s-t\sDUP/) and ($dup_s_pairs eq '') and ($dup_m_pairs eq '') and ($dup_l_pairs eq '')){
         $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_filter\.pl\s-t\sINV/) and ($inv_s_pairs eq '') and ($inv_m_pairs eq '') and ($inv_l_pairs eq '')){
         $line = '#' . $line;
    }
    elsif (($line =~ /merge_SV_calls_multi_tools_filter\.pl\s-t\sINS/) and ($ins_pairs eq '')){
         $line = '#' . $line;
    }

    print OUT $line, "\n";
}
close (FILE);
close (OUT);
