#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;

my $sample_list = '';

my $read_length = 150;

my $non_human = 0;

my $merge_dir = 'Merge_7tools';

my $group = '';

my $help;

GetOptions(
    'sample|s=s' => \$sample_list,
    'read_len|rl=i' => \$read_length,
    'nonhuman|nh=i' => \$non_human,
    'dir|d=s' => \$merge_dir,
    'group|g=s' => \$group,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  merge_SV_vcf.7tools.pl -s <sample-list-file or a sample name> -rl <read_length> (-nh 1 if sample is a non-human species)
  output vcf file name: ${sample_name}.Merge.ALL.vcf
  Tools used: 

  Options:
   --sample or -s <STR>         a file of sample name list (or a sample name can be specified if a single sample is treated)
   --dir or -d <STR>            output directory name, which contains output vcf files [default: Merge_7tools]
   --group or -g <STR>          group ID to be added to output vcf file name, for distingwishing from other sample groups (output vcf file: ${groupID}.${sample_name}.Merge.ALL.vcf) [optional]
   --read_length or -rl <STR>   mean read length in bam files used for SV calling [default: 150]
   --non_human or -nh <INT>     samples are non-human species (human: 0, non-human: 1) [default: 0]
   --help or -h                 output help message
   
=cut

my $script_dir = $Bin;
$script_dir = "$ENV{MOPLINE_DIR}/scripts" if (exists $ENV{MOPLINE_DIR});

my $path = '';
my @ID_list;

my $work_dir = `pwd`;
chomp $work_dir;

my $all_tools = 'STR';

if (-s $sample_list){
    open (FILE, $sample_list) or die "$sample_list is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^$|^#/);
        my @line = split (/\t/, $line);
        my $bam = basename ($line[0]);
        if ($bam =~ /^([^\.\-]+).*\.bam$/){
            my $ID = $1;
            push @ID_list, $ID;
        }
        else{
            push @ID_list, $line[0];
        }
    }
    close (FILE);
}
else{
    push @ID_list, $sample_list;
}

foreach my $ID (@ID_list){
    print STDERR "\n## ID: $ID ##\n";
    system ("mkdir $ID/$merge_dir") if (!-d "$ID/$merge_dir");
    system ("rm -f $ID/$merge_dir/*.vcf");
##    
    my $inGAP_ins_mr1 = 10;
    my $inGAP_ins_mr2 = 18;
    my $inGAP_del_mr1 = 12;
    
    my $inGAP_vcf = "$ID/inGAP/inGAP.$ID.vcf";
    
    my $ins_num = `grep 'SVTYPE=INS' $inGAP_vcf | wc -l`;
    chomp $ins_num;
    
    if ($ins_num > 3200){
        my $mread_set = `$script_dir/stat_inGAP_INS.pl $inGAP_vcf 3100 2700`;
        chomp $mread_set;
        if ($mread_set =~ /(\d+),(\d+)/){
            $inGAP_ins_mr1 = $1;
            $inGAP_ins_mr2 = $2;
        }
    }
    print STDERR "inGAP mr1: $inGAP_ins_mr1 inGAP mr2: $inGAP_ins_mr2\n";
    
    my $del_num = `grep 'SVTYPE=DEL' $inGAP_vcf | wc -l`;
    chomp $del_num;
    
    if ($del_num > 60000){
        $inGAP_del_mr1 = 18;
    }
    elsif ($del_num > 40000){
        $inGAP_del_mr1 = 16;
    }
    elsif ($del_num > 20000){
        $inGAP_del_mr1 = 14;
    }
#INS-start
    my $ins_tools = 'STR';
        
    system ("$script_dir/merge_SV_calls_multi_tools.pl -st INS -rl $read_length -s $ID -d $ID -t $ins_tools -nh $non_human > $ID/$merge_dir/Merge.INS.simple.vcf");

    my $ins_set = 'STR';
   
    system ("$script_dir/merge_SV_calls_multi_tools_filter.pl -t INS -ins $ins_set -v $ID/$merge_dir/Merge.INS.simple.vcf > $ID/$merge_dir/MOP.Merge.INS.vcf");
#DEL-start
    my $del_tools = 'STR';
        
    system ("$script_dir/merge_SV_calls_multi_tools_eachSize.pl -st DEL -rl $read_length -s $ID -d $ID -t $del_tools -nh $non_human > $ID/$merge_dir/Merge.DEL.simple.vcf");
    
    my $del_set_SS = 'STR';
    my $del_set_S = 'STR';
    my $del_set_M = 'STR';
    my $del_set_L = 'STR';
    
    system ("$script_dir/merge_SV_calls_multi_tools_filter.pl -t DEL -del_s $del_set_S -del_m $del_set_M -del_l $del_set_L -v $ID/$merge_dir/Merge.DEL.simple.vcf > $ID/$merge_dir/MOP.Merge.DEL.vcf") if ($del_set_SS eq 'STR') or ($del_set_SS eq '');
    system ("$script_dir/merge_SV_calls_multi_tools_filter.pl -t DEL -del_ss $del_set_SS -del_s $del_set_S -del_m $del_set_M -del_l $del_set_L -v $ID/$merge_dir/Merge.DEL.simple.vcf > $ID/$merge_dir/MOP.Merge.DEL.vcf") if ($del_set_SS ne 'STR') and ($del_set_SS ne '');
#DUP-start    
    my $dup_tools = 'STR';
        
    system ("$script_dir/merge_SV_calls_multi_tools_eachSize.pl -st DUP -rl $read_length -s $ID -d $ID -t $dup_tools -nh $non_human > $ID/$merge_dir/Merge.DUP.simple.vcf");
    
    my $dup_set_S = 'STR';
    my $dup_set_M = 'STR';
    my $dup_set_L = 'STR';
    
    system ("$script_dir/merge_SV_calls_multi_tools_filter.pl -t DUP -dup_s $dup_set_S -dup_m $dup_set_M -dup_l $dup_set_L -v $ID/$merge_dir/Merge.DUP.simple.vcf > $ID/$merge_dir/MOP.Merge.DUP.vcf");
#INV-start        
    my $inv_tools = 'STR';
        
    system ("$script_dir/merge_SV_calls_multi_tools_eachSize.pl -st INV -rl $read_length -s $ID -d $ID -t $inv_tools -nh $non_human > $ID/$merge_dir/Merge.INV.simple.vcf");
    
    my $inv_set_S = 'STR';
    my $inv_set_M = 'STR';
    my $inv_set_L = 'STR';
    
    system ("$script_dir/merge_SV_calls_multi_tools_filter.pl -t INV -inv_s $inv_set_S -inv_m $inv_set_M -inv_l $inv_set_L -v $ID/$merge_dir/Merge.INV.simple.vcf > $ID/$merge_dir/MOP.Merge.INV.vcf");
#Merge-start    
    my $out_vcf = "$ID/$merge_dir/$ID.Merge.ALL.vcf";
    $out_vcf = "$ID/$merge_dir/$group.$ID.Merge.ALL.vcf" if ($group ne '');

    my $merge_command = "$script_dir/merge_SV_calls_ALLtype.pl -rl $read_length -v ";
    $merge_command .= "$ID/$merge_dir/MOP.Merge.DEL.vcf " if (-s "$ID/$merge_dir/MOP.Merge.DEL.vcf");
    $merge_command .= "$ID/$merge_dir/MOP.Merge.DUP.vcf " if (-s "$ID/$merge_dir/MOP.Merge.DUP.vcf");
    $merge_command .= "$ID/$merge_dir/MOP.Merge.INS.vcf " if (-s "$ID/$merge_dir/MOP.Merge.INS.vcf");
    $merge_command .= "$ID/$merge_dir/MOP.Merge.INV.vcf " if (-s "$ID/$merge_dir/MOP.Merge.INV.vcf");
    $merge_command .= "> $out_vcf";
    system ("$merge_command");
    system ("rm $ID/$merge_dir/$group.$ID.Merge.ALL.noGT.vcf") if (-f "$ID/$merge_dir/$group.$ID.Merge.ALL.noGT.vcf");
}
 