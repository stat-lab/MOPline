#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $vcf = '';
my $tool_set = '7tools';
my $tool_custom = '';
my $sample_name = '';
my $sample_dir = '';
my $min_score = 85;
my $help;

my @pre_tool_set = ('CNVnator:DEL:96', 'CNVnator:DUP:93', 'DELLY:DEL:94', 'DELLY:INV:80', 'Lumpy:DEL:92', 'Lumpy:INV:80', 'Manta:DEL:94', 'Manta:INV:80', 'Manta:INS:98', 'MELT:INS:70');
#my @pre_tool_set = ('CNVnator:DEL:96', 'CNVnator:DUP:93', 'Manta:DEL:94', 'Manta:INV:80', 'Manta:INS:98', 'MELT:INS:70');

GetOptions(
    'vcf|v=s' => \$vcf,
    'sample_name|sn=s' => \$sample_name,
    'sample_dir|sd=s' => \$sample_dir,
    'toolset|ts=s' => \$tool_set,
    'min_score|ms=i' => \$min_score,
    'help|h' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  add_genotype_SV_callers_vcf.pl -v <vcf file> -sn <sample name> -sd <sample directory> -ts <toolset>
  (Add genotypes from 6~9 tools used for SV calling to SVs in a vcf file)

  Options:
   --vcf or -v <STR>        SV vcf file path to be genotyped
   --sample_name or -sn <STR>  sample name [mandatory]
   --sample_dir or -sd <STR>   sample directory (optional, sample_dir is regarded as sample_name if not specified)
   --toolset or -ts <STR>    a preset toolset of SV detection tools used (6tools_1, 6tools_2, 7tools, 9tools, or 11tools) or a list file describing tool names in each line [default: 7tools]
   --min_score or -ms <INT> minimum score to accept genotype (if multiple tools call the same genotype at a site, a summed score from these tools is considered) [default: 85]
   --help or -h             output help message
   
=cut

die "SV vcf file to genotype is not specified: \n" if ($vcf eq '');
die "Sample name is not specified: \n" if ($sample_name eq '');

$sample_dir = $sample_name if ($sample_dir eq '');

my @tool_set = ();
my @tools = ();

@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham') if ($tool_set eq '4tools');
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'MELT') if ($tool_set eq '6tools_1');
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'GRIDSS') if ($tool_set eq '6tools_2') or ($tool_set eq '7tools');
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'DELLY', 'MELT', 'Lumpy', 'SoftSV') if ($tool_set eq '9tools');
@tools = ('CNVnator', 'inGAP', 'Manta', 'Wham', 'MATCHCLIP', 'DELLY', 'MELT', 'Lumpy', 'SoftSV', 'forestSV', 'Mobster') if ($tool_set eq '11tools');
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

my $dir = `pwd`;
chomp $dir;
my $corrent_dir = $1 if ($dir =~ /\/([^\/]+)$/);

my %tools;
my %genotype;
my %vcf_files;

foreach my $tool_type (@pre_tool_set){
    my ($tool, $type, $score) = split (/:/, $tool_type);
    ${$tools{$tool}}{$type} = $score;
}
foreach my $tool (@tools){
    next if (!exists $tools{$tool});
    my @vcf_tool = ();
#print STDERR "$tool\n";
    if (-d $sample_dir){
        my @dir = <$sample_dir/*>;
#print STDERR "@dir\n";
        my @dir_tool = grep{$_ =~ /\/$tool$/i} @dir;
#print STDERR "@dir_tool\n";
        if (@dir_tool > 0){
            foreach my $dir (@dir_tool){
                if (-d $dir){
                    my @vcf2 = <$dir/*.vcf>;
                    @vcf_tool = grep{$_ =~ /$dir\/$tool[\.\-_]/i} @vcf2;
                    if (@vcf_tool > 0){
                        $vcf_files{$tool} = $vcf_tool[0];
                    }
                    else{
                        print STDERR "No vcf file prefixed with $sample_name $tool in ../$tool/:\n";
                    }
                }
            }
        }
        else{
            my @dir_tool2 = grep{$_ =~ /\/$tool/i} @dir;
            if (@dir_tool2 > 0){
                foreach my $vcf (@dir_tool2){
                    if ($vcf =~ /$tool\.$sample_name\.vcf$/i){
                        $vcf_files{$tool} = $vcf;
                    }
                }
            }
            else{
                print STDERR "No directory/files with $sample_name $tool name:\n";
            }
        }
    }
    else{
        die "$sample_dir directory is not found:\n";
    }
}

if (scalar keys %vcf_files == 0){
    die "Vcf files, corresponding to the tool set, are not found: \n";
}

foreach my $tool (keys %vcf_files){     # genotyped vcf files
    my $vcf_file = $vcf_files{$tool};
print STDERR "$vcf_file\n";
    open (FILE, $vcf_file) or die "$vcf_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#/);
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        $type = 'INS' if ($type =~ /ALU|LINE1|L1|SVA|HERVK|NUMT|VEI|MEI/);
        my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
        if ($line[7] !~ /GT=([^;]+)/){
#            print STDERR "Warning: No GT tags in $sample_name $vcf_file: \n";
            next;
        }
        my $gt = $1 if ($line[7] =~ /GT=([^;]+)/);
        $gt = './.' if ($gt eq '0/0');
        push @{${${$genotype{$type}}{$chr}}{$pos}}, "$tool==$len==$gt";
    }
    close (FILE);
}

my $min_overlap_ratio = 0.5;
my $ins_sd = 1000;
my $var_sd = 125;

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

open (FILE, $vcf) or die "$vcf is not found: \n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^#/){
        if ($line =~ /fileDate=/){
            $line =~ s/fileDate=\d+/fileDate=$time/;
        }
        elsif ($line =~ /^#CHROM/){
            $line .= "\tFORMAT\t$sample_name";
        }
        print "$line\n";
    }
    next if ($line =~ /^#/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos1 = $line[1];
    my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
    my $len1 = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
    my $end1 = $pos1 + $len1 - 1;
    my %hit;
    foreach my $pos2 (sort {$a <=> $b} keys %{${$genotype{$type}}{$chr}}){
        last if ($end1 + $ins_sd < $pos2);
        foreach my $info (@{${${$genotype{$type}}{$chr}}{$pos2}}){
            my ($tool, $len2, $gt) = split (/==/, $info);
            my $end2 = $pos2 + $len2 - 1;
            next if ($end2 + $ins_sd < $pos1);
            my $flag = 0;
            if ($type eq 'INS'){
                if (abs ($pos1 - $pos2) <= $ins_sd){
                    $flag = 1;
                }
            }
            else{
                if ((abs ($pos1 - $pos2) <= $var_sd) and (abs ($end1 - $end2) <= $var_sd)){
                    $flag = 1;
                }
                if (($pos1 >= $pos2) and ($pos1 <= $end2)){
                    if ($end1 <= $end2){
                        if ($len1 >= $len2 * $min_overlap_ratio){
                            $flag = 1;
                        }
                    }
                    else{
                        if (($end2 - $pos1 >= $len2 * $min_overlap_ratio) and ($end2 - $pos1 >= $len1 * $min_overlap_ratio)){
                            $flag = 1;
                        }
                    }
                }
                elsif (($pos2 >= $pos1) and ($pos2 <= $end1)){
                    if ($end2 <= $end1){
                        if ($len2 >= $len1 * $min_overlap_ratio){
                            $flag = 1;
                        }
                    }
                    else{
                        if (($end1 - $pos2 >= $len2 * $min_overlap_ratio) and ($end1 - $pos2 >= $len1 * $min_overlap_ratio)){
                            $flag = 1;
                        }
                    }
                }
            }
            if ($flag == 1){
                if (!exists $hit{$tool}){
                    $hit{$tool} = "$gt==$pos2";
                }
                else{
                    my ($pre_gt, $pre_pos) = split (/==/, $hit{$tool});
                    if (($pre_gt eq './.') or (($pre_gt ne './.') and ($gt ne './.'))){
                        if (abs ($pos2 - $pos1) < abs ($pre_pos - $pos1)){
                            $hit{$tool} = "$gt==$pos2";
                        }
                    }
                }
            }
        }
    }
    my %GT_score;
    my $top_GT = '';
    my $GT_qual = 0;
    my $top_score = 0;
    if (scalar keys %hit > 0){
        foreach my $tool (keys %hit){
            my ($gt) = split (/==/, $hit{$tool});
            my $GT = '';
            $GT = 'HM' if ($gt eq '1/1');
            $GT = 'HT' if ($gt eq '1/0') or ($gt eq '0/1');
            my $score = 0;
            $score = ${$tools{$tool}}{$type} if (exists ${$tools{$tool}}{$type});
            $GT_score{$GT} += $score if ($GT ne '');
#print STDERR "$pos1\t$len1\t$tool\t$gt\t$score\n" if ($chr eq '1') and ($pos1 == 21399001) and ($type eq 'DUP');
        }
        if (scalar keys %GT_score > 0){
            foreach my $gt2 (sort {$GT_score{$b} <=> $GT_score{$a}} keys %GT_score){
                $top_GT = $gt2;
                $top_score = $GT_score{$gt2};
                last;
            }
        }

    }
    $top_GT = './.' if ($top_GT eq '') or ($top_score < $min_score);
    $top_GT = '0/1' if ($top_GT eq 'HT');
    $top_GT = '1/1' if ($top_GT eq 'HM');
    $line[7] =~ s/SVTYPE=.+?;/SVTYPE=$type;/;
    $line[7] =~ s/SVLEN=/SVLEN=-/ if ($type eq 'DEL') and ($line[7] !~ /SVLEN=-/);
    $line[7] .= ';' . $line[8] if (@line >= 9) and ($line[8] =~ /^TOOLS=/);
    $line[8] = 'GT:GQ:VP:VL';
    $line[9] = "$top_GT:$top_score:$pos1:$len1";
    $line[2] = '.' if ($line[2] =~ /DEL|DUP|INS|INV|TRA/);
    if ($line[7] !~ /END=/){
        $end1 = $pos1 if ($type eq 'INS');
        $line[7] =~ s/;READS/;END=$end1;READS/;
    }   
    print join ("\t", @line), "\n";
}
close (FILE);

