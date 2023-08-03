#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $input_bam = '';
my $ref = '';
my $target_chr = 'ALL';
my $out_prefix = '';
my $cnvnator_path = '';
my $root_dir = '';
my $non_human = 0;
my $build = '37';
my $gap_bed = '-';
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c=s' => \$target_chr,
    'command_path|cp=s' => \$cnvnator_path,
    'root_dir|rd=s' => \$root_dir,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'build=s' => \$build,
    'gap_bed|g=s' => \$gap_bed,
    'mopdir|md=s' => \$MOP_dir,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_CNVnator.pl -b <input bam file> -r <reference directory> -p <prefix of output vcf> -cd <cnvnator directory> -rt <root directory>

  Options:
   --bam or -b <STR>        input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>        reference fasta or reference directory containing chr-separated fasta files [mandatory]
   --target or -c <STR>     target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --command_path or -cp <STR>  path of cnvnator command (e.g., ./CNVnator-0.4.1/cnvnator'. This need not be specified if the directory has been set with $PATH
   --root_dir or -rd <STR>  ROOT package install directory, in which bin and lib directories should be present [mandatory]
   --prefix or -p <STR>     outpout prefix [mandatory]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --build <STR>            reference build (GRCh37, GRCh38, T2T-CHM13) number (37, 38, or T2T) when sample is human [default: 37]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --gap_bed or -g <STR>    gap bed file indicating gap regions in reference (chr start end, separated with tab). 
                            Necessary for non-human species (Data/gap.bed or gap.b38.bed for human)
   --help or -h             output help message
   
=cut

my $cur_dir = `pwd`;
chomp $cur_dir;

my $chr_list = '';

if ($non_human == 0){
    $chr_list = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y';
    $chr_list = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY' if ($build eq '38') or ($build eq 'T2T');
    if ($build eq '38'){
        if ($MOP_dir ne ''){
            $gap_bed = "$MOP_dir/Data/gap.b38.bed";
        }
        else{
            $gap_bed = "$Bin/../../Data/gap.b38.bed";
        }
    }
    elsif ($build eq '37'){
        if ($MOP_dir ne ''){
            $gap_bed = "$MOP_dir/Data/gap.bed";
        }
        else{
            $gap_bed = "$Bin/../../Data/gap.bed";
        }
    }
}
else{
    if (-d $ref){
        my @files = <$ref/*>;
        my @chr;
        foreach my $fasta (sort @files){
            open (FILE, $fasta) or die "$fasta is not found:$!\n";
            while (my $line = <FILE>){
                chomp $line;
                if ($line =~ /^>(\S+)/){
                    push @chr, $1;
                    last;
                }
            }
        }
        $chr_list = join (' ', @chr);
    }
    else{
        my $ref_index = "$ref.fai";
        if (!-f $ref_index){
            system ("ln -s $ref");
            $ref = basename ($ref);
            system ("samtools faidx $ref");
            $ref_index = "$ref.fai";
        }
        my @chr;
        open (FILE, $ref_index) or die "$ref_index is not found:\n";
        while (my $line = <FILE>){
            chomp $line;
            my ($chr) = split (/\t/, $line);
            push @chr, $chr;
        }
        close (FILE);
        $chr_list = join (' ', @chr);
    }
}

print STDERR "gap-bed: $gap_bed\n";

if ($target_chr ne 'ALL'){
    my @tchr = split (/,/, $target_chr);
    my %tchr;
    map{$tchr{$_} = 1} @tchr;
    my @chr = split (/\s+/, $chr_list);
    $chr_list = '';
    foreach (@chr){
        $chr_list .= "$_ " if (exists $tchr{$_});
    }
    $chr_list =~ s/\s$//;
}

$ENV{ROOTSYS} = $root_dir;
$ENV{PATH} = "$root_dir/bin:" . $ENV{PATH};
$ENV{LD_LIBRARY_PATH} = "$root_dir/lib:" . $ENV{LD_LIBRARY_PATH};

$cnvnator_path = 'cnvnator' if ($cnvnator_path eq '');

my $command_extract = "$cnvnator_path -root out.root -chrom $chr_list -tree $input_bam 2>$out_prefix.extract.log";
my $command_histgram = "$cnvnator_path -root out.root -chrom $chr_list -his 1000 -fasta $ref 2>$out_prefix.histgram.log";
$command_histgram = "$cnvnator_path -root out.root -chrom $chr_list -his 1000 -d $ref 2>$out_prefix.histgram.log" if (-d $ref);
my $command_stat = "$cnvnator_path -root out.root -chrom $chr_list -stat 1000 2>$out_prefix.stat.log";
my $command_partition = "$cnvnator_path -root out.root -chrom $chr_list -partition 1000 2>$out_prefix.partition.log";
my $command_call = "$cnvnator_path -root out.root -chrom $chr_list -call 1000 > $out_prefix.out 2>$out_prefix.call.log";

my @out_files;

open (OUT, ">> $out_prefix.command.log");
print OUT "1st command: $command_extract\n";
close (OUT);

system ("$command_extract");

open (OUT, ">> $out_prefix.command.log");
print OUT "2nd command: $command_histgram\n";
close (OUT);

system ("$command_histgram");

open (OUT, ">> $out_prefix.command.log");
print OUT "3rd command: $command_stat\n";
close (OUT);

system ("$command_stat");

open (OUT, ">> $out_prefix.command.log");
print OUT "4th command: $command_partition\n";
close (OUT);

system ("$command_partition");

open (OUT, ">> $out_prefix.command.log");
print OUT "5th command: $command_call\n";
close (OUT);

system ("$command_call");

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_CNVnator_vcf.pl $out_prefix.out $gap_bed > CNVnator.$out_prefix.vcf");
}
else{
    system ("$Bin/convert_CNVnator_vcf.pl $out_prefix.out $gap_bed > CNVnator.$out_prefix.vcf");
}

system ("rm out.root") if (-f 'out.root');

print STDERR "CNVnator run was completed\n";
