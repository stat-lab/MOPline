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
my $inGAP_jar = '';
my $non_human = 0;
my $MOP_dir = '';

my $help;

GetOptions(
    'bam|b=s' => \$input_bam,
    'ref|r=s' => \$ref,
    'target|c=s' => \$target_chr,
    'ingap_jar|ij=s' => \$inGAP_jar,
    'prefix|p=s' => \$out_prefix,
    'non_human|nh=i' => \$non_human,
    'mopdir|md=s' => \$MOP_dir,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_inGAP.pl -b <input bam file> -r <reference fasta> -p <prefix of output vcf> -ij <inGAP-sv jar command path>

  Options:
   --bam or -b <STR>          input bam file, index file (*.bam.bai) should be present in the same directory [mandatory]
   --ref or -r <STR>          reference fasta file, index file (*.fa.fai) should be present in the same directory [mandatory]
   --target or -c <STR>       target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --ingap_jar or -ij <STR>   inGAP-sv jar file (e.g., ..../inGAP_3_1_1/inGAP.jar) [mandatory]
   --non_human or -nh <INT>   samples are non-human species (0: human, 1: non-human) [default: 0]
   --mopdir or -md <STR>      MOPline install directory (optional)
   --prefix or -p <STR>       outpout prefix [mandatory]
   --help or -h               output help message
   
=cut

my %target_chr;
my @tchr = split (/,/, $target_chr) if ($target_chr ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

my @chr_list;

my $ref_index = "$ref.fai";

open (FILE, $ref_index) or die "$ref_index is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr) = split (/\t/, $line);
    if ($non_human == 0){
        if ($target_chr eq 'ALL'){
            push @chr_list, $chr if ($chr =~ /^c*h*r*[\dXY]+$/);
        }
        else{
            push @chr_list, $chr if (exists $target_chr{$chr});
        }
    }
    else{
        if ($target_chr eq 'ALL'){
            push @chr_list, $chr;
        }
        else{
            push @chr_list, $chr if (exists $target_chr{$chr});
        }
    }
}
close (FILE);

my $header = '';
my $seq = '';
my $ref_base = basename ($ref);
$ref_base = $1 if ($ref_base =~ /(.+)\.fa*s*t*a$/);
my %ref_chr;

open (FILE, $ref) or die "$ref is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        if ($seq ne ''){
            if ($non_human == 0){
                if ($header =~ /^c*h*r*[\dXY]+$/){
                    my $ref_chr = "$ref_base.$header.fa";
                    open (OUT, "> $ref_chr");
                    print OUT ">$header\n";
                    print OUT $seq, "\n";
                    close (OUT);
                    $ref_chr{$header} = $ref_chr;
                }
            }
            else{
                my $ref_chr = "$ref_base.$header.fa";
                open (OUT, "> $ref_chr");
                print OUT ">$header\n";
                print OUT $seq, "\n";
                close (OUT);
                $ref_chr{$header} = $ref_chr;
            }
        }
        $seq = '';
        $header = $1;
    }
    else{
        $seq .= $line;
    }
}
close (FILE);
if ($seq ne ''){
    if ($non_human == 0){
        if ($header =~ /^c*h*r*[\dXY]+$/){
            my $ref_chr = "$ref_base.$header.fa";
            open (OUT, "> $ref_chr");
            print OUT ">$header\n";
            print OUT $seq, "\n";
            close (OUT);
            $ref_chr{$header} = $ref_chr;
        }
    }
    else{
        my $ref_chr = "$ref_base.$header.fa";
        open (OUT, "> $ref_chr");
        print OUT ">$header\n";
        print OUT $seq, "\n";
        close (OUT);
        $ref_chr{$header} = $ref_chr;
    }
}


my @out_file;

foreach my $chr (@chr_list){
    my $chr_sam = "$out_prefix.$chr.sam";
    my $out_file = "$out_prefix.$chr.out";
    my $ref_chr = $ref_chr{$chr};
    system ("samtools view -h $input_bam $chr > $chr_sam") if (!-f $chr_sam) or (-z $chr_sam);

    my $command = "java -mx4000m -jar $inGAP_jar SVP -SIZE 1000000 -r $ref_chr -i $chr_sam -o $out_file 2>$out_prefix.inGAP.$chr.log";
    open (OUT2, ">> $out_prefix.command.log");
    print OUT2 "inGAP command: $command\n";
    close (OUT2);
    
    system ("$command");
    push @out_file, "$out_file";

    system ("rm $chr_sam");
}

my $arg = join (' ', @out_file);

if ($MOP_dir ne ''){
    system ("$MOP_dir/scripts/run_SVcallers/convert_inGAP_vcf.pl $arg > inGAP.$out_prefix.vcf");
}
else{
    system ("$Bin/convert_inGAP_vcf.pl $arg > inGAP.$out_prefix.vcf");
}

print STDERR "inGAP run was completed\n";
