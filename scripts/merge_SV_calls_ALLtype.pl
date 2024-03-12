#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);

# merge SV calls from vcf files of multiple types of SV


my @var_file;

my $read_length = 150;
my $var_sd = 125;
my $ins_sd = 200;
my $ins_sd2 = 30;
my $ctx_sd = 20000;

my $min_ctx = 20000;

my $min_DUP_nonDepth_caller = 1000000;
my $min_DEL_nonDepth_caller = 2000000;

my @depth_callers = ('CNVnator');

my $min_overlap_ratio = 0.75;

my $pos_ave = 1;

my $help;

GetOptions(
    'var|v=s{,}' => \@var_file,
    'var_sd|vsd=i' => \$var_sd,
    'read_len|rl=i' => \$read_length,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  merge_SV_calls_ALLtype.pl -v <INS.vcf DEL.vcf DUP.vcf INV.vcf>  > [output vcf]

  Options:
   --var or -v <STR>        list (space-delimited) of vcf files
   --var_sd or -vsd <INT>   standard deviation of break points [default: 125]
   --read_len or -rl <INT>  read length [default: 150]
   --help or -h             output help message
   
=cut

die "vcf files (-v) not specified: \n" if (@var_file == 0);

$var_sd = $read_length;

#print STDERR "< Parameters >\n";
print STDERR "ALL TYPES merging; @var_file\n";

my @svtype = ('INS', 'DUP', 'INV', 'DEL');
my %var;
my %mei_type;
my $numt = 0;
my $vei = 0;

foreach my $var_file (@var_file){
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
    	chomp $line;
    	if ($line =~ /^#/){
    	    next;
    	}
    	my @line = split (/\t/, $line);
    	my $chr = $line[0];
    	my $pos = $line[1];
    	my $type = $line[2];
    	my $subtype = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
        if ($type ne $subtype){
            my $stype = $type;
            $type = $subtype;
            $subtype = $stype;
        }
        if (($type eq 'INS') and ($subtype ne 'INS')){
            if ($subtype eq 'NUMT'){
                $numt = 1;
            }
            elsif ($subtype eq 'VEI'){
                $vei = 1;
            }
            else{
                $mei_type{$subtype} = 1;
            }
        }
    	my $len = 0;
    	$len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
    	my $end = 0;
    	$end = $pos + $len - 1 if ($type ne 'INS');
    	my $chr2 = '';
    	my $pos2 = 0;
    	if ($type eq 'TRA'){
    	    $chr2 = $1 if ($line[7] =~ /CHR2=(.+?);/);
    	    $pos2 = $1 if ($line[7] =~ /POS2=(\d+)/);
    	    $end = 0;
    	    $len = 0;
    	}
    	if (@line > 8){
    	    my $tool = $line[8];
    	    $line[7] .= ';' . $tool;
    	    splice (@line, -1);
    	    $line = join ("\t", @line);
    	}
    	my $flag = 0;
    	if ($type eq 'DEL'){
            if ($len > $min_DEL_nonDepth_caller){
                my $flag = 0;
                foreach my $dpcaller (@depth_callers){
                    if ($line[7] =~ /$dpcaller/i){
                        $flag = 1;
                    }
                }
                if ($flag == 0){
                    next;
                }
            }
        }
	    elsif ($type eq 'DUP'){
            if ($len > $min_DUP_nonDepth_caller){
                my $flag = 0;
                foreach my $dpcaller (@depth_callers){
                    if ($line[7] =~ /$dpcaller/i){
                        $flag = 1;
                    }
                }
                if ($flag == 0){
                    next;
                }
            }
        }
        if ((exists ${$var{$chr}}{$pos}) and ($type eq 'INS')){
            next;
        }
        elsif ((exists ${${$var{$chr}}{$pos}}{'INS'}) and ($type ne 'INS')){
            delete ${${$var{$chr}}{$pos}}{'INS'};
        }
        if (($type eq 'INS') and ($len > 0) and ($len < 40)){
            $line =~ s/SVLEN=$len/SVLEN=0/;
        }
	    ${${$var{$chr}}{$pos}}{$type} = $line;
        my $tool = $1 if ($line[7] =~ /TOOLS=(\S+)/);
    }
}

foreach my $chr (keys %var){        # filter SV with inconsistent lengths called from multiple tools
    foreach my $pos (keys %{$var{$chr}}){
        foreach my $type (keys %{${$var{$chr}}{$pos}}){
            next if ($type eq 'INS');
            my $line = ${${$var{$chr}}{$pos}}{$type};
            my $len = $1 if ($line =~ /SVLEN=-*(\d+)/);
            my $tool = $1 if ($line =~ /TOOLS=(\S+)/);
            my @tools = split (/,/, $tool);
            next if (@tools == 1);
            my @len;
            foreach (@tools){
                next if ($_ =~ /CNVnator/);
                my ($t, $p, $l) = split (/:/, $_);
                push @len, $l;
            }
            next if (@len <= 1);
            my $min_diff = 0;
            my $count1 = 0;
            foreach my $len1 (@len){
                $count1 ++;
                my $count2 = 0;
                foreach my $len2 (@len){
                    $count2 ++;
                    next if ($count2 <= $count1);
                    my $rate = 0;
                    if ($len1 >= $len2){
                        $rate = int ($len2 / $len1 * 100 + 0.5) / 100;
                    }
                    else{
                        $rate = int ($len1 / $len2 * 100 + 0.5) / 100;
                    }
                    if ($rate > $min_diff){
                        $min_diff = $rate;
                    }
                }
            }
            if ($len <= 1000){
                if ($min_diff <= 0.7){
                    delete ${${$var{$chr}}{$pos}}{$type};
                }
            }
            else{
                if ($min_diff <= 0.8){
                    delete ${${$var{$chr}}{$pos}}{$type};
                }
            }
        }
    }
}


my %vcf;
my $pre_chr = '';
my $pre_pos = 0;
my $pre_type = '';
my $pre_subtype = '';
my $pre_len = 0;
my $pre_end = 0;
my $pre_tool_num = 0;
my $pre_line = '';

foreach my $chr (sort keys %var){
    foreach my $pos (sort {$a <=> $b} keys %{$var{$chr}}){
        foreach my $type (sort keys %{${$var{$chr}}{$pos}}){
            my $line = ${${$var{$chr}}{$pos}}{$type};
            my @line = split (/\t/, $line);
            my $subtype = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
            if ($type ne $subtype){
                my $stype = $type;
                $type = $subtype;
                $subtype = $stype;
            }
            my $len = 0;
            $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
            my $end = 0;
            $end = $pos + $len - 1 if ($type ne 'INS');
            $end = $pos if ($type eq 'INS');
            my $chr2 = '';
            my $pos2 = 0;
            if ($type eq 'TRA'){
                $chr2 = $1 if ($line[7] =~ /CHR2=(.+?);/);
                $pos2 = $1 if ($line[7] =~ /POS2=(\d+)/);
                $end = 0;
                $len = 0;
            }
            my $tool = $1 if ($line[7] =~ /TOOLS=(\S+)/);
            my $tool_num = 0;
            $tool_num = $tool =~ s/,/=/g;
            $tool_num ++;
            my $chr02d = $chr;
            $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
            if ($pre_chr ne $chr){
                $pre_pos = 0;
                $pre_end = 0;
            }
	        my $pre_len = $pre_end - $pre_pos + 1;
            if ($pre_end >= $pos){
                if ($type eq 'INS'){
                    if ($pre_type eq 'INS'){
                        if (($subtype ne 'INS') and ($pre_subtype eq 'INS')){
                            delete ${$vcf{$chr02d}}{$pre_pos};
                        }
                        elsif (($subtype eq 'INS') and ($pre_subtype ne 'INS')){
                            next;
                        }
                        else{
                            if ($tool_num > $pre_tool_num){
                                delete ${$vcf{$chr02d}}{$pre_pos};
                            }
                            elsif ($tool_num < $pre_tool_num){
                                next;
                            }
                            else{
                                next;
                            }
                        }
                    }
                    elsif ($pre_type eq 'DEL'){
            			if ($pos - $pre_pos <= $ins_sd2){
            			    delete ${$vcf{$chr02d}}{$pre_pos};
            			}
            			elsif (abs ($pre_end - $pos) <= $ins_sd2){
            			    delete ${$vcf{$chr02d}}{$pre_pos};
            			}
            			else{
            			    next;
            			}
                    }
            		elsif ($pre_type eq 'DUP'){
            			if ($pos - $pre_pos <= $ins_sd2){
            			    next;
            			}
            			else{
            			    if ($tool_num >= $pre_tool_num){
            				    delete ${$vcf{$chr02d}}{$pre_pos} if ($pre_len < 50000);
            			    }
            			    elsif ($tool_num < $pre_tool_num){
            				    next;
            			    }
            			}
        		    }
        		    elsif ($pre_type eq 'INV'){
            			delete ${$vcf{$chr02d}}{$pre_pos};
            		}
                }
                elsif ($type eq 'DEL'){
                    if ($pre_type eq 'INS'){
                        if ($pre_subtype eq 'INS'){
                            delete ${$vcf{$chr02d}}{$pre_pos};
                        }
            			elsif ($pos - $pre_pos <= $ins_sd2){
            			    next;
            			}
            			elsif (abs ($pre_end - $pos) <= $ins_sd2){
            			    next;
            			}
                    }
                    elsif ($pre_type eq 'DUP'){
                        my $overlap = $pre_end - $pos + 1;
                        $overlap = $len if ($end < $pre_end);
                        if (($overlap < $pre_len * $min_overlap_ratio) and ($overlap < $len * $min_overlap_ratio)){
                            delete ${$vcf{$chr02d}}{$pre_pos};
                        }
                    }
                }
                elsif ($type eq 'DUP'){
                    if ($pre_type eq 'DUP'){
                        if ($pre_end >= $end){
                            next;
                        }
                        elsif (($pre_end - $pos >= $pre_len * $min_overlap_ratio) or ($pre_end - $pos >= $len * $min_overlap_ratio)){		# merge pre-DUP + DUP -> new-DUP
                            if ((($pre_end - $pos >= $pre_len * $min_overlap_ratio) and ($pre_end - $pos >= $len * $min_overlap_ratio)) or (($len >= $pre_len) and ($len < $pre_len * 5)) or (($len < $pre_len) and ($pre_len < $len * 5))){
                                my @pre_line = split (/\s+/, $pre_line);
                                my $len2 = $end - $pre_pos + 1;
                                $pre_line[7] =~ s/SVLEN=\d+/SVLEN=$len2/;
                                my $pre_tool = $1 if ($pre_line[7] =~ /TOOLS=(\S+)/);
                                my ($new_tool) = &rm_overlap ($pre_tool, $tool);
                                $pre_line[7] =~ s/TOOLS=[^;]+/TOOLS=$new_tool/;
                                $tool_num = $new_tool =~ s/,/=/g;
                                my $new_line = join ("\t", @pre_line);
                                ${${$vcf{$chr02d}}{$pre_pos}}{$pre_type} = $new_line;
                                $pre_line = $new_line;
                                $pre_len = $len2;
                                $pre_end = $end;
                                $pre_tool_num = $pre_tool_num + $tool_num;
                                next;
                            }
                        }
                    }
                    elsif ($pre_type eq 'DEL'){
                        my $overlap = $pre_end - $pos + 1;
                        $overlap = $len if ($end < $pre_end);
                        if (($overlap < $pre_len * $min_overlap_ratio) and ($overlap < $len * $min_overlap_ratio)){
                            delete ${$vcf{$chr02d}}{$pos};
                        }
                    }
        	        elsif ($pre_type eq 'INS'){
            			if ($pos - $pre_pos <= $ins_sd2){
            			    next;
            			}
            			elsif ($end -$pre_pos <= $ins_sd2){
            			    next;
            			}
            		}
                }
                elsif ($type eq 'INV'){
                    if ($pre_type eq 'INV'){
                        if ($pre_end >= $end){
                            next;
                        }
                        elsif (($pre_end - $pos >= $pre_len * $min_overlap_ratio) or ($pre_end - $pos >= $len * $min_overlap_ratio)){		# merge pre-INV + INV -> new-INV
                            if ((($pre_end - $pos >= $pre_len * $min_overlap_ratio) and ($pre_end - $pos >= $len * $min_overlap_ratio)) or (($len >= $pre_len) and ($len < $pre_len * 5)) or (($len < $pre_len) and ($pre_len < $len * 5))){
                                my @pre_line = split (/\s+/, $pre_line);
                                my $len2 = $end - $pre_pos + 1;
                                $pre_line[7] =~ s/SVLEN=\d+/SVLEN=$len2/;
                                my $pre_tool = $1 if ($pre_line[7] =~ /TOOLS=(\S+)/);
                                my ($new_tool) = &rm_overlap ($pre_tool, $tool);
                                $pre_line[7] =~ s/TOOLS=[^\;]+/TOOLS=$new_tool/;
                                $tool_num = $new_tool =~ s/,/=/g;
                                my $new_line = join ("\t", @pre_line);
                                ${${$vcf{$chr02d}}{$pre_pos}}{$pre_type} = $new_line;
                                $pre_line = $new_line;
                                $pre_len = $len2;
                                $pre_end = $end;
                                $pre_tool_num = $pre_tool_num + $tool_num;
                                next;
                            }
                        }
                    }
                }

            }
            ${${$vcf{$chr02d}}{$pos}}{$type} = $line;
            $pre_line = $line;
            $pre_chr = $chr;
            $pre_pos = $pos;
            $pre_type = $type;
            $pre_subtype = $subtype;
            $pre_len = $len;
            $pre_end = $end;
            $pre_tool_num = $tool_num;
        }
    }
}

my @mei_header;
if (scalar keys %mei_type > 0){
    foreach my $mei (sort keys %mei_type){
        my $line = "##ALT=<ID=INS:ME:$mei,Description=\"Insertion of $mei element\">";
        push @mei_header, $line;
    }
}
if ($numt == 1){
    my $line = "##ALT=<ID=INS:NUMT,Description=\"Insertion of nuclear mitochondrial genome\">";
    push @mei_header, $line;
}
if ($vei == 1){
    my $line = "##ALT=<ID=INS:VEI,Description=\"Insertion of viral genome element\">";
    push @mei_header, $line;
}

my $version = 1.7;
$version = $1 if ($Bin =~ /MOPline[_\-\.](v[\d\.]+)/i);

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

print "##fileformat=VCFv4.0\n";
print "##fileDate=$time\n";
print "##source=MOPline-$version\n";
print "##reference=\n";
print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation (SV)\">\n";
print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles (0 when undefined)\">\n";
print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV\">\n";
print "##INFO=<ID=READS,Number=1,Type=Integer,Description=\"Number of reads supporting the SV allele\">\n";
print "##INFO=<ID=TOOLS,Number=A,Type=String,Description=\"SV calling tool-called position\">\n";
print "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print "##ALT=<ID=DUP,Description=\"Duplication\">\n";
print "##ALT=<ID=INS,Description=\"Insertion of undefined sequence\">\n";
if (@mei_header > 0){
    foreach my $hline (@mei_header){
        print "$hline\n";
    }
}
print "##ALT=<ID=INV:Description=\"Inversion\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

foreach my $chr (sort keys %vcf){
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}}){
        foreach my $type (keys %{${$vcf{$chr}}{$pos}}){
            my @line = split (/\t/, ${${$vcf{$chr}}{$pos}}{$type});
            my $type = $line[2];
            my $subtype = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
            my $subtype2 = '';
            if ($type ne $subtype){
                my $stype = $type;
                $type = $subtype;
                $subtype = $stype;
                if ($subtype =~ /^(.+?),([^,]+)/){
                    $subtype = $1;
                    $subtype2 = $2;
                }
            }
            my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
            $len = 0 - $len if ($type eq 'DEL');
            my $alt = "<$type>" if ($type eq 'DEL') or ($type eq 'DUP') or ($type eq 'INV');
            if ($type eq 'INS'){
                $alt = '<INS' if ($type eq 'INS') and ($subtype eq 'INS');
                $alt = "<INS:ME:$subtype" if ($subtype ne $type) and ($subtype !~ /NUMT|VEI/);
                $alt = "<INS:$subtype" if ($subtype ne $type) and ($subtype =~ /NUMT|VEI/);
                if ($subtype2 ne ''){
                    $alt .= ",INS:MEI:$subtype2" if ($subtype !~ /NUMT|VEI/);
                    $alt .= ",INS:$subtype2" if ($subtype =~ /NUMT|VEI/);
                }
                $alt .= '>';
            }
        
            $line[7] =~ s/SVTYPE=[^\;]+/SVTYPE=$type/;
            $line[7] =~ s/SVLEN=\d+/SVLEN=$len/ if ($type eq 'DEL');
            $line[2] = '.';
            $line[4] = $alt;
            print join ("\t", @line), "\n";
        }
    }
}


sub rm_overlap{
    my ($tool1, $tool2) = @_;
    my $tool3 = $tool1 . ',' . $tool2;
    my @tool3 = split (/,/, $tool3);
    my %find;
    map {$find{$_}++} @tool3;
    my $new_tool = '';
    foreach my $tool (keys %find){
	$new_tool .= $tool . ',';
    }
    $new_tool =~ s/,$//;
    return ($new_tool);
}
