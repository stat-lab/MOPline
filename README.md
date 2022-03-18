# MOPline
Detection and genotyping of structural variants
## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Composition of Required Directories](#composition)
- [Human and Non-human Data](#human)
- [General Usage](general)
	- [\[Step-0\] Run SV detection tools](#step-0)
		- [Notes on SV calling and input bam](#note)
	- [\[Step-1\] Select overlap calls \(high-confidence calls\) from SV call sets(#step-1)
	- [\[Step-2\] Add alignment statistics to SV sites](#step-2)
	- [\[Step-3\] Merge vcf files from multiple samples \(joint-call\)(#step-3)
	- [\[Step-4\] Genotype and SMC(#genotype)
	- [\[Step-5\] Annotate(#annotate)
	- [\[Step-6\] Filter(#filter)
- [Quick Start with Sample Data](#quick)
	- [A: Human SV datasets of 6 samples](#A:)
	- [B: Yeast WGS bam files of 10 samples](#B:)

## Introduction

MOPline accurately and sensitively detects structural variations (SVs) in whole genome sequencing (WGS) data from a single to thousands (and sometinmes tens of thousands) of human or non-human samples. MOPline first selects overlap SV calls (high-quality calls) of each SV type and size-ranges from SV call sets to generate a single vcf file for a sample. The SV call sets are multiple vcf files generated with existing SV detection algorithms (tools). The package is pre-programed to select overlap calls with 6~9 algorithms selected by us, but the user is free to choose any combinations of algorithms. In the second step, alignment statistics on coverage and split reads are added to the SV sites in the vcf files for each sample, providing materials for the next SV genotyping. The third step involves a joint call to merge the vcf files of multiple samples. In the fourth step, genotyping of each SV allele is performed based on multinominal logistic regression using the alignment statistics. In this step, reference alleles are also genotyped with the method described above to recover SVs that were missed due to the selection of SVs with high confidence in the first step. We call this recovery method Supplementing Missing Call (SMC). Finally, SVs are filtered and/or annotated for gene-overlapping SVs.

## Requirements

perl 5 or later versions
[samtools](https://github.com/samtools/samtools)
[vcftools](https://vcftools.github.io/index.html)
java 1.8 or later (for GRIDSS and MELT)

**SV detection tools (algorithms)**
The tools required depend on the tool presets used in MOPline or the tool sets customized by the user. The presets customized in this MOPline package are follows:
**Preset: 7tools (MOPline-7t)**
[CNVnator](https://github.com/abyzovlab/CNVnator)
[GRIDSS](https://github.com/PapenfussLab/grids)
[inGAP](https://sourceforge.net/projects/ingap/files/ingap)
[Manta](https://github.com/Illumina/manta)
[MATCHCLIP](https://github.com/yhwu/matchclips2)
[MELT](https://melt.igs.umaryland.edu)
[Wham](https://github.com/zeeev/wham)
**Preset: 6tools_1 (MOPline-6t-1)**
GRIDSS is excluded from 7tools
**Preset: 6tools_2 (MOPline-6t-2)**
MELT is excluded from 7tools
**Preset: 9tools (MOPline-9t)**
DELLY, Lumpy, SoftSV are added to 6tools_1

**Sample data**
[Sample datasets](http://jenger.riken.jp/en) include human SV call sets from 6 individuals and yeast 10 WGS data. The datasets also include output data created with MOPline.

## Composition of required directories (sample directory and tool directory)

MOPline assumes a directory structure of sample_directory/tool_directory under the working directory, where the sample directory has the sample name or sample ID, Under the sample directory, there are the tool directories with the names of algorithms, such as Manta and inGAP, which is case-sensitive. Under each sample directory, there are also should be bam and its index files. For convenience, when running against a sample (sample name: Sample123), the working directory should contain a Sample123 directory, under which the Sample123.bam and Sample123.bam.bai or their symbolic links should exists. When running with the 7tool preset, seven tool directories (CNVnator, GRIDSS, inGAP, Manta, MATCHCLIP, MELT, and Wham) must exist under the Sample123 directory. In the tool directories, the runs of the corresponding tools are performed.

## Human and Non-human Data

By default, MOPline treats WGS alignment data (bam/cram) and SV call data (vcf) generated based on the human build 37 reference (GRCh37 or GRCh37d5). When using the data based on the human build 38 reference, run several MOPline scripts using the option ‘-build 38’. For non-human species, run several MOPline scripts using the option ‘-nh 1’ and, in some cases, also using other options specifying some reference-specific annotation files for gap, repeat, and gene regions.

## General Usage

### [Step-0] Run SV detection tools

The preset algorithms are a combination of tools that have been evaluated for precision and recall, but the user may use a different combination or use other algorithms not listed here. The final output file from each algorithm must be converted to a MOPline- compatible vcf file using the conversion script in ‘run_SVcallers’ folder of this package. The converted vcf file should contain ‘SVTYPE’, ‘SVLEN’, and ‘READS’ keys in the INFO field, where the READS key represents RSS (reads supporting SV). The conversion scripts, the corresponding algorithms, and their output files are indicated in the table below.

**Table 1.** Conversion scripts for 10 SV detection algorithms
|Algorithm|Command using a conversion script (for a sample name, AB)                      |
| :------ | :---------------------------------------------------------------------------- |
|CNVnator |convert_CNVnator_vcf.pl AB.out gap.bed > CNVnator.AB.vcf                       |
|DELLY    |convert_CNVnator_vcf.pl AB.out > CNVnator.AB.vcf                               |
|GRIDSS   |convert_GRIDSS_vcf.pl AB.vcf.gz > GRIDSS.AB.vcf                                |
|inGAP    |convert_inGAP_vcf.pl AB.chr1.out .... AB.chrY.out > inGAP.AB.vcf               |
|Lumpy    |convert_Lumpy_vcf.pl AB.vcf > Lumpy.AB.vcf                                     |
|Manta    |convert_Manta_vcf.pl results/variants/diploidSV.vcf.gz > Manta.AB.vcf          |
|MATCHCLIP|convert_MATCHCLIP_vcf.pl AB.matchclip.out > MATCHCLIP.AB.vcf                   |
|MELT     |convert_MELT_vcf.pl ALU.final_comp.vcf LINE1.final_comp.vcf SVA.final_comp.vcf <br>HERVK.final_comp.vcf > MELT.AB.vcf|
|SoftSV   |convert_SoftSV_vcf.pl deletions_small.txt insertion_small.txt tandems_small.txt<br> inversions_small.txt deletions.txt tandems.txt inversions.txt > SoftSV.AB.vcf
|Wham     |convert_Wham_vcf.pl AB.vcf > Wham.AB.vcf                                       |

In the case of CNVnator, the second argument of the convert_CNVnator_vcf.pl script must be a gap bed file that indicates the gap regions (a stretch of ‘N’ bases) in the reference genome (gap.bed for human is located in the Data folder in the package, but this file can be omitted for non-human).

For convenience, we have provided scripts to run the selected algorithms in the ‘run_SVcallers’ folder. A run the script with the ‘-h’ option informs us the required arguments and options. We also provide a script for sequential execution of multiple algorithms for a single sample (run_single.pl) and batch scripts for multiple samples using Slurm and LSF job managers (run_batch_slurm.pl, run_batch_LSF.pl). To run these three scripts, specify a configure file that describes the parameters for each algorithm, using a template configure file (config.txt) (an example config file is also provided in the sample data [Sample_data_output/yeast_run]). A tool directory with the name of the algorithm specified in the configure file is automatically created in the sample directories under the working directory, and the algorithm is executed under the tool directory.
