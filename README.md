# MOPline
Detection and Genotyping of Structural Variants
## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Citation](#citation)
- [Install](#install)
- [Composition of Required Directories](#composition)
- [Human and Non-human Data](#hdata)
- [General Usage](#gusage)
	- [\[Step-0\] Preparation of input files](#step0)
		- [(1) Run SV detection tools](#run_sv)
			- [(a) Run user-installed tools using wrapper scripts](#user_install)
			- [(b) Run using Singularity](#singularity)
   			- [(c) Batch runs of SV callers with multiple samples](#batch)
			- [Notes on SV calling and input bam](#notes)
		- [(2) Create coverage files](#create_cov)
	- [\[Step-1\] Select overlap calls \(high-confidence calls\) from SV call sets](#step1)
	- [\[Step-2\] Add alignment statistics to SV sites](#step2)
	- [\[Step-3\] Merge vcf files from multiple samples \(joint-call\)](#step3)
	- [\[Step-4\] Genotype and SMC](#step4)
	- [\[Step-5\] Annotate](#step5)
	- [\[Step-6\] Filter](#step6)
- [Data Required for Non-Human Species](#nhdata)
- [For single sample](#single)
- [Quick Start with Sample Data](#quick)
	- [A: Human SV datasets of 6 samples](#hsample)
	- [B: Yeast WGS bam files of 10 samples](#ysample)

## Introduction

MOPline accurately and sensitively detects structural variations (SVs) in whole genome sequencing (WGS) data from a single to thousands (and sometimes tens of thousands) of human or non-human samples. MOPline first selects overlap SV calls (high-quality calls) of each SV type and size-ranges from SV call sets to generate a single vcf file for a sample. The SV call sets are multiple vcf files generated with existing SV detection algorithms (tools). The package is pre-programed to select overlap calls with 6~9 algorithms selected by us, but the user is free to choose any combinations of algorithms. In the second step, alignment statistics on coverage and split reads are added to the SV sites in the vcf files for each sample, providing materials for the next SV genotyping. The third step involves a joint call to merge the vcf files of multiple samples. In the fourth step, genotyping of each SV allele is performed based on multinominal logistic regression using the alignment statistics. In this step, reference alleles are also genotyped with the method described above to recover SVs that were missed due to the selection of SVs with high confidence in the first step. We call this recovery method Supplementing Missing Call (SMC). Finally, SVs are filtered and/or annotated for gene-overlapping SVs.

## Requirements

perl 5 or later versions  
R 3.0 or later (required library: [nnet](https://cran.r-project.org/web/packages/nnet))  
[samtools](https://github.com/samtools/samtools) (v1.13 or later is required for GRIDSS v2.13-)  
[vcftools](https://vcftools.github.io/index.html) (v0.1.15 or later)  
java 1.8 or later (for GRIDSS and MELT)

#### General input file
- Reference fasta file & its index file  
- Alignment bam/cram file & its index file (only bam file is permitted for several SV detection tools)  
- Gender list file (optional, 1st column: sample name, 2nd column: M|F, separated by tab)

#### SV detection tools (algorithms)  
The tools required depend on the tool presets used in MOPline or the tool sets customized by the user. The necessary SV detection tools must be installed by the user, but it is possible to use tools other than MELT using the singularity definition file provided in this package, as described below.  
The presets customized in this MOPline package are follows:  
- **Preset: 8tools (MOPline-8t)**  
	- [CNVnator](https://github.com/abyzovlab/CNVnator)  
	- [GRIDSS](https://github.com/PapenfussLab/gridss)  
	- [inGAP](https://sourceforge.net/projects/ingap/files/ingap)  
	- [Manta](https://github.com/Illumina/manta)  
	- [MATCHCLIP](https://github.com/yhwu/matchclips2)  
	- [MELT](https://melt.igs.umaryland.edu)  
	- [Wham](https://github.com/zeeev/wham)
 	- [INSurVeyor](https://github.com/kensung-lab/INSurVeyor)  
- **Preset: 7tools (MOPline-7t)**  
	- INSurVeyor is excluded from 8tools
- **Preset: 7tools_1 (MOPline-7t-1)**  
	- GRIDSS is excluded from 8tools  
- **Preset: 7tools_2 (MOPline-7t-2)**  
	- MELT is excluded from 8tools  
- **Preset: 6tools_1 (MOPline-6t-1)**  
	- GRIDSS and INSurVeyor are excluded from 8tools  
- **Preset: 6tools_2 (MOPline-6t-2)**  
	- MELT and INSurVeyor are excluded from 8tools  
- **Preset: 9tools (MOPline-9t)**  
	- DELLY, Lumpy, SoftSV are added to 6tools_1

## Citation

Please cite the following paper when using MOPline for your publication.  

Kosugi et al. Detection of trait-associated structural variations using short-read sequencing. Cell Genomics 3, 100328 (2023). https://doi.org/10.1016/j.xgen.2023.100328

## Install

```
git clone https://github.com/stat-lab/MOPline
```
The Data folder in the MOPline folder contains parameter files, multinomial logistic regression-based model files for genotyping (R.nnet.models), and annotated data files for the human build 37/38/T2T reference. Do not change the name of the files/directories (except config.txt) and the directory structure in the MOPline folder.  
  
MOPline and additional 10 SV detection tools (except for MELT and INSurVeyor) can be executed using [Singularity](#singularity).

### Sample data
[Sample datasets](http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/Sample_data.tar.gz) (or available from https://drive.google.com/drive/folders/1bIEtaaM3xx8POIAf96kV-ImNXTHWPwQQ?usp=sharing) include human SV call sets from 6 individuals and yeast 10 WGS data. The datasets also include output data created with MOPline.

## <a name="composition"></a>Directory structure required (sample directory and tool directory)

MOPline assumes a directory structure of sample_directory/tool_directory under the working directory, where the sample directory has the sample name or sample ID, Under the sample directory, there are the tool directories with the names of algorithms, such as Manta and inGAP, which is case-sensitive. Under each sample directory, there are also should be bam and its index files. For convenience, when running against a sample (sample name: Sample_ID1), the working directory should contain a Sample_ID1 directory, under which the Sample_ID1.bam and Sample_ID1.bam.bai or their symbolic links should exist. When running with the 7tools preset, seven tool directories (CNVnator, GRIDSS, inGAP, Manta, MATCHCLIP, MELT, and Wham) must exist under the Sample_ID1 directory. In the tool directories, the runs of the corresponding tools are performed. In addition, Cov and Merge_7tools folders will be created under each sample directory in Step-0 and Step-1, as described below.  
![image](https://github.com/stat-lab/MOPline/files/10089842/Dir_Struc.pdf)

## <a name="hdata"></a>Human and Non-human Data

By default, MOPline treats WGS alignment data (bam/cram) and SV call data (vcf) generated based on the human build 37 reference (GRCh37 or GRCh37d5). When using the data based on the human build 38 reference or T2T-CHM13.v2.0, run several MOPline scripts using the option ‘--build 38’ or '--build T2T. By specifying the --build option, specific data files associated with the human reference will be automatically selected from the Data folder. For non-human species, run several MOPline scripts using the option ‘-nh 1’ and, in some cases, also using other options specifying some reference-specific annotation files for gap, repeat, and gene regions ([for detail](#nhdata)).

## <a name="gusage"></a>General Usage

### <a name="step0"></a>[Step-0]  Preparation of input files

### <a name="run_sv"></a>(1) Run SV detection tools

The preset algorithms are a combination of tools that have been evaluated for precision and recall, but the user may use a different combination or use other algorithms not listed here. The final output file from each algorithm must be converted to a MOPline-specific vcf file using the conversion script in ‘run_SVcallers’ folder of this package. The converted vcf file should contain ‘SVTYPE’, ‘SVLEN’, and ‘READS’ keys in the INFO field, where the READS key represents RSS (reads supporting SV). Because many short read-based SV detection algorithms, including MOPline, cannot determine the length of INSs, the length of the INS is recorded as 0 or 1 in many INS calls. The conversion scripts, the corresponding algorithms, and their output files are indicated in the table below.

**Table 1.** Conversion scripts for 10 SV detection algorithms
|Algorithm|Command using a conversion script (for a sample name, AB)                      |
| :------ | :---------------------------------------------------------------------------- |
|CNVnator |convert_CNVnator_vcf.pl AB.out gap.bed > CNVnator.AB.vcf                       |
|DELLY    |convert_DELLY_vcf.pl AB.out $nh > DELLY.AB.vcf                               |
|GRIDSS   |convert_GRIDSS_vcf.pl AB.vcf.gz $nh > GRIDSS.AB.vcf                                |
|inGAP    |convert_inGAP_vcf.pl AB.chr1.out .... AB.chrY.out > inGAP.AB.vcf               |
|Lumpy    |convert_Lumpy_vcf.pl AB.vcf $nh > Lumpy.AB.vcf                                     |
|Manta    |convert_Manta_vcf.pl results/variants/diploidSV.vcf.gz $nh > Manta.AB.vcf          |
|MATCHCLIP|convert_MATCHCLIP_vcf.pl AB.matchclip.out $nh > MATCHCLIP.AB.vcf                   |
|MELT     |convert_MELT_vcf.pl ALU.final_comp.vcf LINE1.final_comp.vcf SVA.final_comp.vcf <br>HERVK.final_comp.vcf > MELT.AB.vcf|
|SoftSV   |convert_SoftSV_vcf.pl deletions_small.txt insertion_small.txt tandems_small.txt<br> inversions_small.txt deletions.txt tandems.txt inversions.txt $nh > SoftSV.AB.vcf
|Wham     |convert_Wham_vcf.pl AB.vcf $nh > Wham.AB.vcf                                       |
|INSurVeyor|convert_INSurVeyor_vcf.pl AB.vcf $nh > INSurVeyor.AB.vcf                       |

If the data is for non-human species, specify 1 for $nh. If the data is human, specify 0 for $nh or nothing. In the case of CNVnator, the second argument of the convert_CNVnator_vcf.pl script must be a gap bed file that indicates the gap regions (a stretch of ‘N’ bases) in the reference genome. The gap.bed file for human is located in the Data folder in the package. For non-human species, a gap file can be obtained at [UCSC](https://hgdownload.soe.ucsc.edu/downloads) for some species or created manually, but this file can be omitted.

### <a name="user_install"></a>(a) Run user-installed tools using wrapper scripts

For convenience, we have provided wrapper scripts to run the selected SV detection tools in the ‘run_SVcallers’ folder. The tools are expected to be properly installed by the user. A run of the script with the ‘-h’ option informs us the required arguments and options. We also provide a script for sequential execution of multiple algorithms for a single sample (run_single.pl) and batch scripts for multiple samples using Slurm and LSF job managers (run_batch_slurm.pl, run_batch_LSF.pl). To run these three scripts, specify a configure file that describes the parameters for each algorithm, using a template configure file (MOPline/scripts/run_SVcallers/config.txt) (an example config file is also provided in the sample data [Sample_data_output/yeast_run]). A tool directory with the name of the algorithm specified in the configure file is automatically created in the sample directories under the working directory, and the algorithm is executed under the tool directory.  
  
Example using run_single.pl  
```
run_single.pl -b <input bam> -c <config file> -sd <sample directory>
```
This command creates SV calling results from the SV detection algorithms specified in the config file, and the results are in the tool directories under the sample directory specified with the -sd option.  

### <a name="singularity"></a>(b) Run using Singularity

Other than MELT and INSurVeylor (MELT has license restrictions) the above 9 SV detection algorithms can be run using a Singularity image container built using the Definition file (mopline.def) included in the package. We expect to be run INSurVeylor using Singularity image file (https://github.com/kensung-lab/INSurVeyor/releases) specific to INSurVeylor or installed using conda   
An image file (e.g., mopline.sif) can be obtained using one of the following:  
```
sudo singularity build mopline.sif mopline.def
```
```
singularity build --fakeroot mopline.sif mopline.def
```
Alternatively, a mopline.sif (v1.8.2 1.4 GB) can be obtained at [Jenger site](http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/mopline.sif) :(http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/mopline.sif).  
  
To run CNVnator with singularity  
```
singularity exec <mopline.sif file path> /opt/local/tools/MOPline/scripts/run_SVcallers/run_CNVnator.pl -cp <cnvnator command path> -b <input bam> -r <ref directory> -rd <ROOT install directory> -p <output prefix>   # the path of run_CNVnator.pl is the path in the singularity container.
```
To run the run_single_singularity.pl script for a single sample, install this MOPline package on the appropriate location and run as follows:  
```
~/MOPline_v1.#/scripts/SV_callers/run_single_singularity.pl -s <mopline.sif file path> -b <input bam> -c <config_singularity.txt> -td <temp direcory> -bd <bind directories> -sd <sample directory> 
```
The sif file must be specified as an absolute path. The config_singularity.txt provided in this package can be used for the template, but the lines 'Do not change' indicating the path in the singularity container, must not be changed. Input files (reference.fasta, alignment bam) on the host must be specified with the absolute path. The temp directory corresponds to /tmp on the host or a temporary directory specified with the environmental variable TMPDIR. The bind directories are the paths on the host to add to the singularity container. If the input files used by run_single_singularity.pl are located in a directory other than the working directory, the path to the corresponding directory must be specified with the -bd option (multiple paths separated by commas).

### <a name="batch"></a>(c) Batch runs of SV callers with multiple samples

Users who have access to a job management system that includes Slurm and LSF Job Manager can use batch scripts (run_batch_slurm.pl, run_batch_LSF.pl) for multiple samples.  
If you do not use the job management system, the script run_SVcallers_batch.pl can be used for multiple samples. To use this script, split multiple samples into several batches (e.g., 10 samples x 10 batches for 100 samples), depending on the CPU cores and memory capacity of the computational environment. For each batch、create sample list files showing the sample name and the corresponding bam/cram file path in each line. If the input alignment files are in cram format, the script automatically converts cram files to bam files because some SV callers only accept bam format as input. The --add_MC (-am) option can be used to add MC/MQ tags to bam files, which may be useful for reducing the runtime of the MELT tool. In these cases, a new bam file will be created under the sample directory. A typical usage of this script is as follows:
```
run_SVcallers_batch.pl -sl <sample list file> -r <reference fasta> -sf <mopline.sif singularity file if -c is specified> -c <config file for singularity-based run> -td <temp direcory> -bd <bind directories> -c2 <config file for MELT and/or INSurVeyor> (-am if adding MC/MQ tag in bam file)
(Sample directories must exist under the working directory.)
```
This command executes the SV callers for each sample sequentially. Multiple parallel executions of this command with the split sample list files will reduce the overall runtime. Both the -c and -c2 options can be specified, or just one of them. If the mopline.sif singularity container is not used and a user-installed SV caller is used, define the SV caller parameters in the configuration file specified with -c2.

### <a name="notes"></a>Notes on SV calling and input bam

#### inGAP-sv
inGAP requires sam files and reference fasta files split by chromosome as input files. inGAP is highly sensitive to PCR-duplicates in input sam files. When a sam file containing 10% PCR-duplicates is used, inGAP often generates more than twice as many calls compared with sam without PCR-duplicates.
#### GRIDSS
In GRIDSS, accidental errors can occur, primarily due to the java heap memory size allocations. To solve this problem, it is often set the java maximum heap size (export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -Xmx15g -Xms15g") and the GRIDSS options, ‘15g’ for both, --jvmheap and --otherjvmheap. It is advisable to link the reference fasta and bwa files in the GRIDGSS working directory. GRIDSS is very sensitive to reads with low-quality terminal regions, so it is recommended to trim off low-quality regions of reads before alignment if may fractions of reads contain low-quality regions.
#### Manta
When an input bam is generated with a reference containing many decoy sequences such as GRCh38+decoy, it often takes longer time to complete run. In this case, runtime can be shortened by using the –callRegions option, which specifies a bgzip-compressed bed file indicating only target chromosomes (i.e., chr1, … chrY).
#### MELT
MELT often takes much longer, especially for GRCh38-based alignment data, without MC/MQ tag in the input bam. Also, when using a reference that contains many decoy sequences such as GRCh38+decoy, the -b option can be used a list of chromosomes to exclude (> 1 Mb chr) to reduce execution time. For non-human species, MELT can be run by preparing species-specific mobile element data, as described in the MELT documentation.
#### INSurVeylor
INSurVeylor needs an input BAM/CRAM file containing MC and MQ tags, which can be added to BAM/CRAM using samtools fixmate or picard.jar FixMateInformation command.

#### Input bam
Trimming the read sequence in the fastq file (e.g. using Trimomatic) or removing PCR duplicates in the bam file is recommended, especially if the reads are of low quality. If there are 'N' bases or low quality bases in the terminal regions of the reads, these reads often generate soft clipping alignments in the bam file, which can cause some algorithms to call the wrong breakpoints. Also, when PCR-duplicate reads are present, some algorithms may call wrong SVs. GRIDSS is an SV detection algorithm that is highly sensitive to split-read alignments. if GRIDSS generates tens of thousands of SV calls and takes several days to complete execution (usually completed within a day for 30x human WGS data), it would be an indication that the sequencing data contains many low-quality reads.

The samtools fixmate command can be used to add MC and MQ tags, which is useful for reducing MELT run time, running INSurVeylor, and removing PCR duplicates. Example commands using bwa and samtools fixmate are shown below.
(bwa → bam with MC/MQ tags)
```
bwa mem <ref index> <read_1> <read_2> -M -R “read group header specification” -t <num threads> | samtools view -uShb - | samtools fixmate -um - [out.fm.bam]
```
(sort by coordinate)
```
samtools sort -@ 4 <out.fm.bam> -o [out.sr.bam]
```
(sort by coordinate & remove PCR-duplicates)
```
samtools sort -@ 4 <out.fm.bam> | samtools markdup -ur - [out.rm.sr.bam]
```
### <a name="create_cov"></a>(2) Create coverage files with each sample bam file

In Step-2, MOPline uses coverage files recording read depth and the number of soft-clipped ends for each 50-bp window for each sample.  
Create a coverage file using the create_coverage_file_bam.pl script as follows:
```
mopline create_cov -b <bam_list> -r <reference_fasta> -rl <read_length> -n <num_threads>
```
(-nh 1 if sample is a non-human species)  
**bam_list:** bam/cram list file specifying bam or cram file name per line. The bam_list can also be a list of sample names if the bam file name is `${sample_name}.bam` and exists in the `${sample_name}` directory.  
**read_length:** Mean read length in the bam file

The above command creates a Cov directory under the sample directory, which contains the coverage files for each chromosome (`${sample_name}.chr*.cov.gz`).

For batch jobs, we provide create_coverage_file_bam_single.pl and run_create_cov_batch.pl scripts, the former of which can be used to submit a single bam file job using a job manager such as Slurm and LSF. The run_create_cov_batch.pl script can be used with sample list files, which were split into several batches (e.g., 10 samples x 10 batches for 100 samples), depending on the CPU cores and memory capacity of the computational environment. A sample list file describes sample name in each line. A typical usage of this script is as follows:
```
run_create_cov_batch.pl -s <sample list file> -r <reference fasta file> -rl <read length (bp)>
(Sample directories must exist under the working directory and each sample directory must have a bam file, named ${sample_name}.bam.)
```
This command executes the script create_coverage_file_bam_single.pl sequentially for each sample with a single thread. Multiple parallel executions of this command with the split sample list files will reduce the overall runtime.

### <a name="step1"></a>[Step-1]  Select overlap calls (high-confidence calls) from SV call sets

Once the MOPline-specific vcf files for each algorithm and for each sample are created, selection of overlapping SV calls is first performed. Use the merge_SV_vcf.*tools.pl script in the scripts folder to select overlap calls and high-confidence calls from the SV call sets (vcf files) generated in Step-0 and to merge the selected calls of each SV type and size-ranges for each sample. For high coverage (30x or more) bam files, use mopline subcommands, merge_7tools, merge_6tools_1, merge_6tools_2, or merge_9tools. When using bam files with ~20x coverage, use the subcommand with a ‘_lc’ at the end.
```
mopline merge_7tools -s <sample-list-file or a sample name> -rl <read length> 
```
(-nh 1 if sample is a non-human species)

This command generates a `${sample_name}.Merge.ALL.vcf` file in the ‘Merge_7tools’ folder (default for 7tools preset) created under the sample directory. The INFO field of each SV line displays a TOOLS key indicating which algorithm called the corresponding SV.

### [Using custom algorithm set]

When using a custom algorithm set instead of the preset, users can create their own merge script (corresponding to the merge_SV_vcf.7tools.pl script) using the make_merge_SV_vcf_script.pl script as follows:
```
make_merge_SV_vcf_script.pl -t <algorithm list, comma-separated> -tc <tool-config file> -p <prefix name of output script>
```
(use -h for detailed explanation)  
If the -tc option is not specified in the above command, the tool configuration file (Data/SVtool_pairs_config.txt) is automatically selected. This file specifies the favorable pairs of tools and minimum RSSs of overlap call selection for each of 14 algorithms we have chosen. If additional algorithms are used, the SVtool_pairs_config.txt file can be modified to specify the preferred pairs minimum RSSs for each newly added algorithm for each SV type and size range.

### <a name="step2"></a>[Step-2]  Add alignment statistics to SV sites

In this step, alignment statistics such as DPR, SR, and DPS are added to each SV site in every sample. DPR is the ratio of the depth of the region inside the SV to the adjacent depth, and DPS is the deviation rate of DPR measured in a 50-bp window. SR is the ratio of soft-clipped read ends around the breakpoint to the outside area. To measure these values, a coverage file must first be created for each sample (see [create coverage files](#create_cov)), recording the read depth and the number of soft-clipped ends for each 50-bp window.

Using the coverage files you created, add alignment statistics to each SV site in the `${sample_name}.Merge.ALL.vcf` file created in Step-1. This step also adds reliable genotypes from the genotypes called from several algorithms in Step-0 to some of the SV sites. The command with the add_GT_DPR_vcf.pl script is as follows:
```
mopline add_cov -s <sample_list> -ts <tool_set> -vd <vcf_directory> -n <num_threads> 
```
(--build 38 for human build38, -nh 1 -r <ref.index> -gap <gap_bed> for non-human species)  
**sample_list:** A sample list file showing sample names per line. A sample directory with the same names as the specified sample list must exist under the working directory.  
**tool_set:** Algorithm preset or list file showing algorithm names per line [default: 7tools]  
**vcf_directory:** The name of the directory containing the input vcf files in the sample directories [default: Merge_7tools]

The above command updates a `${sample_name}.Merge.ALL.vcf` file in the vcf_directory and rename the original vcf file as `${sample_name}.Merge.ALL.noAdd.vcf`. Alignment statistics are added to the FORMAT/SAMPLE fields with DR, DS, and SR keys in the updated vcf file.

For batch jobs, we provide a add_GT_DPR_vcf_single.pl script to submit a single vcf file job using a job manager such as Slurm and LSF.

### <a name="step3"></a>[Step-3]  Merge vcf files from multiple samples (joint-call)

Joint calling is performed with the merge_SV_calls_ALLtype.pl script as follows:
```
mopline joint_call -s <sample_list> -md <merge_dir> -od <out_dir> -p <out_prefix>
```
(--build 38 for human build 38, -nh 1 -gap <gap_bed> for non-human species)  
**sample_list:** A sample list file showing sample names per line. A sample directory with the same name as the specified sample list must exist under the working directory. If the sample directory name and the sample name are different, specify them by separating each line with a comma, such as `${sample_directory_name},${sample_name}`. It is assumed that the input vcf files (`${sample_name}.Merge.ALL.vcf`) exist in `${sample_directory}/${merge_dir}`.  
**merge_dir:** Name of the directory containing the input vcf files under the sample directories [default: Merge_7tools]  
**out_dir:** Name of the directory where the output vcf file will be generated [default: the same name as merge_dir]  
**out_prefix:** Prefix name of the output vcf file  
**gap_bed:** A bed file of reference gap regions. For human, the file is automatically selected from the MOPline package. For non-human species, obtain from [UCSC](https://hgdownload.soe.ucsc.edu/downloads) or create manually.

The above command outputs `${out_prefix}.vcf` under `${out_dir}` directory.

### <a name="step4"></a>[Step-4]  Genotype and SMC

In this step, all SV alleles are genotyped based on multinominal logistic regression (model data are in the Data/R.nnet.models folder) and reference alleles are re-genotyped to recover missing SV calls by SMC. This step requires files showing repeat regions in the reference, which are provided for human (Data/simpleRepeat.txt.gz, Data/genomicSuperDups.txt.gz) and are automatically selected.  
Before the run of this step, confirm whether the R nnet package has been installed as follows:  
```
R
>library(nnet)
```
The command using the genotype_SV_SMC_7.4.pl script is as follows:
```
mopline smc -v <input_vcf> -ts <tool_set> -od <out_dir> -p <out_prefix> -n <num_threads>
```
(--build 38 for human build 38, -nh 1 -sr <STR_file> -sd <SD_file> -r <ref_index> for non-human species)  
**input_vcf:** An input vcf file from Step-3  
**tool_set:** Algorithm preset name or a list file showing algorithm names per line [default: 7tools]  
**out_dir:** Name of the directory where the output vcf file will be generated  
**out_prefix:** Prefix name of the output vcf file  
**STR_file:** A simple/short tandem repeat file from [UCSC](https://hgdownload.soe.ucsc.edu/downloads) or [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) output (only for non-human species, can be unspecified)  
**SD_file:** A segmental duplication file from [UCSC](https://hgdownload.soe.ucsc.edu/downloads) (only for non-human species, can be unspecified)  
**ref_index:** A reference fasta index file (mandatory for non-human species)  

The Step-4 corrects the genotypes (given with GT tag) and adds a new tag, MC, to the FORMAT/SAMPLE fields of the output vcf file. The MC tag represents the level of SMC; the lower the MC value, the higher the confidence level. The genotyping step requires a parameter file indicating the minimum SRR for each SV type and algorithm. MOPline provides a default parameter file for 14 pre-selected algorithms, which is automatically selected during this step. If additional algorithms are to be used that are not listed in this parameter file, the user can edit the parameter file by adding parameters of those algorithms (if not edited, all algorithms not present in the parameter file will have a minimum RSS of 3).

This step takes longer to perform as the sample size increases and the genome size increases. For human samples larger than 1,500, it is recommended that this step is performed for each chromosome, which can be done using the -c option (default: ALL).

### <a name="step5"></a>[Step-5]  Annotate (optional)

Step-5 adds gene name/ID and gene region that overlap the SV to the INFO filed (with SVANN key) of the vcf file. The gene region includes exon/CDS (All-exons if the SV completely overlaps all exons), 5’-/3’-UTR, intron, 5’-/3’-flanking regions. Two ranges of the flanking regions are specified by default (5 Kb and 50 Kb); these lengths can be changed with the options, -c5, -c3, -f5, and -f3. These annotations are also added to the FORMAT AN subfield for each sample in an additional output vcf file. For human, the gff3 gene annotation files (Homo_sapiens.GRCh37.87.gff3.gz, Homo_sapiens.GRCh38.104.gff3.gz, or Homo_sapience.T2T-chm13v2.0.ensemble.gff3.gz), downloaded from Ensembl (ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens9), is selected by default. For non-human species, a gff3 annotation file obtained from the Ensembl site must be specified with the -r option. Any input SV vcf file with SVTYPE and SVLEN keys in the INFO field may be used. As input vcf file, an output vcf file from step2, step3, step4, or step6 can be used. The annotate command can be done as follows:
```
mopline annotate -v <input_vcf> -p <out_prefix> -n <num_threads>
```
(--build 38 for human build 38, -nh 1 -r <gff3_file> for non-human species)  
This command generates two output vcf files, `${out_prefix}.annot.vcf` and `${out_prefix}.AS.annot.vcf`. The latter contains annotations for each sample in the FORMAT AN subfield.

### <a name="step6"></a>[Step-6]  Filter (optional)

This step, using the filter_MOPline.pl script, filters out DEL/DUPs with inconsistent DPRs. DUPs associated with gap regions and DUPs overlapping segmental duplications also filtered out based on several criteria (see our paper for detail). For human samples, gap bed and segmental duplication files are automatically selected from the Data directory. For non-human samples, these files can be specified with the -gap and -segdup options. The -ex option can also be used to specify a bed file that indicates regions to exclude SVs. This package provides bed files for human reference builds 37, 38, and T2T-CHM13 that indicate regions where SVs are always indeterminately called due to low quality alignments of short reads. The total lengths of the excluding regions of build 37 and 38 are 26 Mb and 16 Mb, respectively. In addition, SVs within the centromere regions for the human builds 38 and T2T-CHM13 are removed. By default, the bed files of these regions are automatically selected for human but user-provided files can be specified with the -ex and -ec options. If you do not prefer to use these filtering, specify any letter (e.g., -ex 0 or -ec 0) for the options. The input vcf can be from a single sample or from multiple samples but must have the keys DPR, SR, and DPS in the INFO field. As input vcf file, an output vcf file from step2, step3, step4, or step5 can be used.
```
mopline filter -v <input vcf> > [output vcf]
```
(--build 38 for human build 38, -nh 1 -g <gap_bed> -sd <segmental duplication file> for non-human species)

## <a name="nhdata"></a>Data Required for Non-Human Species

For non-human species, the following files are required at specific MOPline steps (for human, these files are provided in this package and automatically selected at each step). Human sample data files are present in the Data directory of this package.

**Table 2.** Files required for non-human species, which must be supplied by users
|Name  |Description                   |Sample file                   |Available site|Required step|
| :--- | :--------------------------- | :--------------------------- | :------ | :------- |
|Gap   |BED file of 'N' regions in the reference fasta|gap.bed     |[UCSC](https://hgdownload.soe.ucsc.edu/downloads)|steps-0,2,3,6
|STR   |Simple tandem repeats (required columns: 2,3,4)|simpleRepeat.txt.gz |[UCSC](https://hgdownload.soe.ucsc.edu/downloads)|step-4
|SegDup|Segmental duplications (required columns: 2,3,4)|genomicSuperDups.txt.gz|[UCSC](https://hgdownload.soe.ucsc.edu/downloads)|steps-4,6
|Cen   |BED file of centromere regions|hg38.centromere.bed|[UCSC](https://hgdownload.soe.ucsc.edu/downloads)|step-6
|Gene  |Gene annotation GFF3 file   |Homo_sapiens.GRCh37.87.gff3.gz|[Ensembl](https://asia.ensembl.org/index.html)|step-5
|Rindex|Reference fasta index       |hs37.fa.fai                   |samtools faidx|steps-1,2

### For data unavailable at the public sites

Gap bed and STR file can be custom created. But MOPline can be executed at any step, even if Gap, STR, SegDup, or Cen file is not specified.  

**Gap bed:** A gap bed file can be created using scripts/gap_fasta_to_bed.pl in this package.  
```
gap_fasta_to_bed.pl [input genome.fasta] [minimum length of gap, default: 1000 (bp)] > [output bed]
```
**STR (SimpleRepeat) file:** STR/SimpleRepeat file can be created using the [TRF](https://github.com/Benson-Genomics-Lab/TRF) and scripts/convert_TRF_data_to_simpleRepeat.pl.  
  
(TRF run)  
```
trf409.linux64 [input genome.fasta] 2 7 7 80 10 50 500 -d -l 1
```
(Convert TRF data to UCSC SimpleRepeat format)  
```
convert_TRF_data_to_simpleRepeat.pl [TRF output *.data file] > [output SimpleRepeat-like file]
```
## <a name="single"></a>For single sample

For MOPline with a single sample, all the steps can be performed, but the required directory structure (sample and tool directories) is required for completion. If SV genotyping is not needed, the steps-3 and -4 can be omitted.

## <a name="quick"></a>Quick Start with Sample Data

The [sample data](http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/Sample_data.tar.gz) (or available from https://drive.google.com/drive/folders/1bIEtaaM3xx8POIAf96kV-ImNXTHWPwQQ?usp=sharing) provided includes two sample datasets: human and yeast (*Saccharomyces cerevisiae*) data. The human data includes SV datasets generated with high coverage WGS datasets of 6 individuals from a 1000 Genomes CEU population (ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data). The yeast data includes WGS bam files for 10 yeast isolates (Peter J. et al., Nature 556, pages 339–344 (2018)) and associated files including reference fasta and annotation files. An [output dataset](http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/Sample_data_output.tar.gz) generated with this dataset using MOPline-7t is also available.

### <a name="hsample"></a>A: Using Human SV datasets of 6 samples

This sample SV data was generated using the 7tools preset (MOPline-7t) with 6 human WGS data (bam files of 30×, 150 bp paired-end sequencing data aligned against the GRCh37 reference), and only SVs corresponding to chromosome 17 were extracted. Each sample directory contains 7 tool directories, and each tool directory contains an SV vcf file generated with the conversion scripts shown in Table 1. Each sample directory contains a Cov directory, which contains a coverage file of chromosome 17 that records alignment statistics from the bam file.
```
cd Human_chr17_b37_6samples
export PATH=$PATH:${MOPline-install-directory-path}
(e.g., export PATH=$PATH:/home/tools/MOPline_v1.7)
```
#### [Step-1] Select overlap calls from SV call sets
```
mopline merge_7tools -s sample_list.txt -rl 150 -d Merge_7tools
```
The command adds alignment statistics and genotype data to `${sample}.Merge.ALL.vcf` for each sample. Check the values given by the DR, DS, and SR tags in the FORMAT/SAMPLE fields of the vcf files.

#### [Step-2] Add alignment statistics to SV sites
Step 2-(1) ‘Create coverage files’ is omitted for this sample data because the Cov directory in each sample director already contains coverage files.
```
mopline add_cov -s sample_list.txt -ts 7tools -vd Merge_7tools -n 4 
```
The command adds alignment statistics and genotype data to `${sample}.Merge.ALL.vcf` for each sample. Check the values given by the DR, DS, and SR tags in the FORMAT/SAMPLE fields of the vcf files.

#### [Step-3] Merge vcf files from multiple samples
```
mopline joint_call -s sample_list.txt -md Merge_7tools -od JointCall -p MOPline
```
The command generates a MOPline.All-samples.vcf in the JointCall directory.

#### [Step-4] Genotype and SMC
```
mopline smc -v JointCall/MOPline.All-samples.vcf -ts 7tools -od MOPline_SMC -p MOPline.smc -n 4
```
A vcf file named MOPline.smc.vcf has been created in the MOPline_SMC directory, with new MC tags in the FORMAT/SAMPLE fields and AF, AC, DPR, SR, and DPS keys in the INFO field.

#### [Step-5] Annotate
```
cd MOPline_SMC
mopline annotate -v MOPline.smc.vcf -p MOPline -n 4
```
This will generate two annotated vcf files in the working directory, MOPline.annot vcf and MOPline.AS.annot.vcf. The latter vcf file contains sample-level gene-overlap annotations in the FORMAT/SAMPLE fields with AN tags.

#### [Step-6] Filter
```
mopline filter -v MOPline.AS.annot.vcf > MOPline.AS.annot.filt.vcf
```

### <a name="ysample"></a>B: Using Yeast WGS bam files of 10 samples

This sample data contains yeast bam alignment files of 30×, 101 bp paired-end WGS data aligned using bwa mem against the S288C *S. cerevisiae* reference for 10 yeast strains. The dataset includes the S288C reference fasta file, *S. cerevisiae* gene annotation gff3 file (ftp://ftp.ensembl.org/), STR repeat file (https://genome.ucsc.edu), andTY1/TY3 retroelement reference files for MELT. The TY1/TY3 MELT reference files were generated according to the [MELT documentation](https://melt.igs.umaryland.edu/manual.php). This tutorial begins with SV calling using the 7tools preset. The commands for all steps have “-nh 1” since the yeast is a non-human species.
```
export PATH=$PATH:${MOPline-install-path}/scripts/run_SVcallers
mkdir yeast_run
cd yeast_run
ln -s ../Yeast_10samples/S288C.fa
ln -s ../Yeast_10samples/S288C.fa.fai
ln -s ../Yeast_10samples/Sc.R64-1-1.104.ensemble.gff3
ln -s ../Yeast_10samples/Sc_simpleRepeat.txt
ln -s ../Yeast_10samples/MELT_lib
ln -s ../Yeast_10samples/bam_list.txt
create_bam_link.pl -b bam_list.txt -bd ../Yeast_10samples
```

#### [Step-0] Run SV detection tools (working directory: yeast_run)

(Single sample mode [sample: CBS457])  
```
mkdir CBS457
cd CBS457
run_single.pl -c ../config.yeast.txt -b CBS457.bam -sn CBS457
```
A sample config file (config.yeast.txt) is contained in the Sample_data_output/yeast_run folder.
A run with this sample config file will generate SV call results of 7 tools in 7 tool directories.
	
(Batch mode using slurm)  
```
run_batch_slurm.pl -c config.yeast.txt -b bam_list.txt -a <account> -p <partition>
```
A run with this sample config file will generate SV call results of 7 tools in 7 tool directories in each of 10 sample directories.

#### [Step-1] Select overlap calls from SV call sets
```
mopline merge_7tools -s bam_list.txt -rl 101 -d Merge_7tools -nh 1 -r S288C.fa.fai
```
The command generates a vcf files named `${sample}.Merge.ALL.vcf` in the Merge_7tools directory generated in each sample directory. For non-human species, a reference index file should be specified with the -r option.

#### [Step-2] Add alignment statistics to SV sites
```
mopline create_cov -b bam_list.txt -r S288C.fa -rl 101 -n 4 -nh 1
```
The command creates a Cov directory under the sample directory, which contains the coverage files for each chromosome. A reference fasta file must be specified with the -r option.
```
mopline add_cov -s bam_list.txt -ts 7tools -vd Merge_7tools -n 4 -nh 1 -r S288C.fa.fai
```
The command adds alignment statistics and genotype data to `${sample}.Merge.ALL.vcf` for each sample. Check the values given in the DR, DS, and SR tags in the FORMAT/SAMPLE fields of the vcf files. For non-human species, a reference index file should be specified with the -r option.

#### [Step-3] Merge vcf files from multiple samples
```
mopline joint_call -s bam_list.txt -md Merge_7tools -od JointCall -p MOPline -nh 1
```
This generates a MOPline.All-samples.vcf in the JointCall directory.

#### [Step-4] Genotype and SMC
```
mopline smc -v JointCall/MOPline.All-samples.vcf -ts 7tools -r S288C.fa.fai -od MOPline_SMC -p MOPline.smc -n 4 -nh 1 -sr Sc_simpleRepeat.txt
```
A vcf file named MOPline.smc.vcf has been created in the MOPline_SMC directory, with new MC tags in the FORMAT/SAMPLE fields and AF, AC, DPR, SR, and DPS keys in the INFO field. For non-human species, the -sr option specifies a simple repeat file.

#### [Step-5] Annotate
```
cd MOPline_SMC
mopline annotate -v MOPline.smc.vcf -p MOPline -n 4 -nh 1 -r Sc.R64-1-1.104.ensemble.gff3
```
This will generate two annotated vcf files in the working directory, MOPline.annot vcf and MOPline.AS.annot.vcf. The latter vcf file contains sample-level gene-overlap annotations with AN tags in the FORMAT/SAMPLE fields. For non-human species, a gff3 annotation file must be specified with the -r option.

#### [Step-6] Filter
```
mopline filter -v MOPline.AS.annot.vcf -nh 1 > MOPline.AS.annot.filt.vcf
```
