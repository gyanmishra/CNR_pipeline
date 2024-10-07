This pipeline processes paired-end raw CUT&RUN/ChIP sequencing reads from the mouse to generate BAM, bedGraph and bigwig files. After processing the reads, it also includes a script to perform peak calling using macs2.

This pipeline is tested to work only in UTSW BioHPC.
\
This pipeline has 2 major components: 

## Processing of Raw reads

Processing and alignment of raw data are in the order.
 
[FASTQC(v0.11.8) >> multiqc(v1.7) >> Trimmomatic(v0.32) >> multiqc((v1.7)) >>bowtie2(v2.4.2) >> samtools(v1.6) (sort, fixmate, markdup) >> bamCoverage (deeptools(v3.5.0))]


The intermediate bam files are removed and only the original alignment BAM file and final duplicate removed BAM files are stored. 

Steps to be followed to run the pipeline. 

**I. Create a Project directory**
It is not compulsory to create the below directories, but this could be one way to manage project folder best.
```
mkdir -p Name_of_project_dir/data # To keep all the required .csv/.tsv files 
mkdir -p Name_of_project_dir/results # This is where the output of this pipeline will be stored
mkdir -p Name_of_project_dir/scripts # This directory can be used to keep the scripts
mkdir -p Name_of_project_dir/Figures # and to store any figures generated from downstream analysis performed
```
e.g
```
mkdir -p MECP2_project/data # To keep all the required .csv/.tsv files 
mkdir -p MECP2_project/results # This is where the output of this pipeline will be stored
mkdir -p MECP2_project/scripts # This directory can be used to keep the scripts
mkdir -p MECP2_project/Figures

```

**II. Create sampleSheet.csv**
*Note: the content of sampleSheet.csv should not contain space*
e.g.
```
SampleID,Read1.fastq.gz,Read2.fastq.gz
7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729_TC240531_A1,TC240531_A1_S1_R1_001.fastq.gz,TC240531_A1_S1_R2_001.fastq.gz
7wk_MECP2_WT_LR_CTX_H3K4me1_ab8895_TC240531_A2,TC240531_A2_S2_R1_001.fastq.gz,TC240531_A2_S2_R2_001.fastq.gz
```
Save the sampleSheet.csv in `MECP2_project/data/`

**III. Download the mm10 bowtie2 index**\
*Note: skip this step if you already have the required index files of the genome*
```
$ mkdir -p MECP2_project/data/bowtie2_mm10_index
$ wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip -P MECP2_project/data/bowtie2_mm10_index
$ unzip MECP2_project/data/bowtie2_mm10_index/mm10.zip -d MECP2_project/data/bowtie2_mm10_index/
```

**IV. Run the pipeline with the following argument as inputs**

Run the below command to process the raw data. 
```
$ git clone https://github.com/gyanmishra/CNR_pipeline.git
$ perl CNR_pipeline/process_CNR.pl \
<result dir> \ # the path of result directory
<raw file directory> \ # the path of directory containing raw files
<bowtie2_index directory> \ # the path of the directory containing bowtie2 index
<sampleSheet.csv> # the path of files with sample Info (e.g sampleSheet.csv)
```
e.g.

```
e.g if Sequencing Run id is "2024_06_11_N2K226_13907_0"

$ perl CNR_pipeline/process_CNR.pl \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/ \
/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/D3aCNR_forGyan/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/bowtie2_mm10_index/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/sampleSheet.csv 
```

*Note : The script `process_CNR.pl` can be run from any directory given the full path of `CNR_Pipeline` directory*

**V. Generate Alignment statistics in tabular format for all the samples**\
*Note: Run the below command from `CNR_pipeline` directory only to generate alignment statistics.*\
see below code
```
$ cd CNR_pipeline
$ sbatch GenerateAlignment_statistics.sh \
<result dir> \ # the path of result directory
<raw file directory> \ # the path of directory containing raw files
<sampleSheet.csv> # the path of files with sample Info (e.g sampleSheet.csv)
```
e.g
```
sbatch GenerateAlignment_statistics.sh \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/ \
/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/D3aCNR_forGyan/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/sampleSheet.csv 
```
open below file to to see the alignment statistics
`/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/stats/Alignment_stats.tsv` 

**Details of output folder:**
1. `MECP2_project/results/2024_06_11_N2K226_13907_0/bam`
    This directory contains the final duplicate removed aligned coordinate sorted bam file and its index.

2. `MECP2_project/results/2024_06_11_N2K226_13907_0/bamCoverage`
    This directory contains RPKM normalized bedGraph and bigwig files with 10bp bin size generated by bamCoverage program of deeptools

3. `MECP2_project/results/2024_06_11_N2K226_13907_0/fastqc_before_trimmomatic`
    This directory contains FASTQC and multiqc output of original fastq files.

4. `MECP2_project/results/2024_06_11_N2K226_13907_0/fastqc`
    This directory contains FASTQC and multiqc output of trimmed fastq files.

5. `MECP2_project/results/2024_06_11_N2K226_13907_0/logs`
    This directory contains slrum script and log file for each sample

6. `MECP2_project/results/2024_06_11_N2K226_13907_0/qc_qualimap`
    This directory contains output directory from Qualimap program which reports alignment statistics.

### Peak calling (macs2)

**I. Prepare tab separated sampleSheet.tsv**
The targetID and controlID should be sampleID used in the sampleSheet.csv file.
```
peakCallingID	targetID	controlID
MECP2_vs_IgG	240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729,240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729	240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729,240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729
```
*save the file in `MECP2_project/data/macs2_sampleSheet.tsv`

**II. Run the peak calling script**\
Run the below command to call peaks. 
```
perl macs2_callPeak.pl --help

Usage: macs2_callPeak.pl \
--resultDir <result dir> \ # the path of result directory
--qvalue \ # use this option to use pvalue
--pvalue \ # use this option to use qvalue (default)
--threshold <set pvalue or qvalue threshold> (default : 0.05)
--sampleSheet <macs2 peak calling sampleSheet>
```

e.g
```
$ perl CNR_pipeline/macs2_callPeak.pl \
--resultDir /work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/ \
--pvalue \
--threshold 0.01 \
--sampleSheet /work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/macs2_sampleSheet.tsv
```
*Note : The script `macs2_callPeak.pl` can be run from any directory given the full path of `CNR_Pipeline` directory** <br >

**Details of macs2 output directory :**

1. `MECP2_project/results/2024_06_11_N2K226_13907_0/macs2_callPeak/MECP2_vs_IgG/` \
    This folder will be created after running `macs2_callPeak.pl` and it will contain following files. 
    - MECP2_vs_IgG_model.r
    - MECP2_vs_IgG_peaks.narrowPeak
    - MECP2_vs_IgG_peaks.xls
    - MECP2_vs_IgG_summits.bed
    - MECP2_vs_Igg_peaks.narrowPeak.noBlacklist.bed
    - MECP2_vs_IgG_summits_noBlacklist.bed


### Motif analysis (HOMER)
```
$ sbatch CNR_pipeline/homer_motif.sh \
<peak file> \ # peak file 
<genome> \ # genome (e.g mm10)
<output directory> # the path of output directory \
<size of the region from center>  # (e.g -50,50)
```