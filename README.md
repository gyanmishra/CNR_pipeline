This pipeline process paired-end raw CUT&RUN/ChIP sequening reads from mouse to generate BAM, bedGraph and bigwig files. After processing the reads, it also includes script to perform peak calling using macs2.

This pipleine is tested to work in UTSW BioHPC
This pipleine has 2 major component : 

## Processing of Raw reads

Processing and alignment of raw data are in the order** 
[FASTQC(v0.11.8) >> multiqc(v1.7) >> Trimmomatic(v0.32) >> multiqc((v1.7)) >>bowtie2(v2.4.2) >> samtools(v1.6) (sort, fixmate, markdup) >> bamCoverage (deeptools(v3.5.0))]


The intermediate bam files are removed and only orginal alignment BAM file and final Aligned files are stored. 

Steps to be followed to run the pipeline. 

**I. Create a Project directory**
It is not compulsory to create below directories, but this could be one way to best manage project folder.
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
save the sampleSheet.csv in `Name_of_project_folder/data/`

**III. Download the mm10 bowtie2 index**\
*Note: skip this step if you already have the reuqired index files of the genome*
```
$ mkdir -p Name_of_project_folder/data/bowtie2_mm10_index
$ wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip -P Name_of_project_folder/data/bowtie2_mm10_index
$ unzip Name_of_project_folder/data/bowtie2_mm10_index/mm10.zip -d Name_of_project_folder/data/bowtie2_mm10_index/
```

**IV. Run the pipeline with following argument as inputs**

Run below command from the CNR_pipeline directory. 
```
$ git clone https://github.com/gyanmishra/CNR_pipeline.git
$ perl CNR_pipeline/process_CNR.pl \
<result dir> \ # the path of result directory
<raw file directory> \ # the path of directory containing raw files
<bowtie2_index directory> \ # the path of directory containing bowtie2 index
<sampleSheet.csv> # the path of files with sample Info (e.g sampleSheet.csv)
```
e.g.

```
Sequencing RUN ID : 2024_06_11_N2K226_13907_0
$ perl CNR_pipeline/process_CNR.pl \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/ \
/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/D3aCNR_forGyan/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/bowtie2_mm10_index/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/data/sampleSheet.csv 
```

*Note : The script `process_CNR.pl` can be run from any directory*

**V. Generate Alignment statisitcs in tabular format for all the samples**
Run below command from the CNR_pipeline directory. 

```
$ sbatch Generate_Alignment_stats.sh \
<result dir> \ # the path of result directory
<raw file directory> \ # the path of directory containing raw files
<sampleSheet.csv> # the path of files with sample Info (e.g sampleSheet.csv)
```
e.g
```
sbatch CNR_pipeline/Generate_Alignment_stats.sh \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/Name_of_project_folder/results/2024_06_11_N2K226_13907_0/ \
/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/D3aCNR_forGyan/ \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/Name_of_project_folder/data/sampleSheet.csv 
```
open below file to to see the alignment statistics
`/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/Name_of_project_folder/results/2024_06_11_N2K226_13907_0/stats/Alignment_stats.tsv` 

**Details of output folder :**
1. `MECP2_project/results/2024_06_11_N2K226_13907_0/bam`
    This directory contains the final duplicate removed alinged coordinate sorted bam file and its index.

2. `MECP2_project/results/2024_06_11_N2K226_13907_0/bamCoverage`
    This directory contains RPKM normalized bedGraph and bigwig files with 10bp bin size generated by bamCoverage program of deeptools

3. `MECP2_project/results/2024_06_11_N2K226_13907_0/fastqc_before_trimmomatic`
    This directory contain FASTQC amd multiqc output of original fastq files.

4. `MECP2_project/results/2024_06_11_N2K226_13907_0/fastqc`
    This directory contain FASTQC and multiqc output of trimmed fastq files.

5. `MECP2_project/results/2024_06_11_N2K226_13907_0/logs`
    This directory contain slrum script and log file for each sample

6. `MECP2_project/results/2024_06_11_N2K226_13907_0/qc_qualimap`
    This directory contain output directory from Qualimap program which reports alignment statistics.

### Peak calling (macs2)

**I. Prepare tab separated sampleSheet.tsv**
The targetID and controlID should be smapleID used in the sampleSheet.csv file.
```
peakCallingID	targetID	controlID
MECP2_vs_IgG	240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729,240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729	240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729,240531TC_A1_7wk_MECP2_WT_LR_CTX_H3K27ac_ab4729
```
*save the file in `Name_of_project_folder/data/macs2_sampleSheet.tsv`

**2. Run the peak calling script**
Run below command from the CNR_pipeline directory. 
```
$ perl CNR_pipeline/macs2_callPeak.pl \
<result directory> \ # the path of result directory
<pvalue or qvalue> \ # pvalue - to find peaks with p-value cut-off of 1e-3; qvalue - to find peaks with qvalue cut-off of 0.05
<macs2_sampleSheet.tsv> # path with name of macs2_sampleSheet.tsv
```

e.g
```
$ perl CNR_pipeline/macs2_callPeak.pl \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/MECP2_project/results/2024_06_11_N2K226_13907_0/ \
pvalue \
/work/OBI/Neuroinformatics_Core/s225347/CNR_pipeline/Name_of_project_folder/data/macs2_sampleSheet.tsv \

```
*Note : The script `macs2_callPeak.pl` can be run from any directory* <br >

**Details of macs2 output directory :**

1. `MECP2_project/results/2024_06_11_N2K226_13907_0/macs2_callPeak/MECP2_vs_IgG/` \
    This folder will be created after running `macs2_callPeak.pl` and it will contain following files. 
    - MECP2_vs_IgG_model.r
    - MECP2_vs_IgG_peaks.narrowPeak
    - MECP2_vs_IgG_peaks.xls
    - MECP2_vs_IgG_summits.bed
    - MECP2_vs_Igg_peaks.narrowPeak.noBlacklist.bed
    - MECP2_vs_IgG_summits_noBlacklist.bed


