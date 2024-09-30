# How to submit peak calling job in BioHPC

create a file and paste the content mentioned below in the file. Edit the file as per the given instruction and save it as something like `run_macs2_peaks.sh`

```
#!/bin/bash

# Specify the job name
#SBATCH --job-name WTBend6-5_vs_WTinput

# replace the number after N to reserve nodes. Max 16 for > 32GB and 64 for 32GB
#SBATCH -N 1

# specify partition among (32GB, 128GB, 256GB, super)
#SBATCH --partition=32GB

# specify the amount of time for running this script (if job completes early before given amount of time then it autamtically gets killed)
#SBATCH -t 0-12:0:0

# The below mentioned file where you want keep the log of run. in case of any error, it will be reported in .err file.

#SBATCH -o WTBend6-5_vs_WTinput_20240909_macs2.out
#SBATCH -e WTBend6-5_vs_WTinput_20240909_macs2.err


# module required to run below macs2 and filter blacklisted regions.
module load macs/2.1.2 bedops/2.4.14 



# peak calling command

macs2 callpeak \
-t test_run/results/240531//bam/WT_M_40wk_BEND6-5_LRCTX_ME240724-A2_N2K263_14062.final.sort.bam \ # provide target bam file
-c test_run/results/240531//bam/WT_M_40wk_input_LRCTX_ME240724-A5_N2K263_14062.final.sort.bam \ # provide control bam file
-f BAMPE \ # change to BAM for single end
-g mm \ # mm is for mouse
-p 1e-3 \ this is p-value cut-off and can be changed or entire argument can be removed, if removed then by default peaks will be reported with qvalue < 0.05
--outdir test_run/results/240531//macs2_callPeak/WTBend6-5_vs_WTinput \ # this is output directory where you want to keep the output of macs2 
--call-summits \
-n WTBend6-5_vs_WTinput # this is the prefix used to name files.

# Now once the peak calling is complete then below command will remove the blacklisted regions. 
# The blacklisted regions are present in mm10 directory inside where this markdown is located.

#remove blacklisted region

bedops -n 1 \
test_run/results/240531//macs2_callPeak/WTBend6-5_vs_WTinput/WTBend6-5_vs_WTinput_peaks.narrowPeak \ #paste here the ouput directory path mentioned in previos command followed by name with suffix _peaks.narrowPeak
/endosome/work/Neuroinformatics_Core/s225347/Pipeline/CNR_pipeline/mm10/mm10.blacklist.bed \ # this is the path of blacklisted bed file.
>test_run/results/240531//macs2_callPeak/WTBend6-5_vs_WTinput/WTBend6-5_vs_WTinput_peaks.narrowPeak.noBlacklist.bed #ouput file name

# similarly for summit file
bedops -n 1 \
test_run/results/240531//macs2_callPeak/WTBend6-5_vs_WTinput/WTBend6-5_vs_WTinput_summits.bed \
/endosome/work/Neuroinformatics_Core/s225347/Pipeline/CNR_pipeline/mm10/mm10.blacklist.bed \
>test_run/results/240531//macs2_callPeak/WTBend6-5_vs_WTinput/WTBend6-5_vs_WTinput_summits_noBlacklist.bed

```

Aftet you have edited the file and provided all the argument then save it and submit the script as job :
`sbatch run_macs2_peaks.sh`

To check the status of job 
`squeue -u sXXXXXX`

Two log files will be created in with suffix .out and .err provided in the option
```
#SBATCH -o WTBend6-5_vs_WTinput_20240909_macs2.out
#SBATCH -e WTBend6-5_vs_WTinput_20240909_macs2.err
```

