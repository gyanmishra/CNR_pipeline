#!/usr/bin/python

import os
import sys
import pandas as pd
import argparse
import subprocess

# list to store files
res = []
store = {}
report_dict = {}


# Iterate directory
def generate_readstats(samplesheet,result_dir,raw_dir):

    sampleSheet = pd.read_csv(args.samplesheet)

    df = pd.DataFrame(columns=['sampleID','Raw_reads(x2)','Trimmed_reads(x2)','Total_Aligned_reads(x2)','Aligned_reads_after_duplicates_removed(x2)'])

    for index, row in sampleSheet.iterrows():

        #print(index)
        Raw_reads = None
        Trimmed_reads = None
        Aligned_reads = None
        Duplicate_removed_aligned_reads = None
        sampleID = row['SampleID']
        readID = row['Read1']

        raw_file_path = os.path.join(raw_dir, row['Read1'])
        raw_reads = "zcat {} | wc -l | awk '{{print $1/4}}'".format(raw_file_path)
        
        trimmed_reads_path = os.path.join(result_dir + 'fastq', row['SampleID'] + '.R1.paired.fastq.gz')
        trimmed_reads = "zcat {} | wc -l | awk '{{print $1/4}}'".format(trimmed_reads_path)

        aligned_bam = os.path.join(result_dir + 'bam', sampleID + '.s.bam')
        aligned_bam_command = "samtools view -c -F 4 {} | awk '{{print $1/2}}'".format(aligned_bam)

        duplicate_removed_bam = os.path.join(result_dir + 'bam', sampleID + '.final.sort.bam')
        duplicate_removed_bam_command = "samtools view -c -F 4 {} | awk '{{print $1/2}}'".format(duplicate_removed_bam)

        #print(duplicate_removed_bam)

        Raw_reads = subprocess.check_output(raw_reads,shell=True).strip()
        #print(Raw_reads)
        Trimmed_reads = subprocess.check_output(trimmed_reads,shell=True).strip()
        #print(Trimmed_reads)
        Aligned_reads = subprocess.check_output(aligned_bam_command,shell=True).strip()
        #print(Aligned_reads)
        Duplicate_removed_aligned_reads = subprocess.check_output(duplicate_removed_bam_command,shell=True).strip()
        #print(Duplicate_removed_aligned_reads)
        df.loc[index] = [sampleID,Raw_reads,Trimmed_reads,Aligned_reads,Duplicate_removed_aligned_reads]

    return(df)



                                   

if __name__ == "__main__":
    
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="A script to generate Alignment statistics")

    # Define arguments with both short and long options
    parser.add_argument('-s', '--samplesheet', type=str, required=True,help='Provide sampleSheet.csv')
    parser.add_argument('-d', '--resultDir', type=str, required=True,help='Provide result directory path')
    parser.add_argument('-r', '--rawDir', type=str, required=True,help='Provide raw directory path')

    args = parser.parse_args()
    #print(args)


    df = generate_readstats(args.samplesheet,args.resultDir,args.rawDir)
    #print(df)
    stats = os.path.join(args.resultDir, 'stats')
    #print(stats)
    if not os.path.exists(stats):
        os.makedirs(stats)

    df.to_csv(stats+"/Alignment_stats.tsv", sep='\t',index=False)