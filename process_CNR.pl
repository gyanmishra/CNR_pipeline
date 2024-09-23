#!/usr/bin/perl

use POSIX qw(strftime);
use FindBin;
use Cwd 'abs_path';

# Get the absolute path of the current script
my $script_path = abs_path($FindBin::Bin);
print "$script_path\n";

# Project Directory
$outputDir = $ARGV[0]; 
print "$outputDir\n";

# Raw data direcetory
$rawFile = $ARGV[1];
print "$rawFile\n";

# Index Directory
$index = $ARGV[2];
print "$index\n";

# results directory
unless (-e "$outputDir/fastqc_before_trimmomatic") {
    `mkdir -p $outputDir/fastqc_before_trimmomatic`;
}
unless (-e "$outputDir/fastq") {
    `mkdir -p $outputDir/fastq`;
}
unless (-e "$outputDir/fastqc") {
    `mkdir -p $outputDir/fastqc`;
}
unless (-e "$outputDir/bam") {
    `mkdir -p $outputDir/bam`;
}
unless (-e "$outputDir/qc_qualimap") {
    `mkdir -p $outputDir/qc_qualimap`;
}
unless (-e "$outputDir/logs") {
    `mkdir -p $outputDir/logs`;
}
unless (-e "$outputDir/bamCoverage") {
    `mkdir -p $outputDir/bamCoverage`;
}


$datestring = strftime "%F", localtime;
$datestring =~ s/-//g ; 


# provide a file having three column 
# 1st column : 1st column as Sample name and 
# 2nd column : Read_1.fastq.gz
# 3rd column : Read_2.fastq.gz

$sampleData  = $ARGV[3];
open my $file, $sampleData or die "Could not open $sampleData: $!";

my $first = 1;
foreach $line (<$file>) 
{
       chomp($line);
       if( $first ) {
            $first = 0;
        }
        else {
              chomp($line);
              @col = split(',',$line);
              chomp($col[0]);
              chomp($col[1]);
              chomp($col[2]);

              $sampleName = $col[0];
              print "$sampleName\n";
              #$sampleName =~ s/\s/_/g;
              $out_filename = $outputDir . "/logs/" . $sampleName . '_' . $datestring .'.sh';

              print "$out_filename\n";
              #print "$sampleName\n";
              open(FH, '>', $out_filename) or die $!;
              print FH "#!/bin/bash\n\n";
              print FH "#SBATCH --job-name CUT_and_RUN\n";
              print FH "#SBATCH -N 1\n";
              print FH "#SBATCH --cpus-per-task=32\n";
              print FH "#SBATCH --partition=32GB\n";
              print FH "#SBATCH -t 0-12:0:0\n";
              print FH "#SBATCH -o $outputDir/logs/$col[0].$datestring.out\n";
              print FH "#SBATCH -o $outputDir/logs/$col[0].$datestring.err\n";
              print FH "#SBATCH --mail-type ALL\n";
              #print FH "#SBATCH --mail-user GyanPrakash.Mishra\@UTSouthwestern.edu\n\n\n";
              print FH "module load fastqc/0.11.8 multiqc/1.7 perl/5.30.1 bowtie2/2.4.2 samtools/1.6 bedtools/2.29.2 UCSC_userApps/v317 Trimmomatic/0.32 deeptools/3.5.0 \n\n";
              #print FH "perl /work/OBI/Neuroinformatics_Core/s225347/CUT_and_RUN/scripts/CNR_pipeline.pl /work/OBI/Neuroinformatics_Core/s225347/CUT_and_RUN/data/CUTandRUN_seq_Data.txt\n\n";
              print FH "\n\n";
              print FH "# COMMAND GROUP 1\n\n";         
              # Merge samples fastq
              #print FH "cat $rawFile/$col[1]*_R1*fastq.gz >$projectDir/results/fastq/$col[0]_R1.fastq.gz\n";
              #print FH "cat $rawFile/$col[1]*_R2*fastq.gz >$projectDir/results/fastq/$col[0]_R2.fastq.gz\n";
              
              ## FASTQC 
              print FH "fastqc $rawFile/$col[1] \\\n-o $outputDir/fastqc_before_trimmomatic\n\n";
              print FH "fastqc $rawFile/$col[2] \\\n-o $outputDir/fastqc_before_trimmomatic\n\n";

              # MultiQC before trimmomatic
              print FH "multiqc -f $outputDir/fastqc_before_trimmomatic -o $outputDir/fastqc_before_trimmomatic\n\n";

              # Trim adapters
              #print FH "cutadapt -a  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -A  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -m 20 -o $projectDir/results/fastq/$col[4]_R1.cutadapt.fastq.gz -p $projectDir/results/fastq/$col[4]_R2.cutadapt.fastq.gz $rawFile/$col[5] $rawFile/$col[6]\n";
              #
              print FH "java -jar \$Trimmomatic PE  \\\n$rawFile/$col[1] \\\n$rawFile/$col[2] \\\n$outputDir/fastq/$sampleName.R1.paired.fastq.gz \\\n$outputDir/fastq/$sampleName.R1.unpaired.fastq.gz \\\n$outputDir/fastq/$sampleName.R2.paired.fastq.gz \\\n$outputDir/fastq/$sampleName.R2.unpaired.fastq.gz \\\nILLUMINACLIP:trimmomatic/adapters/All_TruSeq.fa.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25\n\n";

              ## FASTQC 
              print FH "fastqc $outputDir/fastq/$sampleName.R1.paired.fastq.gz \\\n-o $outputDir/fastqc\n\n";
              print FH "fastqc $outputDir/fastq/$sampleName.R2.paired.fastq.gz \\\n-o $outputDir/fastqc\n\n";

              # MultiQC
              print FH "multiqc -f $outputDir/fastqc -o $outputDir/fastqc\n\n";

              ## Align to mm10 reference genome
              print FH "bowtie2 -p 4 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \\\n-x $index/mm10 \\\n-1 $outputDir/fastq/$sampleName.R1.paired.fastq.gz \\\n-2 $outputDir/fastq/$sampleName.R2.paired.fastq.gz | samtools view -bS -  \\\n>$outputDir/bam/$sampleName.bam\n\n";
              ##
              print FH "samtools sort -@ 8 -n -O BAM \\\n$outputDir/bam/$sampleName.bam \\\n-o $outputDir/bam/$sampleName.s.bam\n\n";
              print FH "samtools fixmate -m \\\n$outputDir/bam/$sampleName.s.bam \\\n$outputDir/bam/$sampleName.s.fm.bam\n\n";
              print FH "samtools sort -@ 8 \\\n$outputDir/bam/$sampleName.s.fm.bam \\\n-o $outputDir/bam/$sampleName.s.fm.sort.bam\n\n";
              ###
              ### Mark Duplicates 
              print FH "samtools markdup -r -s \\\n$outputDir/bam/$sampleName.s.fm.sort.bam \\\n$outputDir/bam/$sampleName.s.fm.sort.md.bam\n\n";
              ##
              print FH "samtools sort -@ 8 \\\n$outputDir/bam/$sampleName.s.fm.sort.md.bam \\\n-o $outputDir/bam/$sampleName.final.sort.bam\n\n";
              ##
              # Generate index
              print FH "samtools index $outputDir/bam/$sampleName.final.sort.bam\n\n";
              
              # Run bamCoverage 
              print FH "bamCoverage --bam $outputDir/bam/$sampleName.final.sort.bam  -o $outputDir/bamCoverage/$sampleName.RPKM.bedGraph \\\n--binSize 10 \\\n--blackListFileName $script_path/mm10/mm10.blacklist_merged.bed \\\n--effectiveGenomeSize 2652783500 \\\n--normalizeUsing RPKM \\\n--outFileFormat \"bedgraph\" \\\n-p \"max/2\"\n\n";
              print FH "bamCoverage --bam $outputDir/bam/$sampleName.final.sort.bam  -o $outputDir/bamCoverage/$sampleName.RPKM.bw \\\n--binSize 10 \\\n--blackListFileName $script_path/mm10/mm10.blacklist_merged.bed \\\n--effectiveGenomeSize 2652783500 \\\n--normalizeUsing RPKM \\\n--outFileFormat \"bigwig\" \\\n-p \"max/2\"\n\n";
              
              print FH "rm $outputDir/bam/$sampleName.bam\n";
              print FH "rm $outputDir/bam/$sampleName.s.fm.bam\n";
              print FH "rm $outputDir/bam/$sampleName.s.fm.sort.bam\n";
              `sbatch $out_filename`;
       }
}
