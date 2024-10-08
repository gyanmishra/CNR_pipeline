#!/usr/bin/perl

use POSIX qw(strftime);
use FindBin;
use Cwd 'abs_path';
use Getopt::Long;

$stats = 'qvalue';
$threshold = 0.05;

# Define options
GetOptions(
    # Project Directory
    'resultDir=s' => \$outputDir,
    'pvalue' => \$pvalue,  # String option
    'qvalue' => \$qvalue,  # String option
    'threshold=f'  => \$threshold,   # Integer option
    'sampleSheet=s' => \$sampleData, 
    'help'   => \$help   # Boolean flag (if present, sets to true)
) or die "Error in command line arguments\n";

if ($help) {
    print "Usage: macs2_callPeak.pl --resultDir <result dir> --qvalue or --pvalue --threshold <set pvalue or qvalue threshold>\n --sampleSheet <macs2 peak calling sampleSheet>\n";
    exit;
}

if ($pvalue){
    $stats = 'pvalue';
}
if ($qvalue){
    $stats = 'qvalue';
}
# Get the absolute path of the current script
my $script_path = abs_path($FindBin::Bin);
#print "$script_path\n";


# provide a tab separated file having three column 
# 1st column : peakCallingID 
# 2nd column : targetFiles 
# 3rd column : controlFiles

$datestring = strftime "%F", localtime;
$datestring =~ s/-//g ;

#$sampleData  = $ARGV[2];
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
            @col = split('\t',$line);
            chomp($col[0]);
            chomp($col[1]);
            chomp($col[2]);

            @target = split(',',$col[1]);
            @control = split(',',$col[2]);
            $t_length = scalar @target;
            $c_length = scalar @control;

            #print "$t_length\n";
            @targets =[];
            @controls =[];

            $check_file = 0;

            for($i=0;$i<$t_length;$i++){
                $targets[$i] = $outputDir . "/bam/" . $target[$i] . '.final.sort.bam';
                if (-e $targets[$i]) {
                    print "$targets[$i] exists!\n";
                } else {
                    #die "$targets[$i] does not exist!\nPlease make sure the correct result folder path has been used\n";
                }

            }
            for($i=0;$i<$c_length;$i++){
                $controls[$i] = $outputDir . "/bam/" . $control[$i] . '.final.sort.bam';
                if (-e $controls[$i]) {
                    print "$controls[$i] exists!\n";
                } else {
                    #die "$controls[$i] does not exist!\nPlease make sure the correct result folder path has been used\n";
                }
            }

            unless (-e "$outputDir/macs2_callPeak") {
                `mkdir -p $outputDir/macs2_callPeak`;
            }
            unless (-e "$outputDir/logs") {
                `mkdir -p $outputDir/logs`;
            }

            #print "@controls\n";
            #print "@targets\n";

            $sampleName = $col[0];
            print "$sampleName\n";
            #$sampleName =~ s/\s/_/g;
            $out_filename = $outputDir . "/logs/" . $sampleName . '_' . $datestring . '_macs2.sh';      
            print "$out_filename\n";
            #print "$sampleName\n";
            open(FH, '>', $out_filename) or die $!;
            print FH "#!/bin/bash\n\n";
            print FH "#SBATCH --job-name $sampleName\n";
            print FH "#SBATCH -N 1\n";
            print FH "#SBATCH --cpus-per-task=32\n";
            print FH "#SBATCH --partition=32GB\n";
            print FH "#SBATCH -t 0-12:0:0\n";
            print FH "#SBATCH -o $outputDir/logs/${sampleName}_${datestring}_macs2.out\n";
            print FH "#SBATCH -e $outputDir/logs/${sampleName}_${datestring}_macs2.err\n";
            print FH "#SBATCH --mail-type ALL\n";
            #print FH "#SBATCH --mail-user GyanPrakash.Mishra\@UTSouthwestern.edu\n\n\n";
            print FH "module load macs/2.1.2 bedops/2.4.14 \n\n";
            #print FH "perl /work/OBI/Neuroinformatics_Core/s225347/CUT_and_RUN/scripts/CNR_pipeline.pl /work/OBI/Neuroinformatics_Core/s225347/CUT_and_RUN/data/CUTandRUN_seq_Data.txt\n\n";
            print FH "\n\n";
            print FH "# peak calling \n";         
            # Merge samples fa          
            if($stats eq 'pvalue'){
                print FH "macs2 callpeak \\\n-t @targets \\\n-c @controls \\\n-f BAMPE \\\n-g mm \\\n-p $threshold \\\n--outdir $outputDir/macs2_callPeak/$sampleName \\\n--call-summits \\\n-n $sampleName\n\n";
                #print FH "#remove blacklisted region\n";
                #print FH "bedops -n 1 $outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak \\\n$script_path/mm10/mm10.blacklist.bed >$outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak.noBlacklist.bed\n\n";
                #print FH "bedops -n 1 $outputDir/macs2_callPeak/$sampleName/${sampleName}_summits.bed \\\n$script_path/mm10/mm10.blacklist.bed >$outputDir/macs2_callPeak/$sampleName/${sampleName}_summits_noBlacklist.bed\n";
            }

            if($stats eq 'qvalue'){
                print FH "macs2 callpeak \\\n-t @targets \\\n-c @controls \\\n-f BAMPE \\\n-g mm \\\n-q $threshold \\\n--outdir $outputDir/macs2_callPeak/$sampleName \\\n--call-summits \\\n-n $sampleName\n\n";
            }    
            
            print FH "#remove blacklisted region\n";
            print FH "bedops -n 1 $outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak \\\n$script_path/mm10/mm10.blacklist.bed >$outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak.noBlacklist.bed\n\n";
            print FH "bedops -n 1 $outputDir/macs2_callPeak/$sampleName/${sampleName}_summits.bed \\\n$script_path/mm10/mm10.blacklist.bed >$outputDir/macs2_callPeak/$sampleName/${sampleName}_summits_noBlacklist.bed\n\n";
            print FH "cut -f 1,2,3 $outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak.noBlacklist.bed | sort-bed - | bedops -m - > $outputDir/macs2_callPeak/$sampleName/${sampleName}_peaks.narrowPeak.noBlacklist_NR.bed";
            `sbatch $out_filename`;  
    }
}
            
