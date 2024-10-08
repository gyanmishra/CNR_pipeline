#### CUT&RUN/ChIPseq data processing

The quality control of FASTQ files was performed using the FASTQC tool (PMID: 24834071). After quality check the reads were trimmed using Trimmomatic (PMID: 24695404) and aligned to the mouse reference genome mm10 using bowtie2 with the following argument (--no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700) (PMID: 22388286). The duplicate reads were removed using samtools (PMID: 19505943). The final bigWig/bedGraph file was generated using bamCoverage (deeptools) (PMID: 24799436) utilizing RPKM methods in 10bp bins.

#### Peak calling

Peak calling was performed using macs2 (PMID: 18798982). mm10 blacklisted regions were removed from the peak files.

#### Motif analysis

The motif analysis were performed using HOMER (PMID: 20513432) with size -50,50 from the center.

#### Correlation of CUT&RUN vs methylation over genic regions. 

Based on the method reported in a previous studies, we correlated [MECP2] CUT&RUN with methylation within the gene body. We quantified the DNA methylation for [CA/CG] dinucleotide using bisulfite sequencing of [10-week-old] mouse cortex for the protein-coding genes with the longest isoform for mm10. The first 3kb of each gene were excluded for mCA/CA calculation due to the presence of a relatively low level of DNA methylation in the promoter region and considered genes with length >4.5kb. To generate smooth-line correlation plots genes were sorted based on the level of CA methylation from highest to lowest. The average [MECP2] RPKM and CA methylation were calculated in a sliding window of bin size of [800] and step sizes of [80].
