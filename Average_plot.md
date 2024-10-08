This document outlines the commands to generate average profile plot and heatmap of CUT&RUN/ChIP-seq signal over genomic regions of interest using [deeptools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) in UTSW BioHPC. 

first reserve a node in BioHPC

```
$ srun --partition=32GB --nodes=1 --pty --time=5:00:00 /bin/bash
```

```
$ module load deeptools/3.5.0
```

**For genomic regions of interest (e.g peak file)**
there are many optional argument for each program that can be added in below command based on the analysis requirement. see [deeptools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html)
first generate matrix using computeMatrix (e.g)
```
computeMatrix reference-point \
-S 7wk_CTX_MECP2WT_Mecp2_CNR_HS031219_N031419.bw \
7wk_CTX_MECP2KO_Mecp2_CNR_HS031219_N031419.bw \
-R  GenomicRegion_of_interest.bed \
--referencePoint center \
       -p "max/2" \
       -b 2000 -a 2000 \
       --skipZeros \
       --missingDataAsZero \
       -o GenomicRegion_of_interest.gz
```

then generate either average plot or heatmap or both

```
plotHeatmap -m GenomicRegion_of_interest.gz \
--missingDataColor 'white' \
-out GenomicRegion_of_interest.pdf \
--colorMap RdYlBu \
--refPointLabel "center" \
--samplesLabel "7wk MECP2 WT" "7wk MECP2 KO" \
--heatmapHeight 3
--whatToShow "plot and heatmap" \ # other option "heatmap only" or "heatmap and colorbar"
```

only Average profile

```
plotProfile -m GenomicRegion_of_interest.gz \
-out GenomicRegion_of_interest.pdf \
--samplesLabel "7wk MECP2 WT" "7wk MECP2 KO" \
--perGroup \
--yAxisLabel "Read Density" \
--xAxisLabel "" \
--plotTitle ""
```

**For genomic regions of interest (e.g genes bed file)**

First download the gtf file for mouse genome 

```
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
```

```
# generate geneBody bed file for gencode annotation
zcat Mus_musculus.GRCm38.102.gtf.gz| awk '{if($3=="gene") print $0}' Mus_musculus.GRCm38.102.gtf | \
grep 'protein_coding' | cut -f 1,4,5,7,9 | \
cut -d ';' -f 1 | sed 's/gene_id "//' - | \
sed 's/"//' | \
awk '{ print "chr"$1,$2,$3,$5,$4}' OFS="\t" >Mus_musculus.GRCm38.102.protein_coding_gene.bed
```

first generate matrix using computeMatrix (e.g)
```
computeMatrix scale-regions \
-S 2wk_S835A_HOM_F_CTX_D3a_vs_Rb_IgG.bw \
2wk_S835A_WT_F_CTX_D3a_vs_Rb_IgG.bw \
-R  Mus_musculus.GRCm38.102.protein_coding_gene.bed \
-p "max/2" \
--regionBodyLength 5000 \
--binSize 10 \
-b 5000 -a 5000 \
--skipZeros \
-p "max/2" \
-o Mus_musculus.GRCm38.102.protein_coding_gene.gz
```

Generate Average profile
```
plotProfile -m Mus_musculus.GRCm38.102.protein_coding_gene.gz \
-out Mus_musculus.GRCm38.102.protein_coding_gene.pdf \
--samplesLabel "2wk S835A HOM CTX D3a" "2wk S835A WT CTX D3a" \
--perGroup \
--yAxisLabel "D3a Read Density" \
--plotTitle ""
```

