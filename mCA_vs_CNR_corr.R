
########################################################################################################################
# Prepare input files for 

# methylation data
mCA_target = "/work/OBI/Neuroinformatics_Core/s225347/DNMT3A_project/data/BSmap/2wk_CTX_Stroud2017_mm10_BSmap_CA.txt.gz"

# CNR bam files (only one sample)
CNR_target = "/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/Eric_New/For_Gyan/210930_H4_S835A_MECP2.final.sort.bam"

# multiple CNR samples
CNR_target = list.files(path = "/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/Eric_New/For_Gyan",
                        pattern = "210930",full.names = TRUE)

# Also can be defined as 
# CNR_target = c("/path1/file.bam","path2/file2.bam")

# Mus musculus GTF file
# Download the file from ensembl
# https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/
gtf = '/work/Neuroinformatics_Core/s225347/Pipeline/CNR_pipeline/Mus_musculus.GRCm38.102.gtf.gz'

########################################################################################################################
# Define bin and step size
binSize <- 399
stepSize <- 40

########################################################################################################################
# Function to check if the required packages is installed or not.
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  } else {
    message(paste(pkg, "is already installed"))
  }
}

# required packages
packages = c("BiocManager","tidyverse","Rsamtools", "GenomicAlignments",
             "GenomicRanges", "SummarizedExperiment","methylKit")

# Apply the function to each package
sapply(packages, install_if_missing)

# Load the required packages
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(methylKit)

# set this path where CNR_pipeline is located.
script_path = "/work/OBI/Neuroinformatics_Core/s225347/Pipeline/CNR_pipeline/"
source(paste0(script_path,"/readBismarkFiles.R"))

# function to calculate average coverage over step and bin size
summarize_logFC_mCA = function(x,binSize = 19,stepSize =4){
  df.bin =list()
  for (i in seq(from=1, to=nrow(x), by=stepSize)){
    print(i)
    start = i
    stop  = i+binSize
    if(stop <= nrow(x)){
      df = x[start:stop,]
      df.bin[[i]] = df  %>% summarize(cov.mean=mean(.[,1]),cov.mean.sd = sd(.[,1]),mCA.mean=mean(.[,2]))
    }
    else{
      return(df.bin)
    }
  }
}

# Function to calculate RPKM from bam file
bamFileRPKM =  function(bamfile){
  sample_name = basename(bamfile)
  gene.GR = protein_coding.noBlacklist
  # Count the number of reads overlapping each gene region
  se <- summarizeOverlaps(features = gene.GR,
                          reads = BamFile(bamfile),
                          mode = "Union",
                          singleEnd = FALSE,  # Set this to FALSE for paired-end data
                          ignore.strand = TRUE)
  
  # Extract the counts
  counts <- assay(se)
  # View the counts
  #counts
  # Calculate gene lengths in kilobases
  gene_lengths_kb <- width(protein_coding.noBlacklist) / 1000
  # Total number of mapped reads (sum of all counts)
  total_mapped_reads <- sum(counts)
  # Total number of mapped reads in millions
  total_mapped_reads_million <- total_mapped_reads / 1e6
  # Calculate RPKM for each gene
  rpkm <- counts / (gene_lengths_kb * total_mapped_reads_million)
  # View the RPKM values
  #rpkm
  # Add RPKM to the metadata of the GRanges object
  mcols(gene.GR)$RPKM <- rpkm
  # View the updated GRanges object
  #protein_coding
  mcols(gene.GR) = mcols(gene.GR)[,'RPKM']
  return(gene.GR)
}

########################################################################################################################

# Prepare the genomic coordinates of protein coding genes and filter out genes less that 4.5kb and also removes 
# first 3kb of region from start of gene.
gtf_data <- import(gtf, format = "gtf")
protein_coding = gtf_data[mcols(gtf_data)$type == "gene"] %>% 
  .[!(is.na(mcols(.)$gene_biotype))] %>% 
  .[mcols(.)$gene_biotype == 'protein_coding'] %>% .[width(.) >4500]
new_width = width(protein_coding) -3000
protein_coding = resize(protein_coding, width = new_width, fix = ifelse(strand(protein_coding) == "+", "start", "end"))
seqlevels(protein_coding) <- paste0("chr", seqlevels(protein_coding))


GR.colnames = c('chr','start','end')
mm10.blacklisted.Regions = read.csv(paste0(script_path,"/mm10/mm10.blacklist.bed"),sep="\t") %>% 
  rename_with(~GR.colnames) %>% GRanges()

protein_coding.noBlacklist = subsetByOverlaps(protein_coding, mm10.blacklisted.Regions, invert = TRUE)
mcols(protein_coding.noBlacklist) <- NULL

########################################################################################################################
# apply the function to calculate the RPKM for each sample in CNR_target object
bamList = list()

for(i in 1:length(CNR_target)){
  rpkm = bamFileRPKM(CNR_target[i])
  bamList[[basename(CNR_target[i])]] = rpkm %>% as.data.frame(check.names = FALSE) %>% 
    mutate(ID = paste0(seqnames,'_',start,"_",end))
  
}
########################################################################################################################

# Read CG or CA methylation BSmap file
meth  = readBismarkCoverage(mCA_target,assembly="unknown",min.cov=5,sample.id = "CA")
meth_at_target = regionCounts(meth,protein_coding.noBlacklist)

# calulate methylation over gene body coordiantes
#head(meth_at_target)
meth_at_target.df = getData(meth_at_target)
meth_at_target.df$mCA = meth_at_target.df$numCs/meth_at_target.df$coverage
meth_at_target.df$ID = paste0(meth_at_target$chr,'_',meth_at_target.df$start,"_",meth_at_target.df$end)


########################################################################################################################
# merge methylation and RPKM values for each gene coordiantes
meth_cnr_merged.list = list()
for(i in 1:length(bamList)){
  meth_cnr = merge(bamList[[i]],meth_at_target.df,by='ID') %>% 
    filter(!is.na(.[,7])) %>% arrange(desc(mCA))
  meth_cnr_merged.list[[names(bamList)[i]]] = meth_cnr[,c(7,15)]
}

meth_cnr_merged.list = lapply(meth_cnr_merged.list, function(x){ x %>% mutate_if(is.character, as.numeric)})

meth_cnr_merged.list = lapply(meth_cnr_merged.list, function(x) x %>% 
                                arrange(desc(mCA)) %>% filter(.[,1] >0))

########################################################################################################################

# calculate average RPKM and mCA of genes in selected bin and step size
tmp = lapply(meth_cnr_merged.list, function(x) 
  summarize_logFC_mCA(x %>% as.data.frame(),binSize = binSize, stepSize = stepSize))
tmp1 = lapply(tmp, function(x) do.call(rbind,x))
tmp2 = plyr::ldply(tmp1,data.frame)
colnames(tmp2)[1] = 'CNR'
meth_cnr_merged.avg = tmp2[,c(2,3,4,1)]


# calculate pearson correlation between mCA and MECP2 on geneBody
# mCA_cnr.corr.df = rbindlist(lapply(
# lapply(meth_cnr_merged.avg, function(x) cor.test(x[,1],x[,2])), tidy),#
# use.names = TRUE, idcol = 'Sample')


########################################################################################################################
# Plot mCA vs CNR  (use below script to customize the plot)
cor.plot = meth_cnr_merged.avg %>% 
  #mutate(CNR= gsub('.final.sort.bam','',CNR)) %>% 
  #mutate(CNR=case_when(grepl('18',CNR) ~ '230405_18_WT_MECP2', grepl('13',CNR)~'230405_13_S835A_MECP2',TRUE ~ CNR)) %>% 
  #filter(grepl('210930',CNR)) %>% 
  #mutate(week = case_when(Age = grepl('210930',CNR)~'2wk cortex',grepl('230405',CNR)~'11wk cortex')) %>% #filter(cov.mean < 2) %>% #head()
  #mutate(CNR = gsub('.*[0-9]_','',CNR)) %>%
  #separate(CNR, c('Genotype', 'Factor'),sep = "_") %>% 
  #mutate(Genotype = factor(Genotype,levels=c("WT","S835A")),
  #        week = factor(week, levels= c('2wk cortex','11wk cortex'))) %>% 
  ggplot(aes(x=mCA.mean,y=cov.mean,color=CNR))+
  #geom_point(size=0.5)+
  geom_line()+
  #facet_wrap(~CNR)+
  #ylim(0,1.5)+
  #coord_cartesian(ylim = c(0, 1.5))+
  geom_ribbon(aes(x=mCA.mean,y=cov.mean,
                  ymin=cov.mean-cov.mean.sd,
                  ymax=cov.mean+cov.mean.sd,fill=CNR),
              alpha=0.1,color='NA') +
  #scale_fill_manual(name = "MECP2 CNR",
  #                  values = c('red','black'),
  #                  aesthetics = c("colour", "fill"),
  #                  labels=c('S835A','WT'))+
  #
  theme_test(20)+
  xlab('mCA/CA')+
  ylab('CUT & RUN\nGene Body coverage ')+
  ggtitle('11wk cortex')
cor.plot

# Add path with filename where you want to save the plot
ggsave("/work/OBI/Neuroinformatics_Core/s225347/DNMT3A_project/Figures/MECP2_CNR_vs_2wk_mCA_210930_WT_S835A.pdf",
       height = 18,width = 25,units = 'cm')
