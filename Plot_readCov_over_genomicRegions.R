library(rtracklayer)

bw_data = list.files("/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/Eric_New/For_Gyan/",pattern = "*.bw",full.names = TRUE)[-3]

target_path = "/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/Eric_New/For_Eric/2wk_S835A_WT_D3a_Dnmt3a_merged_vs_IgG/"

# target regions 
target <- read.csv(paste0(target_path,"2wk_S835A_WT_D3a_Dnmt3a_merged_vs_IgG_peaks.narrowPeak_noBlacklist.bed"),sep = "\t",header = F) %>% 
  dplyr::select(1,2,3) %>% distinct()
colnames(target) = c("chr","start","end")
target.GR = GRanges(target)

read_over_GRanges = function(target, bw_data){
  bw_data = import(bw_data)
  # Find overlaps between the BEDGraph data and the target GRanges
  overlaps <- findOverlaps(target, bw_data)
  
  # Get the overlapping ranges from the BEDGraph data
  bw_hits <- bw_data[subjectHits(overlaps)]
  
  # Get the target GRanges regions corresponding to the overlaps
  target_hits <- target[queryHits(overlaps)]
  
  # Extract scores from the BEDGraph data
  scores <- mcols(bw_hits)$score
  
  # Summarize scores for each target region (sum, mean, etc.)
  # Example: Sum scores for each region
  quantified_scores <- tapply(scores, queryHits(overlaps), sum)
  
  # Add the quantified scores to the original target GRanges
  mcols(target)$score_sum <- NA
  mcols(target)$score_sum[as.numeric(names(quantified_scores))] <- quantified_scores
  return(target)
}

score_over_D3a_peaks = lapply(bw_data, function(x) { read_over_GRanges(target = target.GR, bw_data = x)})
names(score_over_D3a_peaks) = basename(bw_data)
score_over_D3a_peaks.df = lapply(score_over_D3a_peaks,as.data.frame)

# combine all the samples 
rbindlist(score_over_D3a_peaks.df,use.names = TRUE,idcol = "Group1") %>% 
  mutate(Group2 = case_when(grepl('CMVWT',Group1)~'WT',
                            grepl('CMVKO',Group1)~'cKO'),
         Group2 = factor(Group2, levels = c("WT","cKO"))) %>% 
  #reshape2::dcast(Group2 ~  Group1, value.var = "score_sum") %>% head() 
  ggplot(aes(x=Group2,y=log10(score_sum+1)))+
    geom_boxplot() + 
    stat_compare_means() +
    ylab(("log10( Read Density +1)"))+
    xlab("")+
    theme_test(20)
ggsave("/project/OBI/Neuroinformatics_Core/Stroud_lab/shared/Eric_New/For_Eric/Myt1L_over_D3a_sites.pdf",height = 10,width = 8)
          

# 




