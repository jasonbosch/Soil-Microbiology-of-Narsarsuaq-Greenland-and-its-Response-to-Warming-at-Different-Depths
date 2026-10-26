#Jason Bosch

#####SET UP#####

#BioConductor libraries
library("dada2") #Denoising to ASVs

###Directory###
setwd("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/")

#####LOAD DATA#####

# Get joined and trimmed reads
reads_16S <- sort(list.files("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/16S_Cleaned/", pattern="_Cut_Trimmed.fastq", full.names = TRUE))

#Get the summary of previous processing
summary_16S <- read.table("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/16S_Cleaned/16S_Processing_Summary.csv", header = T, sep = ",")

#####ANALYSIS#####

#DADA needs to be run twice, once for each library

for (LIBRARY in c("LibD28","LibD35")) {
  
  # extract sample names
  sample_names <- unique(gsub(paste("^.+",LIBRARY,"_",sep = ""),"",gsub("_R.+$","",reads_16S[grepl(LIBRARY,reads_16S)])))
  
  ##Denoising 16S Sequences into ASVs with DADA
  
  #Estimate the error rates
  err_fwd <- learnErrors(reads_16S[grepl(paste(LIBRARY,".+","R1",sep = "_"),reads_16S)], multithread=7)
  err_rev <- learnErrors(reads_16S[grepl(paste(LIBRARY,".+","R2",sep = "_"),reads_16S)], multithread=7)
  
  #Plot the errors, to verify if everything is ok.
  plot_err_fwd <- plotErrors(err_fwd, nominalQ=TRUE)
  svg(paste("errors_forward_",LIBRARY,".svg",sep = ""), width = 2*84*0.04, height = 2*84*0.04)
  plot_err_fwd
  dev.off()
  plot_err_rev <- plotErrors(err_rev, nominalQ=TRUE)
  svg(paste("errors_reverse_",LIBRARY,".svg",sep = ""), width = 2*84*0.04, height = 2*84*0.04)
  plot_err_rev
  dev.off()
  
  #Perform denoising
  dada_fwds <- dada(reads_16S[grepl(paste(LIBRARY,".+","R1",sep = "_"),reads_16S)], err=err_fwd, multithread=7, pool="pseudo")
  dada_revs <- dada(reads_16S[grepl(paste(LIBRARY,".+","R2",sep = "_"),reads_16S)], err=err_rev, multithread=7, pool="pseudo")
  
  #Merge paired end reads
  mergers <- mergePairs(dada_fwds, reads_16S[grepl(paste(LIBRARY,".+","R1",sep = "_"),reads_16S)], 
                        dada_revs, reads_16S[grepl(paste(LIBRARY,".+","R2",sep = "_"),reads_16S)], 
                        verbose=FALSE, returnRejects=FALSE, minOverlap=12, maxMismatch=0)

  # Create the DADA sequence table
  DADA_16S_Counts_Table <- makeSequenceTable(mergers)
  
  # Remove chimeras
  DADA_16S_Counts_Table_No_Chimera <- removeBimeraDenovo(DADA_16S_Counts_Table, method="consensus", multithread=7, verbose=FALSE)
  
  # get the table that summarises how many reads were discarded for each filtering step
  summary_16S_R <- cbind.data.frame(paste(LIBRARY,sample_names,sep = "_"), rowSums(DADA_16S_Counts_Table), rowSums(DADA_16S_Counts_Table_No_Chimera))
  colnames(summary_16S_R) <- c("Sample", "Merged-Denoised", "Non-chimera")
  
  
  #Merge with the original summary
  summary_16S_Final <- merge.data.frame(summary_16S, summary_16S_R)
  summary_16S_Final$Percent <- signif(summary_16S_Final$`Non-chimera`/summary_16S_Final$Forward*100,4)
  # export filtering information
  write.table(summary_16S_Final, paste("summary_16S_",LIBRARY,"_Final.csv"), col.names = T, row.names = F, sep = ",")
  
  # export the final dada2 table clean from chimeras
  saveRDS(DADA_16S_Counts_Table_No_Chimera, paste("DADA_16S_",LIBRARY,"Counts_Table_No_Chimera.Rds",sep = ""))
  
  # Assign taxonomy to ASVs
  # Using DADA2 silva 138.2 sequences: https://zenodo.org/records/14169026
  set.seed(123456)
  DADA_16S_ASV_Taxonomy <- assignTaxonomy(DADA_16S_Counts_Table_No_Chimera, "~/PostDoc/04_Software/SILVA_138.2_DADA2/silva_nr99_v138.2_toSpecies_trainset.fa.gz", minBoot=70, outputBootstraps=TRUE, tryRC=TRUE, multithread=7)
  
  # Export taxonomy table and bootstrap scores
  saveRDS(DADA_16S_ASV_Taxonomy[[1]], paste("DADA_16S_",LIBRARY,"ASV_Taxonomy.Rds",sep = ""))
  saveRDS(DADA_16S_ASV_Taxonomy[[2]], paste("DADA_16S_",LIBRARY,"ASV_Taxonomy_Scores.Rds",sep = ""))
}

#Session info
# R version 4.6.0 (2026-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Arch Linux
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblas.so.0.3;  LAPACK version 3.12.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Prague
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dada2_1.40.0   Rcpp_1.1.1-1.1
# 
# loaded via a namespace (and not attached):
# [1] generics_0.1.4              SparseArray_1.12.2          bitops_1.0-9                stringi_1.8.7              
# [5] jpeg_0.1-11                 lattice_0.22-9              magrittr_2.0.5              grid_4.6.0                 
# [9] RColorBrewer_1.1-3          plyr_1.8.9                  Matrix_1.7-5                cigarillo_1.2.0            
# [13] scales_1.4.0                Biostrings_2.80.0           codetools_0.2-20            abind_1.4-8                
# [17] cli_3.6.6                   ShortRead_1.70.0            rlang_1.2.0                 crayon_1.5.3               
# [21] XVector_0.52.0              Biobase_2.72.0              DelayedArray_0.38.1         S4Arrays_1.12.0            
# [25] tools_4.6.0                 parallel_4.6.0              reshape2_1.4.5              deldir_2.0-4               
# [29] BiocParallel_1.46.0         dplyr_1.2.1                 interp_1.1-6                ggplot2_4.0.3              
# [33] Rsamtools_2.28.0            hwriter_1.3.2.1             SummarizedExperiment_1.42.0 BiocGenerics_0.58.0        
# [37] vctrs_0.7.3                 R6_2.6.1                    png_0.1-9                   matrixStats_1.5.0          
# [41] stats4_4.6.0                lifecycle_1.0.5             stringr_1.6.0               pwalign_1.8.0              
# [45] Seqinfo_1.2.0               S4Vectors_0.50.0            IRanges_2.46.0              pkgconfig_2.0.3            
# [49] RcppParallel_5.1.11-2       pillar_1.11.1               gtable_0.3.6                glue_1.8.1                 
# [53] tibble_3.3.1                GenomicAlignments_1.48.0    GenomicRanges_1.64.0        tidyselect_1.2.1           
# [57] rstudioapi_0.18.0           MatrixGenerics_1.24.0       latticeExtra_0.6-31         farver_2.1.2               
# [61] compiler_4.6.0              S7_0.2.2      