#Jason Bosch

#####SET UP#####

#BioConductor libraries
library("dada2") #Denoising to ASVs

###Directory###
setwd("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/")

#####LOAD DATA#####

#Import the OTU sequences
OTUs <- read.csv("/run/media/SharedStorage/MBU_Postdoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/ITS_Cleaned/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs.fasta", header=F)
#Get the sequences only
OTU_Sequences <- OTUs[seq(2, nrow(OTUs), by=2), ]

#Get the summary of previous processing
summary_ITS <- read.table("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/ITS_Cleaned/ITS_Processing_Summary.csv", header = T, sep = ",")

#####ANALYSIS#####

# run assignment for OTUs
#Using UNITE 10.0 All_Euk: https://doi.plutof.ut.ee/doi/10.15156/BIO/3301231
Taxa_OTU <- assignTaxonomy(OTU_Sequences, "/home/jason/PostDoc/04_Software/UNITE/sh_general_release_all_19.02.2025/sh_general_release_dynamic_all_19.02.2025.fasta", minBoot=70, outputBootstraps=TRUE, tryRC=TRUE, multithread=6)

# Export taxonomy table and bootstrap scores
saveRDS(Taxa_OTU[[1]], "VSEARCH_ITS_OTU_Taxonomy.Rds")
saveRDS(Taxa_OTU[[2]], "VSEARCH_ITS_OTU_Taxonomy_Scores.Rds")

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