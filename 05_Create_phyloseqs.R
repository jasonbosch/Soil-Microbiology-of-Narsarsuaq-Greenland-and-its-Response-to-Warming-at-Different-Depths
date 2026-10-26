#Jason Bosch

###Libraries###

#CRAN libraries
library("stringr") #To convert text to uppercase

#BioConductor libraries
library("phyloseq") #For phyloseq objects

###Directory###

setwd("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/")

#####LOAD DATA#####

##Get metadata
Metadata <- read.table("~/PostDoc/02_Projects/13_Narsarsuaq/06_Miscellaneous/Metadata_in_use/metadata_narsarsuaq_2020_final.csv",sep="\t",header = T)
Metadata_Final <- sample_data(Metadata)
sample_names(Metadata_Final) <- Metadata_Final$sampleid
# #Extra metadata
# metadata_substrate <- read.csv("~/PostDoc/02_Projects/13_Narsarsuaq/06_Miscellaneous/Metadata_in_use/Narsarsuaq_easy substrate_metadata.csv",sep=",")

##16S data

#LibD28
#Import the ASV table
Counts_16S_LibD28 <- as.data.frame(t(readRDS("DADA_16S_LibD28Counts_Table_No_Chimera.Rds")))
#Rename columns to match sample names
colnames(Counts_16S_LibD28) <- gsub("CAM0+","S",sapply(strsplit(colnames(Counts_16S_LibD28), "_"), "[", 2))
#Import the taxonomy
Taxonomy_16S_LibD28 <- readRDS("DADA_16S_LibD28ASV_Taxonomy.Rds")
Taxonomy_16S_LibD28_final <- tax_table(Taxonomy_16S_LibD28)
#Create phyloseq object
Phyloseq_16S_LibD28 <- phyloseq(otu_table(Counts_16S_LibD28,taxa_are_rows = TRUE), Taxonomy_16S_LibD28_final, Metadata_Final)

#LibD35
#Import the ASV table
Counts_16S_LibD35 <- as.data.frame(t(readRDS("DADA_16S_LibD35Counts_Table_No_Chimera.Rds")))
#Rename columns to match sample names
colnames(Counts_16S_LibD35) <- gsub("CAM0+","S",sapply(strsplit(colnames(Counts_16S_LibD35), "_"), "[", 2))
#Import the taxonomy
Taxonomy_16S_LibD35 <- readRDS("DADA_16S_LibD35ASV_Taxonomy.Rds")
Taxonomy_16S_LibD35_final <- tax_table(Taxonomy_16S_LibD35)
#Create phyloseq object
Phyloseq_16S_LibD35 <- phyloseq(otu_table(Counts_16S_LibD35,taxa_are_rows = TRUE), Taxonomy_16S_LibD35_final, Metadata_Final)

##ITS data

#Import the OTU table
Counts_ITS <- read.table("ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs_Table_Sorted.txt",sep = "\t", header = T)
#Split the two libraries to avoid identical column names
Counts_ITS_LibD28 <- Counts_ITS[,grepl("OTU.ID|LibD28",colnames(Counts_ITS))]
Counts_ITS_LibD37 <- Counts_ITS[,grepl("OTU.ID|LibD37",colnames(Counts_ITS))]
#Rename columns to match sample names
colnames(Counts_ITS_LibD28) <- gsub("ALG0+","S",sapply(strsplit(colnames(Counts_ITS_LibD28), "_"), "[", 2))
colnames(Counts_ITS_LibD28)[1] <- "OTU.ID"
colnames(Counts_ITS_LibD37) <- gsub("ALG0+","S",sapply(strsplit(colnames(Counts_ITS_LibD37), "_"), "[", 2))
colnames(Counts_ITS_LibD37)[1] <- "OTU.ID"
#Join the columns together, summing the shared columns
Counts_ITS_Sum <- Counts_ITS_LibD28
for (SAMPLE in colnames(Counts_ITS_LibD37)[2:ncol(Counts_ITS_LibD37)]) {
  Counts_ITS_Sum[,SAMPLE] <- Counts_ITS_Sum[,SAMPLE] + Counts_ITS_LibD37[,SAMPLE]
}
#Move OTU names to row.names
rownames(Counts_ITS_Sum) <- Counts_ITS_Sum[,"OTU.ID"]
Counts_ITS_Sum <- subset.data.frame(Counts_ITS_Sum, select = -OTU.ID)
#Import the taxonomy
Taxonomy_ITS <- readRDS("VSEARCH_ITS_OTU_Taxonomy.Rds")
#Clean the taxonomy names
for (COL in 1:ncol(Taxonomy_ITS)) {
  Taxonomy_ITS[,COL] <- gsub("[a-z]__","",Taxonomy_ITS[,COL])
}
#Find the names of the sequences
Taxonomy_Names_Raw <- read.table("../02_Cleaned_Data/ITS_Cleaned/ITS_Clusters_0.97_NoSingletons_NoChimeras_OTUs.fasta")
Taxonomy_Names_Table <- data.frame("OTU" = Taxonomy_Names_Raw[seq(1, nrow(Taxonomy_Names_Raw), by=2), ], "Seq" = Taxonomy_Names_Raw[seq(2, nrow(Taxonomy_Names_Raw), by=2), ])
Taxonomy_Names_Table[,"OTU"] <- gsub(">","",gsub(";.+","",Taxonomy_Names_Table[,"OTU"]))
Taxonomy_Names_Table[,"Seq"] <- str_to_upper(Taxonomy_Names_Table[,"Seq"])
#Rename the sequences to the OTU names
Taxonomy_ITS_final <- Taxonomy_ITS
if (identical(Taxonomy_Names_Table[,"Seq"],row.names(Taxonomy_ITS_final))) {
  row.names(Taxonomy_ITS_final) <- Taxonomy_Names_Table[,"OTU"]
}
Taxonomy_ITS_final <- tax_table(Taxonomy_ITS_final)
#Create phyloseq object
Phyloseq_ITS <- phyloseq(otu_table(Counts_ITS_Sum,taxa_are_rows = TRUE), Taxonomy_ITS_final, Metadata_Final)

#####PREPARE THE PHYLOSEQ OBJECTS FOR ANALYSIS#####

#Combine the 16S libraries
Phyloseq_16S <- merge_phyloseq(Phyloseq_16S_LibD28,Phyloseq_16S_LibD35)

#Filter any ASV with fewer than 0.1% of the mean number of reads per sample per library size
Phyloseq_16S_Filtered <- prune_taxa(rowSums(otu_table(Phyloseq_16S)) > round(mean(colSums(otu_table(Phyloseq_16S)))*0.001), Phyloseq_16S)
cat("Removed ",nrow(Phyloseq_16S@otu_table) - nrow(Phyloseq_16S_Filtered@otu_table)," taxa.") #Removed  2456  taxa.
Phyloseq_ITS_Filtered <- prune_taxa(rowSums(otu_table(Phyloseq_ITS)) > round(mean(colSums(otu_table(Phyloseq_ITS)))*0.001), Phyloseq_ITS)
cat("Removed ",nrow(Phyloseq_ITS@otu_table) - nrow(Phyloseq_ITS_Filtered@otu_table)," taxa.") #Removed  1818  taxa.

#Check the taxonomy and remove anything that is incorrectly assigned
#16S
table(Phyloseq_16S_Filtered@tax_table[,"Kingdom"]) #26 Archaea, 2379 Bacteria, 6 Eukaryota
Phyloseq_16S_Filtered <- subset_taxa(Phyloseq_16S_Filtered, Kingdom != "Eukaryota")
table(as.vector(Phyloseq_16S_Filtered@tax_table[,"Order"]=="Chloroplast")) #5
Phyloseq_16S_Filtered <- subset_taxa(Phyloseq_16S_Filtered, Order != "Chloroplast")
table(as.vector(Phyloseq_16S_Filtered@tax_table[,"Family"]=="Mitochondria")) #90
Phyloseq_16S_Filtered <- subset_taxa(Phyloseq_16S_Filtered, Family != "Mitochondria")
#Fungi
table(Phyloseq_ITS_Filtered@tax_table[,"Kingdom"])
#993 Fungi, 1 Ichthyosporia, 1 Rhizaria, 11 Viridiplantae
Phyloseq_ITS_Filtered <- subset_taxa(Phyloseq_ITS_Filtered, tax_table(Phyloseq_ITS_Filtered)[,"Kingdom"] == "Fungi")

#Save the phyloseq objects for later processing
saveRDS(Phyloseq_16S_Filtered, "Phyloseq_16S_Filtered.Rds")
saveRDS(Phyloseq_ITS_Filtered, "Phyloseq_ITS_Filtered.Rds")

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
# [1] phyloseq_1.56.0 stringr_1.6.0  
# 
# loaded via a namespace (and not attached):
# [1] generics_0.1.4      stringi_1.8.7       lattice_0.22-9      digest_0.6.39       magrittr_2.0.5      grid_4.6.0         
# [7] RColorBrewer_1.1-3  iterators_1.0.14    foreach_1.5.2       plyr_1.8.9          jsonlite_2.0.0      Matrix_1.7-5       
# [13] ape_5.8-1           survival_3.8-6      mgcv_1.9-4          scales_1.4.0        permute_0.9-10      Biostrings_2.80.0  
# [19] codetools_0.2-20    ade4_1.7-24         cli_3.6.6           rlang_1.2.0         crayon_1.5.3        XVector_0.52.0     
# [25] Biobase_2.72.0      splines_4.6.0       vegan_2.7-3         tools_4.6.0         parallel_4.6.0      reshape2_1.4.5     
# [31] dplyr_1.2.1         ggplot2_4.0.3       multtest_2.68.0     BiocGenerics_0.58.0 vctrs_0.7.3         R6_2.6.1           
# [37] stats4_4.6.0        lifecycle_1.0.5     Seqinfo_1.2.0       biomformat_1.40.0   S4Vectors_0.50.0    IRanges_2.46.0     
# [43] MASS_7.3-65         cluster_2.1.8.2     pkgconfig_2.0.3     pillar_1.11.1       gtable_0.3.6        data.table_1.18.4  
# [49] glue_1.8.1          Rcpp_1.1.1-1.1      tibble_3.3.1        tidyselect_1.2.1    rstudioapi_0.18.0   farver_2.1.2       
# [55] nlme_3.1-169        igraph_2.3.1        compiler_4.6.0      S7_0.2.2