#Jason Bosch

#####SET UP#####

###Libraries###

#CRAN libraries
library("ggplot2") #Plotting
library("ggrepel") #Plotting assistance
library("cowplot") #Combine multiple plots
library("vegan") #Various functions for microbiome analysis
library("DescTools") #For rounding up values (rares) & differences in means
library("RColorBrewer") #For colour palettes
library("MicrobiomeStat") #Linda differential abundance
library("nlme")#Linear models for soil properties
library("multcompView") #Grouping stats
library("stringr") #Using during guild assignment to correct cases
library("lubridate") #Dates
library("pheatmap") #For drawing the heatmap of differentially abundant taxa
library("emmeans") #Replacement for anova in the richness analysis
library("grid") #Needed when joining the heatmaps

#BioConductor libraries
library("phyloseq") #For phyloseq objects

#Github libraries
library("mirlyn") #For rarefaction (alpha diversity)
library("ggpubfigs") #For friendly palette
library("metagMisc") #Will use for rarefaction (average RA)

###Directory###

setwd("~/PostDoc/02_Projects/13_Narsarsuaq/04_Analysis_Results/Final/")

#####FUNCTIONS#####

#To pull p-values out of a linear model
#https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#####LOAD DATA#####

##Get metadata
Metadata <- read.table("~/PostDoc/02_Projects/13_Narsarsuaq/06_Miscellaneous/Metadata_in_use/metadata_narsarsuaq_2020_final.csv",sep="\t",header = T)
#Set the factors
Metadata$treat <- factor(Metadata$treat, levels = c("C","S","W","SW"),ordered = T)
Metadata$depth <- factor(Metadata$depth, levels = c("0_5","5_10","15_20"),ordered = T)

#Load the phyloseq objects
Physeq_16S <- readRDS("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/Phyloseq_16S_Filtered.Rds")
Physeq_ITS <- readRDS("~/PostDoc/02_Projects/13_Narsarsuaq/02_Cleaned_Data/Phyloseq_ITS_Filtered.Rds")

#Set the factors
for (TAXON in c("16S","ITS")) {
  current_physeq <- get(paste("Physeq",TAXON,sep = "_"))
  current_physeq@sam_data$treat <- factor(current_physeq@sam_data$treat, levels = c("C","S","W","SW"),ordered = T)
  current_physeq@sam_data$depth <- factor(current_physeq@sam_data$depth, levels = c("0_5","5_10","15_20"),ordered = T)
  assign(paste("Physeq",TAXON,sep = "_"),
         current_physeq 
  )
}

#Load fungal guild data
#Import Camelia's altered version of FungalTraits with better identification of ericoid mycorrhiza.
CameliaTraits <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/guilds_from_seed_changed_by_camelia_for_Disko2013Paper.csv",header = F, sep = ",")
for (COL in 1:6) {
  CameliaTraits[,COL] <- gsub("[a-z]__","",CameliaTraits[,COL])
}
#Remove American spelling
CameliaTraits$V8 <- gsub("lichenized","lichenised",CameliaTraits$V8)

#Bring in the temperature and Moisture record
NARS_Temp_Hum <- read.table("~/PostDoc/02_Projects/13_Narsarsuaq/06_Miscellaneous/Metadata_in_use/Snowfence_Narsarsuaq_probe_data_date_time.csv", header = T, sep=",")
#Use the same terminology
NARS_Temp_Hum$Treatment <- gsub("AC","C",NARS_Temp_Hum$Treatment)
NARS_Temp_Hum$Treatment <- gsub("AS","S",NARS_Temp_Hum$Treatment)
NARS_Temp_Hum$Treatment <- gsub("WC","W",NARS_Temp_Hum$Treatment)
NARS_Temp_Hum$Treatment <- gsub("WS","SW",NARS_Temp_Hum$Treatment)

#####ANALYSIS#####

###Constant colours and shapes for figures###

Colours_Fence <- friendly_pal("wong_eight",8)[c(1:5,7)]
names(Colours_Fence) <- c("A","B","C","D","E","F")
Colours_Treatment <- friendly_pal("zesty_four",4)[c(2,4,3,1)]
names(Colours_Treatment) <- c("C","S","W","SW")
Colours_Categories <- friendly_pal("contrast_three",3)
names(Colours_Categories) <- c("Microbial","Enzymatic","Soil")
Shapes_Depth <- c(19,17,15)
names(Shapes_Depth) <- c("0_5","5_10","15_20")
Colours_Depth <- friendly_pal("retro_four",4)[3:1]
names(Colours_Depth) <- c("0_5","5_10","15_20")


###Library sizes###

#16S
summary(colSums(Physeq_16S@otu_table))
hist(colSums(Physeq_16S@otu_table),breaks = seq(0,max(colSums(Physeq_16S@otu_table))+1000,1000))
#ITS
summary(colSums(Physeq_ITS@otu_table))
hist(colSums(Physeq_ITS@otu_table),breaks = seq(0,max(colSums(Physeq_ITS@otu_table))+1000,1000))

###Rarefaction Curve###

# Creation of rarefaction curve data frame
for (MARKER in c("16S","ITS")) {
  #Only need to check in the fractions of library sizes that will get us to the median since we won't plot more
  lib_frac <- median(colSums(get(paste("Physeq",MARKER,sep = "_"))@otu_table))/max(colSums(get(paste("Physeq",MARKER,sep = "_"))@otu_table))
  #Determine the step for 100 points along the curve
  step_size <- lib_frac/100
  #Rarefy at different library sizes
  assign(paste("Physeq_Rarefaction_Curve",MARKER,sep = "_"),
         rarefy_whole_rep(get(paste("Physeq",MARKER,sep = "_")), rep = 50, steps = seq(from = 0.001, to = lib_frac, by = step_size),mc.cores = 6)
  )
}

#Visualization of rarefaction curve
for (MARKER in c("16S","ITS")) {
  assign(paste("Plot_Rarefaction_Curve",MARKER,sep = "_"),
         ggplot(get(paste("Physeq_Rarefaction_Curve",MARKER,sep = "_")), aes(x = LibSize, y = ObsASVCount, group = Sample, colour = sfence)) + 
           geom_line() + 
           geom_vline(xintercept = 3500) +
           scale_y_continuous(limits = c(0,max(get(paste("Physeq_Rarefaction_Curve",MARKER,sep = "_"))$ObsASVCount) * 1.01)) +
           scale_x_continuous(limits = c(0,median(colSums(get(paste("Physeq",MARKER,sep = "_"))@otu_table)))) +
           labs(title = paste(MARKER,"rarefaction curve",sep = " "), x = "Library Size", y = "Richness") +
           scale_colour_manual(values = Colours_Fence) +
           theme_bw() + 
           theme(legend.position = "none")
  )
}

#Going to use minimum library size of 3500.
#We lose 4 ITS samples but it looks much better positioned to capture the richness.
Rare_Lib_Size <- 3500

###Alpha diversity###

#Rarefy communities
#Use mirl to repeatedly rarefy to specified library size 1000 times.
#For some reason this broke after R update!! Seems to be in mclapply, so just using 1 core.
for (MARKER in c("16S","ITS")) {
  assign(paste("Mirl_Object",MARKER,sep = "_"),
         mirl(get(paste("Physeq",MARKER,sep = "_")), 
              libsize = Rare_Lib_Size, 
              rep = 1000, mc.cores = 1L, set.seed = 123456)
  )
}

#Generate a dataframe of alpha-diversity metric from the mirl object
for (MARKER in c("16S","ITS")) {
  assign(paste("AlphaDiv",MARKER,sep = "_"),
         data.frame()
  )
  for (REP in 1:length(get(paste("Mirl_Object",MARKER,sep = "_")))) {
    cat(paste("Calculating ",MARKER," rep ",REP,"\n",sep = ""))
    md <- sample_data(get(paste("Mirl_Object",MARKER,sep = "_"))[[REP]])
    div <- estimate_richness(get(paste("Mirl_Object",MARKER,sep = "_"))[[REP]],measures = c("Observed","Shannon"))
    final <- cbind(md, div)
    assign(paste("AlphaDiv",MARKER,sep = "_"),
           rbind(get(paste("AlphaDiv",MARKER,sep = "_")),final)
    )
  }
}

#Delete the mirl objects as they are quite large and unneeded now
rm(Mirl_Object_16S,Mirl_Object_ITS)

#Visualise the variation in the rarefactions
for (MARKER in c("16S","ITS")) {
  assign(paste("Plot_Observed",MARKER,"Rarefactions",sep = "_"),
         ggplot(get(paste("AlphaDiv",MARKER,sep = "_")), aes(x=treat, y=Observed)) +
           geom_point (aes(colour = sfence),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
           geom_boxplot(fill=NA,size=0.5) +
           theme_bw() +
           scale_y_continuous(limits = c(0,RoundTo(max(get(paste("AlphaDiv",MARKER,sep = "_"))[,"Observed"]),multiple = 50,ceiling)),n.breaks = 6) +
           scale_color_manual(values=Colours_Fence) +
           facet_wrap(~depth) +
           labs(title = paste(MARKER,"richness with 1000X rarefactions",sep = " "), x = NULL, y = "Richness") +
           theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5), 
                 panel.spacing=unit(0,"lines"))
  )
}

#Mean diversity values from the 1000 rarefactions
for (MARKER in c("16S","ITS")) {
  Current_Table <- aggregate(get(paste("AlphaDiv",MARKER,sep = "_"))[,c("sampleid","Observed","Shannon")],.~sampleid,mean)
  Current_Table <- merge.data.frame(Current_Table,Metadata,by="sampleid")
  Current_Table$library_size <-  colSums(get(paste("Physeq",MARKER,sep = "_"))@otu_table)[Current_Table$sampleid]
  Current_Table$treat <- factor(Current_Table$treat, levels = c("C","S","W","SW"),ordered = F)
  Current_Table$depth <- factor(Current_Table$depth, levels = c("0_5","5_10","15_20"),ordered = T)
  assign(paste("ALPHADIV_MEAN",MARKER,sep = "_"),
         Current_Table
  )
}

#Check the correlation with sequencing depth
for (MARKER in c("16S","ITS")) {
  assign(paste("Plot_Diversity_Depth_Correlation",MARKER,sep = "_"),
         ggplot(get(paste("ALPHADIV_MEAN",MARKER,sep = "_")), aes(x=library_size, y=Observed)) +
           geom_point(aes(colour = sfence)) +
           geom_smooth(method = "lm", colour = "black") +
           annotate("text", label = paste("Adj. R^2 = ",formatC(summary(lm(Observed~library_size,get(paste("ALPHADIV_MEAN",MARKER,sep = "_"))))$adj.r.squared, format = "f", digits = 3),"\n",
                                          "p = ",formatC(lmp(lm(Observed~library_size,get(paste("ALPHADIV_MEAN",MARKER,sep = "_")))), format = "e", digits = 1),sep = ""),x = Inf, y = -Inf,hjust=1,vjust=-0.3) +
           theme_bw() +
           scale_y_continuous(n.breaks = 6) +
           scale_x_continuous(n.breaks = 6) +
           labs(x = "Library size", y = "Richness") +
           scale_color_manual(values=Colours_Fence) +
           theme(text=element_text(size=16), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
                 panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  )
}

#Linear models of alpha diversity
Alpha_Summary <- list()
Plots_Alpha <- list()
for (MARKER in c("16S","ITS")) {
  for (METRIC in c("Observed","Shannon")) {

    #Get the data
    current_data <- get(paste("ALPHADIV_MEAN",MARKER,sep = "_"))[, c("treat", "depth", "sfence", METRIC)]
    colnames(current_data) <- c("Treatment", "Depth", "Plot", METRIC)
    #Set reference level
    current_data$Depth <- factor(current_data$Depth, levels = unique(current_data$Depth),ordered = F) #Need to remove the ordering for pairs
    current_data$Depth <- relevel(current_data$Depth, ref="0_5")
    current_data$Treatment <- relevel(current_data$Treatment, ref="C")

    #Perform linear modelling
    set.seed(123456)
    current_linear_model <- lme(as.formula(paste(METRIC,"~Treatment*Depth",sep = )), random=~1|Plot, data=current_data, na.action=na.omit)
    
    #Estimate marginal means
    current_emm_interactions <- emmeans(current_linear_model, ~ Treatment * Depth)
    current_emm_depth <- emmeans(current_linear_model, ~ Depth)
    current_emm_treatment <- emmeans(current_linear_model, ~ Treatment)

    #Dataframe of results
    current_res_interactions <- data.frame(pairs(current_emm_interactions, adjust = "tukey"))
    current_res_depth <- data.frame(pairs(current_emm_depth, adjust = "tukey"))
    current_res_treatment <- data.frame(pairs(current_emm_treatment, adjust = "tukey"))
    
    #Combine results
    current_summary <- rbind(current_res_treatment,current_res_depth,current_res_interactions)

    #Convert interaction table into usable format for plotting
    current_res <- current_res_interactions$p.value
    names(current_res) <- gsub(" ",":",gsub(" - ", "-",current_res_interactions$contrast))
    current_letters <- multcompLetters(current_res)
    current_letters_df <- data.frame("comp" = names(current_letters$Letters), "Letter"=as.vector(current_letters$Letters))
    current_letters_df$treat <- gsub(":.+","",current_letters_df$comp)
    current_letters_df$depth <- gsub("^.+:{1}","",current_letters_df$comp)
    current_letters_df$treat <- factor(current_letters_df$treat, levels = c("C","S","W","SW"),ordered = F)
    current_letters_df$depth <- factor(current_letters_df$depth, levels = c("0_5","5_10","15_20"),ordered = T)
    current_letters_df[,METRIC] <- max(get(paste("ALPHADIV_MEAN",MARKER,sep = "_"))[,METRIC])*1.05

        #Create plot
        current_plot <- ggplot(get(paste("ALPHADIV_MEAN",MARKER,sep = "_")), aes(x=treat, y=!! sym(METRIC))) +
          geom_point (aes(colour = sfence),position = position_jitterdodge(jitter.width = 0.1))  +
          geom_boxplot(fill=NA,size=0.5) +
          stat_summary(geom = "point", fun = mean, shape = 18, size = 3) +
          geom_text(data = current_letters_df, aes(label=Letter)) +
          theme_bw() +
          scale_y_continuous(limits = c(0,RoundTo(max(get(paste("ALPHADIV_MEAN",MARKER,sep = "_"))[,METRIC])*1.05,multiple = ceiling(max(get(paste("ALPHADIV_MEAN",MARKER,sep = "_"))[,METRIC]*0.1)),ceiling)),n.breaks = 6) +
          scale_color_manual(values=Colours_Fence,name="Plot") +
          facet_wrap(~depth, labeller=labeller(depth = c("0_5" = "0-5","5_10"="5-10","15_20"="15-20"))) +
          labs(title = paste(METRIC," ",MARKER,sep = ""), x = NULL, y = ifelse(METRIC=="Observed","Richness",METRIC) ,colour="sfence") +
          theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16),
                panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
        
        #Save results
        Plots_Alpha[[paste(METRIC,MARKER,sep = "_")]] <- current_plot
        Alpha_Summary[[paste(METRIC,MARKER,sep = "_")]] <- current_summary

  }
}

###Beta diversity###

#Distances
#Use recommended Hellinger
for (MARKER in c("16S","ITS")) {
  set.seed(123456)
  assign(paste("Dist_Hellinger",MARKER,sep = "_"),
         as.dist(avgdist(x = as.data.frame(t(get(paste("Physeq",MARKER,sep = "_"))@otu_table)), 
                         sample = Rare_Lib_Size, 
                         distfun = vegdist, dmethod = "hellinger", iterations = 1000))
  )
}

#Check for differences in composition with permanova
#Checks whether the samples have the same centroid. Null hypothesis is same centroid.
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
#R2 shows how much is explained, residuals are unexplained variation.
Permanova_result <- NULL
for (MARKER in c("16S","ITS")) {
  #Prepare the data
  sampledf <- data.frame(sample_data(get(paste("Physeq",MARKER,sep = "_"))))
  #Remove the low library size samples
  sampledf <- sampledf[sampledf$sampleid%in%names(get(paste("Dist_Hellinger",MARKER,sep = "_")))[1][[1]],]
  #Create interaction columns
  sampledf$`treat:depth` <- interaction(sampledf$treat, sampledf$depth)
  #Perform permanova
  set.seed("123456")
  current_adonis <- adonis2(get(paste("Dist_Hellinger",MARKER,sep = "_")) ~ treat*depth, data = sampledf, permutations = 1000, by = "terms", strata = sampledf$sfence)
  #Save result
  Permanova_result[[MARKER]] <- current_adonis
  assign(paste("sampledf",MARKER,sep = "_"),
         sampledf
  )
}

#Create a table of the permanova test results
Beta_Perm_Table <- cbind(as.data.frame(Permanova_result[["16S"]]),as.data.frame(Permanova_result[["ITS"]]))
Beta_Perm_Table[,"Comparison"] <- row.names(Beta_Perm_Table)
Beta_Perm_Table[,"Beta Dispersion (p)"] <- row.names(Beta_Perm_Table)
Beta_Perm_Table[,"Beta Dispersion (p).1"] <- row.names(Beta_Perm_Table)
Beta_Perm_Table[nrow(Beta_Perm_Table)+1,] <- c(rep("Bacteria",5),rep("Fungi",5),"Comparison","Bacteria","Fungi")
Beta_Perm_Table[nrow(Beta_Perm_Table)+1,] <- c(rep(colnames(Beta_Perm_Table[1:5]),2),"Comparison", "Beta Dispersion (p)","Beta Dispersion (p)")
Beta_Perm_Table <- Beta_Perm_Table[c(6,7,1:5),c(11,1:5,12,6:10,13)]
Beta_Perm_Table[c(6,7),c(7,13)] <- NA

#Homogeneity of dispersion tests
for (GROUP in Beta_Perm_Table[3:5,"Comparison"]) {
  beta_b <- betadisper(Dist_Hellinger_16S, sampledf_16S[,GROUP])
  beta_f <- betadisper(Dist_Hellinger_ITS, sampledf_ITS[,GROUP])
  set.seed("123456")
  Beta_Perm_Table[Beta_Perm_Table$Comparison==GROUP,7] <- permutest(beta_b)$tab$`Pr(>F)`[1]
  set.seed("123456")
  Beta_Perm_Table[Beta_Perm_Table$Comparison==GROUP,13] <- permutest(beta_f)$tab$`Pr(>F)`[1]
}

#Adjust betadispers results for multiple testing
Beta_Perm_Table[3:5,7] <- p.adjust(Beta_Perm_Table[3:5,7], method = "fdr")
Beta_Perm_Table[3:5,13] <- p.adjust(Beta_Perm_Table[3:5,13], method = "fdr")

#Neaten the names
Beta_Perm_Table$Comparison <- gsub("treat","Treatment",Beta_Perm_Table$Comparison)
Beta_Perm_Table$Comparison <- gsub("depth","Depth",Beta_Perm_Table$Comparison)

#Calculate the ordinations
for (MARKER in c("16S","ITS")) {
  set.seed(123456)
  assign(paste("Ordination_Hellinger",MARKER,sep = "_"),
         ordinate(get(paste("Physeq",MARKER,sep = "_")), method="NMDS", distance=get(paste("Dist_Hellinger",MARKER,sep = "_")))
  )
}

#Plot the ordinations
Plots_NMDS_Hellinger <- NULL
for (MARKER in c("16S","ITS")) {
  current_plot <- plot_ordination(get(paste("Physeq",MARKER,sep = "_")), get(paste("Ordination_Hellinger",MARKER,sep = "_")), 
                                  color = "treat", shape = "depth") +
    geom_point(size=3) + 
    stat_ellipse(aes(linetype = depth),colour="black") +
    scale_colour_manual(values = Colours_Treatment,name="Treatment") + 
    scale_shape_manual(values = Shapes_Depth, name="Depth", labels = c("0_5" = "0-5","5_10"="5-10","15_20"="15-20")) +
    scale_linetype( labels = c("0_5" = "0-5","5_10"="5-10","15_20"="15-20")) +
    guides(linetype = guide_legend("Depth")) +
    annotate(geom="text", x = Inf, 
             y = Inf, 
             label = paste("Stress: ", round(get(paste("Ordination_Hellinger",MARKER,sep = "_"))$stress, digits = 3)), 
             vjust = "inward", hjust = "inward") +
    theme_bw() +
    labs(title = paste(MARKER,": NMDS (Hellinger)",sep = " ")) +
    theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))
  #Save plot
  Plots_NMDS_Hellinger[[MARKER]] <- current_plot
}

###Taxonomic Composition###

#Get average relative abundance over 1000 rarefactions
#Also has a problem due to the latest version of R not handling parallel computing correctly!
for (MARKER in c("16S","ITS")) {
  set.seed(123456)
  assign(paste("Physeq_Raref_RA",MARKER,sep = "_"),
         phyloseq_mult_raref_avg(get(paste("Physeq",MARKER,sep = "_")),
                                 Rare_Lib_Size,
                                 iter = 1000, parallel = 1)
  )
}

#Create plots and tables of community composition
Taxonomy_RA_plots <- list()
Taxonomy_RA_tables <- list()
Taxonomy_RA_tables_sum <- list()
for (MARKER in c("16S","ITS")) {
  for (LEVEL in c("Phylum","Class","Order","Family","Genus")) {
    cat("Analysing",LEVEL,"for",MARKER,"\n")
    
    #Agglomerate to the correct level and merge rare taxa (<1%)
    current_physeq <- tax_glom(get(paste("Physeq_Raref_RA",MARKER,sep = "_")), taxrank = LEVEL,NArm = FALSE)
    current_physeq@tax_table[taxa_names(filter_taxa(current_physeq, function(x) mean(x) < 0.01, TRUE)),] <- "Rare taxa"
    current_physeq <- merge_taxa(current_physeq,taxa_names(filter_taxa(current_physeq, function(x) mean(x) < 0.01, TRUE)))
    current_physeq@sam_data$treat <- factor(current_physeq@sam_data$treat, levels = c("C","S","W","SW"),ordered = F)
    current_physeq@sam_data$depth <- factor(current_physeq@sam_data$depth, levels = c("0_5","5_10","15_20"),ordered = T)
    #Set colours
    current_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(current_physeq)[,LEVEL])))
    names(current_colours) <- sort(unique(current_physeq@tax_table[,LEVEL]),decreasing = F)
    #Plot relative abundance
    current_plot <-
      plot_bar(current_physeq, fill = LEVEL, x = "sfence") +
      geom_bar(stat = "identity", position = "fill", colour = "black") +
      labs(x= "", y = "Relative Abundance") +
      facet_wrap(~depth*treat, ncol = 4, scales = "free_x", labeller=labeller(depth = c("0_5" = "0-5","5_10"="5-10","15_20"="15-20"))) +
      scale_fill_manual(values = current_colours, labels = names(current_colours), guide = guide_legend(ncol = 1)) +
      theme_bw() +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16), ,axis.text = element_text(size=16))
    
    #Create RA table
    current_RA_table <- as.data.frame(current_physeq@otu_table)
    current_RA_table$Taxon <- as.character(current_physeq@tax_table[row.names(current_RA_table),LEVEL])
    current_RA_table <- current_RA_table[,c(ncol(current_RA_table),1:(ncol(current_RA_table)-1))]

    #Create summary RA table for groupings
    current_RA_table_sum <- as.data.frame(matrix(nrow = nrow(current_RA_table),ncol = 21, dimnames = list(NULL,c("Taxon","Overall",unique(levels(current_physeq@sam_data$treat)),unique(levels(current_physeq@sam_data$depth)),unique(paste(current_physeq@sam_data$treat,current_physeq@sam_data$depth,sep = "-"))))))
    current_RA_table_sum$Taxon <- current_RA_table$Taxon
    current_RA_table_sum$Overall <- rowMeans(current_physeq@otu_table)
    #By treatment
    for (TREATMENT in unique(levels(current_physeq@sam_data$treat))) {
      current_RA_table_sum[,TREATMENT] <- rowMeans(current_physeq@otu_table[,current_physeq@sam_data$treat==TREATMENT])
    }
    #By depth
    for (DEPTH in unique(levels(current_physeq@sam_data$depth))) {
      current_RA_table_sum[,DEPTH] <- rowMeans(current_physeq@otu_table[,current_physeq@sam_data$depth==DEPTH])
    }
    #By combination
    for (COMBINATION in unique(paste(current_physeq@sam_data$treat,current_physeq@sam_data$depth,sep = "-"))) {
      current_RA_table_sum[,COMBINATION] <- rowMeans(current_physeq@otu_table[,current_physeq@sam_data$treat==strsplit(COMBINATION,"-")[[1]][1] &
                                                                                current_physeq@sam_data$depth==strsplit(COMBINATION,"-")[[1]][2]])
    }
    
    #Save to list
    Taxonomy_RA_plots[[paste(LEVEL,"_",MARKER,sep = "")]] <- current_plot
    Taxonomy_RA_tables[[paste(LEVEL,"_",MARKER,sep = "")]] <- current_RA_table
    Taxonomy_RA_tables_sum[[paste(LEVEL,"_",MARKER,sep = "")]] <- current_RA_table_sum
  }
}

###Differential abundance###

#Rarefied abundances were generated during the taxonomy analysis. Can skip commented code block.

# #Get average relative abundance over 1000 rarefactions
# for (MARKER in c("16S","ITS")) {
#   set.seed(123456)
#   assign(paste("Physeq_Raref_RA",MARKER,sep = "_"),
#          phyloseq_mult_raref_avg(get(paste("Physeq",MARKER,sep = "_")),
#                                  Rare_Lib_Size,
#                                  iter = 1000, parallel = 6)
#   )
# }

#Perform differential abundance analysis with Linda
Linda_Outputs <- list()
Linda_Tables <- data.frame()
for (MARKER in c("16S","ITS")) {
  for (LEVEL in c("Phylum", "Class", "Order", "Family", "Genus")) {
    cat("Analysing",LEVEL,"for",MARKER,"\n")
    
    #Agglomerate community to selected level
    current_physeq <- tax_glom(get(paste("Physeq_Raref_RA",MARKER,sep = "_")), taxrank = LEVEL, NArm = FALSE)
    
    #Set factors with a reference for Linda
    current_physeq@sam_data$treat <- factor(current_physeq@sam_data$treat, levels = unique(current_physeq@sam_data$treat), ordered = F)
    current_physeq@sam_data$treat <- relevel(current_physeq@sam_data$treat,ref = "C")
    current_physeq@sam_data$depth <- factor(current_physeq@sam_data$depth, levels = unique(current_physeq@sam_data$depth), ordered = F)
    current_physeq@sam_data$depth <- relevel(current_physeq@sam_data$depth,ref = "5_10")
    
    # set seed
    set.seed(123456)
    # run linda
    # Using prev.filter that removes anything with less than 3 non-zero values. This is because Linda warns they have no statistical power.
    linda_current <- linda(as(otu_table(current_physeq), "matrix"), as(sample_data(current_physeq), "data.frame"), 
                           formula = '~treat+depth+(1|sfence)', alpha = 0.05, prev.filter=3/length(sample_names(current_physeq)), feature.dat.type="proportion", 
                           n.cores=6)
    
    #Save result
    Linda_Outputs[[paste(LEVEL,MARKER,sep = "_")]] <- linda_current
    
    # get results from all the tables in linda_out
    for (TABLE in names(linda_current$output)) {
      
      # get table
      volcano_data <- linda_current$output[[TABLE]]
      
      #Add taxon names and other useful information
      volcano_data[,c("Phylum", "Class", "Order", "Family", "Genus")] <- current_physeq@tax_table[rownames(volcano_data),c("Phylum", "Class", "Order", "Family", "Genus")]
      volcano_data$Marker <- MARKER
      if (grepl("treat", TABLE)) { volcano_data$Comparison <- paste("C vs ",sub("treat","",TABLE),sep = "")}
      if (grepl("depth", TABLE)) { volcano_data$Comparison <- paste("5_10 vs ",sub("depth","",TABLE),sep = "")}
      volcano_data$Level <- LEVEL
      volcano_data$ASV_OTU <- row.names(volcano_data)
      #Add in RA of the base and comparison
      for (ROW in row.names(volcano_data)) {
        if (grepl("treat", TABLE)) { 
          volcano_data[ROW,"Base"] <- mean(current_physeq@otu_table[ROW,current_physeq@sam_data$treat=="C"])
          volcano_data[ROW,"Comp"] <- mean(current_physeq@otu_table[ROW,current_physeq@sam_data$treat==sub("treat","",TABLE)])
        }
        if (grepl("depth", TABLE)) { 
          volcano_data[ROW,"Base"] <- mean(current_physeq@otu_table[ROW,current_physeq@sam_data$depth=="5_10"])
          volcano_data[ROW,"Comp"] <- mean(current_physeq@otu_table[ROW,current_physeq@sam_data$depth==sub("depth","",TABLE)])
        }
      }
      
      #Save Volcano information
      Linda_Tables <- rbind(Linda_Tables,volcano_data)
     }
  }
}

#Adjust p value after multiple testing of Linda
Linda_Tables$padj <- p.adjust(Linda_Tables$pvalue, method = "fdr")
#Apply threshold
Linda_Tables$Threshold <- ifelse(Linda_Tables$log2FoldChange >= 1 & Linda_Tables$padj <= 0.05,"up", ifelse(Linda_Tables$log2FoldChange <= -1  & Linda_Tables$padj <= 0.05, "down", "other"))
#Rearrange columns
Linda_Tables <- Linda_Tables[,c("Marker","Level","Comparison","Threshold","log2FoldChange","pvalue","padj","Phylum","Class","Order","Family","Genus",
                                "Base","Comp","baseMean","lfcSE","stat","reject","df","ASV_OTU")]

#Save a table of the significantly differentially abundant taxa
Linda_Tables_Significant <- Linda_Tables[Linda_Tables$padj<=0.05 & abs(Linda_Tables$log2FoldChange) >= 1,]
Linda_Tables_Significant$Comparison <- gsub("_","-",Linda_Tables_Significant$Comparison)

#Identify the number of unique taxa which differ between at least one pair of samples.
for (MARKER in c("16S","ITS")) {
  for (LEVEL in c("Phylum","Class","Order","Family","Genus")) {
    NUMTREAT <- length(unique(Linda_Tables_Significant[Linda_Tables_Significant$Marker==MARKER & 
                                           Linda_Tables_Significant$Level==LEVEL & 
                                           Linda_Tables_Significant$Comparison%in%c("C vs S","C vs SW","C vs W"),
                                         "ASV_OTU"]))
    NUMDEPTH <- length(unique(Linda_Tables_Significant[Linda_Tables_Significant$Marker==MARKER & 
                                           Linda_Tables_Significant$Level==LEVEL & 
                                           Linda_Tables_Significant$Comparison%in%c("5-10 vs 0-5","5-10 vs 15-20"),
                                         "ASV_OTU"]))
    cat(paste("For ",MARKER," at the ",LEVEL," level ",NUMTREAT, " taxa differ between treatments and ", NUMDEPTH," differ between depths.\n",sep=""))
  }
  cat("\n")
}

#Display the Linda results as a heatmap

#Test at order level
LEVEL <- "Order"
NTAX <- 10

#How many occurrences of significantly differentially abundant taxa have 0 abundance
for (MARKER in c("16S","ITS")) {
  current_sig <- Linda_Tables_Significant[Linda_Tables_Significant$Level==LEVEL & Linda_Tables_Significant$Marker==MARKER,]
  current_physeq <- tax_glom(get(paste("Physeq_Raref_RA",MARKER,sep = "_")), taxrank = LEVEL, NArm = FALSE)
  current_average_RA <- mean(current_physeq@otu_table[current_sig$ASV_OTU,])
  current_percent_zeros <- sum(current_physeq@otu_table[current_sig$ASV_OTU,]==0)/sum(current_physeq@otu_table[current_sig$ASV_OTU,]>=0)*100
  cat("The mean RA for",MARKER,"at",LEVEL,"level is",current_average_RA,"with",current_percent_zeros,"% zero measurements.\n")
}

#Identify the relevant taxa, limited to 10 largest log2fold change per marker.
Heatmap_Taxa <- NULL
for (MARKER in c("16S","ITS")) {
  current_sig <- Linda_Tables_Significant[Linda_Tables_Significant$Level==LEVEL & Linda_Tables_Significant$Marker==MARKER,]
  current_sig <- current_sig[order(-abs(current_sig$log2FoldChange)), ]
  Heatmap_Taxa <- c(Heatmap_Taxa,unique(current_sig$ASV_OTU)[1:NTAX])
}
#Annotate the relevant taxa and create a dataframe of abundances
Heatmap_Taxa_DF <- list()
Heatmap_DF <- list()
for (MARKER in c("16S","ITS")) {
  
  #Agglomerate community to selected level
  current_physeq <- tax_glom(get(paste("Physeq_Raref_RA",MARKER,sep = "_")), taxrank = LEVEL, NArm = FALSE)
  current_heatmap_taxa <- Heatmap_Taxa[Heatmap_Taxa%in%row.names(current_physeq@otu_table)]
  
  #Make a dataframe of the ASV/OTUs and their taxonomy
  current_heatmap_taxa_df <- data.frame(Name = rep("",length(Heatmap_Taxa)), row.names = Heatmap_Taxa)
  for (ASV in current_heatmap_taxa) {
    current_heatmap_taxa_df[ASV,"Marker"] <- MARKER
    if (is.na(tax_table(current_physeq)[ASV,LEVEL][[1]])) {
      current_heatmap_taxa_df[ASV,LEVEL] <- paste("Unidentified",
                                          tax_table(current_physeq)[ASV,][[max(which(!is.na(tax_table(current_physeq)[ASV,])))]])
    } else {
      current_heatmap_taxa_df[ASV,LEVEL] <- tax_table(current_physeq)[ASV,LEVEL][[1]]
    }
  }
  current_heatmap_taxa_df$Name <- NULL
  Heatmap_Taxa_DF[[MARKER]] <- current_heatmap_taxa_df
  
  #Create summary RA table for groupings
  current_RA_table_sum <- as.data.frame(matrix(nrow = NTAX,ncol = 12, dimnames = list(current_heatmap_taxa,unique(paste(current_physeq@sam_data$treat,current_physeq@sam_data$depth,sep = "-")))))
  #By combination
  for (COMBINATION in unique(paste(current_physeq@sam_data$treat,current_physeq@sam_data$depth,sep = "-"))) {
    current_RA_table_sum[,COMBINATION] <- rowMeans(current_physeq@otu_table[,current_physeq@sam_data$treat==strsplit(COMBINATION,"-")[[1]][1] &
                                                                              current_physeq@sam_data$depth==strsplit(COMBINATION,"-")[[1]][2]])[current_heatmap_taxa]
  }
  
  #Matrix of the abundances of differentially abundant taxa
  Heatmap_DF[[MARKER]] <- current_RA_table_sum
}

#Annotate the samples
Heatmap_Sample_DF <- as.data.frame(matrix(ncol = 2,nrow=length(unique(paste(Metadata$treat,Metadata$depth,sep = "-"))),
                                          dimnames = list(unique(paste(Metadata$treat,Metadata$depth,sep = "-")),c("Treatment","Depth"))))
Heatmap_Sample_DF$Treatment <- gsub("-.+","",row.names(Heatmap_Sample_DF))
Heatmap_Sample_DF$Depth <- gsub(".+-","",row.names(Heatmap_Sample_DF))
#Neaten for plotting
Heatmap_Sample_DF$Depth <- gsub("_","-",Heatmap_Sample_DF$Depth)
#Set the factors
Heatmap_Sample_DF$Treatment <- factor(Heatmap_Sample_DF$Treatment, levels = c("C","S","W","SW"),ordered = T)
Heatmap_Sample_DF$Depth <- factor(Heatmap_Sample_DF$Depth, levels = c("0-5","5-10","15-20"),ordered = T)
Heatmap_Sample_DF <- Heatmap_Sample_DF[order(Heatmap_Sample_DF$Depth, Heatmap_Sample_DF$Treatment), ]

#Reorder the heatmap to group by columns
for (MARKER in c("16S","ITS")) {
  Heatmap_DF[[MARKER]] <- Heatmap_DF[[MARKER]][,row.names(Heatmap_Sample_DF)]
}

#Annotate which taxa are differentially abundant in which comparisons
heatmap_linda <- Linda_Tables_Significant[Linda_Tables_Significant$Level==LEVEL & Linda_Tables_Significant$ASV_OTU%in%Heatmap_Taxa,]
Heatmap_Taxa_Annotation <- as.data.frame(matrix(ncol = length(unique(heatmap_linda$Comparison)), nrow = 2*NTAX, dimnames = list(Heatmap_Taxa,rev(unique(heatmap_linda$Comparison)))))
Heatmap_Taxa_Annotation[,1:length(unique(heatmap_linda$Comparison))] <- "Unchanged"
for (ROW in 1:nrow(heatmap_linda)) {
  Heatmap_Taxa_Annotation[heatmap_linda[ROW,"ASV_OTU"],heatmap_linda[ROW,"Comparison"]] <- str_to_sentence(heatmap_linda[ROW,"Threshold"])
}

#Final prep
Colours_Marker <- c("#e9a3c9","#f7f7f7","#a1d76a")
names(Colours_Marker) <- c("Up","Unchanged","Down")
Colours_Heatmap <- Colours_Depth 
names(Colours_Heatmap) <- gsub("_","-",names(Colours_Heatmap))

#Plot the heatmap
Plot_Heatmap_Diff_Abundance <- list()
for (MARKER in c("16S","ITS")) {
  Plot_Heatmap_Diff_Abundance[[MARKER]] <- pheatmap(Heatmap_DF[[MARKER]], 
                                          display_numbers = FALSE,
                                          main = MARKER,
                                          cluster_cols = FALSE,
                                          cluster_rows = TRUE,
                                          gaps_col = c(4,8),
                                          show_colnames = FALSE,
                                          annotation_row = Heatmap_Taxa_Annotation,
                                          annotation_col = Heatmap_Sample_DF,
                                          annotation_colors = list(Treatment = Colours_Treatment, Depth = Colours_Heatmap, '5-10 vs 0-5' = Colours_Marker, '5-10 vs 15-20' = Colours_Marker),
                                          labels_row = Heatmap_Taxa_DF[[MARKER]][rownames(Heatmap_DF[[MARKER]]),LEVEL],
                                          cellwidth = 30,
                                          cellheight = 30,
                                          fontsize = 14,
                                          annotation_legend=ifelse(MARKER=="ITS",TRUE,FALSE))
}

###Fungal Guilds###

#Agglomerate to genus level because that's the level of the database
Physeq_ITS_Guilds <- tax_glom(Physeq_ITS, taxrank="Genus", NArm = FALSE)

#Insert Camelia's FungalTraits data into the phyloseq object using the most common lifestyle for a given genus.
Guilds <- as.data.frame(tax_table(Physeq_ITS_Guilds))
Guilds[,"Lifestyle"] <- "Unidentified"
#Search by genus
#For genus, the lifestyle taken is the most common lifestyle for that genus.
for (GEN in unique(Guilds$Genus)) {
  if (GEN%in%unique(CameliaTraits$V6)) {
    Guilds[which(Guilds$Genus==GEN,"Lifestyle"),"Lifestyle"] <- names(sort(table((CameliaTraits[CameliaTraits$V6 == GEN,"V8"])),decreasing = T)[1])
  }
}
sum(Guilds[,"Lifestyle"]=="Unidentified") #115 unidentified
#Search by family if missing genus
#For family, the lifestyle is taken as the most common lifestyle for that family
for (GEN in unique(Guilds[is.na(Guilds$Genus),"Family"])) {
  if (GEN%in%unique(CameliaTraits$V5)) {
    Guilds[which(Guilds$Family==GEN & is.na(Guilds$Genus)),"Lifestyle"] <- names(sort(table((CameliaTraits[CameliaTraits$V5 == GEN,"V8"])),decreasing = T)[1])
  }
}
sum(Guilds[,"Lifestyle"]=="Unidentified") #77 still unidentified
Guilds$Lifestyle <- str_to_sentence(Guilds$Lifestyle)

#Remove the extraneous data and rebuild the phyloseq object
Guilds <- Guilds[,c("Kingdom","Lifestyle")]
Physeq_ITS_Guilds <- phyloseq(otu_table(Physeq_ITS_Guilds), tax_table(as.matrix(Guilds)), sample_data(Physeq_ITS_Guilds))

#Rarefy for plotting
set.seed(123456)
Physeq_ITS_Guilds_Raref_RA <- phyloseq_mult_raref_avg(Physeq_ITS_Guilds,Rare_Lib_Size,iter = 1000)

#Plot abundance (merge <1%)
Physeq_ITS_Guilds_Raref_RA_Glom <- tax_glom(Physeq_ITS_Guilds_Raref_RA, taxrank = "Lifestyle",NArm = FALSE)
Physeq_ITS_Guilds_Raref_RA_Glom_Dominant <- merge_taxa(Physeq_ITS_Guilds_Raref_RA_Glom,taxa_names(filter_taxa(Physeq_ITS_Guilds_Raref_RA_Glom, function(x) mean(x) < 0.01, TRUE)))
Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@tax_table[is.na(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@tax_table[,"Lifestyle"]),"Lifestyle"] <- "Rare Lifestyles"
Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat <- factor(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat, levels = c("C","S","W","SW"),ordered = F)
Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth <- factor(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth, levels = c("0_5","5_10","15_20"),ordered = T)
#Set colours
colours_guilds <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant)[,"Lifestyle"])))
names(colours_guilds) <- sort(unique(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@tax_table[,"Lifestyle"]),decreasing = F)
#Plot relative abundance
Plot_Guilds <-
  plot_bar(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant, fill = "Lifestyle", x = "sfence") +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x= "", y = "Relative Abundance") +
  facet_wrap(~depth*treat, ncol = 4, scales = "free_x") +
  scale_fill_manual(values = colours_guilds, labels = names(colours_guilds), guide = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16), ,axis.text = element_text(size=16))

#Create RA table without merging the smaller guilds
Guilds_RA_table <- as.data.frame(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@otu_table)
Guilds_RA_table$Guild <- as.character(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@tax_table[row.names(Guilds_RA_table),"Lifestyle"])
Guilds_RA_table <- Guilds_RA_table[,c(ncol(Guilds_RA_table),1:(ncol(Guilds_RA_table)-1))]

#Create summary RA table for groupings
Guilds_RA_table_Sum <- as.data.frame(matrix(nrow = nrow(Guilds_RA_table),ncol = 21, dimnames = list(NULL,c("Guild","Overall",unique(levels(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat)),unique(levels(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth)),unique(paste(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat,Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth,sep = "-"))))))
Guilds_RA_table_Sum$Guild <- Guilds_RA_table$Guild
Guilds_RA_table_Sum$Overall <- rowMeans(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@otu_table)
#By treatment
for (TREATMENT in unique(levels(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat))) {
  Guilds_RA_table_Sum[,TREATMENT] <- rowMeans(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@otu_table[,Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat==TREATMENT])
}
#By depth
for (DEPTH in unique(levels(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth))) {
  Guilds_RA_table_Sum[,DEPTH] <- rowMeans(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@otu_table[,Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth==DEPTH])
}
#By combination
for (COMBINATION in unique(paste(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat,Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth,sep = "-"))) {
  Guilds_RA_table_Sum[,COMBINATION] <- rowMeans(Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@otu_table[,Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$treat==strsplit(COMBINATION,"-")[[1]][1] &
                                                                            Physeq_ITS_Guilds_Raref_RA_Glom_Dominant@sam_data$depth==strsplit(COMBINATION,"-")[[1]][2]])
}

#Perform differential abundance of guilds with Linda

#Perform analysis at Lifestyle level using the object that was already created
#Set factors with a reference for Linda
Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat <- factor(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat, levels = unique(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat), ordered = F)
Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat <- relevel(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat,ref = "C")
Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth <- factor(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth, levels = unique(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth), ordered = F)
Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth <- relevel(Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth,ref = "5_10")

# set seed
set.seed(123456)
# run linda
Linda_Output_Guild <- linda(as(otu_table(Physeq_ITS_Guilds_Raref_RA_Glom), "matrix"), as(sample_data(Physeq_ITS_Guilds_Raref_RA_Glom), "data.frame"), 
                            formula = '~treat+depth+(1|sfence)', alpha = 0.05, prev.filter=3/length(sample_names(Physeq_ITS_Guilds_Raref_RA_Glom)), feature.dat.type="proportion", 
                            n.cores=6)

# get results from all the tables in linda_out
Linda_Tables_Guild <- data.frame()
for (TABLE in names(Linda_Output_Guild$output)) {
  
  # get table
  volcano_data <- Linda_Output_Guild$output[[TABLE]]
  
  #Add taxon names and other useful information
  volcano_data[,c("Lifestyle")] <- Physeq_ITS_Guilds_Raref_RA_Glom@tax_table[rownames(volcano_data),c("Lifestyle")]
  if (grepl("treat", TABLE)) { volcano_data$Comparison <- paste("C vs ",sub("treat","",TABLE),sep = "")}
  if (grepl("depth", TABLE)) { volcano_data$Comparison <- paste("5_10 vs ",sub("depth","",TABLE),sep = "")}
  #Add in RA of the base and comparison
  for (ROW in row.names(volcano_data)) {
    if (grepl("treat", TABLE)) { 
      volcano_data[ROW,"Base"] <- mean(Physeq_ITS_Guilds_Raref_RA_Glom@otu_table[ROW,Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat=="C"])
      volcano_data[ROW,"Comp"] <- mean(Physeq_ITS_Guilds_Raref_RA_Glom@otu_table[ROW,Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$treat==sub("treat","",TABLE)])
    }
    if (grepl("depth", TABLE)) { 
      volcano_data[ROW,"Base"] <- mean(Physeq_ITS_Guilds_Raref_RA_Glom@otu_table[ROW,Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth=="5_10"])
      volcano_data[ROW,"Comp"] <- mean(Physeq_ITS_Guilds_Raref_RA_Glom@otu_table[ROW,Physeq_ITS_Guilds_Raref_RA_Glom@sam_data$depth==sub("depth","",TABLE)])
    }
  }
  
  #Save Volcano information
  Linda_Tables_Guild <- rbind(Linda_Tables_Guild,volcano_data)
}
    
#Adjust p value after multiple testing of Linda
Linda_Tables_Guild$padj <- p.adjust(Linda_Tables_Guild$pvalue, method = "fdr")
#Apply threshold
Linda_Tables_Guild$Threshold <- ifelse(Linda_Tables_Guild$log2FoldChange >= 1 & Linda_Tables_Guild$padj <= 0.05,"up", ifelse(Linda_Tables_Guild$log2FoldChange <= -1  & Linda_Tables_Guild$padj <= 0.05, "down", "other"))
#Rearrange columns
Linda_Tables_Guild <- Linda_Tables_Guild[,c("Comparison","Threshold","log2FoldChange","pvalue","padj","Lifestyle",
                                "Base","Comp","baseMean","lfcSE","stat","reject","df")]

#Save a table of the significantly differentially abundant taxa
Linda_Tables_Guild_Significant <- Linda_Tables_Guild[Linda_Tables_Guild$padj<=0.05 & abs(Linda_Tables_Guild$log2FoldChange) >= 1,]
Linda_Tables_Guild_Significant$Comparison <- gsub("_","-",Linda_Tables_Guild_Significant$Comparison)

###Analyse the soil chem data###

#Create a table to give human readable titles and units for the graphs
Soil_Chem_Key <- as.data.frame(matrix(nrow = 21,ncol = 4, dimnames = list(NULL,c("Variable","Name","Unit","Category"))))
Soil_Chem_Key[,"Variable"] <- colnames(Metadata)[7:27]
Soil_Chem_Key[,"Name"] <- c("pH","Carbon","Nitrogon","Phosphorous","Ergosterol","Laccase","Cellulase","Xylanase","Oxidase",
                            "Peroxidase","Mn perox.","B glucos.","Acid phos.","B Xylos.","Chitinase","Cellobiohyd.",
                            "B galactos.","A Glucos.","Lipase","16S qPCR","18S qPCR")
Soil_Chem_Key[,"Unit"] <- c("","%","%","","ug/g","mU/g","mU/g","mU/g","mU/g",
                            "mU/g","mU/g","U/g","U/g","U/g","U/g","U/g",
                            "U/g","U/g","U/g","10^10 copies/g","10^8copies/g")
Soil_Chem_Key[,"Category"] <- c("Soil","Soil","Soil","Soil","Microbial","Enzymatic","Enzymatic","Enzymatic","Enzymatic",
                                "Enzymatic","Enzymatic","Enzymatic","Enzymatic","Enzymatic","Enzymatic","Enzymatic",
                                "Enzymatic","Enzymatic","Enzymatic","Microbial","Microbial")
Soil_Chem_Key[,"Full"] <- c("pH","Percentage C","Percentage N","Total P","Ergosterol","Laccase","Cellulase","Xylanase","Oxidase",
                            "Peroxidase","Manganese peroxidase","Beta Glucosidase","Acid Phosphatase","Beta Xylosidase","Chitinase","Cellobiohydrolase",
                            "Beta Galactosidase","Alpha Glucosidase","Lipase","16S qPCR","18S qPCR")
Soil_Chem_Key[,"Abbreviation"] <- c("pH","C","N","P","Erg.","Lac.","Cell.","Xyl.","Oxi.",
                            "Per.","Mn per.","B glu.","Acid phos.","B Xyl.","Chi.","Cello.",
                            "B gal.","A Glu.","Lip.","16S","18S")

#Check the soil chemistry data for normality
Plots_Soil_Chem_Normality <- list()
for (VAR in Soil_Chem_Key$Variable) {
  
  #Get the measurements
  current_values <- Metadata[!is.na(Metadata[,VAR]),VAR]
  current_n <- length(current_values)
  
  #Perform statistical tests
  current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
  current_sw <- shapiro.test(current_values)$`p.value`
  
  #Create a normal distribution over the same distance
  norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                           Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                     mean = mean(current_values),
                                     sd = sd(current_values)))
  
  #Adjust to the height of the graph
  scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
  norm_table$Y <- norm_table$Y*scaling_factor
  
  #Plot histogram of the data
  hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                           range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                           length.out = 3))
  current_plot_normality <- ggplot(Metadata[!is.na(Metadata[,VAR]), c("treat", "depth", "sfence", VAR)],aes(!! sym(VAR))) +
    geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
    geom_point(data = norm_table, aes(x = X, y = Y)) +
    annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
    annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
    annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
    labs(title = paste(Soil_Chem_Key[Soil_Chem_Key[,"Variable"]==VAR,"Name"])) +
    scale_y_continuous(breaks = hist_counts) +
    ylab(NULL) +
    xlab(label = paste(Soil_Chem_Key[Soil_Chem_Key[,"Variable"]==VAR,"Unit"])) +
    theme_bw() +
    theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16),
          panel.background = element_rect(fill = alpha(Colours_Categories[Soil_Chem_Key[Soil_Chem_Key$Variable==VAR,"Category"]],alpha = 0.4)))
  
  #Save the plot
  Plots_Soil_Chem_Normality[[VAR]] <- current_plot_normality
}

#Graph the soil properties
Plots_Soil_Chem <- list()
Soil_Chem_Anova <- list()
for (VAR in Soil_Chem_Key$Variable) {
  
  #Perform anova
  set.seed(123456)
  current_aov <- aov(get(VAR) ~ sfence+treat*depth, data=Metadata)
  current_aov_summary <- summary(current_aov)
  
  #Tukey test and prepare CLD
  current_aov_tukey <- TukeyHSD(current_aov,ordered = TRUE)
  current_tukey_table <- as.data.frame(current_aov_tukey$`treat:depth`)
  current_res <- current_tukey_table$`p adj`
  names(current_res) <- row.names(current_tukey_table)
  current_letters <- multcompLetters(current_res)
  current_letters_df <- data.frame("comp" = names(current_letters$Letters), "Letter"=as.vector(current_letters$Letters))
  current_letters_df$treat <- gsub(":.+","",current_letters_df$comp)
  current_letters_df$depth <- gsub("^.+:{1}","",current_letters_df$comp)
  current_letters_df$treat <- factor(current_letters_df$treat, levels = c("C","S","W","SW"),ordered = T)
  current_letters_df$depth <- factor(current_letters_df$depth, levels = c("0_5","5_10","15_20"),ordered = T)
  current_letters_df[,VAR] <- max(Metadata[,VAR],na.rm = T)+(0.05*(max(abs(Metadata[,VAR]),na.rm = T)-min(abs(Metadata[,VAR]),na.rm = T)))
  
  #Perform the plotting
  Metadata_Plot <- Metadata
  if (VAR=="qpcr_bacteria") {
    YLABEL <- bquote(~"10"^{10}~" copies/g")
    Metadata_Plot$qpcr_bacteria <- round(Metadata_Plot$qpcr_bacteria/10000000000,digits = 1)
    current_letters_df$qpcr_bacteria <- round(current_letters_df$qpcr_bacteria/10000000000,digits = 1)
  } else if (VAR=="qpcr_fungi") {
    YLABEL <- bquote(~"10"^{8}~" copies/g")
    Metadata_Plot$qpcr_fungi <- round(Metadata_Plot$qpcr_fungi/100000000,digits = 1)
    current_letters_df$qpcr_fungi <- round(current_letters_df$qpcr_fungi/100000000,digits = 1)
  } else {
    YLABEL <- paste(Soil_Chem_Key[Soil_Chem_Key[,"Variable"]==VAR,"Unit"])
  }
  current_plot <- ggplot(Metadata_Plot, aes(x = depth, y = !! sym(VAR))) +
    geom_boxplot(aes(fill=NULL)) +
    stat_summary(geom = "point", fun = mean, shape = 18, size = 3) +
    geom_point(aes(colour=sfence),position = position_dodge(width = 0.2)) +
    geom_text(data = current_letters_df, aes(label=Letter)) +
    labs(title = paste(Soil_Chem_Key[Soil_Chem_Key[,"Variable"]==VAR,"Full"])) +
    ylab(label = YLABEL) +
    facet_wrap(~treat,ncol=4) +
    scale_color_manual(values = Colours_Fence) +
    guides(colour = guide_legend(title = "Plot")) +
    xlab(label = NULL) +
    theme_bw() +
    theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.spacing=unit(0,"lines"),axis.text = element_text(size=16),
          panel.background = element_rect(fill = alpha(Colours_Categories[Soil_Chem_Key[Soil_Chem_Key$Variable==VAR,"Category"]],alpha = 0.2)))
  
  #Save
  Soil_Chem_Anova[[VAR]] <- current_aov_summary
  Plots_Soil_Chem[[VAR]] <- current_plot
  
}

#Perform linear modelling
Soil_Chem_Linear <- list()
for (VAR in Soil_Chem_Key$Variable) {
  
  #Get the data
  current_data <- Metadata[, c("treat", "depth", "sfence", VAR)]
  colnames(current_data) <- c("treat", "depth", "sfence", "an_env")
  
  #Perform linear modelling
  set.seed(123456)
  an_env_model <- lme(an_env~treat*depth, random=~1|sfence, data=current_data, na.action=na.omit)
  current_summary <- summary(an_env_model)
  current_res <- as.data.frame(current_summary$tTable)
  
  #Save
  Soil_Chem_Linear[[VAR]] <- current_res
  
}

#Make a table of VARs with significant differences according to linear models
Metadata_table <- data.frame()
for (VAR in names(Soil_Chem_Linear)) {
  
  #Create table
  current <- Soil_Chem_Linear[VAR][[1]]
  current_df <- as.data.frame(matrix(nrow = 12,ncol = 7, dimnames = list(NULL,c("Variable","Test","Value","Std.Error","DF","t value","p value"))))
  current_df$Variable <- VAR
  current_df$Test <- str_trim(row.names(current))
  current_df$Value <- current$Value
  current_df$Std.Error <- current$Std.Error
  current_df$DF <- current$DF
  current_df$`t value` <- current$`t-value`
  current_df$`p value` <- current$`p-value`
  
  #Save
  Metadata_table <- rbind(Metadata_table,current_df)
}
#Correct p values for multiple testing
Metadata_table$padj <- p.adjust(Metadata_table$`p value`, method = "fdr")

#Identify the properties which are significant in the linear model
Metadata_table_signif <- Metadata_table
Metadata_table_signif <- Metadata_table_signif[!Metadata_table_signif$Test=="(Intercept)",]
Metadata_table_signif <- Metadata_table_signif[Metadata_table_signif$padj < 0.05,]
#Remove any which are not in the linear relationship
Metadata_table_signif <- Metadata_table_signif[!grepl(".Q",Metadata_table_signif$Test),]
Metadata_table_signif$Test <- gsub(".L","",Metadata_table_signif$Test)
#For interpretation: https://stackoverflow.com/questions/57297771/interpretation-of-l-q-c-4-for-logistic-regression

#Clean up the names
for (VAR in Metadata_table$Variable) {
  Metadata_table[which(Metadata_table$Variable==VAR),"Variable"] <- Soil_Chem_Key[Soil_Chem_Key$Variable==VAR,"Name"]
}
for (VAR in Metadata_table_signif$Variable) {
  Metadata_table_signif[which(Metadata_table_signif$Variable==VAR),"Variable"] <- Soil_Chem_Key[Soil_Chem_Key$Variable==VAR,"Name"]
}
Metadata_table_signif$Test <- str_to_title(Metadata_table_signif$Test)

###Redundancy Analysis (RDA)###

#Rarefied abundances were generated during the taxonomy analysis. Can skip commented code block.

# #Get average relative abundance over 1000 rarefactions
# for (MARKER in c("16S","ITS")) {
#   set.seed(123456)
#   assign(paste("Physeq_Raref_RA",MARKER,sep = "_"),
#          phyloseq_mult_raref_avg(get(paste("Physeq",MARKER,sep = "_")),
#                                  Rare_Lib_Size,
#                                  iter = 1000, parallel = 6)
#   )
# }

#Define lists to store the results from the RDAs
RDA_R_Square_Partial <- list()
RDA_Signif_Test_Model <- list()
RDA_Signif_Test_Terms <- list()
RDA_Figures <- list()
RDA_Arrows <- data.frame()

#Perform RDA at ASV/OTU level
for (MARKER in c("16S","ITS")) {
  
  #Get current phyloseq object
  current_physeq <- get(paste("Physeq_Raref_RA",MARKER,sep = "_"))
  
  #Perform the RDA with the set factors
  set.seed(123456)
  current_RDA <- ordinate(current_physeq, "RDA", formula=as.formula(".~treat+depth+Condition(sfence)"), scale=1)
  
  #Get adjusted R2
  current_r_square_partial <- RsquareAdj(current_RDA)$adj.r.squared
  
  #Test model significance
  set.seed(123456)
  current_RDA_Signif_Test_Model <- anova.cca(current_RDA, step=1000)
  
  #Get significance of terms
  set.seed(123456)
  current_RDA_Signif_Test_Terms <- anova.cca(current_RDA, step=1000, by="terms")
  
  #Create data frame for ggplot plotting
  #Get scores for sites
  current_RDA_df <- as.data.frame(scores(current_RDA, display="sites", scaling = 1)) #Altered to use the same scores settings as plot_ordination
  current_RDA_df$label <- rownames(current_RDA_df)
  current_RDA_df$Treatment <- current_physeq@sam_data[current_RDA_df$label,"treat"]$treat
  current_RDA_df$Depth <- current_physeq@sam_data[current_RDA_df$label,"depth"]$depth
  
  # set factors
  current_RDA_df$Treatment <- factor(current_RDA_df$Treatment, levels = c("C","S","W","SW"),ordered = F)
  current_RDA_df$Depth <- factor(current_RDA_df$Depth, levels= c("0_5","5_10","15_20"),ordered = T)
  
  #Variance explained by axes
  current_RDA_var <- as.data.frame(summary(current_RDA)$cont$importance)
  
  # compute fitting
  set.seed(123456)
  current_var_fitting_enviro <- envfit(current_RDA, as(sample_data(current_physeq), "data.frame")[, colnames(current_physeq@sam_data)[7:27]], perm=1000, na.rm=TRUE, scaling=1)
  
  #Get multiplying factors for the arrows
  #Needs a plot to get the right coordinates
  plot(current_RDA, type="none", scaling=1, main="", display ="sites")
  arrow_factor <- ordiArrowMul(current_var_fitting_enviro, choices = c(1,2))
  # get arrows and scale them
  arrows <- as.data.frame(scores(current_var_fitting_enviro, display="vectors")) * arrow_factor #Scaling means nothing here
  arrows$envvars <- rownames(arrows)
  #Make clean names
  for (ROW in row.names(arrows)) {
    arrows[ROW,"envvars"] <- Soil_Chem_Key[Soil_Chem_Key$Variable==ROW,"Abbreviation"]
  }
  arrows$r2 <- current_var_fitting_enviro$vectors$r
  arrows$Marker <- MARKER
  arrows$p <- current_var_fitting_enviro$vectors$pvals
  # select arrows which have p-value <= 0.05
  sig_arrows <- arrows[names(which(current_var_fitting_enviro$vectors$pvals<=0.05)), ]
    
  #Plot RDA
  current_RDA_plot <- ggplot(current_RDA_df) + #Now matches
    geom_point(mapping=aes(x=RDA1, y=RDA2, colour=Treatment, shape=Depth), size = 4) +
    coord_fixed() +
    scale_shape_manual(values=Shapes_Depth, labels = c("0_5" = "0-5","5_10"="5-10","15_20"="15-20")) +
    scale_colour_manual(values=Colours_Treatment) +
    geom_segment(data=sig_arrows, aes(x=0, xend=RDA1, y=0, yend=RDA2), arrow=arrow(length=unit(0.25, "cm")), colour="black", lwd=0.5) +
    geom_text_repel(data=sig_arrows, aes(x=RDA1, y=RDA2, label=envvars), col="black", alpha = 70,max.overlaps = 20)+
    guides(fill=guide_legend(override.aes=list(shape=22, size=4)), shape=guide_legend(override.aes=list(size=4))) +
    labs(title = paste(MARKER,"pRDA",sep = " "),
         x = paste("RDA1 (",round(current_RDA_var[2,"RDA1"]*100,digits = 2),"%)",sep = ""), 
         y = paste("RDA2 (",round(current_RDA_var[2,"RDA2"]*100,digits = 2),"%)",sep = "")) +
    theme_bw() +
    theme(text=element_text(size=16), aspect.ratio = 1, axis.text = element_text(size=16), plot.title = element_text(hjust = 0.5, size=16))
  
  #Save results to lists
  RDA_R_Square_Partial[[MARKER]] <- current_r_square_partial
  RDA_Signif_Test_Model[[MARKER]] <- current_RDA_Signif_Test_Model
  RDA_Signif_Test_Terms[[MARKER]] <- current_RDA_Signif_Test_Terms
  RDA_Figures[[MARKER]] <- current_RDA_plot
  RDA_Arrows <- rbind(RDA_Arrows, arrows)
  
}

#Find the variables which are brought up in the RDA
#Make a table of the RDA r2 values
RDA_Arrows_signif <- RDA_Arrows[RDA_Arrows$p < 0.05,]
RDA_Arrows_signif <- RDA_Arrows_signif[,c("Marker","r2","envvars","p")]

###Analyse the temperature and Moisture###

#Split temperature and Moisture
NARS_Moisture <- subset(NARS_Temp_Hum,NARS_Temp_Hum$SensorType=="Moist")
NARS_Temperature <- subset(NARS_Temp_Hum,NARS_Temp_Hum$SensorType=="Temp")

#Correct date format
NARS_Temperature$DateTime <- as.POSIXct(paste(NARS_Temperature$Date, NARS_Temperature$Time, sep=" "), format="%Y-%m-%d %H:%M:%S", tz = "UTC")
NARS_Moisture$DateTime <- as.POSIXct(paste(NARS_Moisture$Date, NARS_Moisture$Time, sep=" "), format="%Y-%m-%d %H:%M:%S", tz = "UTC")

#Limit data to the sampling period
NARS_Temperature <- NARS_Temperature[NARS_Temperature$DateTime >= as.POSIXct("2015-07-01", format="%Y-%m-%d", tz = "UTC") & NARS_Temperature$DateTime <= as.POSIXct("2020-08-31", format="%Y-%m-%d", tz = "UTC"),]
NARS_Moisture <- NARS_Moisture[NARS_Moisture$DateTime >= as.POSIXct("2015-07-01", format="%Y-%m-%d", tz = "UTC") & NARS_Moisture$DateTime <= as.POSIXct("2020-08-31", format="%Y-%m-%d", tz = "UTC"),]

#Add factors
NARS_Temperature$Treatment <- factor(NARS_Temperature$Treatment, levels = c("C","S","W","SW"),ordered = T)
NARS_Moisture$Treatment <- factor(NARS_Moisture$Treatment, levels = c("C","S","W","SW"),ordered = T)

#Check control average and variation
summary(NARS_Temperature[NARS_Temperature$Treatment=="C","Value"])
summary(NARS_Moisture[NARS_Moisture$Treatment=="C","Value"])

#Calculate the differences
for (SENSOR in c("Temperature","Moisture")) {
  current_dataframe <- get(paste("NARS",SENSOR,sep = "_"))
  #Calculate the differences for each treatment
  for (TREATMENT in c("C","S","W","SW")) {
    if (identical(current_dataframe[current_dataframe$Treatment==TREATMENT,"DateTime"],
                  current_dataframe[current_dataframe$Treatment=="C","DateTime"])) {
      cat("Match!\n")
      current_dataframe[current_dataframe$Treatment==TREATMENT,"Difference"] <- 
        current_dataframe[current_dataframe$Treatment==TREATMENT,"Value"] - current_dataframe[current_dataframe$Treatment=="C","Value"]
    } else ("No match :( \n")
  }
  #Save
  assign(paste("NARS",SENSOR,sep = "_"),
         current_dataframe
  )
}

#Calculate the average differences over the whole dataset
NARS_Temperature_Diffs_Means <- aggregate(Difference~Treatment,NARS_Temperature[c("Difference","Treatment")],mean)
NARS_Moisture_Diffs_Means <- aggregate(Difference~Treatment,NARS_Moisture[c("Difference","Treatment")],mean)


#Create a short version of the dataset to calculate confidence intervals
#This way we can more easily guarantee the values are paired than in long version
for (SENSOR in c("Temperature","Moisture")) {
  #Set up the table
  current_short <- get(paste("NARS",SENSOR,sep = "_"))[get(paste("NARS",SENSOR,sep = "_"))$Treatment=="C",c("DateTime","Value")]
  colnames(current_short)[2] <- "C"
  #Loop through treatments
  for (TREATMENT in c("S","W","SW")) {
    #Make sure the measurements are in exactly the same order
    if (identical(get(paste("NARS",SENSOR,sep = "_"))[get(paste("NARS",SENSOR,sep = "_"))$Treatment==TREATMENT,"DateTime"],
                  current_short$DateTime)) {
      cat("Match!\n")
      current_short[,TREATMENT] <- get(paste("NARS",SENSOR,sep = "_"))[get(paste("NARS",SENSOR,sep = "_"))$Treatment==TREATMENT,"Value"]
    } else ( stop("No match :( \n") )
  }
  #Save
  assign(paste("NARS_Short",SENSOR,sep = "_"),
         current_short
  )
}

#Calculate the group means and differences from control for the whole dataset
for (SENSOR in c("Temperature","Moisture")) {
  #Set up the table
  current_short <- get(paste("NARS_Short",SENSOR,sep = "_"))
  current_means_diff <- as.data.frame(matrix(nrow = 3,ncol = 4, dimnames = list(NULL,c("Treatment","Difference","Lower","Upper"))))
  current_means_diff$Treatment <- c("S","W","SW")
  #Calculate the difference in means and CI for each row
  for (ROW in 1:nrow(current_means_diff)){
    
    current_means_diff[ROW,c("Difference","Lower","Upper")] <- MeanDiffCI(current_short[!is.na(current_short$C)&!is.na(current_short[,current_means_diff[ROW,"Treatment"]]),current_means_diff[ROW,"Treatment"]],
                                                                          current_short[!is.na(current_short$C)&!is.na(current_short[,current_means_diff[ROW,"Treatment"]]),"C"],
                                                                          method = "basic", na.rm = T, paired = T)
  }
  #Set factor
  current_means_diff$Treatment <- factor(current_means_diff$Treatment, levels = c("S","W","SW"),ordered = T)
  #Save
  assign(paste("Means_Diff_Overall_CI",SENSOR,sep = "_"),
         current_means_diff
  )
}

#Calculate the group means and differences from control for each month
for (SENSOR in c("Temperature","Moisture")) {
  #Set up the table
  current_short <- get(paste("NARS_Short",SENSOR,sep = "_"))
  current_means_diff <- as.data.frame(matrix(nrow = 36,ncol = 5, dimnames = list(NULL,c("Treatment","Month","Difference","Lower","Upper"))))
  current_means_diff$Treatment <- c(rep("S",12),rep("W",12),rep("SW",12))
  current_means_diff$Month <- rep(1:12,3)
  #Calculate the difference in means and CI for each row
  for (ROW in 1:nrow(current_means_diff)){
    
    current_means_diff[ROW,c("Difference","Lower","Upper")] <- MeanDiffCI(current_short[!is.na(current_short$C)&!is.na(current_short[,current_means_diff[ROW,"Treatment"]])&month(current_short$DateTime)==current_means_diff[ROW,"Month"],current_means_diff[ROW,"Treatment"]],
                                                                          current_short[!is.na(current_short$C)&!is.na(current_short[,current_means_diff[ROW,"Treatment"]])&month(current_short$DateTime)==current_means_diff[ROW,"Month"],"C"],
                                                                               method = "basic", na.rm = T, paired = T)
  }
  #Rename months
  current_means_diff$Month <- month.abb[current_means_diff$Month]
  #Set factor
  current_means_diff$Treatment <- factor(current_means_diff$Treatment, levels = c("S","W","SW"),ordered = T)
  current_means_diff$Month <- factor(current_means_diff$Month, levels = month.abb[1:12], ordered = TRUE)
  #Save
  assign(paste("Means_Diff_CI",SENSOR,sep = "_"),
         current_means_diff
  )
}

#Plot the treatment effects
for (SENSOR in c("Temperature","Moisture")) {
  #Plot the data
  if (SENSOR=="Moisture") {
    YLABEL <- bquote("Difference in means ("~"cm"^{3}~".cm"^{-3}~")")
  } else {
    YLABEL <- "Difference in means (°C)"
  }
  assign(paste("Plot_Mean_Diff_CI",SENSOR,sep = "_"),
         ggplot(get(paste("Means_Diff_CI",SENSOR,sep = "_")), aes(x=Treatment, y=Difference)) +
           geom_hline(yintercept = 0,colour = "red") +
           geom_errorbar(aes(ymin = Lower, ymax = Upper, colour = Treatment)) +
           geom_point(aes(colour = Treatment), size = 3) +
           scale_colour_manual(values = Colours_Treatment) +
           theme_bw() +
           facet_wrap(~Month, ncol = 6) +
           scale_y_continuous(limits = c(RoundTo(min(get(paste("Means_Diff_CI",SENSOR,sep = "_"))$Lower)*1.05,multiple = 0.1,floor),
                                         RoundTo(max(get(paste("Means_Diff_CI",SENSOR,sep = "_"))$Upper)*1.05,multiple = 0.1,ceiling)),n.breaks = 6) +
           labs(title = paste("Difference in",SENSOR,"means per month",sep = " "),
                x = "", 
                y = YLABEL) +
           theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 16),
                 panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  )
}

#Plot the actual measurements
Plot_Month_Names <- month.abb
names(Plot_Month_Names) <- 1:12
for (SENSOR in c("Temperature","Moisture")) {
  #Plot the data
  if (SENSOR=="Moisture") {
    YLABEL <- bquote("Soil Moisture ("~"cm"^{3}~".cm"^{-3}~")")
  } else {
    YLABEL <- "Soil Temperature (°C)"
  }
  assign(paste("Plot_Violin",SENSOR,sep = "_"),
         ggplot(get(paste("NARS",SENSOR,sep = "_")), aes(x=Treatment, y=Value)) +
           geom_violin(aes(fill = Treatment)) +
           geom_boxplot(alpha=0) +
           stat_summary(geom = "point", fun = mean, shape = 18, size = 3) +
           scale_fill_manual(values = Colours_Treatment) +
           theme_bw() +
           facet_wrap(~Month, ncol = 6, labeller = as_labeller(Plot_Month_Names)) +
           scale_y_continuous(n.breaks = 6) +
           labs(title = paste("Measurement of soil",str_to_lower(SENSOR),"per month",sep = " "),
                x = "", 
                y = YLABEL) +
           theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 16),
                 panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  )
}


###Publication Exports###

#Main figures

#Alpha Beta Diversity
Alpha_legend <- get_legend(Plots_Alpha$Observed_16S + theme(legend.position = "right")+ guides(colour = guide_legend(ncol = 1)))
Beta_legend <- get_legend(Plots_NMDS_Hellinger$`16S` + theme(legend.position = "right") + guides(colour = guide_legend(ncol = 1)))
Figure_1 <- plot_grid(
    plot_grid(Plots_Alpha$Observed_16S + theme(legend.position = "none") + labs(tag="A"),
              Plots_Alpha$Observed_ITS + theme(legend.position = "none") + labs(tag="B"),
              Alpha_legend,
              ncol = 3, rel_widths = c(1,1,0.3)),
    plot_grid(Plots_NMDS_Hellinger$`16S` + theme(legend.position = "none") + labs(tag="C"),
              Plots_NMDS_Hellinger$ITS + theme(legend.position = "none") + labs(tag="D"),
              Beta_legend,
              ncol = 3, rel_widths = c(1,1,0.3)),
    ncol = 1)
#SVG
svg("Figure_1.svg", width=2*168*0.03937008, height=2*100*0.03937008)
Figure_1
dev.off()
#JPG
jpeg("Figure_1.jpg", width=2*168, height=2*100, units = "mm", res = 300)
Figure_1
dev.off()

#Differential Abundance
plot_heatmap_16S <- grid.grabExpr(print(Plot_Heatmap_Diff_Abundance$`16S`))
plot_heatmap_ITS <- grid.grabExpr(print(Plot_Heatmap_Diff_Abundance$ITS))
Figure_2 <- plot_grid(plot_heatmap_16S, 
                      plot_heatmap_ITS)
#SVG
svg("Figure_2.svg", width=3*168*0.03937008, height=2*100*0.03937008)
Figure_2
dev.off()
#JPG
jpeg("Figure_2.jpg", width=3*168, height=2*100, units = "mm", res = 300)
Figure_2
dev.off()

#pRDA
RDA_legend <- get_legend(RDA_Figures$`16S` + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
Figure_3 <- plot_grid(
  plot_grid(RDA_Figures$`16S` + theme(legend.position = "none") + labs(tag = "A"),
            RDA_Figures$ITS + theme(legend.position = "none") + labs(tag = "B"),
            ncol = 2, greedy = FALSE),
  RDA_legend, ncol = 1, rel_heights = c(10,1),greedy = FALSE)
#SVG
svg("Figure_3.svg", width=2*114*0.03937008, height=2*65*0.03937008)
Figure_3
dev.off()
#JPG
jpeg("Figure_3.jpg", width=2*114, height=2*65, units = "mm", res = 300)
Figure_3
dev.off()

#Supplementary figures

Differences_Legend <- get_legend(Plot_Violin_Temperature + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
Figure_S1 <- plot_grid(
  plot_grid(Plot_Mean_Diff_CI_Temperature + theme(legend.position = "none") + labs(tag="A"),
            Plot_Violin_Temperature + theme(legend.position = "none") + labs(tag="B"),
            Plot_Mean_Diff_CI_Moisture + theme(legend.position = "none") + labs(tag="C"),
            Plot_Violin_Moisture + theme(legend.position = "none") + labs(tag="D"),
            ncol = 2),
  Differences_Legend, ncol = 1, rel_heights = c(10,1))
#SVG
svg("Figure_S1.svg", width=2*228*0.03937008, height=2*114*0.03937008)
Figure_S1
dev.off()
#JPG
jpeg("Figure_S1.jpg", width=2*228, height=2*114, units = "mm", res = 300)
Figure_S1
dev.off()

Environment_Legend <- get_legend(Plots_Soil_Chem$ph + theme(legend.position = "right") + guides(colour = guide_legend(nrow = 1,title = "Plot")))
Figure_S2 <- 
  plot_grid(Plots_Soil_Chem$ph + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_ergo + theme(legend.position = "none"),
            Plots_Soil_Chem$total_per_c + theme(legend.position = "none"),
            Plots_Soil_Chem$qpcr_bacteria + theme(legend.position = "none"),
            Plots_Soil_Chem$total_per_n + theme(legend.position = "none"),
            Plots_Soil_Chem$qpcr_fungi + theme(legend.position = "none"),
            Plots_Soil_Chem$total_p + theme(legend.position = "none"),
            ncol=2, align = "hv")
Figure_S2 <- Figure_S2 + draw_grob(Environment_Legend, x = 0.25, y = -0.375)

#SVG
svg("Figure_S2.svg", width=2*210*0.03937008, height=2*150*0.03937008)
Figure_S2
dev.off()
#JPG
jpeg("Figure_S2.jpg", width=2*210, height=2*150, units = "mm", res = 300)
Figure_S2
dev.off()

Environment_Legend_2 <- get_legend(Plots_Soil_Chem$enz_lac + theme(legend.position = "right") + guides(colour = guide_legend(nrow = 1,title = "Plot")))
Figure_S3 <-
  plot_grid(Plots_Soil_Chem$enz_lac + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_cell + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_xylan + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_oxi + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_peroxi + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_mnp + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_betagluc + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_phos + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_betaxylo + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_chitin + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_hydro + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_betagalac + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_alphagluc + theme(legend.position = "none"),
            Plots_Soil_Chem$enz_lip + theme(legend.position = "none"),
            ncol=3, align = "hv")
Figure_S3 <- Figure_S3 + draw_grob(Environment_Legend_2, x = 0.33, y = -0.4)

#SVG
svg("Figure_S3.svg", width=2*300*0.03937008, height=2*150*0.03937008)
Figure_S3
dev.off()
#JPG
jpeg("Figure_S3.jpg", width=2*300, height=2*150, units = "mm", res = 300)
Figure_S3
dev.off()

Figure_S4 <- plot_grid(
  plot_grid(Plots_Alpha$Shannon_16S + theme(legend.position = "none") + labs(tag="A"),
            Plots_Alpha$Shannon_ITS + theme(legend.position = "none") + labs(tag="B"),
            ncol = 2),
  Alpha_legend, ncol = 2, rel_widths = c(1,0.1))
#SVG
svg("Figure_S4.svg", width=2*168*0.03937008, height=2*50*0.03937008)
Figure_S4
dev.off()
#JPG
jpeg("Figure_S4.jpg", width=2*168, height=2*50, units = "mm", res = 300)
Figure_S4
dev.off()

Figure_S5 <- plot_grid(Taxonomy_RA_plots$Phylum_16S + labs(tag="A") + theme(legend.justification = "left"),
                                        Taxonomy_RA_plots$Order_16S + labs(tag="B") + theme(legend.justification = "left"),
                                        Taxonomy_RA_plots$Genus_16S + labs(tag="C") + theme(legend.justification = "left"),
                                        ncol = 1, align = "hv")
#SVG
svg("Figure_S5.svg", width=2*168*0.03937008, height=2*342*0.03937008)
Figure_S5
dev.off()
#JPG
jpeg("Figure_S5.jpg", width=2*168, height=2*342, units = "mm", res = 300)
Figure_S5
dev.off()

Figure_S6 <- plot_grid(Taxonomy_RA_plots$Phylum_ITS + labs(tag="A") + theme(legend.justification = "left"),
                                        Taxonomy_RA_plots$Order_ITS + labs(tag="B") + theme(legend.justification = "left"),
                                        Taxonomy_RA_plots$Genus_ITS + labs(tag="C") + theme(legend.justification = "left"),
                                        ncol = 1, align = "hv")
#SVG
svg("Figure_S6.svg", width=2*168*0.03937008, height=2*342*0.03937008)
Figure_S6
dev.off()
#JPG
jpeg("Figure_S6.jpg", width=2*168, height=2*342, units = "mm", res = 300)
Figure_S6
dev.off()

#Main tables
write.table(Metadata_table_signif[,c("Variable","Test","Value","Std.Error","padj")], file = "Table_1.csv",quote = T,row.names = F, col.names = T, sep = ",")

#Supplementary tables
write.table(Beta_Perm_Table, file = "Table_S1.csv",quote = T,row.names = F, col.names = F, sep = ",")
write.table(Linda_Tables_Significant,"Table_S2.csv",quote = T,row.names = F, sep = ",")
write.table(Linda_Tables_Guild_Significant,"Table_S3.csv",quote = T,row.names = F, sep = ",")

# sessionInfo()
# R version 4.6.1 (2026-06-24)
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
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] metagMisc_0.6.0.9000 ggpubfigs_1.1.3      mirlyn_1.4.2         phyloseq_1.56.0      emmeans_2.0.4        pheatmap_1.0.13     
# [7] lubridate_1.9.5      stringr_1.6.0        multcompView_0.1-12  nlme_3.1-171         MicrobiomeStat_1.4   RColorBrewer_1.1-3  
# [13] DescTools_0.99.60    vegan_2.7-6          permute_0.9-10       cowplot_1.2.0        ggrepel_0.9.8        ggplot2_4.0.3       
# 
# loaded via a namespace (and not attached):
# [1] jsonlite_2.0.0       rstudioapi_0.19.0    magrittr_2.0.5       TH.data_1.1-5        estimability_2.0.0   farver_2.1.2        
# [7] nloptr_2.2.1         fs_2.1.0             vctrs_0.7.3          multtest_2.68.0      minqa_1.2.8          forcats_1.0.1       
# [13] haven_2.5.5          cellranger_1.1.0     parallelly_1.48.0    mlr3_1.8.0           plyr_1.8.9           sandwich_3.1-3      
# [19] palmerpenguins_0.1.1 rootSolve_1.8.2.4    mlr3tuning_1.7.0     zoo_1.9-1            uuid_1.2-2           igraph_2.3.4        
# [25] lifecycle_1.0.5      iterators_1.0.14     pkgconfig_2.0.3      Matrix_1.7-6         R6_2.6.1             rbibutils_2.4.1     
# [31] future_1.80.0        clue_0.3-68          digest_0.6.39        Exact_3.3            numDeriv_2016.8-1.1  spatial_7.3-19      
# [37] S4Vectors_0.50.3     mlr3misc_0.23.0      labeling_0.4.3       timechange_0.4.0     httr_1.4.9           mgcv_1.9-4          
# [43] compiler_4.6.1       proxy_0.4-29         doParallel_1.0.17    withr_3.0.3          S7_0.2.2             backports_1.5.1     
# [49] spacefillr_0.4.0     MASS_7.3-66          biomformat_1.40.0    fBasics_4052.98      gld_2.6.8            tools_4.6.1         
# [55] ape_5.8-1            glue_1.8.1           lgr_0.5.2            stabledist_0.7-2     checkmate_2.3.4      ade4_1.7-24         
# [61] cluster_2.1.8.3      reshape2_1.4.5       generics_0.1.4       gtable_0.3.6         tzdb_0.5.0           class_7.3-24        
# [67] tidyr_1.3.2          data.table_1.18.6.1  lmom_3.3             hms_1.1.4            XVector_0.52.0       rmutil_1.1.10       
# [73] BiocGenerics_0.58.1  foreach_1.5.2        pillar_1.11.1        bbotk_1.13.1         splines_4.6.1        dplyr_1.2.1         
# [79] lattice_0.23-1       survival_3.8-12      tidyselect_1.2.1     Biostrings_2.80.2    reformulas_0.4.4     mlr3mbo_1.2.1       
# [85] Seqinfo_1.2.0        IRanges_2.46.0       stats4_4.6.1         Biobase_2.72.0       expm_1.0-1           statmod_1.5.2       
# [91] factoextra_2.2.0     timeDate_4052.112    matrixStats_1.5.0    stringi_1.8.9        boot_1.3-32          codetools_0.2-20    
# [97] timeSeries_4052.112  tibble_3.3.1         cli_3.6.6            rpart_4.1.27         xtable_1.8-8         Rdpack_2.6.6        
# [103] Rcpp_1.1.2           readxl_1.5.0.1       globals_0.20.0       stable_1.1.7         parallel_4.6.1       modeest_2.5.0       
# [109] readr_2.2.0          moocore_0.3.2        paradox_1.0.1        lme4_2.0-6           listenv_1.1.0        mvtnorm_1.4-2       
# [115] lmerTest_3.2-1       scales_1.4.0         e1071_1.7-17         purrr_1.2.2          statip_0.2.3         crayon_1.5.3        
# [121] rlang_1.3.0          multcomp_1.4-32     