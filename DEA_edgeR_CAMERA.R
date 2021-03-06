## This code implements a basic RNA-seq data analysis workflow, using edgeR and CAMERA (currently analysing 2 experimental conditions). It is written by
## Asterios Arampatzis and Remo B�ttig and is based on/inspired by the respective package tutorials/guides as well as tutorials/shared code during the course
## "Systems Genomics" in ETH Zurich.

knitr::opts_chunk$set(echo = TRUE)

# Changing my library
myPaths <- .libPaths()
myPaths <- c(myPaths, "C:/RCustom")
.libPaths(myPaths)

# Loading the necessary libraries


rm(list=ls())
library(tidyverse)
library(readr)
library(edgeR)
library(limma)
library(GO.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library('org.Hs.eg.db')
library(KEGGREST)
library(pathview)


## Loading all the read count files and merging them


# Setting the working directory path to the read count text files
setwd(#path)
  
# Getting the full path to all the files
  files <- list.files(#path)
    #all_filenames <- files %>% basename() %>% as.list()
    
# A function that loads and merges all the files
    mergefiles<-function(files){
      #vStore each filename (in a more compact version) so that we can label the columns accordingly
      all_filenames <- files %>% basename() %>% as.list()
      filenames_short<-all_filenames %>% gsub("_readcount_q35_f2.txt", "", .) %>% sub("HKHNCDRXX_", "",.) %>% sub("_S[0-9]", "",.) %>%as.list()
      
      anchor<-as.data.frame(read_table(files[1], col_names = F))
      
      colnames(anchor)<-c(as.character(filenames_short[1]), "ID")
      for (i in c(2:length(files))){
        tmp<-as.data.frame(read_table(files[i], col_names = F))
        colnames(tmp)<-c(as.character(filenames_short[i]), "ID")
        anchor<-merge(anchor, tmp, by = "ID", all = TRUE) #full join
        #colnames(ctrl_merged)<-c("ID","counts_ctrll001","counts_ctrl002")
      }
      #anchor[is.na(anchor)] <- 0 #replace all NA with zero
      anchor
    }
    
merged<-mergefiles(files)
    
    
## Spliting ID column to extract gene symbol
    
    
merged_separated <- separate(data = merged, col = ID, into = c("1", "2", "3", "4", "5", "6", "7", "8"), sep = "\\|")
    
merged_new <- merged_separated[-c(1:5, 7:8)]
colnames(merged_new)[1] <- "Gene Symbol"
    
    
## Creating DGEList object
    
    
# Substitute NA occurences with 0
merged_new[is.na(merged_new)] <- 0
    
# Create the DGEList object
DEA_1_2 <- DGEList(counts=merged_new[,2:ncol(merged)], genes=merged_new[,1]) # Can also remove rows with zero counts using remove.zeros = TRUE
    
    
## Extracting group and condition from column headers
    
    
headers <- grep("L00", names(merged), value = TRUE)
test_group <- sub("^\\D*(\\d+).*$", "\\1", headers)
test_condition <- sub('^[^_]*_(\\d+).*', '\\1', headers)
test_group
test_condition
    
    
## Extracting type of library preparation from column headers
    
    
mRNA_riboZ_headers <- grep("mRNA||riboZ", names(merged), value = TRUE)
mRNA_riboZ_headers <- mRNA_riboZ_headers[mRNA_riboZ_headers!="ID"]
mRNA_riboZ_headers
test_library <- gsub("[^mRN || riboZ]", "", mRNA_riboZ_headers)
test_library
    
    
# A bit of work for the heatmap labeling
    
    
headers_heatmap <- substr(mRNA_riboZ_headers, start = 1, stop = 4)
headers_heatmap_new <- sub("_", ".", headers_heatmap)
headers_heatmap_new <- sub("_", "", headers_heatmap_new)
    
    
## Data exploration and quality assessment
    
    
DEA_1_2_logcounts <- log2(DEA_1_2$counts + 1)
head(DEA_1_2_logcounts)
    
# Example histogram for logcounts (can also do a boxplot)
hist(DEA_1_2_logcounts[,1], main = "Histogram for example sample counts", xlab = "Counts")
    
# Saving
png("1_2_Histogram_before.png")
hist(DEA_1_2_logcounts[,1], main = "Histogram for example sample counts", xlab = "Counts")
dev.off()
    
# Example MA-plot (samples 1_1_mRNA_L001 and 10_1_riboZ_L001)
limma::plotMA(DEA_1_2_logcounts[ ,1:5], xlab = "M", ylab = "A", main = "")
abline(h = 0, col = "red")
    
# MDS plot
plotMDS(DEA_1_2_logcounts, labels=test_condition, cex=1.5, col=test_group)
    
# Saving
png("1_2_MDS.png")
plotMDS(DEA_1_2_logcounts, labels=test_condition, cex=1.5, col=test_group)
dev.off()
    
    
## Filtering & Normalization
    
    
nrows_before_filtering <- nrow(DEA_1_2)
    
# Indices of genes expressed above 1 count per million in at least 50% of the samples
expressed_cpm_index <- rowSums(cpm(DEA_1_2)>1) >= (ncol(merged)-1)/2
    
DEA_1_2 <- DEA_1_2[expressed_cpm_index, , keep.lib.sizes=FALSE]
nrows_after_filtering <- nrow(DEA_1_2)
nrows_before_filtering
nrows_after_filtering
DEA_1_2 <- calcNormFactors(DEA_1_2)
DEA_1_2$samples
#filterByExpr(DEA_1_2)
    
# Data exploration after normalization
    
DEA_1_2_logcounts_after_norm <- log2(DEA_1_2$counts + 1)
#head(DEA_1_2_logcounts)
    
# Example histogram for logcounts (can also do a boxplot)
hist(DEA_1_2_logcounts_after_norm[,1], main = "Histogram for example samle counts", xlab = "Counts")
    
# Saving
png("1_2_Histogram_after.png")
hist(DEA_1_2_logcounts_after_norm[,1], main = "Histogram for example sample counts", xlab = "Counts")
dev.off()
    
# Example MA-plot (two example samples)
limma::plotMA(DEA_1_2_logcounts[ ,1:5], xlab = "M", ylab = "A", main = ""), abline(h = 0, col = "red")
    
    
## Design matrix creation
    
    
Group <- factor(test_group)
Condition <- factor(test_condition)
Library <- factor(test_library)
    
data.frame(Sample=colnames(DEA_1_2), Group, Condition, Library)
    
# Effect of Condition
design_1 <- model.matrix(~Condition)
    
# Effect (additive) of Condition and Library
design_2 <- model.matrix(~Condition + Library)
    
# Effect (nested interaction) of Condition and Library
design_3 <- model.matrix(~Condition + Condition:Library)
    
rownames(design_1) <- colnames(DEA_1_2)
rownames(design_2) <- colnames(DEA_1_2)
rownames(design_3) <- colnames(DEA_1_2)
    
design_1
design_2
design_3
    
    
## Dispersion estimation
    
    
DEA_1_2_design_1 <- estimateDisp(DEA_1_2, design_1, robust=TRUE)
DEA_1_2_design_2 <- estimateDisp(DEA_1_2, design_2, robust=TRUE)
DEA_1_2_design_3 <- estimateDisp(DEA_1_2, design_3, robust=TRUE)
    
DEA_1_2_design_1$common.dispersion
DEA_1_2_design_2$common.dispersion
DEA_1_2_design_3$common.dispersion
    
plotBCV_design_1 <- plotBCV(DEA_1_2_design_1)
plotBCV_design_2 <- plotBCV(DEA_1_2_design_2)
plotBCV_design_3 <- plotBCV(DEA_1_2_design_3)
    
# Saving
png(file="1_2_BCV_design_1.png")
plotBCV(DEA_1_2_design_1)
dev.off()
    
    
## Differential expression analysis
    
    
DEA_1_2_design_1$genes <- select(org.Hs.eg.db,keys=merged_new$`Gene Symbol`,columns="ENTREZID", "SYMBOL")
DEA_1_2_design_2$genes <- select(org.Hs.eg.db,keys=merged_new$`Gene Symbol`,columns="ENTREZID", "SYMBOL")
DEA_1_2_design_3$genes <- select(org.Hs.eg.db,keys=merged_new$`Gene Symbol`,columns="ENTREZID", "SYMBOL")
    
head(DEA_1_2_design_1$genes)
nrow(DEA_1_2_design_1$genes)
    
head(DEA_1_2_design_2$genes)
nrow(DEA_1_2_design_2$genes)
    
head(DEA_1_2_design_3$genes)
nrow(DEA_1_2_design_3$genes)
    
fit_1 <- glmFit(DEA_1_2_design_1, design_1)
fit_2 <- glmFit(DEA_1_2_design_2, design_2)
fit_3 <- glmFit(DEA_1_2_design_2, design_3)
    
colnames(fit_1)
colnames(fit_2)
colnames(fit_3)
    
    
    
lrt_1 <- glmLRT(fit_1)
lrt_2 <- glmLRT(fit_2, coef = 3)
lrt_3 <- glmLRT(fit_3,  coef = 4)
    
topTags(lrt_1)
topTags(lrt_2)
topTags(lrt_3)
    
    
    
summary(decideTests(lrt_1))
summary(decideTests(lrt_2))
summary(decideTests(lrt_3))
    
plotMD(lrt_1)
plotMD(lrt_2)
plotMD(lrt_3)
    
DEA_1_2_lrt_1 <- topTags(lrt_1,n = nrow(lrt_1$table))$table
DEA_1_2_lrt_2 <- topTags(lrt_2,n = nrow(lrt_2$table))$table
DEA_1_2_lrt_3 <- topTags(lrt_3,n = nrow(lrt_3$table))$table
    
# Check for NA occurences
which(is.na(DEA_1_2_lrt_1)==TRUE) 
length(which(is.na(DEA_1_2_lrt_1)==TRUE))  
    
which(is.na(DEA_1_2_lrt_2)==TRUE) 
length(which(is.na(DEA_1_2_lrt_2)==TRUE))  
    
which(is.na(DEA_1_2_lrt_3)==TRUE) 
length(which(is.na(DEA_1_2_lrt_3)==TRUE))  
    
# Remove any NA occurences
DEA_1_2_lrt_1_updated <- na.omit(DEA_1_2_lrt_1)
DEA_1_2_lrt_2_updated <- na.omit(DEA_1_2_lrt_2)
DEA_1_2_lrt_3_updated <- na.omit(DEA_1_2_lrt_3)
    
# Check if there are still NA occurences
which(is.na(DEA_1_2_lrt_1_updated)==TRUE)
which(is.na(DEA_1_2_lrt_2_updated)==TRUE)
which(is.na(DEA_1_2_lrt_3_updated)==TRUE)
    
    
    
# Find genes that are upregulated more than logFC of 3 and have a logCPM of 1.5
length(DEA_1_2_lrt_1_updated$SYMBOL[DEA_1_2_lrt_1_updated$logFC>3 & DEA_1_2_lrt_1_updated$logCPM>1.5])
length(DEA_1_2_lrt_2_updated$SYMBOL[DEA_1_2_lrt_2_updated$logFC>3 & DEA_1_2_lrt_2_updated$logCPM>1.5])
length(DEA_1_2_lrt_3_updated$SYMBOL[DEA_1_2_lrt_3_updated$logFC>3 & DEA_1_2_lrt_3_updated$logCPM>1.5])
    
# Store the Entrez Gene ID's of these upregulated genes
DEA_1_2_up_3logfc_lrt_1 <- DEA_1_2_lrt_1_updated$ENTREZID[which(DEA_1_2_lrt_1_updated$logFC>3 & DEA_1_2_lrt_1_updated$logCPM>1.5)]
DEA_1_2_up_3logfc_lrt_2 <- DEA_1_2_lrt_2_updated$ENTREZID[which(DEA_1_2_lrt_2_updated$logFC>3 & DEA_1_2_lrt_2_updated$logCPM>1.5)]
DEA_1_2_up_3logfc_lrt_3 <- DEA_1_2_lrt_3_updated$ENTREZID[which(DEA_1_2_lrt_3_updated$logFC>3 & DEA_1_2_lrt_3_updated$logCPM>1.5)]
    
# Check if a nth element (e.g. 5th) is an integer
class(DEA_1_2_up_3logfc_lrt_1[5])
class(DEA_1_2_up_3logfc_lrt_2[5])
class(DEA_1_2_up_3logfc_lrt_3[5])
    
# Convert to numeric if not an integer
DEA_1_2_up_3logfc_lrt_1 <- as.numeric(DEA_1_2_up_3logfc_lrt_1)
DEA_1_2_up_3logfc_lrt_2 <- as.numeric(DEA_1_2_up_3logfc_lrt_2)
DEA_1_2_up_3logfc_lrt_3 <- as.numeric(DEA_1_2_up_3logfc_lrt_3)
    
# Check again if a nth element (e.g. 5th) is an integer
class(DEA_1_2_up_3logfc_lrt_1[5])
class(DEA_1_2_up_3logfc_lrt_2[5])
class(DEA_1_2_up_3logfc_lrt_3[5])
    
    
## Gene ontology analysis
    
    
go_DEA_1_2_lrt_1 <- goana(DEA_1_2_up_3logfc_lrt_1, species = "Hs")
go_DEA_1_2_lrt_2 <- goana(DEA_1_2_up_3logfc_lrt_2, species = "Hs")
go_DEA_1_2_lrt_3 <- goana(DEA_1_2_up_3logfc_lrt_3, species = "Hs")
    
    
    
topGO(go_DEA_1_2_lrt_1, ontology = "BP", number = 30)
topGO(go_DEA_1_2_lrt_2, ontology = "BP", number = 30)
topGO(go_DEA_1_2_lrt_3, ontology = "BP", number = 30)
    
    
## Heatmap
    
    
# Calculate log-counts-per-million
DEA_1_2_logcpm <- cpm(DEA_1_2,prior.count = 2,log=TRUE)
    
# Store top 50 genes by p-value (can also do by logFC with the argument sort.by = "logFC")
DEA_1_2_top50_pval_lrt_1 <- topTags(lrt_1, n=50)
entrez_ids_top50_pval_lrt_1 <- rownames(DEA_1_2_top50_pval_lrt_1)
  
DEA_1_2_top50_count_matrix <- matrix(nrow=50,ncol=ncol(DEA_1_2_logcpm))
for(i in 1:length(entrez_ids_top50_pval_lrt_1)){
      DEA_1_2_top50_count_matrix[i,] <- DEA_1_2_logcpm[which(rownames(DEA_1_2_logcpm)==entrez_ids_top50_pval_lrt_1[i]),]
    }
    
rownames(DEA_1_2_top50_count_matrix) <- DEA_1_2_top50_pval_lrt_1$table$SYMBOL[1:50]
colnames(DEA_1_2_top50_count_matrix) <- headers_heatmap_new
    
Heatmap_design_1 <- pheatmap(DEA_1_2_top50_count_matrix, fontsize = 5)
Heatmap_design_1
    
# Saving
pdf("1_2_Heatmap_design_1.pdf")
Heatmap_design_1 <- pheatmap(DEA_1_2_top50_count_matrix, fontsize = 5)
dev.off()
    
    
## KEGG Pathway Analysis
    
    
DEA_1_2_up_3logfc_lrt_1_SYMBOL <- DEA_1_2_lrt_1_updated$SYMBOL[which(DEA_1_2_lrt_1_updated$logFC>3 & DEA_1_2_lrt_1_updated$logCPM>1.5)]
    
DEA_1_2_up_3logfc_lrt_1_SYMBOL
    
# Check which pathways correspond to the upregulated genes
up3_logfc_kegga <- kegga(de=DEA_1_2_up_3logfc_lrt_1,species = "Hs")
up3_logfc_topkeg <- topKEGG(up3_logfc_kegga)
    
pathid <- rownames(up3_logfc_topkeg)[1]
rownames(up3_logfc_topkeg)[1]
    
top_path <- keggGet(pathid)
    
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "#id of top path", species = "hsa", out.suffix = "gse16873")
    
    
## GSEA Analysis
    
    
# Load C5 Gene Set
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata"))
    
# Indexing for CAMERA
index_DEA_1_2_design_1 <- ids2indices(Hs.c5,id=rownames(DEA_1_2_design_1))
index_DEA_1_2_design_2 <- ids2indices(Hs.c5,id=rownames(DEA_1_2_design_2))
index_DEA_1_2_design_3 <- ids2indices(Hs.c5,id=rownames(DEA_1_2_design_3))
    
    
    
# Additional preparation
    
# Calculate log-counts-per-million
DEA_1_2_logcpm <- cpm(DEA_1_2,prior.count = 2,log=TRUE)
    
# Store statistically significant genes by p-value (can also do by logFC with the argument sort.by = "logFC")
DEA_1_2_pval_lrt_1 <- topTags(lrt_1, n=nrow(DEA_1_2), p.value = 0.05)
    
entrez_ids_pval_lrt_1 <- rownames(DEA_1_2_pval_lrt_1)
    
DEA_1_2_count_matrix <- matrix(nrow=nrow(DEA_1_2_pval_lrt_1$table),ncol=ncol(DEA_1_2_logcpm))
for(i in 1:length(entrez_ids_pval_lrt_1)){
      DEA_1_2_count_matrix[i,] <- DEA_1_2_logcpm[which(rownames(DEA_1_2_logcpm)==entrez_ids_pval_lrt_1[i]),]
    }
    
rownames(DEA_1_2_count_matrix) <- DEA_1_2_pval_lrt_1$table$SYMBOL
colnames(DEA_1_2_count_matrix) <- headers_heatmap_new
    
# Use CAMERA
    
camera_DEA_1_2_design_1 <- camera(DEA_1_2_count_matrix, index = index_DEA_1_2_design_1,design = DEA_1_2_design_1$design)
camera_DEA_1_2_design_2 <- camera(DEA_1_2_count_matrix,index = index_DEA_1_2_design_2,design = DEA_1_2_design_2$design)
camera_DEA_1_2_design_3 <- camera(DEA_1_2_count_matrix,index = index_DEA_1_2_design_3,design = DEA_1_2_design_3$design)
    
    
    
    
# Check best matching gene sets
head(camera_DEA_1_2_design_1)
head(camera_DEA_1_2_design_2)
head(camera_DEA_1_2_design_3)
    
# Check worst matching gene sets
tail(camera_DEA_1_2_design_1)
tail(camera_DEA_1_2_design_2)
tail(camera_DEA_1_2_design_3)
    
    
## Volcano plot
    
    
volcano_lrt_1 <- cbind(lrt_1$table$logFC, -log10(lrt_1$table$PValue))
colnames(volcano_lrt_1) <- c("logFC", "negLogPval")
    
volcano_lrt_1_new <- cbind(lrt_1$table$logFC, -log10(lrt_1$table$PValue), lrt_1$genes$SYMBOL)
test_genelabels <- head(volcano_lrt_1_new[, 3], 10)
    
    
volcano_lrt_1_df <- as.data.frame(volcano_lrt_1, c('logFC', 'negLogPval'))
colnames(volcano_lrt_1_df) <- c("logFC", "negLogPval")
    
volcano_lrt_1_ordered <- volcano_lrt_1_df[order(volcano_lrt_1_df$negLogPval, decreasing="TRUE"),]
    
volcano_lrt_1_ordered$genelabels <- ""
volcano_lrt_1_ordered$genelabels[1:10] <- test_genelabels
    
threshold <- lrt_1$table$PValue < 0.05 & abs(lrt_1$table$logFC) > 2
volcano_1_2 <- ggplot(volcano_lrt_1_df) +
  geom_point(aes(x=logFC,y=negLogPval, colour=threshold), show.legend=F) +
  ggtitle("Condition 2 vs Condition 1") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme_bw()
volcano_1_2
    
# Saving
png("1_2_volcano.png")
volcano_1_2
dev.off()
    
    
## Bar plot (with cumulative DEA results across all experiments/designs)
    
    
Results_Bar <- read.table(#import .csv file with cumulative results, header = TRUE, sep = ",")
      
Results_Bar_df <- as.data.frame(Results_Bar)
      
# Saving
pdf("Bar_Results.pdf")
g_bar <- g_bar <- ggplot(Results_Bar_df, aes(x=�..Condition, y=Genes, fill=Direction))
    g_bar + geom_bar(stat = "identity", width = 0.5) +
    theme(axis.text.x = element_text(angle=45,  hjust = 1), panel.background = element_rect(fill = "white", colour = "white", size = 2, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "#BFD5E3"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#BFD5E3"), text = element_text(size=11)) +
    labs(title="Differential Expression Results", 
    subtitle="Number of up/down/nonsignificantly regulated genes", x="", y="Number of genes")
    dev.off()