## This code implements a basic RNA-seq data analysis workflow, using edgeR and CAMERA (currently analysing 2 experimental conditions). It is written by
## Asterios Arampatzis and Remo Bättig and is based on/inspired by the respective package tutorials/guides as well as tutorials during the course
## "Systems Genomics" in ETH Zurich.

knitr::opts_chunk$set(echo = TRUE)

## Changing my library

myPaths <- .libPaths()
myPaths <- c(myPaths, "C:/RCustom")
.libPaths(myPaths)

## Loading the necessary libraries

rm(list=ls())
library(tidyverse)
library(readr)
library(edgeR)
library(limma)
library(GO.db)
library(dplyr)
library(tidyr)

## Loading all the read count files and merging them

# Setting the working directory path to the read count text files
setwd(#path)

# Getting the full path to all the files
files <- list.files(#path)
#all_filenames <- files %>% basename() %>% as.list()

# A function that loads and merges all the files
mergefiles<-function(files){
  #Store each filename (in a more compact version) so that we can label the columns accordingly
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

## Splitting ID column to extract gene symbol

merged_separated <- separate(data = merged, col = ID, into = c("1", "2", "3", "4", "5", "6", "7", "8"), sep = "\\|")
merged_new <- merged_separated[-c(1:5, 7:8)]
colnames(merged_new)[1] <- "Gene Symbol"

## Creating DGEList object

# Substitute NA occurences with 0
merged_new[is.na(merged_new)] <- 0
# Create the DGEList object
DEA_1_2 <- DGEList(counts=merged_new[,2:45], genes=merged_new[,1])

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

## Filtering & Normalization

nrows_before_filtering <- nrow(DEA_1_2)

# Indices of genes expressed above 1 count per million in at least 2 samples
expressed_cpm_index <- rowSums(cpm(DEA_1_2)>1) >= 2

DEA_1_2 <- DEA_1_2[expressed_cpm_index, , keep.lib.sizes=FALSE]
nrows_after_filtering <- nrow(DEA_1_2)
nrows_before_filtering
nrows_after_filtering
DEA_1_2 <- calcNormFactors(DEA_1_2)
DEA_1_2$samples
# <- filterByExpr(DEA_1_2)
#DEA_1_2 <- DEA_1_2[keep, , keep.lib.sizes=FALSE]

## Data exploration

#plotMDS(DEA_1_2)

## Design matrix creation

Group <- factor(test_group)
Condition <- factor(test_condition)
#Library <- factor(test_library)

data.frame(Sample=colnames(DEA_1_2), Group, Condition)
design_1 <- model.matrix(~Group + Condition)
#design_2 <- model.matrix(~Library + Library:Condition)
rownames(design_1) <- colnames(DEA_1_2)
#rownames(design_2) <- colnames(DEA_1_2)
design_1
#design_2

## Dispersion estimation

DEA_1_2_design_1 <- estimateDisp(DEA_1_2, design_1, robust=TRUE)
#DEA_1_2_design_2 <- estimateDisp(DEA_1_2, design_2, robust=TRUE)
DEA_1_2_design_1$common.dispersion
#DEA_1_2_design_2$common.dispersion
plotBCV_design_1 <- plotBCV(DEA_1_2_design_1)
#plotBCV_design_2 <- plotBCV(DEA_1_2_design_2)

## Differential expression analysis

library('org.Hs.eg.db')
DEA_1_2_design_1$genes <- select(org.Hs.eg.db,keys=merged_new$`Gene Symbol`,columns="ENTREZID", "SYMBOL")

head(DEA_1_2$genes)
nrow(DEA_1_2$genes)

fit_1 <- glmFit(DEA_1_2_design_1, design_1)
#fit_2 <- glmFit(DEA_1_2_design_2, design_2)

lrt_1 <- glmLRT(fit_1)
#lrt_2 <- glmLRT(fit_2)
topTags(lrt_1)
#topTags(lrt_2)
summary(decideTests(lrt_1))
#summary(decideTests(lrt_2))
#plotMD(lrt_1)
#plotMD(lrt_2)

DEA_1_2_list <- topTags(lrt_1,n = nrow(lrt_1$table))$table

# Check for NA occurences
which(is.na(DEA_1_2_list)==TRUE) 
length(which(is.na(DEA_1_2_list)==TRUE))  

## Remove any NA occurences
DEA_1_2_list_updated <- na.omit(DEA_1_2_list)

### Chere if there are still NA occurences
which(is.na(DEA_1_2_list_updated)==TRUE)

# Find genes that are upregulated more than logFC of 3 and have a logCPM of 1.5
length(DEA_1_2_list_updated$genes[DEA_1_2_list_updated$logFC>3 & DEA_1_2_list_updated$logCPM>1.5])

# Store the Entrez Gene ID's of these upregulated genes
DEA_1_2_up_3logfc <- DEA_1_2_list_updated$ENTREZID[which(DEA_1_2_list_updated$logFC>3 & DEA_1_2_list_updated$logCPM>1.5)]

# Check if a nth element (e.g. 5th) is an integer
class(DEA_1_2_up_3logfc[5])

# Convert to numeric if not an integer
DEA_1_2_up_3logfc <- as.numeric(DEA_1_2_up_3logfc)

# Check again if a nth element (e.g. 5th) is an integer
class(DEA_1_2_up_3logfc[5])


## Gene ontology analysis

go_DEA_1_2 <- goana(DEA_1_2_up_3logfc)
topGO(go_DEA_1_2, ontology = "BP", number = 30)

## GSEA Analysis

# Load C5 Gene Set
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata"))

# Indexing for CAMERA
index_DEA_1_2 <- ids2indices(Hs.c5,id=rownames(DEA_1_2_design_1))

# Use CAMERA
camera_DEA_1_2 <- camera(DEA_1_2_design_1,index = index_DEA_1_2,design = DEA_1_2_design_1$design_1)

# Check best matching gene sets
head(camera_DEA_1_2)

# Check worst matching gene sets
tail(camera_DEA_1_2)

## Visualization

library(pheatmap)
# Calculate log-counts-per-million
DEA_1_2_logcpm <- cpm(DEA_1_2,prior.count = 2,log=TRUE)