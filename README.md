# RNA-seq data analysis workflow with edgeR and CAMERA

A pipeline for RNA-seq data analysis that includes all the usual steps (DEA, gene set enrichment and functional annotation) using edgeR and CAMERA.

**Note:**  This is a group project during the lab course "Next-Generation Sequencing" in ETH Zurich.

A summary of the workflow employed for the data analysis is illustrated below:

![RNA-seq data analysis workflow](https://github.com/Ariwor/RNA-seq-edgeR-CAMERA/blob/master/Workflow.png)

### Quality control & adapter trimming

The quality of the raw reads was assessed with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html##fastqc), adapters were trimmed using [trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic) (parameters: PE -phred33, TRAILING: 30, MINLEN: 30).

### Transcriptome alignment

After adapter trimming, the paired reads were aligned against a reference transcriptome (gencode.v33.transcripts.fa) with [Bowtie2 v2.4.1](https://www.nature.com/articles/nmeth.1923).

### Differential expression analysis, annotation & interpretation

The resulting BAM files of the Bowtie2 transcriptome alignment were filtered with [samtools v1.10](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) (parameters -p 35 -f 2: phred score ≥ 35 and each segment properly aligned according to the aligner) and converted into read counts. Differential expression analysis (DEA) and annotation were performed with the [EdgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458) and [limma](https://academic.oup.com/nar/article/43/7/e47/2414268) packages respectively in the [Rstudio](https://rstudio.com/) environment(version 3.6.1). Specifically:

* The data were filtered for transcripts with > 1 count per million (cpm) in at least 50% of the samples.
* `calcNormFactors` was used to account for compositional differences between the libraries.
* Dispersion was estimated with `estimateDisp`.
* The [org.Hs.eg.db package](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) was used to obtain the corresponding gene-EntrezID’s.
* Three design matrices were defined according to experimental design, specifically (i) Design 1: Condition, (ii) Design 2: Condition + Library (_additive_ effect), (iii) Design 3: Condition + Condition:Library (_nested_ interaction). For the scope of this report, only Design 1 is employed in most of the downstream analysis.
* `glmFit` and `glmLRT` were used to fit a model and conduct likelihood ratio tests.
* `topTags` was used to extract a table of the top differentially expressed genes. From this table, all genes
with a log-fold change (logFC) > 3 and log-counts per million (logCPM) > 1.5 were selected for gene ontology (GO) annotation.
* `goana` was used to test for over-representation of GO terms / KEGG pathways. The top GO terms were extracted with `topGO`.
* Human C5 GO gene sets were [obtained](http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata) and after determining the read counts per gene (by summing up the counts of corresponding transcripts/isoforms), we used `camera` for competitive gene set testing (filtering for genes that had a p-value below 0.05).
* For visualization, volcano-/bar-plots and heatmaps were generated using the [ggplot2](https://www.springer.com/gp/book/9783319242750) and [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) packages respectively.
