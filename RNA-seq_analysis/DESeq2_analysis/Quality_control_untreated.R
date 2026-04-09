#!/usr/bin/env RScript

#command-line executable:  chmod +x Quality_control_untreated.R
# run as: ./Quality_control_untreated.R

# quality control; visualizing sample relatedness on regularized log transformed count data for HAP1 knockout RNA-seq data from Schick et al., 2019

library(DESeq2)
library("vsn")
library("biomaRt")
library("pheatmap")
library("dplyr")
library("RColorBrewer")
library(cluster)
library(factoextra)
library(stats)
library(clValid)
library(readr)
library("RNAseqQC")
library("ensembldb")
library("ggplot2")
library("purrr")
library("tidyr")
library("tibble")
library("magrittr")
library('dbx')
library("RMySQL")
library("biomaRt")
library(org.Hs.eg.db)
library(clusterProfiler)

save_heatmap <- function(x, filename, width=7, height=7) {
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

############## Upload counts data ##############
# counts: table containing count data for each gene (as rownames) and sample (as column names)
counts <- read_csv("baf_complex.expression_counts.gene_level.csv")

colnames(counts) <- gsub("RNA-seq_HAP1_", "", colnames(counts))

counts<- counts[ , -which(names(counts) %in% c("WT_r1_GFP","WT_r2_GFP", "ARID1A_r1_618_10","ARID1A_r2_618_10","ARID1A_r3_618_10","SMARCA4_r1_2878_4","SMARCA4_r2_2878_4" ,"SMARCA4_r3_2878_4"  ))]

colnames(counts) <- stringr::str_extract(colnames(counts), "[^_]*_[^_]*")

counts<- counts[ , -which(names(counts) %in% c("ACTB_r1","ACTB_r2","ACTB_r3", "BCL11A_r1",  "BCL11A_r2" , "BCL11A_r3","BCL11B_r1" , "BCL11B_r2" , "BCL11B_r3" , "BCL7A_r1" ,  "BCL7A_r2" ,  "BCL7A_r3"  , "BCL7B_r1" ,  "BCL7B_r2" ,
                                               "BCL7B_r3" ,  "DPF1_r1" ,   "DPF1_r2"  ,  "DPF1_r3" ,   "DPF2_r1", "DPF2_r2" ,   "DPF2_r3"   , "DPF3_r1"  ,  "DPF3_r2"  ,  "DPF3_r3"  , "SMARCA2_r1", "SMARCA2_r2", "SMARCA2_r3" ,"SMARCC2_r1" ,"SMARCC2_r2", "SMARCC2_r3" ,
                                               "SMARCD2_r1", "SMARCD2_r2", "SMARCD2_r3", "SMARCD3_r1" ,"SMARCD3_r2", "SMARCD3_r3"   ))]

counts <- counts %>% column_to_rownames(var = "gene_name")


##############  Upload contrast table ##############
#contrastTable: table describing the DESeq results that will be compared for each contrast

contrastTable <- read.delim("contrasts.txt", header = TRUE)

############## Upload sample table & create necessary reference levels ##############
#sampleTable: table describing the condition of each sample, including cell line, treatment, batch number, replicate number

sampleTable <- read.delim("targets.txt", header = TRUE,stringsAsFactors = TRUE)

#filter for cell lines used in Incucyte
sampleTable <- subset(sampleTable, genotype %in% c("C631","ARID1A_KO", "ARID1B_KO","ARID2_KO","BRD7_KO","PHF10_KO","PBRM1_KO","BRD9_KO","SMARCD1_KO","SMARCC1_KO","SMARCA4_KO"))

rownames(sampleTable) <- sampleTable$sample

# create base level for genotype as C631 (WT)
sampleTable$genotype <- as.factor(sampleTable$genotype)

sampleTable$genotype <- relevel(sampleTable$genotype, ref ="C631")

# get rownames of sampleTable in same order as how samples are listed in colnames of counts
sampleTable <- sampleTable[colnames(counts), ]

############## Create input data for DESeq ##############
# create model matrix using design formula
ml <- model.matrix(~ genotype, sampleTable)

# remove comparisons that don't exist (i.e. BRD9 KO treated with HU)
all.zero <- apply(ml, 2, function(x) all(x==0))
idx <- which(all.zero)
ml <- ml[,-idx]

# create input data for DESeq
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleTable,
                              design = ml)

############## regularized log transformation ################

# apply regularized log transformation to count data
rlogdf <- rlog(dds, blind=TRUE)

############## Heatmap of sample-to-sample euclidean distances ################

# heatmap of distance matrix provides overview over similarities and dissimilarities btwn samples
sampleDists <- dist(t(assay(rlogdf)),method="euclidean")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlogdf$genotype)
colnames(sampleDistMatrix) <- paste(rlogdf$genotype)#NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dist_heatmap <- pheatmap(sampleDistMatrix,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         col=colors,
                         fontsize_col = 7,
                         fontsize_row = 7,
                         fontsize = 11,
                         main = "Heatmap of sample-to-sample distances")


NewFilename = sprintf("euclidean_dist_rlog_heatmap.pdf")


save_heatmap(dist_heatmap, NewFilename)
