#!/usr/bin/env RScript

#command-line executable:  chmod +x DE_analysis_untreated.R
# run as: ./DE_analysis_untreated.R

# Script for calculating log2 fold change of gene expression KO treat/dmso vs WT treat/dmso for the HAP1 BAF knockout RNA-seq dataset from Schick et al., 2019

library(DESeq2)
library(stringr)
library("dplyr")
library("plyr")
library(tidyverse)
library(stats)

############## Upload counts data ##############
# counts: table containing count data for each gene (as rownames) and sample (as column names)

counts <- read_csv("baf_complex.expression_counts.gene_level.csv")

colnames(counts) <- gsub("RNA-seq_HAP1_", "", colnames(counts))

# only keep same 10 cell lines investigated in the other datasets (i.e. Incucyte)
counts<- counts[ , -which(names(counts) %in% c("WT_r1_GFP","WT_r2_GFP", "ARID1A_r1_618_10","ARID1A_r2_618_10","ARID1A_r3_618_10","SMARCA4_r1_2878_4","SMARCA4_r2_2878_4" ,"SMARCA4_r3_2878_4"))]

colnames(counts) <- stringr::str_extract(colnames(counts), "[^_]*_[^_]*")

counts <- counts %>% column_to_rownames(var = "gene_name")

##############  Upload contrast table ##############
#contrastTable: table describing the DESeq results that will be compared for each contrast

contrastTable <- read.delim("contrasts.txt", header = TRUE)

############## Upload sample table & create necessary reference levels ##############
#sampleTable: table describing the condition of each sample, including cell line, treatment, batch number, replicate number

sampleTable <- read.delim("untreated_targets.txt", header = TRUE,stringsAsFactors = TRUE)

rownames(sampleTable) <- sampleTable$sample

genotypeList <- unique(sampleTable$genotype)
genotypeList <- genotypeList[genotypeList != "C631"]

for(genotype in genotypeList){

  sampleTablesubset <- sampleTable[sampleTable$genotype %in% c("C631",genotype), ]

  # create base level for treatment as DMSO
  sampleTablesubset$genotype <- as.factor(sampleTablesubset$genotype)
  sampleTablesubset$genotype <- relevel(sampleTablesubset$genotype, ref ="C631")

  # get rownames of sampleTable in same order as how samples are listed in colnames of counts
  countsSubset <- counts[,rownames(sampleTablesubset) ]

  ############## Create input data for DESeq ##############
  # create model matrix using design formula
  ml <- model.matrix(~ genotype, sampleTablesubset)

  # remove comparisons that don't exist (i.e. BRD9 KO treated with HU)
  all.zero <- apply(ml, 2, function(x) all(x==0))
  idx <- which(all.zero)
  ml <- ml[,-idx]

  # create input data for DESeq
  dds <- DESeqDataSetFromMatrix(countData = countsSubset,
                                colData = sampleTablesubset,
                                design = ml)

  ############## Run DE analysis ##############
  # run DE analysis
  dds2 <- DESeq(dds)

  resultsNames(dds2)

  res<- lfcShrink(dds2, coef=2, type="ashr",format="DataFrame")

  contrastTablesubset <- contrastTable[gsub("-","",str_extract(contrastTable$contrast,"^(.*)[-]"))==genotype,]

  res$cell_line <- genotype

  NewFilename = sprintf("%s.txt",contrastTablesubset$contrast.name)

  write.table(res, file =NewFilename, sep = "\t",row.names = TRUE, col.names = NA)
}
