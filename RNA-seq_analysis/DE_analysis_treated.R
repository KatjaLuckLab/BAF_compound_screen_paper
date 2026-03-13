#!/usr/bin/env RScript

#command-line executable:  chmod +x DE_analysis_treated.R
# run as: ./DE_analysis_treated.R

# Script computes log2 fold change in gene expreesion KO treat/dmso vs WT treat/dmso for each tested condition from the treated RNA-seq data from Spang et al., 2026

library(DESeq2)
library(stringr)
library("dplyr")
library("plyr")
library(tidyverse)

############## Upload counts data ##############
# counts: table containing count data for each gene (as rownames) and sample (as column names)

# load counts data
counts <- read.delim("/subread_counts.txt", header = TRUE,row.names = 1)

colnames(counts) <- gsub("sfb_schick_2025_01_", "", colnames(counts))

colnames(counts) <- gsub("HAP1_", "", colnames(counts))

colnames(counts) <- gsub(".cutadapt", "", colnames(counts))


##############  Upload contrast table ##############
#contrastTable: table describing the DESeq results that will be used for each contrast

#load interaction table
contrastTable <- read.delim("interactions.txt", header = TRUE)

############## Upload sample table & create necessary reference levels ##############
#sampleTable: table describing the condition of each sample, including cell line, treatment, batch number, replicate number

# load targets table
sampleTable <- read.delim("treated_targets.txt", header = TRUE)

# set sample names as the rownames of dataframe
rownames(sampleTable) <- sampleTable$sample

# loop through the different compound and knockout genotype combinations
sampleTableKOsubset <- subset(sampleTable, genotype != "C631")

for(genotype in unique(sampleTableKOsubset$genotype)){

  sampleTableKOsubset2 <- sampleTableKOsubset[sampleTableKOsubset$genotype == genotype, ]

  treatmentList <- unique(sampleTableKOsubset2$treatment)
  treatmentList <- treatmentList[treatmentList != "DMSO"]

  for(treatment in treatmentList){

    sampleTablecmpndsubset <- sampleTable[sampleTable$genotype == genotype & sampleTable$treatment == treatment,]

    sampleTablesubset <- sampleTable[sampleTable$genotype %in% c("C631",genotype) & sampleTable$treatment %in% c("DMSO",treatment) & sampleTable$batch ==unique(sampleTablecmpndsubset$batch), ]


    # create base level for genotype as C631 (WT)
    sampleTablesubset$genotype <- as.factor(sampleTablesubset$genotype)
    sampleTablesubset$genotype <- relevel(sampleTablesubset$genotype, ref ="C631")


    # create base level for treatment as DMSO
    sampleTablesubset$treatment <- as.factor(sampleTablesubset$treatment)
    sampleTablesubset$treatment <- relevel(sampleTablesubset$treatment, ref ="DMSO")

    # get rownames of sampleTable in same order as how samples are listed in colnames of counts
    countsSubset <- counts[,rownames(sampleTablesubset) ]


    ############## Create input data for DESeq ##############
    # create model matrix using design formula
    ml <- model.matrix(~ genotype + treatment + genotype:treatment, sampleTablesubset)

    # create input data for DESeq
    dds <- DESeqDataSetFromMatrix(countData = countsSubset,
                                  colData = sampleTablesubset,
                                  design = ml)

    ############## Run DE analysis ##############
    # run DE analysis
    dds2 <- DESeq(dds)

    resultsNames(dds2)

    res<- lfcShrink(dds2, coef=4, type="ashr",format="DataFrame")

    contrastTablesubset <- contrastTable[contrastTable$interaction==resultsNames(dds2)[4],]

    res$cell_line <- gsub("__","",str_extract(contrastTablesubset$interaction.name,"^(.*)[_]{2}"))

    res$treatment <- gsub("__","",str_extract(contrastTablesubset$interaction.name,"[_]{2}(.*)"))

    NewFilename = sprintf( "%s.txt",contrastTablesubset$interaction.name)

    write.table(res, file =NewFilename, sep = "\t",row.names = TRUE, col.names = NA)
  }

}
