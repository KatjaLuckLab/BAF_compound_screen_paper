#!/usr/bin/env RScript

#command-line executable:  chmod +x Reactome_DEG_untreated.R
# run as: ./Reactome_DEG_untreated.R

# this script is for visualizing DEGs included in the general pathways downloaded from the Reactome database for HAP1 BAF knockout RNA-seq data from Schick et al., 2019
################################################## Loading Necessary Libraries and Source Code #########################################################################
setwd("~/GitHub_projects/IMB/BAF_Project/BAF_RNA-Seq_Analysis")
source("./R/clustering_functions.R")


library("dplyr")
library("plyr")
library("readr")
library("biomaRt")
library(cluster)
library(factoextra)
library(stats)
library(clValid)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(ReactomePA)
library(NbClust)
library(fpc)
library(EnhancedVolcano)
library(stringr)
library(openxlsx)
library(ggvenn)
library(openxlsx)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()}


#RNA seq dataset:
dataset_path<-" "

#DE analysis comparison performed
comparison <- "KO_vs_WT_DEGs"
filename_pattern <-"C631*\\.txt"

SubDir <- sprintf("/%s/%s/Genes/Protein-coding/", dataset_path,comparison)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Exists!")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Created!")
}

setwd(sprintf("/%s/%s/Genes/Protein-coding/",dataset_path,comparison))

########## upload all gene lists - entire pathways and subpathways ###########
pathways_list <- c("Cell_Cycle_Checkpoints","DNA_Damage_Repair","Growth_Signaling")

pathways_df_list <- list()

subpathways_files_list <- list()

for(n in seq_along(pathways_list)){

  pathway_df <- list.files(path = sprintf("/Reactome_Gene_Lists/%s/",pathways_list[n]), pattern = "*.xlsx", full.names = TRUE) %>%
          lapply(read.xlsx) %>%
          bind_rows

  pathway_df <- pathway_df %>% mutate(Identifier = factor(Identifier)) %>% distinct(Identifier,.keep_all = TRUE)

  pathway_df<- clusterProfiler::bitr(geneID = pathway_df$Identifier, fromType = "UNIPROT", toType = c("ENSEMBL","GENENAME","SYMBOL"), OrgDb = "org.Hs.eg.db")


  pathways_df_list[[n]] <- pathway_df

  subpathways_files_list[[n]] <- list.files(path = sprintf("/Reactome_Gene_Lists/%s/",pathways_list[n]),pattern = "\\.xlsx$")

}
names(pathways_df_list) <- pathways_list

names(subpathways_files_list) <-pathways_list


########## create metadata containing sub-pathway gene subsets ###########

subpathways_df_list_list <- list()

for(n in seq_along(subpathways_files_list)){
  subpathways_df_list <- list()

  for(i in seq_along(subpathways_files_list[[n]])){

    subpathway_df <- read.xlsx(gsub(" ","",paste("/Reactome_Gene_Lists/",names(subpathways_files_list)[n],"/",subpathways_files_list[[n]][i])))

    subpathway_df <- subpathway_df %>% mutate(Identifier = factor(Identifier)) %>% distinct(Identifier,.keep_all = TRUE)

    subpathway_df<- bitr(geneID = subpathway_df$Identifier, fromType = "UNIPROT", toType = c("ENSEMBL","GENENAME","SYMBOL"), OrgDb = "org.Hs.eg.db")

    subpathway_df <- subpathway_df %>% mutate(pathway = names(subpathways_files_list)[n],
                                              subpathway = gsub(".xlsx","",subpathways_files_list[[n]][i]))

    subpathways_df_list[[i]] <- subpathway_df

  }

  names(subpathways_df_list) <- gsub(".xlsx","",subpathways_files_list[[n]])

  subpathways_df_list <- bind_rows(subpathways_df_list)

  subpathways_df_list_list[[n]] <- subpathways_df_list
}

names(subpathways_df_list_list) <- pathways_list

metadata <- bind_rows(subpathways_df_list_list)

########## upload all DE data all KOs & all Compounds & merge all data - filter for genes in Reactome database curated lists ##########

  de_data <- list.files(path = sprintf("/%s/%s/",dataset_path,comparison),
                        pattern = filename_pattern, full.names = TRUE) %>%
    lapply(read.delim) %>%
    bind_rows

  # remove NA
  de_data <- de_data[!is.na(de_data$padj),]

  de_data <- de_data %>% dplyr::rename(SYMBOL = X)

  #filter for cell lines used in Incucyte
  de_data <- subset(de_data, cell_line %in% c("ARID1A_KO", "ARID1B_KO","ARID2_KO","BRD7_KO","PHF10_KO","PBRM1_KO","BRD9_KO","SMARCD1_KO","SMARCC1_KO","SMARCA4_KO"))

################ Filter Pathway Lists for Significantly DEGs in at least one condition ##################

DEG_pathways_df_ensembl_list <- list()

signif_DEG_pathways_df_ensembl_list <- list()

DEG_pathways_df_hgnc_list <- list()

for(n in seq_along(pathways_df_list)){

  pathway_df_list_subset <- pathways_df_list[[n]]

  de_data_subset <- subset(de_data, SYMBOL %in% c(pathways_df_list[[n]]$SYMBOL))

  # converting significance values into stars
  de_data_subset <- de_data_subset %>% mutate(significance =case_when(padj <0.001 ~ "***",
                                                                      padj < 0.01 ~ "**",
                                                                      padj < 0.05 ~ "*",
                                                                      TRUE ~ " "))


  de_data_subset_hgnc <- de_data_subset%>% dplyr::left_join(pathway_df_list_subset[,c(2,3,4)])

  # converting significance values into stars
  de_data_subset_hgnc <- de_data_subset_hgnc %>% mutate(significance =case_when(padj <0.001 ~ "***",
                                                                                padj < 0.01 ~ "**",
                                                                                padj < 0.05 ~ "*",
                                                                                TRUE ~ " "))

  DEG_pathways_df_ensembl_list[[n]] <- de_data_subset

  DEG_pathways_df_hgnc_list[[n]] <- de_data_subset_hgnc

  de_data_subset2 <- de_data_subset %>%
    group_by(SYMBOL, cell_line) %>%
    dplyr::filter(padj < 0.05)

  signif_DEG_pathways_df_ensembl_list[[n]] <- de_data_subset2

}

names(DEG_pathways_df_ensembl_list) <- names(pathways_df_list)

names(DEG_pathways_df_hgnc_list) <- names(pathways_df_list)

names(signif_DEG_pathways_df_ensembl_list) <- names(pathways_df_list)


###### loop through all pathways to plot heatmap ##########

for(l in seq_along(DEG_pathways_df_hgnc_list)){

DEG_pathways_df_hgnc_list_subset <- DEG_pathways_df_hgnc_list[[l]]

signif_DEG_pathways_df_ensembl_list_subset<- signif_DEG_pathways_df_ensembl_list[[l]]

metadata_subset <-subset(metadata, pathway==names(DEG_pathways_df_hgnc_list)[l])

print(names(DEG_pathways_df_hgnc_list)[l])

# filter for gene symbols with DE value for more than one unique condition (cell line + compound) - want to capture DEG across conditions
de_data_subset <- DEG_pathways_df_hgnc_list_subset %>% group_by(SYMBOL) %>% dplyr::filter(n_distinct(cell_line)>1)

de_data_subset <- de_data_subset %>% distinct(SYMBOL,log2FoldChange, cell_line,.keep_all = TRUE)

##########  filter for ENSEMBL IDs of DEGs with p-value < 0.01 in at least one condition   ##########
de_data_subset_path <- subset(de_data_subset, SYMBOL %in% signif_DEG_pathways_df_ensembl_list_subset$SYMBOL)

# filter out possible duplicate entries
de_data_subset_path <- de_data_subset_path %>% distinct(SYMBOL,log2FoldChange, cell_line,.keep_all = TRUE)


if(isTRUE(names(DEG_pathways_df_hgnc_list)[l]=="DNA_Damage_Repair" & filter_DEG_subset == TRUE)==TRUE){

  # taking subset of DDR genes that have p-value < 0.001 for at least one condition
  de_data_subset_path2 <- de_data_subset_path %>% group_by(SYMBOL) %>% filter(padj < 0.001) %>% filter(n()>0)

  de_data_subset_path <- subset(de_data_subset_path, SYMBOL %in% de_data_subset_path2$SYMBOL)

  #subset gene symbols in metadata
  metadata_subset <- subset(metadata_subset, rownames(metadata_subset)%in%de_data_subset_path$SYMBOL)

}


if(isTRUE(names(DEG_pathways_df_hgnc_list)[l]=="Growth_Signaling" & filter_DEG_subset == TRUE)==TRUE){

  #taking subset of growth signaling genes that have similar functions pertaining to the pathways highlighted in manuscript
  de_data_subset_path <- subset(de_data_subset_path, SYMBOL %in% c("TGFB3","TGFBR3","BMP2","GREM2","CHRDL1","BAMBI","MET","KIT","PDGFRA","PDGFRB","IGF1R","FLT4","AXL",
                                                                   "FZD10","SFRP1","ROR1","TLE1","COL1A1","COL3A1","COL5A2","COL11A1","COL4A3",
                                                                   "FN1","THBS1","THBS2","LAMA2","LAMA4","MMP2","ARHGAP6","ARHGEF5","DLC1","TIAM1","TIAM2",
                                                                   "ACTA2","ACTN2","WASL","ADORA1","ADRA2C","CNR1","S1PR3","PTGER3"))

  #subset gene symbols in metadata
  metadata_subset <- subset(metadata_subset, rownames(metadata_subset)%in%de_data_subset_path$SYMBOL)
}

metadata_subset <- metadata_subset %>% mutate(subpathway = gsub("_"," ",subpathway))

metadata_subset <- metadata_subset %>% distinct(SYMBOL,.keep_all = TRUE)

rownames(metadata_subset) <- metadata_subset[,4]

metadata_subset <- metadata_subset[c("subpathway")]

####### create HeatmapDF object ########

de_data_subset_path <- de_data_subset_path%>%
  mutate(cell_line = gsub("_", " ", cell_line))

HeatmapDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = length(unique(de_data_subset_path$SYMBOL)), ncol = length(unique(de_data_subset_path$cell_line))))

de_data_subset_path <- de_data_subset_path%>%
  mutate(cell_line = gsub("_", " ", cell_line))

colnames(HeatmapDF) <- unique(de_data_subset_path$cell_line)

rownames(HeatmapDF) <- unique(de_data_subset_path$SYMBOL)


CL_list <- unique(de_data_subset_path$cell_line)

gene_list <- unique(de_data_subset_path$SYMBOL)

for(k in 1:length(CL_list)){

  de_data_cl_subset <- subset(de_data_subset_path, cell_line==CL_list[k])

  for(j in 1:length(gene_list)){


    de_data_cl_gene_subset <- subset(de_data_cl_subset, SYMBOL ==gene_list[j])

    if (nrow(de_data_cl_gene_subset) == 0){
      next} else

        HeatmapDF[j,k] <- de_data_cl_gene_subset$log2FoldChange
  }
}

HeatmapDF <- as.matrix(HeatmapDF)

HeatmapDF <-HeatmapDF[complete.cases(HeatmapDF), ]


# do same thing for p-values

PValueDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = length(unique(de_data_subset_path$SYMBOL)), ncol = length(unique(de_data_subset_path$cell_line))))

colnames(PValueDF) <-  unique(de_data_subset_path$cell_line)

rownames(PValueDF) <- unique(de_data_subset_path$SYMBOL)

CL_list <- unique(de_data_subset_path$cell_line)

gene_list <- unique(de_data_subset_path$SYMBOL)

for(k in 1:length(CL_list)){

  de_data_cl_subset <- subset(de_data_subset_path, cell_line==CL_list[k])

  for(j in 1:length(gene_list)){


    de_data_cl_gene_subset <- subset(de_data_cl_subset, SYMBOL ==gene_list[j])

    if (nrow(de_data_cl_gene_subset) == 0){
      next} else

        PValueDF[j,k] <- de_data_cl_gene_subset$significance
  }
}

PValueDF <- as.matrix(PValueDF)

PValueDF <-PValueDF[complete.cases(PValueDF), ]


########### Run different clustering methods & K number of clusters ###########

compound <- "All"
rowklist <- 2:10

UnsupervisedResultsDF <- data.frame(matrix(ncol=8,nrow=0))

t_HeatmapDF <- t(HeatmapDF)

# run multiple clustering methods and cluster sizes on genes

for(K in rowklist){
  # run K-means clustering
  rowInput <- runKMeans(HeatmapDF,K,compound)

  rowInput <- c("row",rowInput)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, rowInput)

  #run Hierarchical clustering
  rowInput <- runHierarchicalClustering(HeatmapDF,K,compound)

  rowInput <- c("row",rowInput)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, rowInput)
}

colnames(UnsupervisedResultsDF) <- c("type","Compound","Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score")


# run multiple clustering methods and cluster sizes on conditions tested

rowklist<-c(2:3)
for(K in rowklist){
  # run K-means clustering
  rowInput <- runKMeans(t_HeatmapDF,K,compound)

  rowInput <- c("col",rowInput)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, rowInput)

  #run Hierarchical clustering
  rowInput <- runHierarchicalClustering(t_HeatmapDF,K,compound)

  rowInput <- c("col",rowInput)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, rowInput)
}

colnames(UnsupervisedResultsDF) <- c("type","Compound","Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score")


########### Compare performance of the different clustering methods and sizes ###########

bestClustMethodDF <- data.frame(matrix(ncol=9,nrow=0))

ListOfUnsupervisedMethodComps <- runUnsupervisedAlgorithmComparison(UnsupervisedResultsDF)

bestClustMethodDF <- ListOfUnsupervisedMethodComps$bestClustMethod

# save plots for comparing all methods

SubDir <- sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/")

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Saving clustering method comparisons")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Saving clustering method comparisons")
}

SubDir <- sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/")

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Saving clustering method comparisons")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Saving clustering method comparisons")
}


# for row-specific clustering

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/silhouette2_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/dunn_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/connectivity_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/connectivity_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/ch_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering/wss_row.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodsrow$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


# for column-specific clustering

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/silhouette2_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/dunn_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/connectivity_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/connectivity_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/ch_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/column_clustering/wss_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


# setting specific color range for plots
scaling <- "scaled"
genetype <-gsub("_"," ",names(DEG_pathways_df_hgnc_list)[l])
compound <-"No"
columnk <- 0

# for case when columns should be visualized in specific order
col_order <- c("ARID1A KO", "ARID1B KO","ARID2 KO","BRD7 KO","PHF10 KO","PBRM1 KO","BRD9 KO","SMARCD1 KO","SMARCC1 KO","SMARCA4 KO")

####### get number of DEGs that are up or downregulated for each KO   #######
de_data3 <- DEG_pathways_df_hgnc_list_subset %>% distinct(SYMBOL, cell_line,log2FoldChange,.keep_all = TRUE)

de_data3 <- de_data3 %>% mutate(cell_line = gsub("_KO", "", cell_line))

de_data_subset_DEGs <- de_data3 %>% mutate(regulation = case_when(log2FoldChange > 0 & padj < 0.05 ~ "Upregulation",
                                                                 log2FoldChange < 0 &  padj < 0.05 ~ "Downregulation",
                                                                        TRUE ~ "not significantly DE"))


de_data_subset_DEGs <- de_data_subset_DEGs %>% group_by(cell_line, regulation) %>% dplyr::summarise(Number_of_DEGs = n())


de_data_subset_DEGs <- subset(de_data_subset_DEGs,regulation!="not significantly DE" )


de_data_subset_DEGs <- de_data_subset_DEGs %>% mutate(Number_of_DEGs = case_when(regulation== "Downregulation" ~ -Number_of_DEGs,
                                                               TRUE ~ Number_of_DEGs))

de_data_subset_DEGs <- de_data_subset_DEGs %>%  mutate(regulation = factor(
  regulation,
  levels=c("Upregulation", "Downregulation")))

  de_data_subset_DEGs <- de_data_subset_DEGs %>%
  mutate(cell_line = factor(cell_line, levels = c("ARID1A", "ARID1B", "ARID2","BRD7","PHF10","PBRM1","BRD9","SMARCD1","SMARCC1","SMARCA4")))

up_down_plot <- ggplot(de_data_subset_DEGs, aes(x=cell_line, y=Number_of_DEGs, fill=regulation)) +
  geom_bar(stat="identity", position="stack") +theme_classic() +
  labs(
    col='Cell line modification') +
  labs(x=paste("Knockout\n"),
       y=paste("Number of differentially\n expressed genes"))+ scale_y_continuous(labels = abs) + theme(legend.title = element_blank(),
                                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
theme(axis.text=element_text(size=63),
      legend.text = element_text(size = 54),
      legend.title = element_text(size = 54),
      axis.text.x = element_text(size = 63),
      axis.text.y = element_text(size = 63),
      axis.title=element_text(size=63,face="bold"),
      plot.margin = margin(1.25,0.90,0.90,0.90, "cm"),
      legend.position=c(.219,.236),
      legend.key.size=unit(3,"lines"),
      legend.background=element_blank())

NewFilename <- sprintf("/%s/%s/Genes/%s/%s_DEGs_up_down_plot.pdf",dataset_path, comparison,names(DEG_pathways_df_hgnc_list)[l],names(DEG_pathways_df_hgnc_list)[l])

ggsave(filename=NewFilename, plot=up_down_plot,width = 4.3, height = 4, dpi = 150, units = "in", device='pdf')

if(isTRUE(genetype=="Cell Cycle Checkpoints")==TRUE){
  minColorRange <- -1.5
  maxColorRange <- 1.5

  bestClustMethodDF[which(bestClustMethodDF$type=="row"),4] <- 3
  bestClustMethodDF[which(bestClustMethodDF$type=="col"),4] <- 0
  bestClustMethodDF[which(bestClustMethodDF$Method=="Hierarchical"),3] <- "K-means"
} else if(isTRUE(genetype=="DNA Damage Repair")==TRUE){
  minColorRange <- -2
  maxColorRange <- 2

  bestClustMethodDF[which(bestClustMethodDF$type=="row"),4] <- 2
  bestClustMethodDF[which(bestClustMethodDF$type=="col"),4] <- 0
  bestClustMethodDF[which(bestClustMethodDF$Method=="Hierarchical"),3] <- "K-means"
} else{
  minColorRange <- -2
  maxColorRange <- 2

  bestClustMethodDF[which(bestClustMethodDF$type=="row"),4] <- 2 #3 when filter_DEG_subset == FALSE
  bestClustMethodDF[which(bestClustMethodDF$type=="col"),4] <- 0
  bestClustMethodDF[which(bestClustMethodDF$Method=="Hierarchical"),3] <- "K-means"
}

# apply desired order to the input data
HeatmapDF <- HeatmapDF[ , col_order]

# apply desired order to the input data
PValueDF <- PValueDF[ , col_order]

print(paste("Number of DE growth sign. genes: ",nrow(HeatmapDF)))

plotlist<- getClusteredHeatmap(HeatmapDF,bestClustMethodDF,columnk,scaling,minColorRange,maxColorRange,metadata=metadata_subset,PValueDF = PValueDF)

#creating necessary directories
SubDir <- sprintf("./Unsupervised_Clustering/Heatmap/")

if(isTRUE(file.exists(SubDir))==TRUE){
  print(paste("Plotting",bestClustMethodDF$Method,"clustered heatmap with" ,bestClustMethodDF$Cluster_Num, "clusters"))
} else {
  dir.create(SubDir, recursive = TRUE)
  print(paste("Plotting",bestClustMethodDF$Method,"clustered heatmap with" ,bestClustMethodDF$Cluster_Num, "clusters"))
}

# save clustered heatmap
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_%s_heatmap_%s.pdf",bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1], bestClustMethodDF$Cluster_Num[2],compound)

save_pheatmap_pdf(plotlist$Heatmap, NewFilename)

#save average silhouette plot - row clustering
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_avg_sil_row_%s.pdf", bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1],compound)

ggsave(filename=NewFilename, plot=plotlist$AvgSilrow,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


#save average silhouette plot - col clustering
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_avg_sil_col_%s.pdf", bestClustMethodDF$Method[2],
                      bestClustMethodDF$Cluster_Num[2],compound)

ggsave(filename=NewFilename, plot=plotlist$AvgSilcol,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

# save scatter plot
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_scatter_plot_%s.pdf", bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1],compound)

ggsave(filename=NewFilename, plot=plotlist$ClusterViz,width = 15, height = 10, dpi = 150, units = "in", device='pdf')
}
