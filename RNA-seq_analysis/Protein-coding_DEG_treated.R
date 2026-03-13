#!/usr/bin/env RScript

#command-line executable:  chmod +x Protein-coding_DEG_treated.R
# run as: ./Protein-coding_DEG_treated.R

# Script for analyzing protein-coding genes for the treated RNA-seq data from Spang et al., 2026 and visualizing in a clustered heatmap

################################################## Loading Necessary Libraries and Source Code #########################################################################
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
library(enrichplot)
library(ggtree)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#RNA seq dataset:
dataset_path<-" "

#DE analysis comparison performed
comparison <- "Compound_vs_DMSO_KO_vs_WT_DEGs"
filename_pattern <-"KO.*\\.txt"


SubDir <- sprintf("/%s/%s/Genes/Protein-coding/", dataset_path,comparison)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Exists!")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Created!")
}

setwd(sprintf("/%s/%s/Genes/Protein-coding/",dataset_path,comparison))


  ########## upload all DE data all KOs & all Compounds & merge all data ##########

  de_data <- list.files(path = sprintf("/%s/%s/",dataset_path,comparison),
                        pattern = filename_pattern, full.names = TRUE) %>%
    lapply(read.delim) %>%
    bind_rows

  # remove NA
  de_data <- de_data[!is.na(de_data$padj),]


  # rename column containing ensembl gene IDs and remove decimal values to match gene IDs in ensembl for downstreatm analysis
  de_data <- de_data %>% dplyr::rename(ensembl_gene_id = X)
  de_data <- de_data %>% mutate(ensembl_gene_id = gsub("\\..*","",ensembl_gene_id))

  # filter for protein-coding genes

  metadata<- clusterProfiler::bitr(geneID = de_data$ensembl_gene_id, fromType = "ENSEMBL", toType = c("ENSEMBL","GENETYPE","SYMBOL"), OrgDb = "org.Hs.eg.db")


  metadata <- subset(metadata, GENETYPE =="protein-coding")

  # left join in order to get gene names
  metadata <- metadata %>% dplyr::rename(ensembl_gene_id = ENSEMBL, hgnc_symbol = SYMBOL)

  # only keep protein-coding genes
  de_data <- subset(de_data, ensembl_gene_id %in% metadata$ENSEMBL)

  ##########  filter for genes that are significant in at least one BAF KO for each treatment ##########

  de_data2 <- de_data %>%
    group_by(ensembl_gene_id, treatment) %>%
    dplyr::filter(padj < 0.001)

  ##########  Merge Compound & Cell Line Names together ##########

  de_data <- de_data %>% mutate(condition = paste(treatment,cell_line,sep="_"))

  de_data2 <- de_data2 %>% mutate(condition = paste(treatment,cell_line,sep="_"))


  de_data <- de_data %>% dplyr::rename(ENSEMBL = ensembl_gene_id)

  de_data2 <- de_data2 %>% dplyr::rename(ENSEMBL = ensembl_gene_id)


  ####### create HeatmapDF object ########

  de_data_subset <- subset(de_data, ENSEMBL %in% de_data2$ENSEMBL)

  de_data_subset <- de_data_subset %>% distinct(ENSEMBL,log2FoldChange, condition,.keep_all = TRUE)

  de_data_subset <- de_data_subset%>%
    mutate(condition = gsub("_", " ", condition))

  HeatmapDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = length(unique(de_data_subset$ENSEMBL)), ncol = length(unique(de_data_subset$condition))))

  de_data_subset <- de_data_subset%>%
    mutate(condition = gsub("_", " ", condition))

  colnames(HeatmapDF) <- unique(de_data_subset$condition)

  rownames(HeatmapDF) <- unique(de_data_subset$ENSEMBL)


  CL_list <- unique(de_data_subset$condition)

  gene_list <- unique(de_data_subset$ENSEMBL)

  for(k in 1:length(CL_list)){

    de_data_cl_subset <- subset(de_data_subset, condition==CL_list[k])

    for(j in 1:length(gene_list)){


      de_data_cl_gene_subset <- subset(de_data_cl_subset, ENSEMBL ==gene_list[j])

      if (nrow(de_data_cl_gene_subset) == 0){
        next} else

          HeatmapDF[j,k] <- de_data_cl_gene_subset$log2FoldChange
    }
  }


HeatmapDF <- as.matrix(HeatmapDF)

HeatmapDF <-HeatmapDF[complete.cases(HeatmapDF), ]

t_HeatmapDF <- t(HeatmapDF)


########### Run different clustering methods & K number of clusters ###########

compound <- "all"
rowklist <- 2:20

UnsupervisedResultsDF <- data.frame(matrix(ncol=8,nrow=0))

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

rowklist <- 2:2

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

SubDir <- sprintf("./Unsupervised_Clustering/Algorithm_Comparison/row_clustering")

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Saving clustering method comparisons")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Saving clustering method comparisons")
}

SubDir <- sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering")

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Saving clustering method comparisons")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Saving clustering method comparisons")
}

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


##### saving column clustering metrics results

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/silhouette2_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/dunn_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/connectivity_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/connectivity_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/ch_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

NewFilename = sprintf("./Unsupervised_Clustering/Algorithm_Comparison/col_clustering/wss_col.pdf")

ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethodscol$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


# setting specific color range for plots
scaling <- "scaled"
genetype <-"Protein-coding"
compound <-"All"
minColorRange <- -1
maxColorRange <- 1
columnk <- 0

# prefer columns to be visualized in specific order
col_order <- c("Ada ARID1A KO", "Ada ARID1B KO", "Ada ARID2 KO","Ada BRD9 KO","Ada SMARCA4 KO","CPT ARID1A KO", "CPT ARID1B KO", "CPT ARID2 KO","CPT BRD9 KO","CPT SMARCA4 KO",
               "HU ARID1A KO", "HU ARID1B KO", "HU ARID2 KO","HU SMARCA4 KO","MLN ARID1A KO", "MLN ARID1B KO", "MLN ARID2 KO","MLN SMARCA4 KO",
               "MMS ARID1A KO", "MMS ARID1B KO", "MMS ARID2 KO","MMS BRD9 KO","MMS SMARCA4 KO","Pbc ARID1A KO", "Pbc ARID1B KO", "Pbc ARID2 KO","Pbc BRD9 KO","Pbc SMARCA4 KO")

# updating number of clusters
bestClustMethodDF[which(bestClustMethodDF$type=="row"),4] <- 5
bestClustMethodDF[which(bestClustMethodDF$type=="col"),4] <- 0
bestClustMethodDF[which(bestClustMethodDF$Method=="Hierarchical"),3] <- "K-means"

# apply desired order to the input data
HeatmapDF <- HeatmapDF[ , col_order]

print(paste("Number of DE growth sign. genes: ",nrow(HeatmapDF)))

plotlist<- getClusteredHeatmap(HeatmapDF,bestClustMethodDF,columnk,scaling,minColorRange,maxColorRange)

#creating necessary directories
SubDir <- sprintf("./Unsupervised_Clustering/Heatmap/")

if(isTRUE(file.exists(SubDir))==TRUE){
  print(paste("Plotting",bestClustMethodDF$Method,"clustered heatmap with" ,bestClustMethodDF$Cluster_Num, "clusters"))
} else {
  dir.create(SubDir, recursive = TRUE)
  print(paste("Plotting",bestClustMethodDF$Method,"clustered heatmap with" ,bestClustMethodDF$Cluster_Num, "clusters"))
}

# save clustered heatmap
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_%s_%s_%s_heatmap.pdf",bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1], bestClustMethodDF$Cluster_Num[2],genetype,comparison)

save_pheatmap_pdf(plotlist$Heatmap, NewFilename)

#save average silhouette plot - row clustering
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_avg_sil_row_%s_%s.pdf", bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1],genetype,comparison)

ggsave(filename=NewFilename, plot=plotlist$AvgSilrow,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

#save average silhouette plot - col clustering
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_avg_sil_col_%s_%s.pdf", bestClustMethodDF$Method[2],
                      bestClustMethodDF$Cluster_Num[2],genetype,comparison)

ggsave(filename=NewFilename, plot=plotlist$AvgSilcol,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

# save scatter plot
NewFilename = sprintf("./Unsupervised_Clustering/Heatmap/%s_%s_scatter_plot_%s_%s.pdf", bestClustMethodDF$Method[1],
                      bestClustMethodDF$Cluster_Num[1],genetype,comparison)

ggsave(filename=NewFilename, plot=plotlist$ClusterViz,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


###### Annotating k-means clusters & running GO enrichment #######

SubDir <- sprintf("/%s/%s/Genes/Protein-coding/Unsupervised_Clustering/Heatmap/GO_enrichment/",dataset_path,comparison)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("Exists!")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Created!")
}


setwd(sprintf("/%s/%s/Genes/Protein-coding/Unsupervised_Clustering/Heatmap/GO_enrichment/",dataset_path,comparison))

# for(k in k_list){

set.seed(123)
k2m_data <- kmeans(HeatmapDF, as.numeric(bestClustMethodDF$Cluster_Num[1]), nstart=25)

# Row annotation with cluster assignment
row_annot <- data.frame(Cluster = factor(k2m_data$cluster))
rownames(row_annot) <- rownames(HeatmapDF)

# Order genes by their assigned cluster
cluster_order <- order(k2m_data$cluster)
deviations_ordered <- HeatmapDF[cluster_order, ]

# Also reorder the annotation to match
row_annot_ordered <- row_annot[cluster_order, , drop = FALSE]

# Assume `row_annot_ordered` has a column called "Cluster"
levels <- unique(row_annot_ordered$Cluster)

# Create a list of genes per cluster
gene_clusters <- split(rownames(HeatmapDF), k2m_data$cluster)

for (i in seq_along(gene_clusters)) {
  write.table(gene_clusters[[i]],
              file = paste0("Cluster_", i, "_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)}

# Extract gene-to-cluster mapping
gene_to_cluster <- k2m_data$cluster
gene_clusters <- split(names(gene_to_cluster), gene_to_cluster)

# Loop and write full table per cluster
for (i in seq_along(gene_clusters)) {
  genes <- gene_clusters[[i]]

  # Subset the heatmap matrix for this cluster
  cluster_data <- deviations_ordered[genes, , drop = FALSE]

  # Create a data frame with gene names and values
  df <- data.frame(
    Gene = rownames(cluster_data),
    Cluster = i,
    cluster_data,
    row.names = NULL)

  # Write to file
  write.csv(df,
            file = paste0("Cluster_", i, "_genes_with_values.csv"), row.names = TRUE)}

###### Clusters GO Enrichment ######

# Create a list to store results
enrichment_results <- list()

# Loop through each cluster
for (i in seq_along(gene_clusters)) {
  cat("\nRunning enrichment for Cluster", i, "\n")

  # Get genes for the current cluster
  genes <- gene_clusters[[i]]

  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Run GO enrichment for Biological Process (BP)
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

  # Store the results for each cluster
  enrichment_results[[i]] <- ego}

# Plotting only top 10 terms that are significant (padjust < 0.05) & saving these in a list
barplots <-c()
geneplots <- c()
treeplots <- c()
mapplots <- c()
n<-0
m<-0
s<-0

# Loop over enrichment results and plot if results exist
for (i in seq_along(enrichment_results)) {

  # Set a number of top terms to display (e.g., top 10)
  top_terms <- 10

  ego <- enrichment_results[[i]]

  # Check if the enrichment result exists and contains significant rows
  if (!is.null(ego) && "result" %in% slotNames(ego)) {
    res_df <- ego@result
    if (nrow(res_df) > 0 && any(res_df$p.adjust < 0.05)) {
      cat("Plotting Cluster", i, "\n")

      # Create barplot
      p <- barplot(ego, showCategory = top_terms,
                   title = paste("Cluster", i, "GO Enrichment"))

      p2 <- heatplot(ego, showCategory=top_terms)+ labs(title=paste("Cluster", i, "GO Enrichment")) + theme(text = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1))

      if(isTRUE(nrow(ego)>5)==TRUE){
        ego2 <- pairwise_termsim(ego)
        p3 <- treeplot(ego2, cladelab_offset=8, tiplab_offset=.3, fontsize_cladelab =4) +
          hexpand(.2)

        p4 <- emapplot(ego2)}

      #create list of plots
      if(isTRUE(!(is.null(p)))==TRUE){
        n<- 1+n
        barplots[[n]] <- p
        geneplots[[n]] <- p2
        if(isTRUE(nrow(ego)>2)==TRUE){
          if(isTRUE(!(is.null(p3)))==TRUE){
            m<- 1+m
            treeplots[[m]] <- p3}
          if(isTRUE(!(is.null(p4)))==TRUE){
            s<- 1+s
            mapplots[[s]] <- p4}}
      }

      # Save the plot
      ggsave(paste0("Cluster_", i, "_GO_Enrichment.png"),
             plot = p, width = 8, height = 10)

      # Save the plot
      ggsave(paste0("Cluster_", i, "_GO_Enrichment_genes.png"),
             plot = p2, width = 10, height = 6)

      if(isTRUE(nrow(ego)>2)==TRUE){
        # Save the plot
        ggsave(paste0("Cluster_", i, "_GO_Enrichment_treeplot.png"),
               plot = p3, width = 10, height = 6)

        # Save the plot
        ggsave(paste0("Cluster_", i, "_GO_Enrichment_mapplot.png"),
               plot = p4, width = 10, height = 6)}

    } else {
      cat("No significant GO terms for Cluster", i, "- skipping plot\n")
    }
  } else {
    cat("No enrichment result for Cluster", i, "- skipping plot\n")
  }
}

# Save the list of plots
if(isTRUE(length(barplots)>1)==TRUE){
  ggsave(paste0("All_Clusters_GO_Enrichment.png"),
         plot = gridExtra::grid.arrange(grobs = barplots),width = 10, height = 14)
  ggsave(paste0("All_Clusters_GO_Enrichment_genes.png"),
         plot = gridExtra::grid.arrange(grobs = geneplots),width = 10, height = 12)
  ggsave(paste0("All_Clusters_GO_Enrichment_treeplot.png"),
         plot = gridExtra::grid.arrange(grobs = treeplots))
  ggsave(paste0("All_Clusters_GO_Enrichment_mapplot.png"),
         plot = gridExtra::grid.arrange(grobs = mapplots))

}

# Save individual CSV files for each cluster
for (i in seq_along(enrichment_results)) {
  ego <- enrichment_results[[i]]

  # Check if enrichment result exists and contains significant terms
  if (!is.null(ego) && "result" %in% slotNames(ego)) {
    res_df <- ego@result
    if (nrow(res_df) > 0) {
      # Save to CSV
      write.csv(res_df, file = paste0("Cluster_", i, "_GO_Enrichment.csv"), row.names = FALSE)
      cat("Saved GO enrichment results for Cluster", i, "\n")
    } else {
      cat("No significant results for Cluster", i, " - skipping CSV\n")
    }
  } else {
    cat("No enrichment result for Cluster", i, " - skipping CSV\n")
  }
}


