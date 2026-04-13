############################################
# HAP1 KO Transcriptome Deviation Analysis (PES Analysis)
# Author: Eleni Aretaki
# Date: 2026-04-01
#
# Purpose:
# 1. Filter protein-coding genes
# 2. Compute KO vs all preferential expression scores (PES)
# 3. Cluster genes using k-means
# 4. Generate Heatmap visualization
# 5. Perform GO enrichment per cluster
############################################


###### Setup environment ######

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ComplexHeatmap"
), ask = FALSE, update = FALSE)

library(edgeR)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

library(ggplot2)


###### Load and prepare expression matrix ######

rownames(protein_coding_gene_counts_hap1) <- protein_coding_gene_counts_hap1$gene_name
protein_coding_gene_counts_hap1$gene_name <- NULL

expr_matrix <- protein_coding_gene_counts_hap1


###### Keep protein-coding genes only ######

protein_coding_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

expr_matrix <- expr_matrix[
  rownames(expr_matrix) %in% protein_coding_symbols,
]


###### Define KO panel ######

selected_KOs <- c(
  "ARID1A", "ARID1B", "SMARCA4", "ARID2", "SMARCD1",
  "SMARCC1", "PHF10", "PBRM1", "BRD7", "BRD9", "WT"
)


###### Clone exclusion patterns ######

exclude_patterns <- c(
  "ARID1A_r[0-9]+_618_10",
  "SMARCA4_r[0-9]+_2878_4",
  "WT_r[0-9]+_GFP"
)


###### Helper: select KO columns ######

get_ko_columns <- function(ko_gene, expr_matrix, exclude_patterns) {
  
  ko_cols <- grep(ko_gene, colnames(expr_matrix), value = TRUE)
  
  bad_cols <- grep(
    paste(exclude_patterns, collapse = "|"),
    ko_cols,
    value = TRUE
  )
  
  setdiff(ko_cols, bad_cols)
}


###### Build filtered expression matrix ######

selected_columns <- unlist(
  lapply(
    selected_KOs,
    get_ko_columns,
    expr_matrix = expr_matrix,
    exclude_patterns = exclude_patterns
  )
)

expr_matrix_filtered <- expr_matrix[, selected_columns]


###### Convert counts to log-CPM ######

dge <- DGEList(counts = expr_matrix_filtered)
cpm_matrix <- cpm(dge, log = TRUE, prior.count = 1)


###### Compute KO deviation (PES) scores ######

compute_ko_deviation <- function(ko_gene, expr_matrix, exclude_patterns) {
  
  cat("Processing KO:", ko_gene, "\n")
  
  ko_columns <- get_ko_columns(ko_gene, expr_matrix, exclude_patterns)
  
  ko_mean <- rowMeans(expr_matrix[, ko_columns, drop = FALSE])
  
  other_columns <- setdiff(colnames(expr_matrix), ko_columns)
  
  other_mean <- rowMeans(expr_matrix[, other_columns, drop = FALSE])
  
  iqr_values <- apply(
    expr_matrix[, other_columns, drop = FALSE],
    1,
    function(x) IQR(x, na.rm = TRUE)
  )
  
  (ko_mean - other_mean) / iqr_values
}


###### Calculate deviations (PES scores) ######

all_deviations <- lapply(
  selected_KOs,
  compute_ko_deviation,
  expr_matrix = cpm_matrix,
  exclude_patterns = exclude_patterns
)

deviations_matrix <- do.call(cbind, all_deviations)
colnames(deviations_matrix) <- selected_KOs

write.csv(deviations_matrix, "1_vs_all_deviations_matrix.csv")


###### Filter significant genes for clustering ######

deviation_filter <- apply(
  deviations_matrix,
  1,
  function(x) any(x > 1 | x < -1)
)

deviations_filtered <- deviations_matrix[deviation_filter, ]
data_for_clustering <- na.omit(as.matrix(deviations_filtered))


###### K-means clustering ######

wss <- sapply(1:15, function(k){
  kmeans(data_for_clustering, centers = k, nstart = 10)$tot.withinss
})

plot(
  1:15, wss,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Total within-cluster sum of squares",
  main = "Elbow Plot for K-means Clustering"
)

k <- 9
set.seed(123)

kmeans_result <- kmeans(data_for_clustering, centers = k)


###### Order rows by cluster ######

row_annot <- data.frame(Cluster = factor(kmeans_result$cluster))
rownames(row_annot) <- rownames(data_for_clustering)

cluster_order <- order(kmeans_result$cluster)

deviations_ordered <- data_for_clustering[cluster_order, ]
row_annot_ordered <- row_annot[cluster_order, , drop = FALSE]


###### Column order ######

cell_line_order <- c(
  "WT", "ARID1A","ARID1B",
  "ARID2", "BRD7", "PHF10","PBRM1",
  "BRD9", "SMARCD1",
  "SMARCC1","SMARCA4"
)

deviations_ordered <- deviations_ordered[, cell_line_order]


###### Cluster colors ######

cluster_levels <- levels(row_annot_ordered$Cluster)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Paired"))(length(cluster_levels)),
  cluster_levels
)


###### ComplexHeatmap palette ######

my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

my_palette_CH <- colorRamp2(
  seq(-5, 5, length.out = 100),
  my_palette
)


###### Row annotation ######

row_annot_df <- data.frame(Cluster = row_annot_ordered$Cluster)

row_ha <- rowAnnotation(
  Cluster = row_annot_df$Cluster,
  col = list(
    Cluster = cluster_colors
  ),
  simple_anno_size = unit(2, "mm"),
  show_annotation_name = FALSE,
  show_legend = FALSE
)


###### ComplexHeatmap plot ######

pdf("HAP1_TiP_heatmap.pdf",
    width = 3,
    height = 4
)

Heatmap(
  deviations_ordered,
  name = "Deviation",
  col = my_palette_CH,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  left_annotation = row_ha,
  width = unit(5, "cm"),
  heatmap_legend_param = list(
    title = "",
    legend_height = unit(2, "cm"),
    grid_width = unit(2, "mm"),
    title_gp = gpar(fontsize = 6),
    direction = "vertical",
    labels_gp = gpar(fontsize = 6)
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

dev.off()


write.csv(
  deviations_ordered,
  "1_vs_all_deviations_filtered_matrix_9.csv"
)


###### Export clusters ######

gene_to_cluster <- kmeans_result$cluster
gene_clusters <- split(names(gene_to_cluster), gene_to_cluster)

for (i in seq_along(gene_clusters)) {
  
  genes <- gene_clusters[[i]]
  
  cluster_data <- deviations_ordered[genes, , drop = FALSE]
  
  df <- data.frame(
    Gene = rownames(cluster_data),
    Cluster = i,
    cluster_data,
    row.names = NULL
  )
  
  write.csv(
    df,
    file = paste0("Cluster_", i, "_genes_with_values.csv"),
    row.names = TRUE
  )
}


###### GO enrichment ######

enrichment_results <- list()

for (i in seq_along(gene_clusters)) {
  
  cat("\nRunning enrichment for Cluster", i, "\n")
  
  genes <- gene_clusters[[i]]
  
  entrez_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  ego <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  enrichment_results[[i]] <- ego
}


###### GO plots ######

for (i in seq_along(enrichment_results)) {
  
  ego <- enrichment_results[[i]]
  
  if (!is.null(ego) && "result" %in% slotNames(ego)) {
    
    res_df <- ego@result
    
    if (nrow(res_df) > 0 && any(res_df$p.adjust < 0.05)) {
      
      p <- barplot(
        ego,
        showCategory = 10,
        title = paste("Cluster", i, "GO Enrichment")
      )
      
      ggsave(
        paste0("Cluster_", i, "_GO_Enrichment.png"),
        plot = p,
        width = 8,
        height = 6
      )
    }
  }
}


###### Save GO tables ######

for (i in seq_along(enrichment_results)) {
  
  ego <- enrichment_results[[i]]
  
  if (!is.null(ego) && "result" %in% slotNames(ego)) {
    
    res_df <- ego@result
    
    if (nrow(res_df) > 0) {
      
      write.csv(
        res_df,
        file = paste0("Cluster_", i, "_GO_Enrichment.csv"),
        row.names = FALSE
      )
    }
  }
}

