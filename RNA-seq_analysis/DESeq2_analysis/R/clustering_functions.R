#' Run Hierarchical Clustering

#' This function runs hierarchical clustering on the log2 fold-change in gene expression and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the log2 fold-change in gene expression computed from DESeq2.
#' @param K Numerical value representing the number of clusters to test.
#' @param compound Variable defining the treatment type of the samples

#' @return
#'\item{A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score).}
#' @author Caroline Barry
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export


runHierarchicalClustering <- function(InputDF2,K,compound){


  ########### Hierachical Clustering - Average Linkage Method ############

  dist_man <- dist(InputDF2, method="manhattan")

  hc_m2 <- hclust(d=dist_man, method="average")

  coph_m2 <- cophenetic(hc_m2)
  cor(dist_man,coph_m2)

  groupward6 <- cutree(hc_m2, k = K)

  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Silhouette Analysis ########################

  average <- eclust(InputDF2, "hclust", k = K, hc_metric = "manhattan",hc_method = "average", graph = F)


  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Dunn Index ########################

  averagea <- cluster.stats(dist(InputDF2), average$cluster)


  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Connectivity ########################

  ConnectivityScore <-connectivity(distance = NULL, average$cluster, Data = InputDF2, neighbSize = 20,
                                   method = "euclidean")


  ########## save performance metrics to dataframe ########

  return(c(compound,"Hierarchical", K, averagea$avg.silwidth, averagea$dunn,averagea$within.cluster.ss,averagea$ch, ConnectivityScore))

}

#########################################################################################################################################################################################################################################################################

#' Run K-Means Clustering

#' This function runs k-means clustering on the log2 fold-change in gene expression and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the log2 fold-change in gene expression computed from DESeq2.
#' @param K Numerical value representing the number of clusters to test.
#' @param compound Variable defining the treatment type of the samples

#' @return
#'\item{A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score).}
#' @author Caroline Barry
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export
#########################################################################################################################################################################################################################################################################
runKMeans <- function(InputDF2,K,compound){

  k2m_data <- kmeans(InputDF2, K, nstart=25)

  ######################## cluster validation for K-means - Silhouette Analysis ########################

  k2m_data <- factoextra::eclust(InputDF2, "kmeans", k = K, nstart = 25, graph = F)


  ######################## cluster validation for K-means - Dunn Index ########################
  km_stats <- cluster.stats(dist(InputDF2), k2m_data$cluster)

  ######################## cluster validation for K-means - CH Index ########################
  # km_stats$ch

  ######################## cluster validation for K-means - Connectivity ########################
  ConnectivityScore <- connectivity(distance = NULL, k2m_data$cluster, Data = InputDF2, neighbSize = 20,
                                    method = "euclidean")

  ########## save performance metrics to dataframe ########

  return(c(compound,"K-means", K, km_stats$avg.silwidth, km_stats$dunn,km_stats$within.cluster.ss,km_stats$ch, ConnectivityScore))

}

#########################################################################################################################################################################################################################################################################

#' Comparing the Performance of Clustering Algorithms

#' This function plots the performance of various clustering algorithms on the log2 fold-change in gene expression and returns them in a list. The best scoring clustering method and number of clusters is also returned.
#' @param UnsupervisedResultsDF A dataframe containing the performance of k-means  hierarchical clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score).

#' @return
#'\item{A list of plots comparing the performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and a dataframe containing the best clustering method and number of clusters across the genes (rows) and conditions (columns).}
#' @author Caroline Barry
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @export

runUnsupervisedAlgorithmComparison <- function(UnsupervisedResultsDF){

  ################ Determining the Best Clustering Algorithm Based on Clustering Analysis Results #####################
  UnsupervisedResultsDF%<>% mutate(Cluster_Size= str_pad(Cluster_Num, width = 2,side ="left",'0'))

  UnsupervisedResultsDF%<>% mutate(Method_Clust = paste(Cluster_Size,Method, sep="_"))

  UnsupervisedResultsDF_col <- subset(UnsupervisedResultsDF, type =="col")

  UnsupervisedResultsDF <- subset(UnsupervisedResultsDF, type =="row")

  # plot the average silhouette scores for all clustering algorithms
  AllAvgSilPlotObject<-ggplot(UnsupervisedResultsDF, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score") +
    theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),legend.title = element_text(size=12),
          plot.title=element_text(size=25),plot.subtitle=element_text(size=18))


  # plot the Dunn Indexes for all clustering algorithms
  AllDunnPlotObject<-ggplot(UnsupervisedResultsDF, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  # plot the connectivity scores for all clustering algorithms
  AllConnPlotObject<-ggplot(UnsupervisedResultsDF, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the WSS for all clustering algorithms
  AllWSSPlotObject<-ggplot(UnsupervisedResultsDF, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the CH index for all clustering algorithms
  AllCHPlotObject<-ggplot(UnsupervisedResultsDF, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  #now doing same for column clustering metrics

  # plot the average silhouette scores for all clustering algorithms
  AllAvgSilPlotObjectcol<-ggplot(UnsupervisedResultsDF_col, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score") +
    theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),legend.title = element_text(size=12),
          plot.title=element_text(size=25),plot.subtitle=element_text(size=18))

  # plot the Dunn Indexes for all clustering algorithms
  AllDunnPlotObjectcol<-ggplot(UnsupervisedResultsDF_col, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  # plot the connectivity scores for all clustering algorithms
  AllConnPlotObjectcol<-ggplot(UnsupervisedResultsDF_col, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the WSS for all clustering algorithms
  AllWSSPlotObjectcol<-ggplot(UnsupervisedResultsDF_col, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the CH index for all clustering algorithms
  AllCHPlotObjectcol<-ggplot(UnsupervisedResultsDF_col, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  bestClustMethod <- UnsupervisedResultsDF %>%
    arrange(desc(as.numeric(as.character(CH_Index))))

  bestClustMethod <- bestClustMethod[1,]

  if(isTRUE(!(is.null(dim(UnsupervisedResultsDF_col))))==TRUE){
    bestClustMethod_col <- UnsupervisedResultsDF_col %>%
      arrange(desc(as.numeric(as.character(CH_Index))))

    bestClustMethod_col <- bestClustMethod_col[1,]

    bestClustMethod <- rbind(bestClustMethod,bestClustMethod_col)
  } else{ print("No clustering of columns")}

  return(list(AllMethodsrow = list(AvgSil = AllAvgSilPlotObject, Dunn = AllDunnPlotObject, Connectivity = AllConnPlotObject,
                                   WSS = AllWSSPlotObject, CH = AllCHPlotObject),
              AllMethodscol = list(AvgSil = AllAvgSilPlotObjectcol, Dunn = AllDunnPlotObjectcol, Connectivity = AllConnPlotObjectcol,
                                   WSS = AllWSSPlotObjectcol, CH = AllCHPlotObjectcol),
              bestClustMethod = bestClustMethod))

}


#########################################################################################################################################################################################################################################################################
#' Plot Clustered Heatmap of Differentially Expressed Genes

#' This function plots a clustered heatmap of the Differentially Expressed Genes (DEGs). These are clustered using either K-means or Hierarchical clustering based on the user-defined method and number of clusters. The clustered heatmap, average sihouette
#' values plotted for column and row clustering, and PCA plot of the conditions are returned.
#' @param HeatmapDF Matrix containing log2fold-change gene expression values, in format of rows are genes and columns are individual conditions.
#' @param bestClustMethod Row of values containing clustering method and number of clusters.
#' @param columnk Numerical variable represneting the user-specified k number of clusters for the columns of the heatmap.
#' @param scaling Character variable that is user-specified where min and max values for color scale are specified ("scaled") or are simply the max and min values of the dataset ("unscaled"; default).
#' @param minColorRange Numerical variable represneting the minimum value of the color scale. Default is set to NULL.
#' @param maxColorRange Numerical variable represneting the maximum value of the color scale. Default is set to NULL.
#' @param PValueDF Martrix containing significance asterisks (representing the adjusted p-values from DESeq2) for each gene and condition. Default is set to NULL.
#' @param metadata Dataframe containing pathway and subpathway information for the genes. These labels will be used for annotating the heatmap. Default is set to NULL.
#' @author Caroline Barry
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export


getClusteredHeatmap <- function(HeatmapDF,bestClustMethod,columnk,scaling="unscaled",minColorRange=NULL,maxColorRange=NULL,PValueDF=NULL,metadata=NULL){

  bestClustMethodcol <- subset(bestClustMethod, type =="col")

  bestClustMethod <- subset(bestClustMethod, type =="row")

  # in case PValueDF == NULL, create empty matrix
  if(isTRUE(is.matrix(PValueDF))==FALSE){

  PValueDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = nrow(HeatmapDF), ncol = ncol(HeatmapDF)))

  colnames(PValueDF) <- colnames(HeatmapDF)

  rownames(PValueDF) <- rownames(HeatmapDF)

  CL_list <- colnames(PValueDF)

  gene_list <- rownames(PValueDF)

  for(k in 1:length(CL_list)){

    for(j in 1:length(gene_list)){

      PValueDF[j,k] <- " "
    }
  }

  PValueDF <- as.matrix(PValueDF)
  } else


  ###### if best clustering method is K-means ########
  if(isTRUE(bestClustMethod$Method=="K-means")==TRUE){

    set.seed(123)
    k2m_data <- kmeans(HeatmapDF, as.numeric(as.character(bestClustMethod$Cluster_Num)), nstart=25)

    # Row annotation with cluster assignment
    row_annot <- data.frame(Cluster = factor(k2m_data$cluster))
    rownames(row_annot) <- rownames(HeatmapDF)

    # Order genes by their assigned cluster
    cluster_order <- order(k2m_data$cluster)
    deviations_ordered <- HeatmapDF[cluster_order, ]

    # Also reorder the annotation to match
    row_annot_ordered <- row_annot[cluster_order, , drop = FALSE]


    if (isTRUE(genetype!="Protein-coding"&genetype!="BAF subunit")==TRUE){

      row_annot_ordered <- merge(row_annot_ordered,metadata, by=0, all=TRUE)

      rownames(row_annot_ordered)=row_annot_ordered$Row.names
      row_annot_ordered<-row_annot_ordered[,-1]


      row_annot_ordered <- row_annot_ordered %>% dplyr::rename(row_cluster = Cluster)

      row_annot_ordered <-row_annot_ordered[complete.cases(row_annot_ordered), ]


      # Assume `row_annot_ordered` has a column called "row_cluster"
      levels <- unique(row_annot_ordered$row_cluster)

      # Generate a palette (e.g., Set1)

      colors1 <- colorRampPalette(brewer.pal(9, "Set1"))(as.numeric(as.character(bestClustMethod$Cluster_Num)))

      colors <- setNames(colors1, levels)

      # Assume `row_annot_ordered` has a column called "subpathway"
      subpathway_levels <- unique(row_annot_ordered$subpathway)

      # Generate a palette (e.g., Set1)

      colorssp <- colorRampPalette(brewer.pal(9, "Set1"))(length(subpathway_levels))

      subpath_colors <- setNames(colorssp, subpathway_levels)

      # Pass to pheatmap
      annotation_colors <- list(
        row_cluster = colors,
        subpathway = subpath_colors)


    } else{row_annot_ordered <- row_annot_ordered %>% dplyr::rename(row_cluster = Cluster)

    # Assume `row_annot_ordered` has a column called "Cluster"
    levels <- unique(row_annot_ordered$row_cluster)

    # Generate a palette (e.g., Set1)
    colors1 <- colorRampPalette(brewer.pal(9, "Set1"))(as.numeric(as.character(bestClustMethod$Cluster_Num)))

    colors <- setNames(colors1, levels)

    # Pass to pheatmap
    annotation_colors <- list(
      row_cluster = colors)}

    #order Pvalue DF based on clustering order
    pvalues_ordered <- PValueDF[cluster_order,]

    # Define color palette (e.g., 100 steps from blue to white to red)
    col_length <- 11 #5 #7
    my_palette <- colorRampPalette(c("blue", "white", "red"))(col_length)

    if(isTRUE(columnk!=0)==TRUE){

      # if(isTRUE(nrow(bestClustMethodcol)!=0)==TRUE){

      t_HeatmapDF <- t(HeatmapDF)

      col_k2m_data <- kmeans(t_HeatmapDF, as.numeric(as.character(bestClustMethodcol$Cluster_Num)), nstart=25)

      # Row annotation with cluster assignment
      col_annot <- data.frame(Cluster = factor(col_k2m_data$cluster))
      rownames(col_annot) <- rownames(t_HeatmapDF)

      # Order genes by their assigned cluster
      col_cluster_order <- order(col_k2m_data$cluster)
      col_deviations_ordered <- t_HeatmapDF[col_cluster_order, ]

      deviations_ordered_col <- deviations_ordered[,rownames(col_deviations_ordered)]

      pvalues_ordered_col <- pvalues_ordered[,rownames(col_deviations_ordered)]

      # Also reorder the annotation to match
      col_annot_ordered <- col_annot[col_cluster_order, , drop = FALSE]

      # Assume `col_annot_ordered` has a column called "Cluster"
      col_annot_ordered <- col_annot_ordered %>% dplyr::rename(col_cluster = Cluster)

      col_levels <- unique(col_annot_ordered$col_cluster)

      # Generate a palette (e.g., Set1)
      colors2 <- colorRampPalette(brewer.pal(9, "Set1"))(as.numeric(as.character(bestClustMethodcol$Cluster_Num)))

      col_colors <- setNames(colors2, col_levels)

      # Pass to pheatmap
      if (isTRUE(is.data.frame(metadata))==TRUE){
        annotation_colors <- list(
          row_cluster = colors,
          subpathway = subpath_colors,
          col_cluster = col_colors
        )} else{annotation_colors <- list(
          row_cluster = colors,
          col_cluster = col_colors)}


      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        # Define custom breaks from -2 to +2 (centered at 0)
        my_breaks <- seq(minval, maxval, length.out = col_length + 1)  # must be length(colors) + 1

        comprehensive_heatmap <- pheatmap::pheatmap(deviations_ordered_col,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = FALSE,   # still cluster KOs
                                                    annotation_row = row_annot_ordered,
                                                    annotation_col = col_annot_ordered,
                                                    annotation_colors =annotation_colors,
                                                    color = my_palette,
                                                    breaks = my_breaks,
                                                    legend_breaks = c(minval, 0, maxval),
                                                    fontsize_col = 7,
                                                    fontsize_row = 5,
                                                    show_rownames = FALSE,
                                                    fontsize=6,
                                                    main = paste0(compound,  " treatment - ", genetype, " genes"))
      } else{

        comprehensive_heatmap <- pheatmap::pheatmap(deviations_ordered_col,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = FALSE,   # still cluster KOs
                                                    annotation_row = row_annot_ordered,
                                                    annotation_col = col_annot_ordered,
                                                    annotation_colors =annotation_colors,
                                                    show_rownames = FALSE,
                                                    color = my_palette,
                                                    breaks = my_breaks,
                                                    legend_breaks = c(minval, 0, maxval),
                                                    fontsize_col = 7,
                                                    fontsize_row = 3,
                                                    main = paste0(compound,  " treatment - ", genetype, " genes"))

      }


      #now produce silhouette plot for column clustering
      k2m_datacol <- factoextra::eclust(t(InputDF), "kmeans", k = as.numeric(as.character(bestClustMethodcol$Cluster_Num)), nstart = 25, graph = F)
      AvgSilPlotObjectcol <-factoextra::fviz_silhouette(k2m_datacol, palette = colors2,
                                                        ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                                          axis.title=element_text(size=14,face="bold"),
                                                                                          legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                                          plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                                          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                                          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


    } else{
      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        # Define custom breaks from -2 to +2 (centered at 0)
        my_breaks <- seq(minval, maxval,  length.out = col_length + 1)  # must be length(colors) + 1

        comprehensive_heatmap <- pheatmap::pheatmap(deviations_ordered,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = FALSE,   # still cluster KOs
                                                    annotation_row = row_annot_ordered,
                                                    annotation_colors =annotation_colors,
                                                    color = my_palette,
                                                    breaks = my_breaks,
                                                    legend_breaks = c(minval, 0, maxval),
                                                    fontsize_col = 5.5,
                                                    fontsize_row = 3,
                                                    show_rownames = TRUE,
                                                    cellwidth = 7,
                                                    cellheight = 2.7,
                                                    fontsize=5.5,
                                                    angle_col = 90,
                                                    display_numbers = pvalues_ordered,
                                                    main = paste0(compound,  " treatment - ", genetype, " genes"))


      } else{

        comprehensive_heatmap <- pheatmap::pheatmap(deviations_ordered,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = TRUE,   # still cluster KOs
                                                    annotation_row = row_annot_ordered,
                                                    annotation_colors =annotation_colors,
                                                    show_rownames = TRUE,
                                                    color = my_palette,
                                                    breaks = my_breaks,
                                                    fontsize_col = 10,
                                                    main = paste0(compound,  " treatment - ", genetype, " genes\n"))

      }
    }

    ClusterPlotObject<-factoextra::fviz_cluster(k2m_data, data = InputDF,
                                                ellipse.type = "convex",
                                                star.plot = TRUE,
                                                geom="point",
                                                ggtheme = theme_minimal(),
                                                main = "K-means Clustering"
    )  + theme(axis.text.y = element_text(size = 13),axis.text.x = element_text(size = 13),
               axis.title=element_text(size=14,face="bold"),
               legend.text = element_text(size=13),legend.title = element_text(size=14),
               plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
               axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
               axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))



    k2m_data <- factoextra::eclust(HeatmapDF, "kmeans", k = as.numeric(as.character(bestClustMethod$Cluster_Num)), nstart = 25, graph = F)
    AvgSilPlotObjectrow <-factoextra::fviz_silhouette(k2m_data, palette = colors1,
                                                      ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                                        axis.title=element_text(size=14,face="bold"),
                                                                                        legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                                        plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                                        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                                        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

  }


  ###### if best clustering method is hierarchical ########
  if(isTRUE(bestClustMethod$Method=="Hierarchical")==TRUE){

    dist_man <- dist(HeatmapDF, method="manhattan")

    hc_m2 <- hclust(d=dist_man, method="average")

    groupward6 <- cutree(hc_m2, k = as.numeric(as.character(bestClustMethod$Cluster_Num)))

    # Define color palette (e.g., 100 steps from blue to white to red)
    col_length <- 11
    my_palette <- colorRampPalette(c("blue", "white", "red"))(col_length)

    if(columnk==0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        # Define custom breaks from -2 to +2 (centered at 0)
        my_breaks <- seq(minval, maxval, length.out = col_length + 1)  # must be length(colors) + 1

        comprehensive_heatmap <- pheatmap::pheatmap(InputDF,
                                                    cluster_rows = TRUE,  # preserve k-means order
                                                    cluster_cols = FALSE,   # still cluster KOs
                                                    cutree_rows = as.numeric(as.character(bestClustMethod$Cluster_Num)),
                                                    annotation_row = row_annot_ordered,
                                                    annotation_colors = annotation_colors,
                                                    show_rownames = FALSE,
                                                    color = my_palette,
                                                    breaks = my_breaks,
                                                    fontsize_col = 10,
                                                    main = paste0(compound, " treatment\n"))

      } else{
        comprehensive_heatmap <- pheatmap::pheatmap(InputDF,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = TRUE,   # still cluster KOs
                                                    cutree_rows = as.numeric(as.character(bestClustMethod$Cluster_Num)),
                                                    annotation_row = row_annot_ordered,
                                                    annotation_colors = annotation_colors,
                                                    show_rownames = FALSE,
                                                    color = my_palette,
                                                    fontsize_col = 10,
                                                    main = paste0(compound, " treatment\n"))

      }
    }

    if(isTRUE(columnk!=0)==TRUE){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        # Define custom breaks from -2 to +2 (centered at 0)
        col_length <- 11
        my_breaks <- seq(minval, maxval, length.out = col_length + 1)  # must be length(colors) + 1

        comprehensive_heatmap <- pheatmap::pheatmap( InputDF,
                                                     show_rownames = FALSE,
                                                     cutree_rows = as.numeric(as.character(bestClustMethod$Cluster_Num)),
                                                     cutree_cols = as.numeric(as.character(bestClustMethodcol$Cluster_Num)),
                                                     color = my_palette,
                                                     breaks = my_breaks,
                                                     main = paste0(compound, " treatment\n"))

      } else{

        comprehensive_heatmap <- pheatmap::pheatmap(deviations_ordered,
                                                    cluster_rows = FALSE,  # preserve k-means order
                                                    cluster_cols = FALSE,   # still cluster KOs
                                                    annotation_row = row_annot_ordered,
                                                    annotation_colors = annotation_colors,
                                                    show_rownames = FALSE,
                                                    color = my_palette,
                                                    fontsize_col = 10,
                                                    main = paste0(compound, " treatment\n"))
      }
    }

    ClusterPlotObject<-fviz_dend(hc_m2, k = as.numeric(as.character(bestClustMethod$Cluster_Num)),
                                 cex = 0.5,
                                 color_labels_by_k = TRUE,
                                 rect = TRUE ,
                                 main="Hierarchical Clustering")


    average <- eclust(HeatmapDF, "hclust", k = as.numeric(as.character(bestClustMethod$Cluster_Num)), hc_metric = "manhattan",hc_method = "average", graph = F)
    AvgSilPlotObject<-fviz_silhouette(average, palette = "jco",
                                      ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                        axis.title=element_text(size=14,face="bold"),
                                                                        legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                        plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

    #now produce silhouette plot for column clustering
    averagecol <- eclust(HeatmapDF, "hclust", k = as.numeric(as.character(bestClustMethodcol$Cluster_Num)), hc_metric = "manhattan",hc_method = "average", graph = F)
    AvgSilPlotObjectcol <-factoextra::fviz_silhouette(averagecol, palette = "jco",
                                                      ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                                        axis.title=element_text(size=14,face="bold"),
                                                                                        legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                                        plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                                        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                                        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
  }


  if(isTRUE(columnk!=0)==TRUE){return(list(Heatmap = comprehensive_heatmap, ClusterViz = ClusterPlotObject,AvgSilrow = AvgSilPlotObjectrow,AvgSilcol = AvgSilPlotObjectcol))
  } else{return(list(Heatmap = comprehensive_heatmap, ClusterViz = ClusterPlotObject,AvgSilrow = AvgSilPlotObjectrow))}
}
