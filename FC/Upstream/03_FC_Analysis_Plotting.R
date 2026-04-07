############################################
# Flow Cytometry Analysis & Plotting 
# Author: Eleni Aretaki
# Date: 2026-04-08 (last edited)
############################################

############################################
# Phenotype-specific cell cycle distribution
############################################

compute_cellcycle_distribution <- function(
    gs,
    phenotype_pops,
    phenotype_label,
    combine_fun = c("OR", "AND")
) {
  
  combine_fun <- match.arg(combine_fun)
  
  results <- lapply(sampleNames(gs), function(sample) {
    
    gh <- gs[[sample]]
    
    pop_indices <- lapply(phenotype_pops, function(pop) {
      gh_pop_get_indices(gh, pop)
    })
    
    phenotype_idx <- Reduce(
      switch(combine_fun,
             OR  = `|`,
             AND = `&`),
      pop_indices
    )
    
    g1 <- gh_pop_get_indices(gh, "G1")
    s  <- gh_pop_get_indices(gh, "S")
    g2 <- gh_pop_get_indices(gh, "G2")
    
    total_pheno <- sum(phenotype_idx)
    
    if (total_pheno == 0) {
      return(data.frame(
        Sample = sample,
        Total = 0,
        G1 = NA,
        S  = NA,
        G2 = NA
      ))
    }
    
    data.frame(
      Sample = sample,
      Total = total_pheno,
      G1 = sum(phenotype_idx & g1) / total_pheno * 100,
      S  = sum(phenotype_idx & s)  / total_pheno * 100,
      G2 = sum(phenotype_idx & g2) / total_pheno * 100
    )
  })
  
  df <- do.call(rbind, results)
  colnames(df)[2] <- paste0(phenotype_label, "_total")
  
  return(df)
}



############################################
# Negative population distribution
############################################

compute_negative_distribution <- function(
    gs,
    positive_pops,
    phenotype_label,
    parent_pop = "Singlets"
) {
  
  results <- lapply(sampleNames(gs), function(sample) {
    
    gh <- gs[[sample]]
    
    parent_idx <- gh_pop_get_indices(gh, parent_pop)
    
    pos_indices <- lapply(positive_pops, function(pop) {
      gh_pop_get_indices(gh, pop)
    })
    
    pos_combined <- Reduce(`|`, pos_indices)
    
    phenotype_idx <- parent_idx & !pos_combined
    
    g1 <- gh_pop_get_indices(gh, "G1")
    s  <- gh_pop_get_indices(gh, "S")
    g2 <- gh_pop_get_indices(gh, "G2")
    
    total_pheno <- sum(phenotype_idx)
    
    if (total_pheno == 0) {
      return(data.frame(
        Sample = sample,
        Total = 0,
        G1 = NA,
        S  = NA,
        G2 = NA
      ))
    }
    
    data.frame(
      Sample = sample,
      Total = total_pheno,
      G1 = sum(phenotype_idx & g1) / total_pheno * 100,
      S  = sum(phenotype_idx & s)  / total_pheno * 100,
      G2 = sum(phenotype_idx & g2) / total_pheno * 100
    )
  })
  
  df <- do.call(rbind, results)
  colnames(df)[2] <- paste0(phenotype_label, "_total")
  
  return(df)
}



############################################
# Gate overview plotting
############################################

plot_gate_overview <- function(
    gs,
    x_channel,
    y_channel,
    subset_pop,
    gates_to_draw,
    gates_to_annotate = gates_to_draw,
    annotation_df = NULL,
    well_column = "well_ID",
    cellline_column = "modification",
    treatment_column = "treatment",
    bins = 200,
    ncol = 8,
    output_file,
    width = 16,
    height = 10
) {
  
  p <- ggcyto(
    gs,
    aes(
      x = !!sym(x_channel),
      y = !!sym(y_channel)
    ),
    subset = subset_pop
  ) +
    geom_hex(bins = bins) +
    ggcyto_par_set(limits = "data") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  for (g in gates_to_draw) {
    p <- p + geom_gate(g, size = 0.2)
  }
  
  for (g in gates_to_annotate) {
    p <- p + geom_stats(
      gate = g,
      type = c("gate_name", "percent"),
      adjust = 0.4,
      size = 2,
      fill = NA
    )
  }
  
  if (!is.null(annotation_df)) {
    
    annotation_df <- annotation_df %>%
      mutate(
        custom_label = paste(
          .data[[cellline_column]],
          .data[[treatment_column]],
          sep = " - "
        )
      )
    
    label_map <- setNames(
      annotation_df$custom_label,
      annotation_df[[well_column]]
    )
    
    p <- p +
      facet_wrap(
        ~well,
        ncol = ncol,
        labeller = labeller(well = label_map)
      )
    
  } else {
    
    p <- p + facet_wrap(~well, ncol = ncol)
    
  }
  
  ggsave(
    output_file,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
}