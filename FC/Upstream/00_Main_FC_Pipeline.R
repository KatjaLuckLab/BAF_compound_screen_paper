############################################
# Main Flow Cytometry Analysis Workflow
# Author: Eleni Aretaki
# Date: 2026-04-08 (last edited)
############################################

#' Main Flow Cytometry Analysis Workflow
#'
#' This function orchestrates the complete workflow for analyzing flow cytometry data.
#' It includes environment setup, data loading, gating, plotting, and statistics generation.
#'
#' @details 
#' The function automates processing for multiple folders within a specified working directory. 
#' For each subfolder, it performs the following steps:
#' - Loads and preprocesses .fcs files.
#' - Applies single-cell, DNA damage, and cell cycle gates dynamically.
#' - Computes key statistics, including Hoechst.A peak values, EdU statistics, and cell counts.
#' - Generates summary CSV files and a variety of plots, such as cell cycle and DNA damage plots.
#'
#' @section Outputs:
#' - CSV files summarizing counts, and counts of extra phenotypes.
#' - PDF files containing cells, single cells, cell cycle plots, and DNA damage plots.
#' 
#' @section Requirements:
#' - R packages: `flowCore`, `flowStats`, `flowViz`, `ggcyto`, `openCyto`, `tidyverse`, `gridExtra`, `vioplot`, `gtools`.
#' - Subfolders should contain .fcs files and be named according to replicates or experimental groups.
#'
#' @return
#' This function does not return any objects. It generates and saves results to the specified directories.
#'
#' @seealso 
#' `setup_environment`, `load_and_preprocess_data`, `add_cell_cycle_gates_dynamic`, 
#' `generate_summary`, etc
#'
main <- function() {
  
  ############################################
  # 1. Setup Environment
  ############################################
  
  cat("Setting up environment...\n")
  setup_environment()
  
  ############################################
  # 2. Define paths and parameters
  ############################################
  
  subfolders <- list.dirs(working_dir, recursive = FALSE)
  annotation_df <- read.csv(annotation_file)
  
  ############################################
  # 3. Iterate through replicates
  ############################################
  
  cat("Loading and preprocessing data...\n")
  
  # Iterate through each subfolder/replicate
  for (subfolder in subfolders) {
    
    # Extract unique identifier from subfolder name
    folder_name <- basename(subfolder)
    
    # Message for tracking progress
    message("Processing folder: ", folder_name)
    
    # Retrieve the PARP.A and yH2AX.A gate based on the folder name
    gate_setting1 <- gate_settings1[[folder_name]]
    gate_setting2 <- gate_settings2[[folder_name]]
    gate_setting3 <- gate_settings3[[folder_name]]
    
    # # If the gate_setting is NA (i.e., no setting for the folder), warn and skip
    # if (is.na(gate_setting)) {
    #   warning("No gate setting specified for folder: ", folder_name)
    #   next
    # }
    
    ############################################
    # 4. Load and preprocess FCS data
    ############################################
    
    fs <- load_and_preprocess_data(subfolder)
    gs <- transform_to_gatingset(fs)
    
    # Get all sample names from the gating set
    all_samples <- sampleNames(gs)
    
    ############################################
    # 5. Apply gating hierarchy
    ############################################
    
    cat("Applying cell gate...\n")
    gs <- add_cells_gate(gs)
    
    cat("Applying single cell gate...\n")
    gs <- add_single_cell_gates(gs)
    
    cat("Applying DNA-damage gates...\n")
    gs <- add_dna_damage_gates(gs, parp_gate1 = gate_setting1, parp_gate2 = gate_setting2, yh2ax_gate = gate_setting3)
    
    cat("Applying EdU-positive gate...\n")
    gs <- add_edu_positive_gate(gs)
    
    ############################################
    # 6. Calculate midpoints for dynamic gating
    ############################################
    
    # cat("Calculating midpoints for EdU-positive cells...\n")
    midpoint_file <- file.path(subfolder, paste0(folder_name, "_midpoints.csv"))
    midpoints_df <- calculate_midpoints_for_all_samples(gs, all_samples, gating_path = "S", output_file = midpoint_file)
    message("Midpoints saved for folder: ", folder_name)
    
    ############################################
    # 7. Apply dynamic cell cycle gates
    ############################################
    
    cat("Applying dynamic cell cycle gates...\n")
    for (i in 1:nrow(midpoints_df)) {
      sample_name <- midpoints_df$Sample[i]
      midpoint <- midpoints_df$Midpoint[i]
      
      # Apply dynamic cell cycle gates for each sample
      gs[[sample_name]] <- add_cell_cycle_gates_dynamic(gs[[sample_name]], midpoint)
    }
    
    message("Gates applied for folder: ", folder_name)
    
    ############################################
    # 8. Generating statistics
    ############################################
    
    cat("Generating statistics...\n")

    generate_summary(gs,
                    file.path(subfolder, paste0(folder_name, "_countdata.csv"))
    )
    
    ############################################
    # 9. Gate overview plots
    ############################################
    
    plot_gate_overview(
      gs,
      x_channel = "FSC.A",
      y_channel = "SSC.A",
      subset_pop = "root",
      gates_to_draw = c("Cells"),
      annotation_df = annotation_df,
      output_file = file.path(subfolder,
                              paste0(folder_name, "_cells_gate.pdf"))
    )
    
    plot_gate_overview(
      gs,
      x_channel = "Hoechst.A",
      y_channel = "Hoechst.H",
      subset_pop = "Cells",
      gates_to_draw = c("Singlets"),
      annotation_df = annotation_df,
      output_file = file.path(subfolder, paste0(folder_name, "_singlets_gate.pdf"))
    )
    
    plot_gate_overview(
      gs,
      x_channel = "PARP.A",
      y_channel = "yH2AX.A",
      subset_pop = "Singlets",
      gates_to_draw = c("Negative",
                        "Apoptosis-Positive",
                        "DSB-Positive",
                        "Double-Positive"),
      annotation_df = annotation_df,
      output_file = file.path(subfolder, paste0(folder_name, "_dna_damage_gate.pdf"))
    )
    
    plot_gate_overview(
      gs,
      x_channel = "Hoechst.A",
      y_channel = "EdU.A",
      subset_pop = "Singlets",
      gates_to_draw = c("G1", "S", "G2"),
      annotation_df = annotation_df,
      output_file = file.path(subfolder,
                              paste0(folder_name, "_cell_cycle_gate.pdf"))
    )

    
    ############################################
    # 10. Phenotype-specific cell cycle distributions
    ############################################
    
    cat("Computing phenotype-specific cell cycle distributions...\n")

    parp_pos  <- compute_cellcycle_distribution(gs,
                                                c("Apoptosis-Positive", "Double-Positive"),
                                                "PARP_pos")

    yh2ax_pos <- compute_cellcycle_distribution(gs,
                                                c("DSB-Positive", "Double-Positive"),
                                                "yH2AX_pos")

    parp_neg  <- compute_negative_distribution(gs,
                                               c("Apoptosis-Positive", "Double-Positive"),
                                               "PARP_neg")

    yh2ax_neg <- compute_negative_distribution(gs,
                                               c("DSB-Positive", "Double-Positive"),
                                               "yH2AX_neg")
    
    
    ############################################
    # 11. Export phenotype distributions
    ############################################
    
    write.csv(parp_pos,
              file.path(subfolder, paste0(folder_name, "_PARP_pos_cellcycle.csv")),
              row.names = FALSE)

    write.csv(yh2ax_pos,
              file.path(subfolder, paste0(folder_name, "_yH2AX_pos_cellcycle.csv")),
              row.names = FALSE)

    write.csv(parp_neg,
              file.path(subfolder, paste0(folder_name, "_PARP_neg_cellcycle.csv")),
              row.names = FALSE)

    write.csv(yh2ax_neg,
              file.path(subfolder, paste0(folder_name, "_yH2AX_neg_cellcycle.csv")),
              row.names = FALSE)
    
  }
  
  ############################################
  # 12. Finish
  ############################################
  
  cat("\nAll folders processed. Analysis complete!\n")
}



# Execute the script
source("FC_Parameters_Exp1.R")
source("FC_Parameters_Exp2.R")
main()

