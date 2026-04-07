############################################
# Flow Cytometry Gating Functions
# Author: Eleni Aretaki
# Date: 2026-04-07 (last edited)
############################################

############################################
# Cells gate
############################################

add_cells_gate <- function(gs) {
  
  g.cells <- polygonGate(
    filterId = "Cells",
    "FSC.A" = c(0, 0, 25e4, 25e4, 25e3),
    "SSC.A" = c(5e4, 25e4, 25e4, 0, 0)
  )
  
  gs_pop_add(gs, g.cells)
  recompute(gs)
  
  return(gs)
}



############################################
# Single cell gate
############################################

add_single_cell_gates <- function(gs) {
  
  g.singlets <- polygonGate(
    filterId = "Singlets",
    "Hoechst.A" = c(2e4, 20e4, 20e4, 2e4),
    "Hoechst.H" = c(1e4, 14e4, 20e4, 6e4)
  )
  
  gs_pop_add(gs, g.singlets, parent = "Cells")
  recompute(gs)
  
  return(gs)
}



############################################
# DNA damage gates
############################################

add_dna_damage_gates <- function(gs, parp_gate1, parp_gate2, yh2ax_gate) {
  
  gate_double_neg <- polygonGate(
    filterId = "Negative",
    PARP.A = c(0, parp_gate1, parp_gate2, 0),
    yH2AX.A = c(0, 0, yh2ax_gate, yh2ax_gate)
  )
  
  gate_apoptosis_pos <- polygonGate(
    filterId = "Apoptosis-Positive",
    PARP.A = c(parp_gate1, 4, 4, parp_gate2),
    yH2AX.A = c(0, 0, yh2ax_gate, yh2ax_gate)
  )
  
  gate_dsb_pos <- polygonGate(
    filterId = "DSB-Positive",
    PARP.A = c(0, parp_gate2, parp_gate2, 0),
    yH2AX.A = c(yh2ax_gate, yh2ax_gate, 5, 5)
  )
  
  gate_double_pos <- polygonGate(
    filterId = "Double-Positive",
    PARP.A = c(parp_gate2, 4, 4, parp_gate2),
    yH2AX.A = c(yh2ax_gate, yh2ax_gate, 5, 5)
  )
  
  gs_pop_add(gs, gate_double_neg, parent = "Singlets")
  gs_pop_add(gs, gate_apoptosis_pos, parent = "Singlets")
  gs_pop_add(gs, gate_dsb_pos, parent = "Singlets")
  gs_pop_add(gs, gate_double_pos, parent = "Singlets")
  
  recompute(gs)
  
  return(gs)
}



############################################
# EdU positive gate
############################################

add_edu_positive_gate <- function(gs) {
  
  gate_edu_pos <- polygonGate(
    filterId = "S",
    list(
      Hoechst.A = c(0, 17e4, 17e4, 0),
      EdU.A = c(1.8, 2.4, Inf, Inf)
    )
  )
  
  gs_pop_add(gs, gate_edu_pos, parent = "Singlets")
  recompute(gs)
  
  return(gs)
}



############################################
# Calculate EdU-positive midpoints
############################################

calculate_midpoints_for_all_samples <- function(
    gs,
    selected_samples,
    gating_path = "S",
    output_file = NULL
) {
  
  midpoints <- list()
  
  for (sample in selected_samples) {
    
    cytoset_data <- gs_pop_get_data(gs[[sample]], gating_path)
    flowframe_data <- cytoset_data[[1]]
    
    data <- as.data.frame(exprs(flowframe_data))
    hoechst_values <- data[["Hoechst.A"]]
    
    p5  <- quantile(hoechst_values, 0.05, na.rm = TRUE)
    p95 <- quantile(hoechst_values, 0.95, na.rm = TRUE)
    
    midpoint <- (p5 + p95) / 2
    
    midpoints[[sample]] <- midpoint
  }
  
  midpoints_df <- data.frame(
    Sample = names(midpoints),
    Midpoint = unlist(midpoints)
  )
  
  if (!is.null(output_file)) {
    write_csv(midpoints_df, output_file)
    message("Sample-specific midpoints saved to: ", output_file)
  }
  
  return(midpoints_df)
}



############################################
# Dynamic G1 / G2 gates
############################################

add_cell_cycle_gates_dynamic <- function(gs, midpoint) {
  
  slope <- 0.6 / 170000
  edu_boundary <- 1.8 + slope * midpoint
  
  gate_g1_like <- polygonGate(
    filterId = "G1",
    EdU.A = c(0, 0, edu_boundary, 1.8),
    Hoechst.A = c(0, midpoint, midpoint, 0)
  )
  
  gate_g2_like <- polygonGate(
    filterId = "G2",
    EdU.A = c(0, 0, 2.4, edu_boundary),
    Hoechst.A = c(midpoint, 17e4, 17e4, midpoint)
  )
  
  gs_pop_add(gs, gate_g1_like, parent = "Singlets")
  gs_pop_add(gs, gate_g2_like, parent = "Singlets")
  
  recompute(gs)
  
  return(gs)
}