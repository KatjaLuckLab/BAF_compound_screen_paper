############################################
# Flow Cytometry Set-Up & Preprocessing Functions
# Author: Eleni Aretaki
# Date: 2026-04-07 (last edited)
############################################

############################################
# Flow Cytometry – Environment Setup
############################################

setup_environment <- function() {
  
  # Install BiocManager if missing
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # Required packages
  required_packages <- c(
    "flowCore",
    "flowStats",
    "flowViz",
    "ggcyto",
    "openCyto",
    "tidyverse",
    "gridExtra",
    "vioplot"
  )
  
  
  # Install + load packages
  for (pkg in required_packages) {
    
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg)
    
    library(pkg, character.only = TRUE)
  }
  
  
  # Additional CRAN package
  if (!requireNamespace("gtools", quietly = TRUE))
    install.packages("gtools")
  
  library(gtools)
}



############################################
# Load and preprocess FCS data
############################################

load_and_preprocess_data <- function(subfolder) {
  
  fs <- read.flowSet(
    path = subfolder,
    pattern = ".fcs",
    alter.names = TRUE
  )
  
  pData(fs)$well <- gsub(
    ".*_(.*)_.*_.*.fcs",
    "\\1",
    sampleNames(fs)
  )
  
  
  ############################################
  # Rename channels
  ############################################
  
  colnames(fs)[colnames(fs) == "UV355nm.450.50.A"] <- "Hoechst.A"
  colnames(fs)[colnames(fs) == "UV355nm.450.50.H"] <- "Hoechst.H"
  colnames(fs)[colnames(fs) == "UV355nm.450.50.W"] <- "Hoechst.W"
  
  colnames(fs)[colnames(fs) == "BL488nm.530.30.A"] <- "EdU.A"
  colnames(fs)[colnames(fs) == "BL488nm.530.30.H"] <- "EdU.H"
  
  colnames(fs)[colnames(fs) == "YG561nm.586.15.A"] <- "PARP.A"
  colnames(fs)[colnames(fs) == "YG561nm.586.15.H"] <- "PARP.H"
  
  colnames(fs)[colnames(fs) == "RL640nm.670.30.A"] <- "yH2AX.A"
  colnames(fs)[colnames(fs) == "RL640nm.670.30.H"] <- "yH2AX.H"
  
  return(fs)
}



############################################
# Transform and convert to GatingSet
############################################

transform_to_gatingset <- function(fs) {
  
  trans <- transformList(
    c("EdU.A", "PARP.A", "yH2AX.A"),
    logicleTransform()
  )
  
  fs_trans <- transform(fs, trans)
  
  gs <- GatingSet(fs_trans)
  
  return(gs)
}