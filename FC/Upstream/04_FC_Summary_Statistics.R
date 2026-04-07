############################################
# Flow Cytometry Summary Statistics
# Author: Eleni Aretaki
# Date: 2026-04-08 (last edited)
############################################

############################################
# Summary statistics export
############################################

generate_summary <- function(gs, output_file) {
  
  ps <- gs_pop_get_count_with_meta(gs) %>%
    mutate(
      percent_of_parent = Count / ParentCount
    ) %>%
    select(
      sampleName,
      well,
      Population,
      Count,
      ParentCount,
      percent_of_parent
    ) %>%
    arrange(
      factor(
        well,
        levels = str_sort(unique(well), numeric = TRUE)
      )
    )
  
  write_csv(ps, output_file)
  
  cat("Summary statistics saved to:", output_file, "\n")
}