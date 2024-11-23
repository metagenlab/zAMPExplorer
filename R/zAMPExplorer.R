#' @title Launch the Shiny App
#' @description This function launches the Shiny app and installs any missing dependencies.
#' @export
zAMPExplorer <- function() {
  
  # Function to check and install missing packages
  check_and_install <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing missing package:", pkg))
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(pkg, ask = FALSE)
      }
    }
  }
  
  # Install critical Bioconductor packages
  critical_packages <- c("phyloseq", "microbiome", "Biostrings", "Maaslin2", "microbiomeMarker", "DirichletMultinomial",
"microbiome", "ComplexHeatmap", "MicrobiotaProcess", "microbiomeutilities", "DirichletMultinomial", "InteractiveComplexHeatmap",
 "TreeSummarizedExperiment")
  check_and_install(critical_packages)
  
  # Launch the Shiny app
  library(shiny)
  shinyApp(ui = ui, server = server)
}
