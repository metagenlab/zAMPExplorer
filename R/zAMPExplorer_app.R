#' @title Launch the Shiny App
#' @description This function launches the Shiny app and installs any missing dependencies.
#' @export
zAMPExplorer_app <- function() {
  # List of required packages
  packages <- c(
    "shiny", "phyloseq", "microViz", "microbiome", "mia", "miaViz",
    "MicEco", "metagMisc", "phyloseq.extended", "radiant.data",
    "ComplexHeatmap", "MicrobiotaProcess", "microbiomeMarker",
    "microbiomeutilities", "DirichletMultinomial",
    "InteractiveComplexHeatmap", "TreeSummarizedExperiment"
  )
  
  # Function to install missing packages
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      if (pkg %in% c("phyloseq", "microbiome", "ComplexHeatmap", "TreeSummarizedExperiment",
                     "InteractiveComplexHeatmap", "DirichletMultinomial", "MicrobiotaProcess", "mia", "miaViz")) {
        # Install Bioconductor packages
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        # Install CRAN packages
        install.packages(pkg)
      }
    }
  }

  # Install any missing packages
  invisible(lapply(packages, install_if_missing))

  # Locate the Shiny app directory within the package and run the app
  appDir <- system.file("shiny", package = "zAMPExplorer")
  if (appDir == "") {
    stop("Could not find the Shiny app directory. Try re-installing the `zAMPExplorer` package.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
