# Dependencies.R

# Function to check if a package is installed and install it if not
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

# List of required packages
required_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "phyloseq", "microViz",
  "microbiome", "ggplot2", "cowplot", "RColorBrewer", "circlize",
  "merTools", "DirichletMultinomial", "reshape2", "dplyr", "VennDiagram",
  "microbiomeutilities", "vegan", "gridExtra", "ggpubr", "microbiomeMarker",
  "mia", "miaViz", "TreeSummarizedExperiment", "metagMisc", "MicEco",
  "ComplexHeatmap", "plotly", "ggplotify", "InteractiveComplexHeatmap",
  "DT", "htmlwidgets", "webshot2"
)

# Install missing packages and load all
install_if_missing(required_packages)
