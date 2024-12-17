source("renv/activate.R")


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap", "ggplot2"), update = FALSE)

if (!requireNamespace("microViz", quietly = TRUE)) {
  remotes::install_github("david-barnett/microViz", ref = "0.12.5")
}
