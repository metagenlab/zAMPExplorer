# install_dependencies.R
# This script installs all required dependencies for the zAMPExplorer package.

# Function to install a CRAN package if not already installed
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste("CRAN package already installed:", pkg))
  }
}

# Function to install a Bioconductor package if not already installed
install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message(paste("Bioconductor package already installed:", pkg))
  }
}

# Function to install a GitHub package if not already installed
install_if_missing_github <- function(repo, pkg_name) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message(paste("Installing GitHub package:", pkg_name, "from", repo))
    remotes::install_github(repo, dependencies = TRUE, quiet = TRUE)
  } else {
    message(paste("GitHub package already installed:", pkg_name))
  }
}

# Function to install microViz from R-Universe
install_microViz <- function() {
  if (!requireNamespace("microViz", quietly = TRUE)) {
    message("Installing microViz from R-Universe...")
    install.packages(
      "microViz",
      repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
    )
  } else {
    message("microViz is already installed.")
  }
}

# List of packages to install
cran_packages <- c("UpSetR", "writexl", "ggplot2")
bioc_packages <- c("Maaslin2", "MicrobiotaProcess", "phyloseq",
                   "InteractiveComplexHeatmap", "DirichletMultinomial", "microbiome")
github_packages <- list(
  list(repo = "yanlinlin82/ggvenn", pkg = "ggvenn"),
  list(repo = "ropensci/plotly", pkg = "plotly")
)

# Install CRAN packages
message("Installing CRAN packages...")
lapply(cran_packages, install_if_missing_cran)

# Install Bioconductor packages
message("Installing Bioconductor packages...")
lapply(bioc_packages, install_if_missing_bioc)

# Install GitHub packages
message("Installing GitHub packages...")
lapply(github_packages, function(pkg) {
  install_if_missing_github(pkg$repo, pkg$pkg)
})

# Install microViz from R-Universe
install_microViz()

message("All required dependencies have been installed successfully!")
