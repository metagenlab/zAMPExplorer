
# zAMPExplorer

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe
Status](https://metagenlab.r-universe.dev/badges/zAMPExplorer)](https://metagenlab.r-universe.dev)
<!-- badges: end -->

zAMPExplorer is an interactive Shiny application designed to facilitate
the downstream analysis of outputs generated by the [zAMP
pipeline](https://zamp.readthedocs.io/en/latest/) . zAMP itself is a
comprehensive DADA2-based bioinformatics pipeline for processing 16S
rRNA gene amplicon metagenomic sequences, offering convenient,
reproducible, and scalable analysis from raw sequencing reads to
taxonomic profiles. The output of zAMP is a phyloseq object, which
serves as the input for zAMPExplorer.

A typical `phyloseq` object contains the following components:

``` r
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]
```

For more information on the [phyloseq
package](https://rdrr.io/bioc/phyloseq/man/phyloseq-package.html),
please visit the official [phyloseq
documentation](https://rdrr.io/bioc/phyloseq/man/phyloseq.html).

zAMPExplorer enables users to perform a wide range of microbiota and
statistical analyses, including compositional barplot, relative
abundance heatmap, community diversity (alpha diversity), community
similarity through unsupervised (NMDS/PCoA) and supervised (RDA)
ordinations, differential abundance testing using MaAsLin, and community
typing (or clustering) of microbial profiles using Dirichlet Multinomial
Mixtures (DMM). All of these analyses are made accessible through an
intuitive graphical interface, bridging the gap between complex
command-line bioinformatics processing and user-friendly data
exploration.

## Prerequisites and installation

### Prerequisites

- **Operating system**: Windows, macOS, or Linux
- **R**: Version 4.3.2 or later
- **RStudio**: Recommended for running the Shiny app
- **zAMP**: zAMPExplorer is designed to work with output generated from
- **microViz**: Some functionalities of zAMPExplorer depend on [microViz](https://david-barnett.github.io/microViz/). Since it is not available on CRAN, please install it before running the app.
  

## Installation

### Method 1: Install using R-universe

``` r

# Install zAMPExplorer from R-universe
install.packages('zAMPExplorer', repos = c('https://metagenlab.r-universe.dev', 'https://cloud.r-project.org'))
```

### Method 2: Install from source (GitHub)

``` r

# If not already installed 
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

#Install zAMPExplorer
remotes::install_github("metagenlab/zAMPExplorer", dependencies = TRUE)
```

### Method 3: Clone the repository and install zAMPExplorer locally

``` r
# Clone repository (in bash)
git clone https://github.com/metagenlab/zAMPExplorer.git

# Navigate to the app's directory
cd zAMPExplorer
```

``` r

# Install zAMPExplorer
install.packages(".", repos = NULL, type = "source")
```

\#Then load the library and run the app

``` r
# Load zAMPExplorer
library(zAMPExplorer)

# Launch the application
zAMPExplorer::zAMPExplorer_app()
```

## Overview of the interface

zAMPExplorer is divided into several tabs, each dedicated to a specific
type of analysis:

![Upload Data Tab](inst/figures/1-2.png) \### Check Phyloseq Components
Upload your `phyloseq` object here. This object is the output of the
zAMP pipeline and serves as the input for all downstream analyses in the
app.

<figure>
<img src="inst/figures/2.png" alt="Phyloseq Component Tab" />
<figcaption aria-hidden="true">Phyloseq Component Tab</figcaption>
</figure>

### Reads QC

Visualize the total number of reads per sample and other read quality
metrics.

<figure>
<img src="inst/figures/3.png" alt="Reads QC Tab" />
<figcaption aria-hidden="true">Reads QC Tab</figcaption>
</figure>

### Taxa Overview

Explore the estimated number of organisms in each sample at different
taxonomic levels.

![Taxa Overview Tab](inst/figures/4.png) \### Compositional Barplot
Create interactive barplots to visualize the relative abundance of taxa
within your samples.

<figure>
<img src="inst/figures/5.png" alt="Compositional Barplot Tab" />
<figcaption aria-hidden="true">Compositional Barplot Tab</figcaption>
</figure>

### Heatmap

Generate heatmaps to visualize the relative abundance of taxa across
samples and sample groups.

### Alpha Diversity

Analyze and visualize alpha diversity metrics, comparing diversity
within groups.

![Alpha Diversity Tab](inst/figures/6.png) \### Beta Diversity Explore
beta diversity using different distance matrices to assess
similarities/differences in microbial communities between samples.

### Differential abundance testing (MAaslin2)

Determining associations between microbial features (e.g., taxa) and
metadata.

![DAT Tab](inst/figures/7.png) \### Community Typing (DMM) Perform
community typing using Dirichlet Multinomial Mixture models to infer the
optimal number of community types inside the dataset.

<figure>
<img src="inst/figures/8.png" alt="DMM Tab" />
<figcaption aria-hidden="true">DMM Tab</figcaption>
</figure>

### RDA Plot

Perform redundancy analysis (RDA) to explore the association between
your samples and explanatory variables.

## Troubleshooting and FAQ

### Common issues

# Issue

App fails to launch. - Solution: Ensure all dependencies are installed
by running `source("Dependencies.R")`.

# Issue

Plots are not displaying correctly. - Solution: Verify your R version is
4.0 or later and that all necessary libraries are installed.

### FAQ

- Q: What file types can I upload?
  - A: zAMPExplorer supports `.rds` files containing `phyloseq` objects.
- Q: How do I update the app with new features?
  - A: Pull the latest updates from the GitHub repository, or fork the
    project and make changes to your fork.

## References

The following R packages are integral to the functionality of zAMP
Explorer. We highly recommend consulting their respective documentation:

## References

The following R packages are integral to the functionality of
`zAMPExplorer`. For detailed information about each package, visit the
provided links:

- **[shiny](https://CRAN.R-project.org/package=shiny)**: Shiny Web
  Application Framework for R.
- **[shinydashboard](https://CRAN.R-project.org/package=shinydashboard)**:
  Create Dashboards with ‘Shiny’.
- **[shinyWidgets](https://github.com/dreamRs/shinyWidgets)**: Custom
  Input Widgets for Shiny.
- **[webshot2](https://CRAN.R-project.org/package=webshot2)**: Taking
  Screenshots of Web Pages.
- **[htmlwidgets](https://github.com/ramnathv/htmlwidgets)**: HTML
  Widgets for R.
- **[dplyr](https://dplyr.tidyverse.org)**: A grammar of data
  manipulation.
- **[phyloseq](https://rdrr.io/bioc/phyloseq/man/phyloseq-package.html)**:
  Handling and Analysis of High-Throughput Microbial Community Data.
- **[DT](https://cran.r-project.org/web/packages/DT/index.html)**: An R
  interface to the JavaScript library DataTables.
- **[ggplot2](https://ggplot2.tidyverse.org/)**: Create Elegant Data
  Visualizations Using the Grammar of Graphics.
- **[ggpubr](https://rpkgs.datanovia.com/ggpubr/)**: Publication-Ready
  Plots Based on ggplot2.
- **[plotly](https://CRAN.R-project.org/package=plotly)**: Create
  Interactive Web Graphics via ‘plotly.js’.
- **[vegan](https://CRAN.R-project.org/package=vegan)**: Community
  Ecology Package.
- **[RColorBrewer](https://renenyffenegger.ch/notes/development/languages/R/packages/RColorBrewer/index)**:
  ColorBrewer Palettes.
- **[ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)**:
  Make Complex Heatmaps.
- **[InteractiveComplexHeatmap](https://academic.oup.com/bioinformatics/article/38/5/1460/6448211?login=false)**:
  Make Interactive Complex Heatmaps in R.
- **[microViz](https://github.com/david-barnett/microViz)**: Tools to
  make microbial community data more accessible.
- **[microbiome](https://microbiome.github.io/tutorials/)**: Microbiome
  Analytics.
- **[DirichletMultinomial](https://bioconductor.org/packages/release/bioc/html/DirichletMultinomial.html)**:
  Dirichlet-Multinomial Mixture Model Machine Learning for Microbiome
  Data.

For detailed information about each package, visit the provided links.

## Contributing

Thank you for using zAMPExplorer. We hope it helps you gain deeper
insights into your microbiome data (16S amplicons). Please feel free to
contribute, suggest features, or report any issues you encounter.

## Acknowledgments

We would like to thank the developers of the R packages and tools
integrated into zAMPExplorer. Please make sure to acknowledge their
contributions in any publications or projects using this tool.
