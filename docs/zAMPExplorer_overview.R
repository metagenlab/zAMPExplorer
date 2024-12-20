## ----include = FALSE, options(webshot.app.timeout = 120)----------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "vignettes/figures/",
  out.width = "100%"
)

## -----------------------------------------------------------------------------
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
## sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5 tips and 4 internal nodes ]
## refseq()      DNAStringSet:      [ 5 reference sequences ]

## ----install-dependencies, message=FALSE, warning=FALSE, eval=TRUE------------

# Optional: You can install dependencies using the install_dependencies.R script or install them independently

# Download the install_dependencies.R script to install required dependencies.
# download.file("https://github.com/metagenlab/zAMPExplorer/tree/latest/R/install_dependencies.R",
#               "install_dependencies.R")

# Run the script to install all dependencies
# source("install_dependencies.R")


# If not already installed 
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }

#Install zAMPExplorer
remotes::install_github("metagenlab/zAMPExplorer", ref = "latest", dependencies = TRUE)


## ----Install from R-universe, message=FALSE, warning=FALSE, eval=TRUE---------

# Install zAMPExplorer from R-universe
# install.packages('zAMPExplorer', repos = c('https://metagenlab.r-universe.dev', 'https://cloud.r-project.org'))


## ----Install locally, message=FALSE, warning=FALSE, eval=TRUE-----------------
# Clone repository (in bash)
#git clone https://github.com/metagenlab/zAMPExplorer.git

# Navigate to the app's directory
# cd zAMPExplorer


## ----install-dependencies locally, message=FALSE, warning=FALSE, eval=TRUE----
#Install dependencies in R
#source("install_dependencies.R")

# Install zAMPExplorer
#install.packages(".", repos = NULL, type = "source")


## ----load the library, message=FALSE, warning=FALSE, eval=TRUE----------------
# Load zAMPExplorer
library(zAMPExplorer)

# Launch the application
zAMPExplorer::zAMPExplorer_app()

## ----Docker installation1, message=FALSE, warning=FALSE, eval=TRUE------------

#To verify that Docker is installed, run (in bash):
#docker --version


## ----Docker installation2, message=FALSE, warning=FALSE, eval=TRUE------------

#Option1: Pull the pre-Built image. 
#You can pull the pre-built Docker image directly from a container registry (replace with your registry info if applicable):

#docker pull <your-dockerhub-username>/zampexplorer:latest

#Option 2: Build the image locally
#If cloned the zAMPExplorer source code from GitHub, navigate to the directory containing the Dockerfile and run:

#docker build --platform linux/x86_64 -t zampexplorer:latest .


## ----Docker installation3, message=FALSE, warning=FALSE, eval=TRUE------------

#docker run --rm -p 3838:3838 zampexplorer:latest

#--rm: Removes the container after it stops.
#-p 3838:3838: Maps the container’s Shiny app port 3838 to your local port 3838.


## ----Docker installation4, message=FALSE, warning=FALSE, eval=TRUE------------

#Open a web browser and navigate to:
#http://localhost:3838


## ----session_info, echo=FALSE-------------------------------------------------

sessionInfo()

