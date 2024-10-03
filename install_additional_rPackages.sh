#!/bin/bash
# install additional packages not available through conda
Rscript -e 'devtools::install_github("vmikk/metagMisc")'
Rscript -e 'remotes::install_github("david-barnett/microViz")'
Rscript -e 'devtools::install_github("microsud/microbiomeutilities")'
Rscript -e 'remotes::install_github("Russel88/MicEco")'
Rscript -e 'install.packages("VennDiagram", repos = "http://cran.us.r-project.org")'
Rscript -e 'remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")' 
