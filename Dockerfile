FROM mambaorg/micromamba
LABEL org.opencontainers.image.source=https://github.com/metagenlab/zAMPExplorer
LABEL org.opencontainers.image.description="Shiny app for amplicon data visualization and stats"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /pkg
RUN micromamba config set extract_threads 1 && \
    micromamba install -n base -y -f /pkg/env.yml && \
    micromamba clean -afy

ARG MAMBA_DOCKERFILE_ACTIVATE=1 
RUN R -e "install.packages('zAMPExplorer', repos = c('https://metagenlab.r-universe.dev', 'https://cloud.r-project.org'))"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh","R","-e","options(browser = 'false');library(zAMPExplorer);zAMPExplorer::zAMPExplorer_app()"]
